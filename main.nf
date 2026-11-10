#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ══════════════════════════════════════════════════════════════════════════════
//  ASAP Parkinson's snRNA-seq Pipeline
//
//  Step 1 — PREPARE_SAMPLES    : filter samplesheet, derive paths
//  Step 2 — ambient correction : SOUPX | CELLBENDER | skip
//  Step 3 — QC_PER_SAMPLE      : scrublet, MT%, adaptive gene filter
//  Step 4 — VALIDATE_SEX       : XIST / chrY heatmap + exclusion
//  Step 5 — integration        : INTEGRATE_HARMONY | INTEGRATE_SCVI
//  Step 6 — CELLSNP + COMPARE_GENOTYPES (runs in parallel, optional)
//           cellsnp-lite pseudo-bulk genotyping + cross-sample concordance
//
//  Key parameters
//    --ambient_tool     soupx (default) | cellbender | skip
//    --integration_tool harmony (default) | scvi
//    --sex_filter       all (default) | M | F
//    --extra_filter     "col1:val1,val2;col2:val3"  (empty = none)
//    --genotype_check   true (default) | false
// ══════════════════════════════════════════════════════════════════════════════

log.info """
╔══════════════════════════════════════════════════════════════════╗
║   ASAP snRNA-seq Pipeline  —  mLOY / Parkinson's               ║
╠══════════════════════════════════════════════════════════════════╣
║  samplesheet      : ${params.samplesheet}
║  outdir           : ${params.outdir}
║  ambient_tool     : ${params.ambient_tool}
║  integration_tool : ${params.integration_tool}
║  sex_filter       : ${params.sex_filter}
║  extra_filter     : ${params.extra_filter ?: '(none)'}
║  genotype_check   : ${params.genotype_check}
║  skip_genotyping  : ${params.skip_genotyping}
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 1 — Samplesheet validation + path derivation
//  Applies hypoxia filter + any extra_filter.
//  Derives raw_matrix_path (for CellBender) from matrix_path.
// ══════════════════════════════════════════════════════════════════════════════
process PREPARE_SAMPLES {
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path samplesheet

    output:
    path "samples_valid.csv",   emit: samples_csv
    path "samples_missing.txt", emit: missing

    script:
    def ef = params.extra_filter ? "--extra_filter '${params.extra_filter}'" : ""
    """
    prepare_samples.py \\
        --samplesheet ${samplesheet} \\
        --out_csv     samples_valid.csv \\
        --out_missing samples_missing.txt \\
        ${ef}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 2a — SoupX (R, CPU)
//  Reads outs/ via load10X(), estimates contamination with autoEstCont().
//  Output: per-sample MTX directory.
// ══════════════════════════════════════════════════════════════════════════════
process SOUPX {
    tag "${sample_id}"
    maxForks params.max_parallel
    publishDir "${params.outdir}/ambient_reports", mode: 'copy', pattern: "*_report.csv"

    input:
    tuple val(sample_id), val(meta)

    output:
    tuple val(sample_id), val(meta), path("${sample_id}_soupx"), emit: corrected
    path "${sample_id}_soupx_report.csv",                        emit: report

    script:
    """
    soupx.R \\
        --sample_id   ${sample_id} \\
        --matrix_path ${meta.matrix_path} \\
        --outdir      ${sample_id}_soupx \\
        --report      ${sample_id}_soupx_report.csv
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 2b — CellBender 0.3.0 (GPU)
//  Variational autoencoder for ambient RNA removal (Pyro + PyTorch).
//  Input : raw_feature_bc_matrix (.h5 preferred, directory fallback).
//  Output: ${sample_id}_cellbender_filtered.h5  (cells only, corrected counts).
//
//  Tune cellbender_total_droplets in nextflow.config (default 25 000).
//  Set cellbender_expected_cells = 0 to let CellBender auto-detect.
// ══════════════════════════════════════════════════════════════════════════════
process CELLBENDER {
    tag "${sample_id}"
    maxForks params.max_parallel
    publishDir "${params.outdir}/ambient_reports", mode: 'copy', pattern: "*.pdf"
    publishDir "${params.outdir}/ambient_reports", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), val(meta)

    output:
    tuple val(sample_id), val(meta),
          path("${sample_id}_cellbender_filtered.h5"), emit: corrected
    path "${sample_id}_cellbender.pdf",  optional: true, emit: report
    path "${sample_id}_cellbender.log",  optional: true, emit: log

    script:
    def raw = meta.raw_matrix_path ?: ""
    if (!raw)
        error "raw_matrix_path is empty for ${sample_id} — required for CellBender. " +
              "Check that prepare_samples.py could find raw_feature_bc_matrix."

    def cuda_flag  = params.cellbender_cuda ? "--cuda" : ""
    def cells_flag = (params.cellbender_expected_cells as int) > 0
                     ? "--expected-cells ${params.cellbender_expected_cells}" : ""
    """
    # Patch torch.save/load to use cloudpickle — fixes "cannot pickle weakref.ReferenceType"
    # bug in CellBender 0.3+ with PyTorch >= 2.0 (Pyro model objects contain weak refs).
    # cloudpickle is always available in the cellbender container (pyro-ppl dependency).
    cat > cb_wrapper.py << 'PYEOF'
import sys, torch, cloudpickle as cp

_orig_save = torch.save
def _patched_save(obj, f, **kw):
    try:
        _orig_save(obj, f, **kw)
    except TypeError:
        fh = open(f, 'wb') if isinstance(f, str) else f
        cp.dump(obj, fh)
        if isinstance(f, str):
            fh.close()

_orig_load = torch.load
def _patched_load(f, **kw):
    try:
        return _orig_load(f, **kw)
    except Exception:
        fh = open(f, 'rb') if isinstance(f, str) else f
        result = cp.load(fh)
        if isinstance(f, str):
            fh.close()
        return result

torch.save = _patched_save
torch.load = _patched_load

sys.argv = ['cellbender'] + sys.argv[1:]
from cellbender.base_cli import main
sys.exit(main())
PYEOF

    python3 cb_wrapper.py remove-background \\
        --input   ${raw} \\
        --output  ${sample_id}_cellbender.h5 \\
        --total-droplets-included ${params.cellbender_total_droplets} \\
        --fpr     ${params.cellbender_fpr} \\
        ${cells_flag} \\
        ${cuda_flag}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 3 — Per-sample QC (scanpy + scrublet)
//  --corrected accepts: MTX directory (SoupX), .h5 file (CellBender), or
//  the filtered_feature_bc_matrix dir itself (skip mode).
//  Detection is automatic inside qc_per_sample.py.
// ══════════════════════════════════════════════════════════════════════════════
process QC_PER_SAMPLE {
    tag "${sample_id}"
    maxForks params.max_parallel
    publishDir "${params.outdir}/per_sample_qc", mode: 'copy', pattern: "*.h5ad"
    publishDir "${params.outdir}/per_sample_qc", mode: 'copy', pattern: "*.EXCLUDED"

    input:
    tuple val(sample_id), val(meta), path(corrected_output)

    output:
    tuple val(sample_id), path("${sample_id}.h5ad"), optional: true, emit: h5ad
    path "*.EXCLUDED",                               optional: true, emit: excluded

    script:
    def meta_json = groovy.json.JsonOutput.toJson(meta)
    """
    cat > meta.json << 'METAEOF'
${meta_json}
METAEOF

    qc_per_sample.py \\
        --sample_id  ${sample_id} \\
        --corrected  ${corrected_output} \\
        --filtered   ${meta.matrix_path} \\
        --meta_json  meta.json
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 4 — Sex annotation validation
// ══════════════════════════════════════════════════════════════════════════════
process VALIDATE_SEX {
    publishDir "${params.outdir}/sex_validation", mode: 'copy'

    input:
    path h5ad_files
    path samplesheet

    output:
    path "validated_samples.csv",              emit: valid_csv
    path "excluded_samples.csv",               emit: excluded_csv
    path "sex_corrections.csv",                emit: corrections
    path "ASAP_samplesheet_sexcorrected.xlsx", emit: corrected_xlsx, optional: true
    path "heatmap_chrY_XIST_per_patient.png",  emit: heatmap

    script:
    """
    validate_sex.py \\
        --h5ad_dir    . \\
        --samplesheet ${samplesheet} \\
        --outdir      .
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 5a — Harmony integration (CPU)
// ══════════════════════════════════════════════════════════════════════════════
process INTEGRATE_HARMONY {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path h5ad_files
    val  sex_filter

    output:
    path "ASAP_integrated.h5ad", emit: integrated
    path "plots/*",              emit: plots
    path "tables/*",             emit: tables

    script:
    """
    integrate.py \\
        --h5ad_dir   . \\
        --outdir     . \\
        --sex_filter ${sex_filter}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 6b — scVI-tools integration (GPU)
// ══════════════════════════════════════════════════════════════════════════════
process INTEGRATE_SCVI {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path h5ad_files
    val  sex_filter

    output:
    path "ASAP_integrated.h5ad", emit: integrated
    path "scvi_model/",          emit: model
    path "plots/*",              emit: plots
    path "tables/*",             emit: tables

    script:
    """
    integrate_scvi.py \\
        --h5ad_dir   . \\
        --outdir     . \\
        --sex_filter ${sex_filter}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 6a — cellsnp-lite pseudo-bulk genotyping (per sample)
//  Piles up reads from all cells in the BAM at known SNP positions (1000G VCF).
//  Outputs a pseudo-bulk VCF with AD/DP/OTH counts (not per-cell).
//  BAM and barcodes paths are derived automatically from meta.matrix_path.
// ══════════════════════════════════════════════════════════════════════════════
process CELLSNP {
    tag "${sample_id}"
    maxForks params.max_parallel
    publishDir "${params.outdir}/genotype_check", mode: 'copy', pattern: "*.vcf.gz"

    input:
    tuple val(sample_id), val(meta)

    output:
    tuple val(sample_id), path("cellSNP.${sample_id}.vcf.gz"), emit: vcf

    script:
    // Derive BAM and barcodes paths from matrix_path
    // matrix_path  = .../outs/filtered_feature_bc_matrix
    // bam_path     = .../outs/possorted_genome_bam.bam
    // barcodes     = .../outs/filtered_feature_bc_matrix/barcodes.tsv.gz
    def outs_dir = meta.matrix_path.replaceAll(/\/filtered_feature_bc_matrix.*$/, '')
    def bam_path = "${outs_dir}/possorted_genome_bam.bam"
    def barcodes = "${meta.matrix_path}/barcodes.tsv.gz"
    """
    mkdir -p snp_out

    cellsnp-lite \\
        -s  ${bam_path} \\
        -b  ${barcodes} \\
        -T  ${params.cellsnp_vcf} \\
        -O  snp_out \\
        -p  ${task.cpus} \\
        --minCOUNT ${params.cellsnp_min_count} \\
        --cellTAG CB \\
        --UMItag UB \\
        --gzip

    mv snp_out/cellSNP.base.vcf.gz cellSNP.${sample_id}.vcf.gz
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 6b — Cross-sample genotype concordance
//  Compares all within-patient pairs; flags swaps / mislabeling.
//  Uses pseudo-bulk VCFs from CELLSNP.
// ══════════════════════════════════════════════════════════════════════════════
process COMPARE_GENOTYPES {
    publishDir "${params.outdir}/genotype_check", mode: 'copy'

    input:
    path vcf_files
    path samplesheet

    output:
    path "concordance_all_pairs.tsv",        emit: all_pairs
    path "concordance_same_patient.tsv",     emit: same_patient
    path "concordance_flags.tsv",            emit: flags
    path "genotype_concordance_heatmap.pdf", emit: heatmap

    script:
    """
    compare_genotypes.py \\
        --vcf_dir     . \\
        --samplesheet ${samplesheet} \\
        --min_dp      ${params.cellsnp_min_count} \\
        --min_called  ${params.cellsnp_min_called} \\
        --outdir      .
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 6b-alt — Regenerate concordance PDF without re-running cellsnp
//  Used when skip_genotyping = true.
//  Requires concordance_all_pairs.tsv to already exist in outdir/genotype_check.
// ══════════════════════════════════════════════════════════════════════════════
process REPLOT_GENOTYPES {
    publishDir "${params.outdir}/genotype_check", mode: 'copy'

    input:
    path tsv          // concordance_all_pairs.tsv staged into work dir
    path samplesheet

    output:
    path "genotype_concordance_heatmap.pdf", emit: heatmap

    script:
    """
    compare_genotypes.py \\
        --skip_genotyping \\
        --samplesheet ${samplesheet} \\
        --min_dp      ${params.cellsnp_min_count} \\
        --min_called  ${params.cellsnp_min_called} \\
        --outdir      .
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  WORKFLOW
// ══════════════════════════════════════════════════════════════════════════════
workflow {

    // ── 1. Parse samplesheet ──────────────────────────────────────────────────
    PREPARE_SAMPLES(file(params.samplesheet))

    samples_ch = PREPARE_SAMPLES.out.samples_csv
        .splitCsv(header: true)
        .map { row ->
            def meta = row.collectEntries { k, v -> [k, v ?: ""] }
            tuple(row.sample_id, meta)
        }

    // ── 2. Ambient RNA correction — conditional on ambient_tool ───────────────
    def ambient_ch

    if (params.ambient_tool == "cellbender") {
        CELLBENDER(samples_ch)
        ambient_ch = CELLBENDER.out.corrected
    } else if (params.ambient_tool == "soupx") {
        SOUPX(samples_ch)
        ambient_ch = SOUPX.out.corrected
    } else {
        // "skip" — pass the filtered matrix directly; QC reads it as-is
        ambient_ch = samples_ch.map { sid, meta ->
            tuple(sid, meta, file(meta.matrix_path))
        }
    }

    // ── 3. QC per sample ──────────────────────────────────────────────────────
    QC_PER_SAMPLE(ambient_ch)

    // ── 4. Sex validation ─────────────────────────────────────────────────────
    h5ads_all = QC_PER_SAMPLE.out.h5ad
        .map { sid, h5ad -> h5ad }
        .collect()

    VALIDATE_SEX(h5ads_all, file(params.samplesheet))

    // ── 5. Filter h5ads to sex-validated samples ─────────────────────────────
    verified_h5ads = QC_PER_SAMPLE.out.h5ad
        .combine(VALIDATE_SEX.out.valid_csv)
        .filter { sid, h5ad, csv ->
            def ids = csv.text.split('\n').drop(1).collect { it.split(',')[0].trim() }
            sid in ids
        }
        .map { sid, h5ad, _csv -> h5ad }
        .collect()

    // ── 6. Integration ────────────────────────────────────────────────────────
    if (params.integration_tool == "scvi") {
        INTEGRATE_SCVI(verified_h5ads, params.sex_filter)
    } else {
        INTEGRATE_HARMONY(verified_h5ads, params.sex_filter)
    }

    // ── 7. Genotype check (parallel — runs alongside ambient correction) ──────
    //  CELLSNP reads directly from the BAM, no dependency on QC or ambient steps.
    if (params.genotype_check) {
        if (params.skip_genotyping) {
            // Reuse existing concordance TSV — just regenerate the PDF
            existing_tsv = file("${params.outdir}/genotype_check/concordance_all_pairs.tsv")
            if (!existing_tsv.exists()) {
                error "skip_genotyping = true but ${existing_tsv} does not exist.\n" +
                      "Run the pipeline with skip_genotyping = false first."
            }
            REPLOT_GENOTYPES(existing_tsv, file(params.samplesheet))
        } else {
            CELLSNP(samples_ch)

            vcfs_all = CELLSNP.out.vcf
                .map { sid, vcf -> vcf }
                .collect()

            COMPARE_GENOTYPES(vcfs_all, file(params.samplesheet))
        }
    }
}


workflow.onComplete {
    log.info """
    Pipeline ${workflow.success ? 'COMPLETED ✓' : 'FAILED ✗'}
      ambient    : ${params.ambient_tool}
      integrate  : ${params.integration_tool}
      sex        : ${params.sex_filter}
      genotype   : ${params.genotype_check}  skip_genotyping=${params.skip_genotyping}
      duration   : ${workflow.duration}
      output     : ${params.outdir}
    """.stripIndent()
}
