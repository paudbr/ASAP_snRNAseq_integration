#!/usr/bin/env python3
"""
create_samplesheet.py — Build ASAP_samplesheet.xlsx from all PMDBS team metadata.

Works with BOTH local directories and GCS paths (gs://...).
When pmdbs_v3_dir is a GCS path, reads CSVs directly via google-cloud-storage —
no local download needed. Used as Step 0 of the Nextflow pipeline on Google Batch.

Usage (manual / local):
  python create_samplesheet.py \\
      --pmdbs_v3_dir  ~/workspace/01_PMDBS/pmdbs-sc-rnaseq-v3 \\
      --pmdbs_v3_gcs  gs://asap-crn-pmdbs-sc-rnaseq-collection-v3 \\
      --output        ASAP_samplesheet.xlsx

Usage (inside pipeline — pmdbs_v3_dir IS the GCS bucket):
  create_samplesheet.py \\
      --pmdbs_v3_dir  gs://asap-crn-pmdbs-sc-rnaseq-collection-v3 \\
      --pmdbs_v3_gcs  gs://asap-crn-pmdbs-sc-rnaseq-collection-v3 \\
      --output        ASAP_samplesheet.xlsx

Output columns:
  IDENTITY        : ASAP_sample_id | sample_id | dataset_id | team | batch
  ANATOMY         : brain_region | pool_id | source_sample_id
  CLINICAL        : condition | sex | age_at_death | age_at_onset | disease_duration
  QUALITY         : pmi | post_mortem_ph
  NEUROPATHOLOGY  : braak_asyn | braak_nft | hoehn_yahr | sn_neuronal_loss
                    apoe_e4 | cognitive_status | path_ad_level
  GENETICS        : known_pd_mutation | pd_mutation
  DATA            : bam_path | matrix_path
"""

import argparse
import glob
import io
import os
import sys

import pandas as pd


# ── GCS helpers ───────────────────────────────────────────────────────────────

def _is_gcs(path: str) -> bool:
    return str(path).startswith("gs://")


def _gcs_parts(gcs_path: str) -> tuple[str, str]:
    """Parse gs://bucket/prefix/... → (bucket_name, prefix_with_trailing_slash)."""
    without = gcs_path.removeprefix("gs://")
    parts   = without.split("/", 1)
    bucket  = parts[0]
    prefix  = (parts[1].rstrip("/") + "/") if len(parts) > 1 and parts[1] else ""
    return bucket, prefix


def _gcs_client():
    try:
        from google.cloud import storage
        return storage.Client()
    except ImportError:
        sys.exit("ERROR: google-cloud-storage not installed.\n"
                 "  pip install google-cloud-storage")


def list_gcs_subdirs(gcs_root: str) -> list[tuple[str, str]]:
    """Return [(gcs_dir_path, name), ...] for immediate subdirectories."""
    bucket_name, prefix = _gcs_parts(gcs_root)
    client = _gcs_client()
    blobs  = client.list_blobs(bucket_name, prefix=prefix, delimiter="/")
    list(blobs)  # must exhaust iterator to populate .prefixes
    result = []
    for p in blobs.prefixes:
        name = p.rstrip("/").rsplit("/", 1)[-1]
        result.append((f"gs://{bucket_name}/{p}", name))
    return sorted(result, key=lambda x: x[1])


def find_gcs_csv(gcs_prefix: str, keyword: str) -> str | None:
    """First .csv/.tsv blob under gcs_prefix whose filename contains keyword."""
    bucket_name, prefix = _gcs_parts(gcs_prefix)
    client = _gcs_client()
    for blob in client.list_blobs(bucket_name, prefix=prefix):
        fname = blob.name.rsplit("/", 1)[-1]
        if keyword.lower() in fname.lower() and fname.endswith((".csv", ".tsv")):
            return f"gs://{bucket_name}/{blob.name}"
    return None


def read_gcs_table(gcs_path: str) -> pd.DataFrame:
    """Stream a CSV/TSV from GCS into a DataFrame without touching disk."""
    bucket_name, blob_name = _gcs_parts(gcs_path)
    blob_name = blob_name.rstrip("/")
    client = _gcs_client()
    content = client.bucket(bucket_name).blob(blob_name).download_as_text()
    sep = "\t" if gcs_path.endswith(".tsv") else ","
    return pd.read_csv(io.StringIO(content), sep=sep, dtype=str).rename(columns=str.strip)


# ── local helpers ─────────────────────────────────────────────────────────────

def find_local_csv(directory: str, keyword: str) -> str | None:
    for ext in ("*.csv", "*.tsv"):
        for path in sorted(glob.glob(os.path.join(directory, ext))):
            if keyword.lower() in os.path.basename(path).lower():
                return path
    return None


def read_local_table(path: str) -> pd.DataFrame:
    sep = "\t" if path.endswith(".tsv") else ","
    return pd.read_csv(path, sep=sep, dtype=str).rename(columns=str.strip)


# ── dispatch wrappers (GCS or local transparently) ────────────────────────────

def find_csv(base: str, keyword: str) -> str | None:
    return find_gcs_csv(base, keyword) if _is_gcs(base) else find_local_csv(base, keyword)


def read_table(path: str) -> pd.DataFrame:
    return read_gcs_table(path) if _is_gcs(path) else read_local_table(path)


def list_subdirs(root: str) -> list[tuple[str, str]]:
    """Return [(path, name), ...] for immediate subdirs of root."""
    if _is_gcs(root):
        return list_gcs_subdirs(root)
    root = os.path.expanduser(root)
    dirs = sorted(d for d in glob.glob(os.path.join(root, "*")) if os.path.isdir(d))
    return [(d, os.path.basename(d)) for d in dirs]


# ── normalisation maps ────────────────────────────────────────────────────────

CONDITION_MAP = {
    "idiopathic pd":        "PD",
    "pd":                   "PD",
    "parkinson":            "PD",
    "parkinson's disease":  "PD",
    "control":              "Control",
    "healthy control":      "Control",
    "normal":               "Control",
    "lbd":                  "LBD",
    "lewy body dementia":   "LBD",
    "msa":                  "MSA",
}

REGION_MAP = {
    "substantia nigra":          "SN",
    "sn":                        "SN",
    "cingulate":                 "CG",
    "cingulate cortex":          "CG",
    "anterior cingulate cortex": "CG",
    "prefrontal cortex":         "PFCTX",
    "pfc":                       "PFCTX",
    "amygdala":                  "AMG",
    "amg":                       "AMG",
    "putamen":                   "PUT",
    "put":                       "PUT",
    "middle temporal gyrus":     "MTG",
    "mtg":                       "MTG",
    "entorhinal cortex":         "ENT",
}


def norm_condition(s: pd.Series) -> pd.Series:
    return s.str.strip().str.lower().map(CONDITION_MAP).fillna(s.str.strip())


def norm_region(s: pd.Series) -> pd.Series:
    return s.str.strip().str.lower().map(REGION_MAP).fillna(s.str.strip())


def col(df: pd.DataFrame, name: str) -> pd.Series:
    if name in df.columns:
        return df[name].astype(str).str.strip().replace("nan", pd.NA)
    return pd.Series(pd.NA, index=df.index, dtype=str)


def coalesce(*series) -> pd.Series:
    result = series[0].copy()
    for s in series[1:]:
        result = result.where(result.notna() & (result.str.strip() != ""), s)
    return result


# ── per-team loader ───────────────────────────────────────────────────────────

def load_team(team_dir: str, team_name: str, pmdbs_v3_gcs: str) -> pd.DataFrame:
    sample_path   = find_csv(team_dir, "SAMPLE")
    clinpath_path = find_csv(team_dir, "CLINPATH")
    subject_path  = find_csv(team_dir, "SUBJECT")
    data_path     = find_csv(team_dir, "DATA")

    if not sample_path:
        print(f"  [{team_name}] no SAMPLE CSV — skipping")
        return pd.DataFrame()

    sample = read_table(sample_path)
    n = len(sample)

    if clinpath_path and "ASAP_subject_id" in sample.columns:
        clinpath = read_table(clinpath_path)
        keep = [c for c in [
            # identity
            "ASAP_subject_id",
            # clinical
            "age_at_death", "age_at_onset", "disease_duration_in_years",
            "duration_pmi", "primary_diagnosis", "cause_of_death",
            # neuropathology — Braak/staging
            "path_braak_asyn", "path_braak_nft", "hoehn_yahr_stage",
            "sn_neuronal_loss",
            # comorbidity / genetics
            "apoe_e4_status", "cognitive_status",
            "pathology_ad_level", "path_ad_level",   # name varies by team
            "known_pathogenic_mutation", "pd_pathogenic_mutation",
        ] if c in clinpath.columns]
        sample = sample.merge(clinpath[keep].drop_duplicates("ASAP_subject_id"),
                              on="ASAP_subject_id", how="left")

    if subject_path and "ASAP_subject_id" in sample.columns:
        subject = read_table(subject_path)
        keep = [c for c in [
            "ASAP_subject_id", "sex", "gp2_phenotype", "source_subject_id", "race",
        ] if c in subject.columns]
        sample = sample.merge(subject[keep].drop_duplicates("ASAP_subject_id"),
                              on="ASAP_subject_id", how="left")

    bam_map: dict[str, str] = {}
    if data_path:
        data   = read_table(data_path)
        ft_col = next((c for c in data.columns if "file_type" in c.lower()), None)
        bam_rows = (data[data[ft_col].str.lower().str.strip() == "bam"]
                    if ft_col else data[data.apply(
                        lambda r: any("bam" in str(v).lower() for v in r), axis=1)])
        for _, row in bam_rows.iterrows():
            uri = str(row.get("gcp_uri", "")).strip()
            sid = str(row.get("ASAP_sample_id", "")).strip()
            if uri and uri.lower() not in ("na", "nan", "") and sid:
                bam_map[sid] = uri

    out = pd.DataFrame()

    # ── identity ──────────────────────────────────────────────────────────────
    out["ASAP_sample_id"] = col(sample, "ASAP_sample_id")
    out["sample_id"]      = col(sample, "sample_id")
    ds = col(sample, "ASAP_dataset_id")
    out["dataset_id"] = ds if not ds.isna().all() else pd.Series(team_name, index=sample.index)
    out["team"]             = team_name
    out["batch"]            = col(sample, "batch")

    # ── anatomy ───────────────────────────────────────────────────────────────
    region_raw = coalesce(
        col(sample, "region_level_2"), col(sample, "region_level_1"),
        col(sample, "region"),         col(sample, "brain_region"),
    )
    out["brain_region"]     = norm_region(region_raw.fillna(""))
    out["pool_id"]          = col(sample, "pool_id")
    out["source_sample_id"] = col(sample, "source_sample_id")

    # ── clinical ──────────────────────────────────────────────────────────────
    cond_raw = coalesce(
        col(sample, "primary_diagnosis"), col(sample, "condition_id"),
        col(sample, "condition"),         col(sample, "gp2_phenotype"),
    )
    out["condition"]          = norm_condition(cond_raw.fillna(""))
    out["sex"]                = col(sample, "sex").str.upper()
    out["age_at_death"]       = col(sample, "age_at_death")
    out["age_at_onset"]       = col(sample, "age_at_onset")
    out["disease_duration"]   = col(sample, "disease_duration_in_years")
    out["cause_of_death"]     = col(sample, "cause_of_death")

    # ── sample quality ────────────────────────────────────────────────────────
    out["pmi"]                = col(sample, "duration_pmi")
    out["post_mortem_ph"]     = col(sample, "post_mortem_ph")

    # ── neuropathology & staging ──────────────────────────────────────────────
    out["braak_asyn"]         = col(sample, "path_braak_asyn")
    out["braak_nft"]          = col(sample, "path_braak_nft")
    out["hoehn_yahr"]         = col(sample, "hoehn_yahr_stage")
    out["sn_neuronal_loss"]   = col(sample, "sn_neuronal_loss")
    out["apoe_e4"]            = col(sample, "apoe_e4_status")
    out["cognitive_status"]   = col(sample, "cognitive_status")
    # pathology_ad_level column name varies slightly across teams
    out["path_ad_level"]      = coalesce(
        col(sample, "pathology_ad_level"), col(sample, "path_ad_level")
    )

    # ── genetics ──────────────────────────────────────────────────────────────
    out["known_pd_mutation"]  = col(sample, "known_pathogenic_mutation")
    out["pd_mutation"]        = col(sample, "pd_pathogenic_mutation")

    # ── data paths ────────────────────────────────────────────────────────────
    out["bam_path"]           = out["ASAP_sample_id"].map(bam_map)
    # v3 collection ships CellBender-corrected files — use them directly
    # and run the pipeline with ambient_tool = "skip"
    out["matrix_path"]        = (
        pmdbs_v3_gcs.rstrip("/") + "/" + team_name + "/"
        + out["ASAP_sample_id"] + "_1.cellbender_filtered.h5"
    )

    print(f"  [{team_name}]  {n} samples"
          + (f"  |  {len(bam_map)} BAM URIs" if bam_map else "  |  no BAMs")
          + ("  |  clinpath ✓" if clinpath_path else "")
          + ("  |  subject ✓"  if subject_path  else ""))
    return out


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build ASAP PMDBS samplesheet from metadata CSVs (local or GCS)"
    )
    parser.add_argument("--pmdbs_v3_dir", required=True,
        help="Root dir with team subdirs — local path or gs:// URI")
    parser.add_argument("--pmdbs_v3_gcs", required=True,
        help="GCS root used to construct matrix_path (gs://bucket/...)")
    parser.add_argument("--teams", default=None,
        help="Comma-separated team subdirs to include (default: all)")
    parser.add_argument("--output", default="ASAP_samplesheet.xlsx")
    args = parser.parse_args()

    root = args.pmdbs_v3_dir
    if not _is_gcs(root):
        root = os.path.expanduser(root)
        if not os.path.isdir(root):
            sys.exit(f"ERROR: directory not found: {root}")

    if args.teams:
        team_names = [t.strip() for t in args.teams.split(",")]
        sep = "/" if _is_gcs(root) else os.sep
        team_entries = [(root.rstrip("/") + "/" + t, t) for t in team_names]
    else:
        team_entries = list_subdirs(root)
        team_names   = [n for _, n in team_entries]

    print(f"\nFound {len(team_entries)} team directories under {root}:")
    for name in team_names:
        print(f"  · {name}")
    print()

    frames = []
    for tdir, tname in team_entries:
        df = load_team(tdir, tname, args.pmdbs_v3_gcs)
        if not df.empty:
            frames.append(df)

    if not frames:
        sys.exit("ERROR: no samples assembled — check input path")

    ss = pd.concat(frames, ignore_index=True)

    col_order = [
        # identity
        "ASAP_sample_id", "sample_id", "dataset_id", "team", "batch",
        # anatomy / library
        "brain_region", "pool_id", "source_sample_id",
        # clinical
        "condition", "sex", "age_at_death", "age_at_onset",
        "disease_duration", "cause_of_death",
        # sample quality
        "pmi", "post_mortem_ph",
        # neuropathology & staging
        "braak_asyn", "braak_nft", "hoehn_yahr", "sn_neuronal_loss",
        "apoe_e4", "cognitive_status", "path_ad_level",
        # genetics
        "known_pd_mutation", "pd_mutation",
        # data
        "bam_path", "matrix_path",
    ]
    for c in col_order:
        if c not in ss.columns:
            ss[c] = pd.NA
    ss = ss[col_order]

    print(f"\n{'─'*62}")
    print(f"Total samples  : {len(ss)}")
    print(f"\nBy team:")
    for t, n in ss["team"].value_counts().items():
        print(f"  {t:<42} {n:>4}")
    print(f"\nBrain regions  : {dict(ss['brain_region'].value_counts())}")
    print(f"Conditions     : {dict(ss['condition'].value_counts())}")
    print(f"Sex            : {dict(ss['sex'].value_counts())}")
    # neuropathology coverage
    for col_name in ("braak_asyn", "hoehn_yahr", "apoe_e4", "known_pd_mutation"):
        n_filled = ss[col_name].notna().sum()
        print(f"{col_name:<22}: {n_filled}/{len(ss)} filled")
    nb = ss["bam_path"].isna().sum()
    if nb:
        print(f"\n⚠  {nb} samples without BAM path")
    print(f"{'─'*62}")

    print(f"\nWriting → {args.output}")
    ss.to_excel(args.output, index=False, engine="openpyxl")
    print("Done ✓")


if __name__ == "__main__":
    main()
