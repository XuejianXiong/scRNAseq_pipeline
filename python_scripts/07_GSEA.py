#!/usr/bin/env python3
"""
Step 7: GSEA Analysis + Visualization + HTML Report for scRNA-seq data.

This step:
- Runs GSEA for each group (from DE results)
- Generates heatmaps of leading edge genes
- Generates dot plots and bar plots of top pathways
- Produces an HTML summary report
"""

import sys
import warnings
from pathlib import Path
import datetime
import os

import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
from logzero import logger
import yaml

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action="ignore", category=FutureWarning)
sns.set(style="whitegrid")

# -----------------------------
# Load Config and Parameters
# -----------------------------
if len(sys.argv) < 3:
    print("Usage: step7_gsea.py <project_name> <config.yaml>")
    sys.exit(1)

PROJ_NAME = sys.argv[1]
CONFIG_FILE = Path(sys.argv[2])

if not CONFIG_FILE.exists():
    print(f"Config file not found: {CONFIG_FILE}")
    sys.exit(1)

with open(CONFIG_FILE, "r") as f:
    cfg = yaml.safe_load(f)

try:
    GSEA_PARAMS = cfg["gsea"]
    GENE_SET_LIBRARY = GSEA_PARAMS["gene_set_library"]
    MIN_GENESET_SIZE = GSEA_PARAMS["min_geneset_size"]
    MAX_GENESET_SIZE = GSEA_PARAMS["max_geneset_size"]
    PERMUTATION_NUM = GSEA_PARAMS["permutation_num"]
    SEED = GSEA_PARAMS["seed"]
except KeyError as e:
    print(f"Missing GSEA parameter in config.yaml: {e}")
    sys.exit(1)

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("07_log.txt", PROJ_NAME)
GSEA_DIR = Path(f"{RESULTS_DIR}/GSEA")
GSEA_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = RESULTS_DIR / "06_de_data.h5ad"
HTML_REPORT = RESULTS_DIR / "07_report.html"

sc.settings.verbosity = 1
sc.settings.logfile = LOG_FILE


# -----------------------------
# Helper Functions
# -----------------------------
def plot_leading_edge_heatmap(
    adata: sc.AnnData, lead_genes_str: str, groupby: str, outdir: Path, grp: str
) -> None:
    """Plot a heatmap of leading edge genes for a given group."""
    lead_genes = lead_genes_str.split(";")
    genes_in_data = [g for g in lead_genes if g in adata.var_names]
    if not genes_in_data:
        logger.warning(f"No leading edge genes found in adata for group {grp}")
        return
    expr = adata[:, genes_in_data].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    df_expr = pd.DataFrame(expr, columns=genes_in_data, index=adata.obs_names)
    df_expr["group"] = adata.obs[groupby]
    df_expr = df_expr.sort_values("group")
    df_expr_numeric = df_expr.drop(columns="group")

    plt.figure(figsize=(10, 6))
    sns.heatmap(df_expr_numeric.T, cmap="vlag", center=0, xticklabels=False)
    plt.title(f"Leading Edge Genes Heatmap: {groupby} = {grp}")
    plt.ylabel("Genes")
    plt.xlabel("Cells (sorted by group)")
    plt.tight_layout()
    outfile = outdir / f"heatmap_leading_edge_{groupby}_{grp}.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    logger.info(f"Saved heatmap: {outfile}")


def plot_dotplot_top_pathways(
    df_gsea: pd.DataFrame, outdir: Path, groupby_col: str, grp: str
) -> None:
    """Create a dot plot of the top enriched pathways."""
    df_top = df_gsea.sort_values("FDR q-val").head(10).copy()
    sizes = (
        df_top["Gene %"]
        .str.split("/")
        .apply(lambda x: int(x[0]) / int(x[1]) if len(x) == 2 else 0)
    )
    sizes = sizes * 1000

    plt.figure(figsize=(8, 5))
    scatter = plt.scatter(
        x=df_top["NES"],
        y=df_top["Term"],
        s=sizes,
        c=-np.log10(df_top["FDR q-val"] + 1e-10),
        cmap="viridis",
        edgecolor="black",
    )
    plt.colorbar(scatter, label="-log10(FDR q-val)")
    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.title(f"Top 10 Enriched Pathways Dot Plot\n{groupby_col} = {grp}")
    plt.tight_layout()
    outfile = outdir / f"dotplot_top_pathways_{groupby_col}_{grp}.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    logger.info(f"Saved dot plot: {outfile}")


def plot_barplot_nes(
    df_gsea: pd.DataFrame, outdir: Path, groupby_col: str, grp: str
) -> None:
    """Create a bar plot of normalized enrichment scores (NES)."""
    df_top = df_gsea.sort_values("FDR q-val").head(10).iloc[::-1]
    plt.figure(figsize=(8, 5))
    sns.barplot(x="NES", y="Term", data=df_top, palette="coolwarm")
    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.title(f"Top 10 Enriched Pathways Bar Plot\n{groupby_col} = {grp}")
    plt.tight_layout()
    outfile = outdir / f"barplot_nes_{groupby_col}_{grp}.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    logger.info(f"Saved bar plot: {outfile}")


def run_gsea_for_group(
    adata: sc.AnnData, groupby_col: str, grp: str, gsea_dir: Path
) -> dict:
    """Run GSEA for a single group and generate plots."""
    key_added = f"rank_genes_{groupby_col}"
    if key_added not in adata.uns:
        logger.warning(f"Key '{key_added}' not found — skipping {groupby_col} = {grp}")
        return {}

    df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)
    ranked_gene_list = df.set_index("names")["logfoldchanges"].sort_values(
        ascending=False
    )

    outdir = gsea_dir / f"{groupby_col}_{grp}"
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        gsea_res = gp.prerank(
            rnk=ranked_gene_list,
            gene_sets=GENE_SET_LIBRARY,
            outdir=str(outdir),
            format="png",
            min_size=MIN_GENESET_SIZE,
            max_size=MAX_GENESET_SIZE,
            permutation_num=PERMUTATION_NUM,
            seed=SEED,
            verbose=True,
            no_plot=False,
        )
    except Exception as e:
        logger.error(f"GSEA failed for {groupby_col} = {grp}: {e}")
        return {}

    res_csv = outdir / "gseapy.gene_set.prerank.report.csv"
    if res_csv.exists():
        df_enrich = pd.read_csv(res_csv).sort_values("FDR q-val")
        if not df_enrich.empty:
            lead_genes_str = df_enrich.iloc[0].get("Lead_genes", "")
            if isinstance(lead_genes_str, str) and lead_genes_str:
                plot_leading_edge_heatmap(
                    adata, lead_genes_str, groupby_col, outdir, grp
                )
            plot_dotplot_top_pathways(df_enrich, outdir, groupby_col, grp)
            plot_barplot_nes(df_enrich, outdir, groupby_col, grp)
            logger.info(f"Saved GSEA results for {groupby_col} = {grp}")
            return {"csv": res_csv, "outdir": outdir}
    return {}


def generate_html_report(
    report_data: dict, output_file: Path, results_dir: Path
) -> None:
    """Generate an HTML report summarizing GSEA results."""
    with open(output_file, "w") as f:
        f.write(f"<html><head><title>Step 7: GSEA Report</title></head><body>\n")
        f.write(f"<h1>Step 7 Report: GSEA and Visualization</h1>\n")
        f.write(
            f"<p>Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n"
        )

        for groupby_col, grp_data in report_data.items():
            f.write(f"<h2>Grouping: {groupby_col}</h2>\n")
            for grp, info in grp_data.items():
                f.write(f"<h3>Group: {grp}</h3>\n")
                df_preview = pd.read_csv(info["csv"]).head(10)
                f.write(df_preview.to_html(index=False))
                for plot_name in [
                    "heatmap_leading_edge",
                    "dotplot_top_pathways",
                    "barplot_nes",
                ]:
                    img_path = (
                        Path(info["outdir"]) / f"{plot_name}_{groupby_col}_{grp}.png"
                    )
                    rel_img_path = os.path.relpath(img_path, start=results_dir)
                    if img_path.exists():
                        f.write(f"<h4>{plot_name.replace('_',' ').title()}</h4>\n")
                        f.write(f'<img src="{rel_img_path}" width="700"><br><br>\n')

        f.write("</body></html>\n")
    logger.info(f"Generated HTML report at {output_file}")


# -----------------------------
# Main Pipeline
# -----------------------------
def main() -> None:
    logger.info("Step 7 started: GSEA and Visualization")
    adata = sc.read(INPUT_FILE)
    logger.info(f"Loaded DE data: {adata.n_obs} cells × {adata.n_vars} genes")

    if "treatment" in adata.obs.columns:
        GROUPINGS = ["treatment", "leiden"]
    elif "sample" in adata.obs.columns:
        GROUPINGS = ["sample", "leiden"]
    else:
        logger.error("No 'treatment' or 'sample' column found — skipping GSEA")
        sys.exit(1)

    report_data = {}
    for groupby_col in GROUPINGS:
        report_data[groupby_col] = {}
        for grp in sorted(adata.obs[groupby_col].unique()):
            result = run_gsea_for_group(adata, groupby_col, grp, GSEA_DIR)
            if result:
                report_data[groupby_col][grp] = result

    generate_html_report(report_data, HTML_REPORT, RESULTS_DIR)

    logger.info("✅ Step 7 complete: GSEA analysis and report generated.")


if __name__ == "__main__":
    main()
