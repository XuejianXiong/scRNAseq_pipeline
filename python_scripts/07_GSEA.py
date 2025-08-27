# -----------------------------
# Step 7: GSEA + Visualization + HTML Report
# -----------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import os
import datetime
from logzero import logger
from pathlib import Path
import warnings
import sys

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# Plot Functions
# -----------------------------
def plot_leading_edge_heatmap(
        adata: sc.AnnData, 
        lead_genes_str: str, 
        groupby: str, 
        outdir: str, 
        grp: str
    ) -> None:
    lead_genes = lead_genes_str.split(';')
    genes_in_data = [g for g in lead_genes if g in adata.var_names]
    if not genes_in_data:
        logger.warning(f"No leading edge genes found in adata for group {grp}")
        return
    expr = adata[:, genes_in_data].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    df_expr = pd.DataFrame(expr, columns=genes_in_data, index=adata.obs_names)
    df_expr['group'] = adata.obs[groupby]
    df_expr = df_expr.sort_values('group')
    df_expr_numeric = df_expr.drop(columns='group')

    plt.figure(figsize=(10, 6))
    sns.heatmap(df_expr_numeric.T, cmap='vlag', center=0, xticklabels=False)
    plt.title(f"Leading Edge Genes Heatmap: {groupby} = {grp}")
    plt.ylabel("Genes")
    plt.xlabel("Cells (sorted by group)")
    plt.tight_layout()
    outfile = f"{outdir}/heatmap_leading_edge_{groupby}_{grp}.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    logger.info(f"Saved heatmap: {outfile}")

def plot_dotplot_top_pathways(
        df_gsea: pd.DataFrame, outdir: str, groupby_col: str, grp: str
    ) -> None:
    df_top = df_gsea.sort_values('FDR q-val').head(10).copy()
    sizes = df_top['Gene %'].str.split('/').apply(lambda x: int(x[0])/int(x[1]) if len(x) == 2 else 0)
    sizes = sizes * 1000

    plt.figure(figsize=(8, 5))
    scatter = plt.scatter(
        x=df_top['NES'], y=df_top['Term'],
        s=sizes,
        c=-np.log10(df_top['FDR q-val'] + 1e-10),
        cmap='viridis',
        edgecolor='black'
    )
    plt.colorbar(scatter, label='-log10(FDR q-val)')
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.title(f"Top 10 Enriched Pathways Dot Plot\n{groupby_col} = {grp}")
    plt.tight_layout()
    outfile = f"{outdir}/dotplot_top_pathways_{groupby_col}_{grp}.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    logger.info(f"Saved dot plot: {outfile}")

def plot_barplot_nes(
        df_gsea: pd.DataFrame, outdir: str, groupby_col: str, grp: str
    ) -> None:
    df_top = df_gsea.sort_values('FDR q-val').head(10).iloc[::-1]
    plt.figure(figsize=(8, 5))
    sns.barplot(x='NES', y='Term', data=df_top, palette='coolwarm')
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.title(f"Top 10 Enriched Pathways Bar Plot\n{groupby_col} = {grp}")
    plt.tight_layout()
    outfile = f"{outdir}/barplot_nes_{groupby_col}_{grp}.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    logger.info(f"Saved bar plot: {outfile}")


# -----------------------------
# User-Adjustable Parameters
# -----------------------------
PROJ_NAME = "cropseq"  
#PROJ_NAME = "retina"  

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("07_log.txt", PROJ_NAME)
GSEA_DIR = Path(f"{RESULTS_DIR}/GSEA")
GSEA_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = RESULTS_DIR / "06_de_data.h5ad"
HTML_REPORT = RESULTS_DIR / "07_report.html"

GENE_SET_LIBRARY = "KEGG_2016"
MIN_GENESET_SIZE = 15
MAX_GENESET_SIZE = 500
PERMUTATION_NUM = 100
SEED = 42

# -----------------------------


logger.info("Step 7 started: GSEA and Visualization")
# -----------------------------
# Load data
# -----------------------------
adata = sc.read(INPUT_FILE)
logger.info(f"Loaded DE data: {adata.n_obs} cells × {adata.n_vars} genes")

if "treatment" in adata.obs.columns:
    GROUPINGS = ['treatment', 'leiden']  # Grouping columns for GSEA
elif "sample" in adata.obs.columns:
    GROUPINGS = ['sample', 'leiden']  # Grouping columns for GSEA
else:
    logger.error(
        "Metadata 'treatment' and 'sample' not found in adata.obs, skipping GSEA plot."
    )
    sys.exit(1)  # Stop execution immediately

# -----------------------------
# Run GSEA
# -----------------------------
report_data = {}

for groupby_col in GROUPINGS:
    key_added = f"rank_genes_{groupby_col}"
    if key_added not in adata.uns:
        logger.warning(f"Key '{key_added}' not found in AnnData.uns — skipping {groupby_col}")
        continue

    groups = sorted(adata.obs[groupby_col].unique())
    report_data[groupby_col] = {}

    for grp in groups:
        logger.info(f"Running GSEA for {groupby_col} = {grp}")
        df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)
        ranked_gene_list = df.set_index('names')['logfoldchanges'].sort_values(ascending=False)

        outdir = f"{GSEA_DIR}/{groupby_col}_{grp}"
        Path(outdir).mkdir(parents=True, exist_ok=True)

        try:
            gsea_res = gp.prerank(
                rnk=ranked_gene_list,
                gene_sets=GENE_SET_LIBRARY,
                outdir=outdir,
                format='png',
                min_size=MIN_GENESET_SIZE,
                max_size=MAX_GENESET_SIZE,
                permutation_num=PERMUTATION_NUM,
                seed=SEED,
                verbose=True,
                no_plot=False
            )
        except Exception as e:
            logger.error(f"GSEA failed for {groupby_col} = {grp}: {e}")
            continue

        if gsea_res.res2d.empty:
            logger.warning(f"No significant enrichment for {groupby_col} = {grp}")
            continue

        res_csv = f"{outdir}/gseapy.gene_set.prerank.report.csv"
        if os.path.exists(res_csv):
            df_enrich = pd.read_csv(res_csv).sort_values('FDR q-val')
            if not df_enrich.empty:
                lead_genes_str = df_enrich.iloc[0].get('Lead_genes', "")
                if isinstance(lead_genes_str, str) and lead_genes_str:
                    plot_leading_edge_heatmap(adata, lead_genes_str, groupby_col, outdir, grp)
                else:
                    logger.warning(f"No leading edge genes info for {groupby_col} = {grp}")

                plot_dotplot_top_pathways(df_enrich, outdir, groupby_col, grp)
                plot_barplot_nes(df_enrich, outdir, groupby_col, grp)

                report_data[groupby_col][grp] = {
                    "csv": res_csv,
                    "outdir": outdir
                }
                logger.info(f"Saved GSEA results for {groupby_col} = {grp}")
            else:
                logger.warning(f"Empty enrichment DataFrame for {groupby_col} = {grp}")
        else:
            logger.warning(f"GSEA result CSV missing for {groupby_col} = {grp}")

# -----------------------------
# Generate HTML Report
# -----------------------------
with open(HTML_REPORT, "w") as f:
    f.write(f"""
<html><head><title>Step 7: GSEA Report</title></head><body>
<h1>Step 7 Report: GSEA and Visualization</h1>
<p>Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
""")

    for groupby_col, grp_data in report_data.items():
        f.write(f"<h2>Grouping: {groupby_col}</h2>\n")
        for grp, info in grp_data.items():
            f.write(f"<h3>Group: {grp}</h3>\n")
            df_preview = pd.read_csv(info['csv']).head(10)
            f.write(df_preview.to_html(index=False))
            for plot_name in ["heatmap_leading_edge", "dotplot_top_pathways", "barplot_nes"]:
                img_path = f"{info['outdir']}/{plot_name}_{groupby_col}_{grp}.png"
                rel_img_path = os.path.relpath(img_path, start=RESULTS_DIR)
                if os.path.exists(img_path):
                    f.write(f'<h4>{plot_name.replace("_", " ").title()}</h4>\n')
                    f.write(f'<img src="{rel_img_path}" width="700"><br><br>\n')

    f.write("</body></html>")
logger.info(f"Generated HTML report at {HTML_REPORT}")


print("✅ Step 7 complete: GSEA analysis and report generated.")

