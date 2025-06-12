import scanpy as sc
import pandas as pd
import numpy as np
import os
import datetime
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# Step 7: GSEA (Gene Set Enrichment Analysis) + Visualization + HTML Report
# -----------------------------

# -----------------------------
# Define paths
# -----------------------------
input_file = "results/06_de_data.h5ad"
output_dir = "results"
gsea_dir = os.path.join(output_dir, "GSEA")
html_report = os.path.join(output_dir, "07_report.html")
log_file = os.path.join(output_dir, "07_log.txt")

# Create output directories
os.makedirs(gsea_dir, exist_ok=True)

# -----------------------------
# Initialize log
# -----------------------------
def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{timestamp}] {msg}"
    print(full_msg)
    with open(log_file, "a") as f:
        f.write(full_msg + "\n")

# -----------------------------
# Plotting functions
# -----------------------------
def plot_leading_edge_heatmap(adata, lead_genes_str, groupby, outdir, grp):
    lead_genes = lead_genes_str.split(';')
    genes_in_data = [g for g in lead_genes if g in adata.var_names]
    if not genes_in_data:
        log(f"No leading edge genes found in adata for group {grp}")
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
    outfile = os.path.join(outdir, f"heatmap_leading_edge_{groupby}_{grp}.png")
    plt.savefig(outfile, dpi=150)
    plt.close()
    log(f"Heatmap saved to {outfile}")

def plot_dotplot_top_pathways(df_gsea, outdir, groupby_col, grp):
    df_top = df_gsea.sort_values('FDR q-val').head(10).copy()
    # Convert 'Gene %' like '12/100' to numeric size
    sizes = df_top['Gene %'].str.split('/').apply(lambda x: int(x[0])/int(x[1]) if len(x)==2 else 0)
    sizes = sizes * 1000  # scale for visibility

    plt.figure(figsize=(8, 5))
    scatter = plt.scatter(
        x=df_top['NES'], y=df_top['Term'],
        s=sizes,
        c=-np.log10(df_top['FDR q-val']+1e-10), cmap='viridis', edgecolor='black'
    )
    plt.colorbar(scatter, label='-log10(FDR q-val)')
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.title(f"Top 10 Enriched Pathways Dot Plot\n{groupby_col} = {grp}")
    plt.grid(True, axis='x')
    plt.tight_layout()

    outfile = os.path.join(outdir, f"dotplot_top_pathways_{groupby_col}_{grp}.png")
    plt.savefig(outfile, dpi=150)
    plt.close()
    log(f"Dot plot saved to {outfile}")

def plot_barplot_nes(df_gsea, outdir, groupby_col, grp):
    df_top = df_gsea.sort_values('FDR q-val').head(10).copy()
    df_top = df_top[::-1]  # reverse for horizontal bar plot

    plt.figure(figsize=(8, 5))
    sns.barplot(x='NES', y='Term', data=df_top, palette='coolwarm')
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.title(f"Top 10 Enriched Pathways Bar Plot\n{groupby_col} = {grp}")
    plt.tight_layout()

    outfile = os.path.join(outdir, f"barplot_nes_{groupby_col}_{grp}.png")
    plt.savefig(outfile, dpi=150)
    plt.close()
    log(f"Bar plot saved to {outfile}")

# -----------------------------
# Load AnnData
# -----------------------------
adata = sc.read(input_file)
log(f"Loaded DE data: {adata.shape[0]} cells × {adata.shape[1]} genes")

# -----------------------------
# Settings for GSEA
# -----------------------------
groupings = ['treatment', 'leiden']
min_geneset_size = 15
max_geneset_size = 500
permutation_num = 100
seed = 42
gene_set_library = "KEGG_2016"

# Data structure to store report info
report_data = {}

# -----------------------------
# Run GSEA + Save results + Plots
# -----------------------------
for groupby_col in groupings:
    key_added = f"rank_genes_{groupby_col}"
    if key_added not in adata.uns:
        log(f"⚠️ Key '{key_added}' not found in AnnData.uns — skipping {groupby_col}")
        continue

    log(f"Starting GSEA for grouping: '{groupby_col}'")

    groups = sorted(adata.obs[groupby_col].unique())
    report_data[groupby_col] = {}

    for grp in groups:
        log(f"  - GSEA for group: '{grp}'")

        df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)
        ranked_gene_list = df.set_index('names')['logfoldchanges'].sort_values(ascending=False)

        outdir = os.path.join(gsea_dir, f"{groupby_col}_{grp}")
        os.makedirs(outdir, exist_ok=True)
        try:
            gsea_res = gp.prerank(
                rnk=ranked_gene_list,
                gene_sets=gene_set_library,
                outdir=outdir,
                format='png',
                min_size=min_geneset_size,
                max_size=max_geneset_size,
                permutation_num=permutation_num,
                seed=seed,
                verbose=True,
                no_plot=False
            )

            if gsea_res.res2d.empty:
                log(f"    ⚠️ No significant enrichment results for {groupby_col}={grp}")
            else:
                log(f"    Found {len(gsea_res.res2d)} enrichment results for {groupby_col}={grp}")

            if not gsea_res.res2d.empty:
                top_terms = gsea_res.res2d.sort_values('FDR q-val').head(3)['Term'].tolist()
                for term in top_terms:
                    try:
                        gsea_res.plot(terms=[term], ofname=os.path.join(outdir, f"manual_plot_{term}.png"))
                        log(f"    ✅ Manually generated plot for {term}")
                    except Exception as e:
                        log(f"    ❌ Failed to plot {term}: {e}")

        except Exception as e:
            log(f"    ❌ GSEA failed for {groupby_col}={grp}: {e}")
            continue

        # Read enrichment results CSV for the report and extra plots
        res_csv = os.path.join(outdir, "gseapy.gene_set.prerank.report.csv")
        if os.path.exists(res_csv):
            df_enrich = pd.read_csv(res_csv).sort_values('FDR q-val')
            if not df_enrich.empty:
                # Call new plots
                lead_genes_str = df_enrich.iloc[0]['Lead_genes'] if 'Lead_genes' in df_enrich.columns else None
                if lead_genes_str and isinstance(lead_genes_str, str) and lead_genes_str != "":
                    plot_leading_edge_heatmap(adata, lead_genes_str, groupby_col, outdir, grp)
                else:
                    log(f"    ⚠️ No leading edge genes info for {groupby_col}={grp}")

                plot_dotplot_top_pathways(df_enrich, outdir, groupby_col, grp)
                plot_barplot_nes(df_enrich, outdir, groupby_col, grp)
            else:
                log(f"    ⚠️ Empty enrichment DataFrame for {groupby_col}={grp}")

            report_data[groupby_col][grp] = {
                "csv": res_csv,
                "outdir": outdir
            }
            log(f"    ✅ GSEA results and plots saved in {outdir}")
        else:
            log(f"    ⚠️ GSEA results CSV missing in {outdir}")

# -----------------------------
# Generate HTML report
# -----------------------------
with open(html_report, "w") as f:
    f.write(f"""
<html><head><title>GSEA Enrichment Report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; }}
h1, h2, h3 {{ color: #2C3E50; }}
table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
th {{ background-color: #f2f2f2; }}
tr:nth-child(even) {{ background-color: #f9f9f9; }}
a {{ color: #2980B9; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
</style>
</head><body>
<h1>GSEA Enrichment Report</h1>
<p>Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
""")

    for groupby_col, grp_data in report_data.items():
        f.write(f"<h2>Grouping: {groupby_col}</h2>")
        for grp, info in grp_data.items():
            f.write(f"<h3>Group: {grp}</h3>")
            csv_path = info['csv']
            outdir = info['outdir']
            # Show enrichment table preview
            df_preview = pd.read_csv(csv_path).head(10)
            f.write(df_preview.to_html(index=False))
            # Embed plots
            for plotname in ["heatmap_leading_edge", "dotplot_top_pathways", "barplot_nes"]:
                img_path = os.path.join(outdir, f"{plotname}_{groupby_col}_{grp}.png")
                
                # make img_path relative to the HTML file location
                rel_img_path = os.path.relpath(img_path, start=output_dir)
                if os.path.exists(img_path):
                    f.write(f'<h4>{plotname.replace("_", " ").title()}</h4>')
                    f.write(f'<img src="{rel_img_path}" style="max-width:800px; max-height:600px;"><br><br>')

    f.write("</body></html>")

log(f"HTML report saved to {html_report}")
