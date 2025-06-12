import scanpy as sc
import pandas as pd
import numpy as np
import os
import datetime
import gseapy as gp

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

            print(gsea_res.res2d.columns)

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

        # Read enrichment results CSV for the report
        res_csv = os.path.join(outdir, "gseapy.gene_set.prerank.report.csv")
        if os.path.exists(res_csv):
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
.img-container {{ margin-bottom: 30px; }}
</style>
</head><body>
<h1>Gene Set Enrichment Analysis (GSEA) Report</h1>
<p><b>Date:</b> {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
""")

    for groupby_col, groups_dict in report_data.items():
        f.write(f"<h2>GSEA Results by {groupby_col.capitalize()}</h2>\n")
        for grp, paths in groups_dict.items():
            csv_rel = os.path.relpath(paths["csv"], start=os.path.dirname(html_report))
            outdir_rel = os.path.relpath(paths["outdir"], start=os.path.dirname(html_report))

            # Load top 10 enriched gene sets for preview
            df_enrich = pd.read_csv(paths["csv"]).sort_values('FDR q-val').head(10)
            html_table = df_enrich.to_html(index=False, classes="table table-striped", border=0)

            # Try to find a representative enrichment plot image
            # The typical gseapy plot files are named like "gseapy.prerank.<geneset_name>.png"
            import glob
            img_files = glob.glob(os.path.join(paths["outdir"], "gseapy.prerank.*.png"))
            # Pick top 3 plots to display (if any)
            top_imgs = img_files[:3]

            f.write(f"""
            <h3>{groupby_col.capitalize()} = {grp}</h3>
            <p>Top 10 enriched gene sets (FDR q-value sorted):</p>
            {html_table}
            <p>Full CSV: <a href="{csv_rel}">{csv_rel}</a></p>
            """)

            if top_imgs:
                f.write("<div class='img-container'><p>Representative enrichment plots:</p>")
                for img_path in top_imgs:
                    img_rel = os.path.relpath(img_path, start=os.path.dirname(html_report))
                    f.write(f'<img src="{img_rel}" alt="Enrichment plot" width="400" style="margin-right:20px;">')
                f.write("</div>")
            else:
                f.write("<p><i>No enrichment plots found.</i></p>")

    f.write("</body></html>")

log(f"Generated HTML report at {html_report}")

print("✅ Step 7 complete: GSEA enrichment analysis + visualization + report done and saved.")
