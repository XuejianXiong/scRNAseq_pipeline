import scanpy as sc
import pandas as pd
import numpy as np
import os
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# Step 6: Differential Expression Analysis + Report
# -----------------------------

# -----------------------------
# User-Adjustable Parameters / Paths
# -----------------------------
INPUT_FILE = "results/05_clustered_data.h5ad"
OUTPUT_DIR = "results"
CSV_DIR = os.path.join(OUTPUT_DIR, "DE_csv")
FIG_DIR = "figures"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "06_de_data.h5ad")
HTML_REPORT = os.path.join(OUTPUT_DIR, "06_report.html")
LOG_FILE = os.path.join(OUTPUT_DIR, "06_log.txt")

# Create output directories if not exist
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# -----------------------------
# Logging function
# -----------------------------
def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{timestamp}] {msg}"
    print(full_msg)
    with open(LOG_FILE, "a") as f:
        f.write(full_msg + "\n")

# -----------------------------
# Volcano plot helper function
# -----------------------------
def plot_volcano(df, title, save_path):
    df['pvals_adj'] = df['pvals_adj'].fillna(1)
    df['log10_padj'] = -np.log10(df['pvals_adj'].replace(0, 1e-300))
    df['significant'] = (df['pvals_adj'] < 0.05) & (df['logfoldchanges'].abs() > 1)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x='logfoldchanges',
        y='log10_padj',
        hue='significant',
        palette={True: 'red', False: 'grey'},
        legend=False,
        alpha=0.7
    )
    plt.axhline(-np.log10(0.05), linestyle='--', color='blue')
    plt.axvline(1, linestyle='--', color='blue')
    plt.axvline(-1, linestyle='--', color='blue')
    plt.title(title)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-Value')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()

# -----------------------------
# Load clustered AnnData object
# -----------------------------
adata = sc.read(INPUT_FILE)
log(f"Loaded clustered data: {adata.shape[0]} cells × {adata.shape[1]} genes")

# -----------------------------
# Differential Expression Analysis
# -----------------------------
GROUPINGS = ['treatment', 'leiden']  # Group by these obs columns

for groupby_col in GROUPINGS:
    log(f"Starting DE analysis grouped by '{groupby_col}'")
    key_added = f"rank_genes_{groupby_col}"
    sc.tl.rank_genes_groups(adata, groupby=groupby_col, method='wilcoxon', key_added=key_added, pts=True)

    unique_groups = sorted(adata.obs[groupby_col].unique())
    for grp in unique_groups:
        df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)

        # Save CSV file
        csv_path = os.path.join(CSV_DIR, f"06_DE_{groupby_col}_{grp}.csv")
        df.to_csv(csv_path, index=False)
        log(f"Saved DE CSV for {groupby_col}={grp} at {csv_path}")

        # Save volcano plot
        fig_path = os.path.join(FIG_DIR, f"06_volcano_{groupby_col}_{grp}.png")
        plot_volcano(df, f"Volcano Plot: {groupby_col}={grp}", fig_path)
        log(f"Saved volcano plot for {groupby_col}={grp} at {fig_path}")

# -----------------------------
# Generate HTML Report
# -----------------------------
with open(HTML_REPORT, "w") as f:
    f.write(f"""
<html><head><title>Differential Expression Report</title></head><body>
<h1>Differential Expression Analysis Report</h1>
<p><b>Date:</b> {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
""")
    for groupby_col in GROUPINGS:
        f.write(f"<h2>DE Results by {groupby_col.capitalize()}</h2>\n")
        for grp in sorted(adata.obs[groupby_col].unique()):
            csv_name = f"DE_csv/06_DE_{groupby_col}_{grp}.csv"
            fig_name = f"../figures/06_volcano_{groupby_col}_{grp}.png"

            df_preview = pd.read_csv(os.path.join(CSV_DIR, f"06_DE_{groupby_col}_{grp}.csv")).head(10)
            html_table = df_preview.to_html(index=False, classes="table table-striped", border=0)

            f.write(f"""
<h3>{groupby_col.capitalize()}: {grp}</h3>
<p>Top 10 DE genes:</p>
{html_table}
<p>Full CSV: <a href="{csv_name}">{csv_name}</a></p>
<img src="{fig_name}" alt="Volcano plot {groupby_col} {grp}" width="600"><br><br>
""")
    f.write("</body></html>\n")

log(f"Generated HTML report at {HTML_REPORT}")

# -----------------------------
# Save updated AnnData object
# -----------------------------
adata.write(OUTPUT_FILE)
log(f"Saved DE results to {OUTPUT_FILE}")

print("✅ Step 6 complete: Differential expression analysis done and saved.")
