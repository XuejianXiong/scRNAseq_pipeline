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
# Define paths
# -----------------------------
input_file = "results/05_clustered_data.h5ad"
output_dir = "results"
csv_dir = "results/DE_csv"
fig_dir = "figures"

# Create output directories
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)
output_file = os.path.join(output_dir, "06_de_data.h5ad")
html_report = os.path.join(output_dir, "06_report.html")
log_file = os.path.join(output_dir, "06_log.txt")


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
# Define plot_volcano function
# -----------------------------
def plot_volcano(df, title, save_path):
    # Prepare data for volcano plot
    df['pvals_adj'] = df['pvals_adj'].fillna(1)
    df['log10_padj'] = -np.log10(df['pvals_adj'].replace(0, 1e-300))
    df['significant'] = (df['pvals_adj'] < 0.05) & (df['logfoldchanges'].abs() > 1)

    plt.figure(figsize=(8,6))
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
    
    return None


# -----------------------------
# Load clustered AnnData
# -----------------------------
adata = sc.read(input_file)
log(f"Loaded integrated data: {adata.shape[0]} cells × {adata.shape[1]} genes")

# -----------------------------
# Differential Expression
# -----------------------------
# Columns for DE grouping
#groupings = ['treatment', 'leiden', 'cluster']
groupings = ['treatment', 'leiden']

for groupby_col in groupings:
    log(f"Starting DE analysis grouped by '{groupby_col}'")
    key_added = f"rank_genes_{groupby_col}"
    sc.tl.rank_genes_groups(adata, groupby=groupby_col, method='wilcoxon', key_added=key_added, pts=True)

    unique_groups = sorted(adata.obs[groupby_col].unique())
    for grp in unique_groups:
        # Get DE results for this group
        df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)
        
        # Save CSV
        csv_path = os.path.join(csv_dir, f"06_DE_{groupby_col}_{grp}.csv")
        df.to_csv(csv_path, index=False)
        log(f"Saved DE CSV for {groupby_col}={grp} at {csv_path}")

        # Plot volcano
        fig_path = os.path.join(fig_dir, f"06_volcano_{groupby_col}_{grp}.png")
        plot_volcano(df, f"Volcano Plot: {groupby_col}={grp}", fig_path)
        log(f"Saved volcano plot for {groupby_col}={grp} at {fig_path}")

#print((adata.obs["leiden"] == adata.obs["cluster"]).all())


# -----------------------------
# Generate HTML Report
# -----------------------------
with open(html_report, "w") as f:
    f.write(f"""
    <html><head><title>Differential Expression Report</title></head><body>
    <h1>Differential Expression Analysis Report</h1>
    <p><b>Date:</b> {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    """)

    for groupby_col in groupings:
        f.write(f"<h2>DE Results by {groupby_col.capitalize()}</h2>\n")
        for grp in sorted(adata.obs[groupby_col].unique()):
            csv_name = f"DE_csv/06_DE_{groupby_col}_{grp}.csv"
            fig_name = f"../figures/06_volcano_{groupby_col}_{grp}.png"
            
            # Read top 10 for inline table preview
            df_preview = pd.read_csv(os.path.join(csv_dir, f"06_DE_{groupby_col}_{grp}.csv")).head(10)
            html_table = df_preview.to_html(index=False, classes="table table-striped", border=0)
            
            f.write(f"""
            <h3>{groupby_col.capitalize()}: {grp}</h3>
            <p>Top 10 DE genes:</p>
            {html_table}
            <p>Full CSV: <a href="{csv_name}">{csv_name}</a></p>
            <img src="{fig_name}" alt="Volcano plot {groupby_col} {grp}" width="600"><br><br>
            """)
    f.write("</body></html>\n")

log(f"Generated HTML report at {html_report}")

# -----------------------------
# Save updated AnnData object
# -----------------------------
adata.write(output_file)
log("Saved DE results to results/06_de_data.h5ad")

print("✅ Step 6 complete: Differential expression analysis done and saved.")

