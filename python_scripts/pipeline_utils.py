from pathlib import Path
import logging
import scanpy as sc

def setup_dirs_logs(log_name: str, project: str = None) -> tuple[Path, Path, Path]:
    """
    Initialize project directories, configure logging, and set Scanpy plotting parameters.

    This utility function ensures reproducibility and organization in bioinformatics pipelines 
    by creating dedicated results and figures directories, setting up logging to a timestamped 
    log file, and configuring Scanpy plotting behavior.

    Parameters
    ----------
    log_name : str
        Name of the log file (e.g., '02_preprocessing.log').
    project : str, optional
        Project or dataset name (e.g., 'cropseq', 'retina'). Directories will be created under
        'results/{project}' and 'figures/{project}'. If None, top-level directories are used.

    Returns
    -------
    RESULTS_DIR : Path
        Path to the results directory.
    FIGURES_DIR : Path
        Path to the figures directory.
    LOG_FILE : Path
        Path to the log file.
    """
    # Define directories
    RESULTS_DIR = Path("results") / project if project else Path("results")
    FIGURES_DIR = Path("figures") / project if project else Path("figures")
    
    # Create directories if they do not exist
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Define log file path
    LOG_FILE = RESULTS_DIR / log_name

    # Reset any existing logging handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Configure logging
    logging.basicConfig(
        filename=LOG_FILE,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    # Configure Scanpy plotting
    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURES_DIR)
    sc.settings.verbosity = 2

    logging.info(f"Directories initialized: results -> {RESULTS_DIR}, figures -> {FIGURES_DIR}")
    logging.info(f"Logging configured: {LOG_FILE}")

    return RESULTS_DIR, FIGURES_DIR, LOG_FILE
