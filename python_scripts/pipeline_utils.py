from pathlib import Path
from logzero import logger, logfile
import scanpy as sc


def setup_dirs_logs(log_name: str, project: str) -> tuple[Path, Path, Path]:
    """
    Initialize project directories, configure logging, and set Scanpy plotting parameters.

    Ensures reproducibility and organization by creating dedicated data, results, and figures
    directories, configuring logging to a project-specific log file, and setting Scanpy plotting behavior.

    Parameters
    ----------
    log_name : str
        Name of the log file (e.g., '02_preprocessing.log').
    project : str
        Project or dataset name (e.g., 'cropseq', 'retina'). Directories will be created under
        'data/{project}', 'results/{project}', and 'figures/{project}'.

    Returns
    -------
    RESULTS_DIR : Path
        Path to the results directory.
    FIGURES_DIR : Path
        Path to the figures directory.
    LOG_FILE : Path
        Path to the log file.

    Raises
    ------
    ValueError
        If project is None or an empty string.
    """
    # Validate project
    if not project or not project.strip():
        raise ValueError("The 'project' parameter is required and cannot be empty.")

    # Define directories
    RESULTS_DIR = Path("results") / project
    FIGURES_DIR = Path("figures") / project

    # Create directories if they do not exist
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Define log file path
    LOG_FILE = RESULTS_DIR / log_name

    # Configure logzero to log to both console and file
    logfile(LOG_FILE, mode="w")  # overwrite log file each run

    # Log directory setup
    logger.info(
        f"Directories initialized: results -> {RESULTS_DIR}, figures -> {FIGURES_DIR}"
    )
    logger.info(f"Logging configured to write both console and file: {LOG_FILE}")

    # Configure Scanpy plotting
    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURES_DIR)
    sc.settings.verbosity = 2

    return RESULTS_DIR, FIGURES_DIR, LOG_FILE
