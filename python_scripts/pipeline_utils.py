from pathlib import Path
import logging
import scanpy as sc

def setup_dirs_logs(log_name: str) -> tuple[Path, Path, str]:
    """
    Set up pipeline directories, logging, and Scanpy settings.

    This function creates the 'results' and 'figures' directories (if they 
    do not exist), configures logging to write to the specified log file 
    within the results directory, and sets Scanpy's plotting parameters.

    Parameters
    ----------
    log_name : str
        The filename for the log file (e.g., '02_log.txt').

    Returns
    -------
    tuple of Path
        A tuple containing:
        - RESULTS_DIR : Path
            Path to the 'results' directory.
        - FIGURE_DIR : Path
            Path to the 'figures' directory.
        - LOG_FILE : str
            Path of the log file.
    """
    RESULTS_DIR = Path("results")
    FIGURE_DIR = Path("figures")
    LOG_FILE = RESULTS_DIR / log_name

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=LOG_FILE,
        level=logging.INFO,
        format="%(asctime)s %(message)s"
    )

    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURE_DIR)
    sc.settings.verbosity = 2

    return RESULTS_DIR, FIGURE_DIR, LOG_FILE
