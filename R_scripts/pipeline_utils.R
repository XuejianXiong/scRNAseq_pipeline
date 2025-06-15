# -----------------------------
# Functions
# -----------------------------

# Function to set the log message
log_msg <- function(msg) {
  
  # Generate current timestamp in YYYY-MM-DD HH:MM:SS format
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Construct log line with timestamp prefix
  line <- sprintf("[%s] %s", timestamp, msg)
  
  # Print log message to console
  cat(line, "\n")
  
  # Append log message to the specified log file
  # Assumes 'log_file' is a global variable set elsewhere in the script
  write(line, file = log_file, append = TRUE)
}

