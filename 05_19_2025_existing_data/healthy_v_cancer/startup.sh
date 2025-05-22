#!/bin/bash

# Startup script for a mock analysis

# Check if the script is being sourced
(return 0 2>/dev/null)
if [ $? -ne 0 ]; then
  echo "ERROR: This script must be sourced, not executed."
  echo "Use: source setup.sh"
  exit 1
fi

# Load required modules
module load R/4.4.2
module load sra-tools
module load gcc


module list
