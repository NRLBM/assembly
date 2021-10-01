#!/bin/bash

# Set bash in strict mode etc.
set -euo pipefail

# Check if a single argument is provided
if [[ $# -ne 1 ]]
then
  echo "This script needs a single argument (R1 read file) and prints result to STDOUT. Exiting"
  exit 1
else
  # Sketch R1 read file and parse using awk and printf
  # Based on ~200 test isolates, this is on average an underestimate of 4% but never deviated more than -10% or +10% for bacterial WGS data
  mash sketch -o tmpfile -k 32 -r -m 3 "$1" 2>&1 1>/dev/null | awk 'NR == 1 {print $NF}' | xargs -i printf "%.0f" {}
  rm tmpfile.msh
fi
