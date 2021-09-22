#!/bin/bash

set -euo pipefail

if [[ $# -ne 1 ]]
then
  echo "This script needs a single argument (R1 read file) and prints result to STDOUT. Exiting"
  exit 1
else
  mash sketch -o tmpfile -k 32 -r -m 3 "$1" 2>&1 1>/dev/null | awk 'NR == 1 {print $NF}' | xargs -i printf "%.0f" {}
  rm tmpfile.msh
fi
