#!/bin/bash

set -euo pipefail

if [ $# -ne 2 ]
then
  echo "Usage: bash $0 [CLUSTER LIST] [SKA FASTA OUTPUT DIR]"
  exit 1
fi

CLUSTER_LIST="$1"
SKA_OUTPUT_DIR="$2"

mkdir -p "$SKA_OUTPUT_DIR"

while read NAME FILE
do
  ska fasta -o "$SKA_OUTPUT_DIR"/"$NAME" "$FILE"
done < "$CLUSTER_LIST"
