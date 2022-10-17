#!/bin/bash

set -euo pipefail

DB_PATH=${1%/*}
DB_FILE=${1##*/}

mkdir -p "$DB_PATH"

wget --spider "$2"

wget -O "$1" "$2"

cd "$DB_PATH"

makeblastdb -in "$DB_FILE" -title emmtyper_db -dbtype nucl -hash_index
