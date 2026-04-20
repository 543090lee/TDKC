#!/bin/bash

set -u 
set -e 

NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
DATA_DIR="data"
RSYNC_SERVER="rsync://$NCBI_SERVER"
FTP_SERVER="ftp://$NCBI_SERVER"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

function download_file() {
  file="$1"
  if [ -n "${TDKC_USE_FTP:-}" ]; then
    wget -q "${FTP_SERVER}${file}"
  else
    rsync --no-motd "${RSYNC_SERVER}${file}" .
  fi
}

for subsection in gb wgs; do
  fname="nucl_${subsection}.accession2taxid.gz"
  if [ ! -e "${fname%.gz}" ] && [ ! -e "$fname" ]; then
    download_file "/pub/taxonomy/accession2taxid/nucl_${subsection}.accession2taxid.gz"
  else
    1>&2 echo "accession2taxid files are already present, gonna skip."
  fi
done

if ls ./*.accession2taxid.gz 2>/dev/null | grep -q '.'; then
  gunzip ./*.accession2taxid.gz
  1>&2 echo " done."
fi

if [ ! -e "nodes.dmp" ]; then
  download_file "/pub/taxonomy/taxdump.tar.gz"
  tar -zxf taxdump.tar.gz nodes.dmp
  rm taxdump.tar.gz
else
  1>&2 echo "nodes.dmp is already present, gonna skip."
fi

1>&2 echo "Taxonomy files are all downloaded."