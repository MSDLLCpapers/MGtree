#!/usr/bin/env bash

# test MGtree scripts
# John M. Gaspar, Scott Norton, Samantha Sholes
# Version 1.0 2024-10-03

HOME_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# input test-data/ref files
input_nwk0="${HOME_DIR}/test-data/testNoV_CDC_reference.nwk"
input_nwk1="${HOME_DIR}/test-data/testNoV_CDC_reference_id.nwk"
input_nwk2="${HOME_DIR}/test-data/testNoV_CDC_reference_nodesnamed.nwk"
input_rename="${HOME_DIR}/test-data/testNoV_CDC_rename_nodes.csv"
input_sam="${HOME_DIR}/test-data/testSRR11028620_mapped.sam"
input_tsv="${HOME_DIR}/test-data/testSRR11028620_parsed.tsv"
ref_gen="${HOME_DIR}/test-data/testSRR11028620_typed.tsv"
ref_leaf="${HOME_DIR}/test-data/testSRR11028620_leafnodes.tsv"
ref_json="${HOME_DIR}/test-data/testSRR11028620_qname.json"

# output test files
output_nwk1="${HOME_DIR}/test-data/output_testNoV_CDC_reference_id.nwk"
output_nwk2="${HOME_DIR}/test-data/output_testNoV_CDC_reference_nodesnamed.nwk"
output_tsv="${HOME_DIR}/test-data/output_testSRR11028620_parsed.tsv"
output_gen="${HOME_DIR}/test-data/output_testSRR11028620_typed.tsv"
output_leaf="${HOME_DIR}/test-data/output_testSRR11028620_leafnodes.tsv"
output_json="${HOME_DIR}/test-data/output_testSRR11028620_qname.json"

# remove previous output files
rm -f "$output_nwk1" "$output_nwk2" "$output_tsv" \
  "$output_gen" "$output_leaf" "$output_json"
nsucceeded=0
nfailed=0

# test updateNewick.py naming nodes
python "${HOME_DIR}/bin/updateNewick.py" \
  -i "$input_nwk0" \
  -o "$output_nwk1"
diff "$input_nwk1" "$output_nwk1"
code=$?
if [ $code -eq 0 ]; then
  ((++nsucceeded))
else
  echo "updateNewick.py naming nodes FAILED" 1>&2
  ((++nfailed))
fi

# test updateNewick.py genotyping nodes
python "${HOME_DIR}/bin/updateNewick.py" \
  -i "$input_nwk1" \
  -n "$input_rename" \
  -o "$output_nwk2"
diff "$input_nwk2" "$output_nwk2"
code=$?
if [ $code -eq 0 ]; then
  ((++nsucceeded))
else
  echo "updateNewick.py genotyping nodes FAILED" 1>&2
  ((++nfailed))
fi

# test parseSAM.py
python "${HOME_DIR}/bin/parseSAM.py" \
  -i "$input_sam" \
  -o "$output_tsv" \
  -S
diff "$input_tsv" "$output_tsv"
code=$?
if [ $code -eq 0 ]; then
  ((++nsucceeded))
else
  echo "parseSAM.py FAILED" 1>&2
  ((++nfailed))
fi

# test MGtree.py
python "${HOME_DIR}/bin/MGtree.py" \
  -i "$input_tsv" \
  -n "$input_nwk2" \
  -o "$output_gen" \
  -t "$output_leaf" \
  -q "$output_json" \
  -g

diff "$ref_gen" "$output_gen"
code=$?
if [ $code -eq 0 ]; then
  ((++nsucceeded))
else
  echo "MGtree.py genotype output FAILED" 1>&2
  ((++nfailed))
fi

diff "$ref_leaf" "$output_leaf"
code=$?
if [ $code -eq 0 ]; then
  ((++nsucceeded))
else
  echo "MGtree.py leaf output FAILED" 1>&2
  ((++nfailed))
fi

diff "$ref_json" "$output_json"
code=$?
if [ $code -eq 0 ]; then
  ((++nsucceeded))
else
  echo "MGtree.py json output FAILED" 1>&2
  ((++nfailed))
fi

echo "Summary: $nsucceeded successes, $nfailed failures" 1>&2
exit $nfailed
