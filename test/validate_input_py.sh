#!/bin/bash

conda activate parallel_kofam_scan

mkdir -p ~/sdb/parallel_kofam_scan_test_outs/validate_input_test_outs
# this script creates the outputdir
rm -rf ~/sdb/parallel_kofam_scan_test_outs/validate_input_test_outs && \
    python /home/navi/My-Github-Projects/Parallel-kofamscan/parallel_kofam_scan/workflow/scripts/validate_input.py \
        -i /home/navi/My-Github-Projects/Parallel-kofamscan/test/data/Sample-1.fna \
        -o ~/sdb/parallel_kofam_scan_test_outs/validate_input_test_outs

rm -rf ~/sdb/parallel_kofam_scan_test_outs

