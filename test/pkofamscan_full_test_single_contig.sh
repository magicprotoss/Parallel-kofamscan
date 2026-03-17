#!/bin/bash

# not working, manually run the activate conda env first
# source activate parallel_kofam_scan

# during dev
# mkdir -p ~/sdb/parallel_kofam_scan_test_outs/main
# mkdir -p ~/sdb/parallel_kofam_scan_test_outs/wd

# rm -rf ~/sdb/parallel_kofam_scan_test_outs/* && \
#     python /home/navi/My-Github-Projects/Parallel-kofamscan/parallel_kofam_scan/pkofamscan \
#         -i /home/navi/My-Github-Projects/Parallel-kofamscan/test/data/SM1510S13_test.faa \
#         -o ~/sdb/parallel_kofam_scan_test_outs/main \
#         --working_dir ~/sdb/parallel_kofam_scan_test_outs/wd \
#         -db /mnt/kofamscan-db/ \
#         -p 56

# final test before release

mkdir -p ~/sdb/parallel_kofam_scan_test_outs
rm -rf ~/sdb/parallel_kofam_scan_test_outs/* && \
    python /home/navi/My-Github-Projects/Parallel-kofamscan/parallel_kofam_scan/pkofamscan \
        -i /home/navi/My-Github-Projects/Parallel-kofamscan/test/data/SM1510S13_test.faa \
        -o ~/sdb/parallel_kofam_scan_test_outs \
        -db /mnt/kofamscan-db/ \
        -p 56

# clean up
# rm -rf ~/sdb/parallel_kofam_scan_test_outs