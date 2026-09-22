#!/usr/bin/env bash
# VM-only continuation of the 2026-09-21 validation run. Launch inside tmux.
set -Eeuo pipefail

repo=/home/navi/pfamscan_tests/Parallel-kofamscan
stock=/home/navi/pfamscan_tests/repo_simple_test/stock-20260921
real=/home/navi/pfamscan_tests/700genome_test/run-20260921
scratch=/mnt/sdb/pfamscan-tests
export PATH=/home/navi/miniconda3/envs/parallel_kofam_scan/bin:/home/navi/miniconda3/bin:$PATH
mkdir -p "$real"
trap 'result=$?; printf "%s\n" "$result" > "$real/controller.exit-code"' EXIT

wait_for_job() {
    local status_file=$1 session=$2
    while [[ ! -f "$status_file" ]]; do
        if ! tmux has-session -t "$session" 2>/dev/null; then
            # Recheck after observing session exit to avoid a completion race.
            [[ -f "$status_file" ]] || { echo "Job $session disappeared without a status file"; return 1; }
        fi
        sleep 10
    done
    [[ $(< "$status_file") == 0 ]]
}

printf 'waiting-for-stock\n' > "$real/stage"
wait_for_job "$stock/exit-code" kofam-stock-20260921
export TMPDIR="$stock/tmp"
printf 'auditing-stock\n' > "$real/stage"
python "$repo/test/verify_vm_results.py" "$stock/KEGG_annotations/snakemake_config.yaml" \
    "$stock/verification.json" > "$stock/verification.log" 2>&1

printf 'preparing-scratch\n' > "$real/stage"
wait_for_job /home/navi/pfamscan_tests/scratch-setup.exit-code kofam-scratch-setup
mountpoint -q /mnt/sdb
[[ $(findmnt -n -o SOURCE /mnt/sdb) == /dev/sdb ]]
mkdir -p "$scratch/database/profiles" "$scratch/work-20260921" "$scratch/tmp-20260921"
export TMPDIR="$scratch/tmp-20260921"

# Cache the exact supplied database to reduce dependence on the network share.
printf 'caching-database\n' > "$real/stage"
rsync -a --partial /mnt/kofamscan-db/ko_list "$scratch/database/"
rsync -a --partial /mnt/kofamscan-db/profiles/ "$scratch/database/profiles/"
sha256sum /mnt/kofamscan-db/ko_list "$scratch/database/ko_list" > "$real/ko-list-sha256.txt"
conda list -n parallel_kofam_scan --explicit > "$real/environment-main.txt"
conda list -n parallel-kofamscan.dependency.kofamscan --explicit > "$real/environment-kofam.txt"
find "$repo/parallel_kofam_scan" "$repo/test" -type f \
    \( -name '*.py' -o -name '*.smk' -o -name '*.yaml' -o -name pkofamscan \) \
    -print0 | sort -z | xargs -0 sha256sum > "$real/source-sha256.txt"
mapfile -d '' inputs < <(find /home/navi/pfamscan_tests/700genome_test/input_bakta_faas \
    -maxdepth 1 -type f -name '*.faa' -print0 | sort -z)
[[ ${#inputs[@]} -gt 0 ]]
printf '%s\n' "${inputs[@]}" > "$real/input-manifest.txt"
printf 'Starting real workflow with %s input FAA files\n' "${#inputs[@]}"

printf 'real-workflow\n' > "$real/stage"
/usr/bin/time -v -o "$real/workflow-time.txt" \
    python -u "$repo/test/vm_workflow.py" -i "${inputs[@]}" -o "$real" \
    -d "$scratch/work-20260921" -db "$scratch/database" -p 24 -t 4 -c 100000 \
    > "$real/workflow.log" 2>&1

printf 'auditing-real-results\n' > "$real/stage"
/usr/bin/time -v -o "$real/audit-time.txt" \
    python "$repo/test/verify_vm_results.py" "$real/KEGG_annotations/snakemake_config.yaml" \
    "$real/verification.json" > "$real/verification.log" 2>&1
printf 'complete\n' > "$real/stage"
echo 'Stock and real-data validation completed successfully.'
