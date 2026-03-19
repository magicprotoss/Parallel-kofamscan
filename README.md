# Parallel-kofamscan

A snakemake workflow to split large unigenes.faa into chucks in order to prevent [kofam_scan](https://github.com/takaram/kofam_scan) from stuck at formatting results using ruby. This workflow also filters the results down to top-hits and saves it to a excel sheet for the user.

## Installation

1.  clone the repo

``` bash
git clone https://github.com/magicprotoss/Parallel-kofamscan
```

2.  install create 2 conda envs for the main shell as well as dependencies

``` bash
# main env, named parallel_kofam_scan
conda env create -f Parallel-kofamscan/parallel_kofam_scan/envs/main.yaml
# dependency env, named parallel-kofamscan.dependency.kofamscan
conda env create -f Parallel-kofamscan/parallel_kofam_scan/envs/kofamscan.yaml
```

3.  move the repo folder under the main env, soft link the scripts to the bin sub-dir, and add exec permission

``` bash
mv Parallel-kofamscan <your-path-to-miniconda3>/envs/parallel-kofamscan/
cd <your-path-to-miniconda3>/envs/parallel-kofamscan/Parallel-kofamscan
ln -s parallel_kofam_scan/pkofamscan ../bin/ && chmod u+x ../bin/pkofamscan
ln -s parallel_kofam_scan/workflow ../bin && chmod u+x ../bin/workflow/scripts/*.py
```

To use it, simply activate the conda env and run pkofamscan

``` bash
conda activate parallel_kofam_scan

pkofamscan --help
usage: Parallel-kofamscan [-h] -i INPUT_PATH [INPUT_PATH ...] -o OUTPUT_PATH -db PATH_TO_KEGG_DATABASE
                          [-d WORKING_DIR] [-p WORKERS] [-t THREADS_PER_WORKER] [-c CHUNK_SIZE] [-f FORCE]
                          [-e MIN_E_VALUE] [-s MIN_SCORE]

A snakemake workflow to split unigenes into chucks in order to reduce kofamscan's runtime

options:
  -h, --help            show this help message and exit
  -i INPUT_PATH [INPUT_PATH ...], --input_path INPUT_PATH [INPUT_PATH ...]
                        Unigenes to annotate, filenames will be converted to sample-ids, gzipped files are
                        supported as well
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Directory to save the annotation results, the output would be stored in an sub-dir
                        called KEGG_annotations
  -db PATH_TO_KEGG_DATABASE, --path_to_KEGG_database PATH_TO_KEGG_DATABASE
                        Path to the kofamscan database, which contains a 'ko_list' file and a 'profiles'
                        directory
  -d WORKING_DIR, --working_dir WORKING_DIR
                        Directory to store intermediate results, default is
                        /tmp/parallel_kofamscan/<your_uid>/<hash_of_output_dir>_<hash_of_current_datetime>,
                        if set, the temp dir would be <working_dir>/parallel_kofamscan/<your_uid>/<hash_of_ou
                        tput_dir>_<hash_of_current_datetime>
  -p WORKERS, --workers WORKERS
                        Number of cores to use for when running the workflow, default is 0, which indicates
                        the workflow will use the num of cores on the machine - 1, or the num of chuncks *
                        16, whichever is smaller
  -t THREADS_PER_WORKER, --threads_per_worker THREADS_PER_WORKER
                        Number of theards to use when running the exec_annotation cmd, default is 8
  -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Number of unigenes per chunk, default is 100000
  -f FORCE, --force FORCE
                        whether to overwrite existing result, default is no
  -e MIN_E_VALUE, --min_e_value MIN_E_VALUE
                        Minimum E-value threshold to retain a KO hit during top-hit filtering, default is
                        0.001
  -s MIN_SCORE, --min_score MIN_SCORE
                        Minimum hmmsearch score threshold to retain a KO hit during top-hit filtering,
                        default is 100
```
