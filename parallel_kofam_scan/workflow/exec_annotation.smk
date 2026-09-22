#############
### Setup ###
#############
ruleorder: validate_input > concat_shuffle_split > exec_annotation > collect_results_run_top_hit
from snakemake.utils import min_version
min_version("9.16.0")

import glob
import os
os.umask(0o002)

# how do I pass the parameters from the helper script pkofamscan and invoke the validate_input rule in parallel?
# here's an example opts object created by the helper script pkofamscan
# {'input_path': ['/home/navi/My-Github-Projects/Parallel-kofamscan/test/data/Sample-1.fna'], 'output_path': '/home/navi/sdb/parallel_kofam_scan_test_outs/main/KEGG_annotations', 'ko_list': '/mnt/kofamscan-db/ko_list', 'profiles': '/mnt/kofamscan-db/profiles', 'working_dir': '/tmp/parallel_kofamscan/1000/30350551ec0c37f0544002c7053032c3_0f4d922f384f104498c71327af6ae664', 'workers': 56, 'threads': 8, 'chunck_size': 100000}
# basically, given a list of input_path, the working_dir, and the num workers:
# snakemake will execute the validate_input rule for each input_path in parallel, using a maximum of num workers, 
# and the output of each validate_input rule will be stored in the working_dir with a unique sample_id for each input_path \
# --this part is handeled by the script validate_input.py itself, it will compute fps using the input fp

working_dir = config['working_dir']
output_dir = config['output_path']
log_dir = os.path.join(output_dir, "logs")
os.makedirs(log_dir, exist_ok=True)
scripts_dir = os.path.join(workflow.basedir, 'scripts')

# Helper function to get sample ID from filename
# need to consider .gz as well
# it's dumb bt it works TT
def get_sample_id(filepath):
    sample_id = os.path.basename(filepath)
    if sample_id.endswith('.fna'):
        sample_id = sample_id[:-4]  # Remove .fna extension
    elif sample_id.endswith('.fa'):
        sample_id = sample_id[:-3]  # Remove .fa extension
    elif sample_id.endswith('.fasta'):
        sample_id = sample_id[:-6]  # Remove .fasta extension
    elif sample_id.endswith('.faa'):
        sample_id = sample_id[:-4]  # Remove .faa extension
    elif sample_id.endswith('.fna.gz'):
        sample_id = sample_id[:-7]  # Remove .fna.gz extension
    elif sample_id.endswith('.fa.gz'):
        sample_id = sample_id[:-6]  # Remove .fa.gz extension
    elif sample_id.endswith('.fasta.gz'):
        sample_id = sample_id[:-9]  # Remove .fasta.gz extension
    elif sample_id.endswith('.faa.gz'):
        sample_id = sample_id[:-7]  # Remove .faa.gz extension
    

    return sample_id


# Map sample IDs to input paths
# Check for duplicate sample IDs
sample_ids = [get_sample_id(fp) for fp in config['input_path']]
if len(sample_ids) != len(set(sample_ids)):
    from collections import Counter
    duplicates = [item for item, count in Counter(sample_ids).items() if count > 1]
    raise ValueError(f"Duplicate sample IDs found: {', '.join(duplicates)}. Please check your input files.")

SAMPLES = {sid: fp for sid, fp in zip(sample_ids, config['input_path'])}

def get_annotation_targets(wildcards):
    checkpoint_output = checkpoints.concat_shuffle_split.get(**wildcards).output[0]
    split_dir = os.path.join(os.path.dirname(checkpoint_output), "all_hashed_unigenes_no_dup_split")
    files = glob.glob(os.path.join(split_dir, "stdin.part_*.fasta"))
    if not files:
        return []
    part_ids = [os.path.basename(f).split('.')[1] for f in files]
    return [os.path.join(working_dir, "exec_annotation", pid, "annotations.tsv") for pid in part_ids]


#####################################
### Rule all checks final outputs ###
#####################################


rule all:
    input:
        os.path.join(output_dir, "final_results", "top-hits.xlsx"),
        os.path.join(output_dir, "final_results", "hmm-hits.parquet"),
        os.path.join(output_dir, "final_results", "sequence-metadata-map.parquet")

#######################
### Validate inputs ###
#######################

# You need to specify the output files in the Snakemake rule for two main reasons, even if the python script itself doesn't explicitly take them as arguments (but rather infers them from the output directory):

# Dependency Graph: Snakemake builds a DAG (Directed Acyclic Graph) of jobs. It needs to know what files a rule produces so it can determine if downstream rules (which will use these files as input) can be run. If you don't list the outputs, Snakemake won't know that validate_input produces header_map.parquet, unigene.faa, etc., and thus won't know how to create them if another rule asks for them.
# Tracking Completion: Snakemake uses the existence and modification timestamps of output files to determine if a rule has completed successfully and if it needs to be re-run. If you don't list the outputs, Snakemake might run the rule, but it won't "know" the job achieved anything, and might re-run it unnecessarily or fail to trigger subsequent steps.
# In your case, validate_input.py takes an output directory (-o), but Snakemake tracks files. Listing the specific files inside that directory tells Snakemake "when this rule is done, expect these specific files to exist in that directory".

rule validate_input:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        os.path.join(working_dir, "validate_input", "{sample}", "header_map.parquet"),
        os.path.join(working_dir, "validate_input", "{sample}", "unigene_hash.faa"),
        os.path.join(working_dir, "validate_input", "{sample}", "done")
    log:
        os.path.join(log_dir, "validate_input", "{sample}.log")
    params:
        script = scripts_dir + "/validate_input.py",
        out_dir = lambda wildcards: os.path.join(working_dir, "validate_input", wildcards.sample)
    threads: 1
    shell:
        "python3 {params.script:q} -i {input:q} -o {params.out_dir:q} > {log:q} 2>&1"

##########################################################
### Concat unigenes, shuffle them and split to chuncks ###
##########################################################

checkpoint concat_shuffle_split:
    input:
        expand(os.path.join(working_dir, "validate_input", "{sample}", "done"), sample=SAMPLES.keys())
    output:
        os.path.join(working_dir, "concat_shuffle_split", "done")
    log:
        os.path.join(log_dir, "concat_shuffle_split.log")
    params:
        script = scripts_dir + "/concat_shuffle_split.py",
        validate_input_dir = os.path.join(working_dir, "validate_input"),
        out_dir = os.path.join(working_dir, "concat_shuffle_split"),
        chunk_size = config['chunk_size']
    threads: config['workers']
    shell:
        "python3 {params.script:q} -i {params.validate_input_dir:q} -o {params.out_dir:q} -c {params.chunk_size} -t {threads} > {log:q} 2>&1"

###########################
### Run exec_annotation ###
###########################

rule exec_annotation:
    input:
        os.path.join(working_dir, "concat_shuffle_split", "all_hashed_unigenes_no_dup_split", "stdin.{part_id}.fasta")
    output:
        annotations = os.path.join(working_dir, "exec_annotation", "{part_id}", "annotations.tsv")
    log:
        os.path.join(log_dir, "exec_annotation", "kofam_scan_{part_id}.log")
    params:
        ko_list = config['ko_list'],
        profiles = config['profiles'],
        tmp_dir = lambda wildcards: os.path.join(working_dir, "exec_annotation", wildcards.part_id, "tmp"),
    threads: config['threads']
    shell:
        "conda run -n parallel-kofamscan.dependency.kofamscan exec_annotation -o {output.annotations:q} {input:q} -k {params.ko_list:q} -p {params.profiles:q} --cpu {threads} -f detail-tsv --tmp-dir {params.tmp_dir:q} > {log:q} 2>&1"

###########################################
### Collect results and perform top-hit ###
###########################################

rule collect_results_run_top_hit:
    input:
        annotations = get_annotation_targets,
        header_maps = expand(os.path.join(working_dir, "validate_input", "{sample}", "header_map.parquet"), sample=SAMPLES.keys())
    output:
        os.path.join(output_dir, "final_results", "top-hits.xlsx"),
        os.path.join(output_dir, "final_results", "hmm-hits.parquet"),
        os.path.join(output_dir, "final_results", "sequence-metadata-map.parquet")
    log:
        os.path.join(log_dir, "collect_results.log")
    benchmark:
        os.path.join(log_dir, "collect_results.benchmark.tsv")
    params:
        script = scripts_dir + "/collect_results.py",
        working_dir = working_dir,
        output_dir = os.path.join(output_dir, "final_results"),
        min_e_value = float(config['min_e_value']),
        min_score = int(config['min_score'])
    threads: 1
    shell:
        "python3 {params.script:q} -wd {params.working_dir:q} -o {params.output_dir:q} -t {threads} -e {params.min_e_value} -s {params.min_score} > {log:q} 2>&1"

onsuccess:
    import shutil
    print(f"Workflow finished successfully. Removing temporary working directory: {working_dir}")
    if os.path.exists(working_dir):
        shutil.rmtree(working_dir)
