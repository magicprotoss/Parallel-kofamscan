"""Run the normal workflow on the test VM while retaining intermediates for audit.

Accepts the same arguments as pkofamscan. Use a fresh output directory and run
inside tmux. This harness deliberately disables the workflow's cleanup hook.
"""

import importlib.machinery
import importlib.util
from pathlib import Path
import subprocess

import yaml


def main():
    repo = Path(__file__).resolve().parents[1]
    cli_path = repo / "parallel_kofam_scan" / "pkofamscan"
    loader = importlib.machinery.SourceFileLoader("pkofamscan_cli", str(cli_path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    cli = importlib.util.module_from_spec(spec)
    loader.exec_module(cli)
    opts = cli.args_to_opts(cli.parse_args())
    output = Path(opts["output_path"])
    config = output / "snakemake_config.yaml"
    config.write_text(yaml.safe_dump(opts), encoding="utf-8")
    print(f"Retaining intermediates at {opts['working_dir']}", flush=True)
    subprocess.run([
        "snakemake", "--snakefile", str(repo / "parallel_kofam_scan/workflow/exec_annotation.smk"),
        "--configfile", str(config), "--cores", str(opts["workers"]),
        "--rerun-incomplete", "--printshellcmds", "--directory", str(output),
        "--no-hooks",
    ], check=True)


if __name__ == "__main__":
    main()
