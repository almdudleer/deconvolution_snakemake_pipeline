"""
python3 run_linseedv2.py -l \
  --snakemake_path "out/21_08_2022/snakemake" \
  --data_path "out/datasets" \
  --inits_path "out/21_08_2022/inits" \
  --results_path "out/21_08_2022/results" \
  --reports_path "out/21_08_2022/reports" \
  --blocks_path "out/21_08_2022/blocks.csv" \
  --num_inits 10 --min_ct 6 --max_ct 10\
  --min_mad 10 --filter_genes 2000 \
  --dataset "HNSC_HK_CC_MALE" \
  --analysis_name "HNSC_HK_CC_MALE_MAD10_GENES_2000" \
  --knn_filter
"""

from optparse import OptionParser
from subprocess import Popen
import datetime as dt
from pathlib import Path
import yaml
import shutil

parser = OptionParser()
parser.add_option("-d", "--dataset")
parser.add_option("-a", "--analysis_name")
parser.add_option("--data_path")
parser.add_option("--snakemake_path")
parser.add_option("--inits_path")
parser.add_option("--results_path")
parser.add_option("--reports_path")
parser.add_option("--blocks_path")
parser.add_option("--snakefile_path", default="Snakefile")
parser.add_option("--dt")
parser.add_option("--init_strategy")
parser.add_option("--num_inits", type=int, default=5)
parser.add_option("--top_mad", type=int)
parser.add_option("--min_mad", type=int)
parser.add_option("--min_median", type=int)
parser.add_option("--filter_genes", type=int, default=0)
parser.add_option("--knn_filter", action="store_true", dest="knn_filter", default=False)
parser.add_option("--filter_samples", type=int, default=0)
parser.add_option("--scale_iterations", type=int, default=20)
parser.add_option("--min_ct", type=int)
parser.add_option("--max_ct", type=int)
parser.add_option("-l", action="store_true", dest="local", default=True)
parser.add_option("-b", action="store_false", dest="local")
parser.add_option("--apply_filters", action="store_true", dest="apply_filters", default=False)

(options, args) = parser.parse_args()

DT_STAMP = options.dt
if DT_STAMP is None:
    DT_STAMP = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

snakemake_path = Path(options.snakemake_path)
reports_path = Path(options.reports_path)
reports_path.mkdir(parents=True, exist_ok=True)

processes = []
with open(reports_path / f"{DT_STAMP}.html", "w+") as report:
    report.writelines("<html>")
    report.writelines(f"<head><title>{options.analysis_name}</title></head>")
    report.writelines("<body>")
    report.writelines(f"Dataset: {options.dataset}<br>")
    report.writelines(f"Analysis: {options.analysis_name}<br>")
    report.writelines(f"Report ID: {DT_STAMP}<br><br>")

    for ct in range(options.min_ct, options.max_ct + 1):
        WORK_DIR = snakemake_path / f"ct{ct}"

        dirs_symlinks = {
            WORK_DIR / "resources": None,
            WORK_DIR / "config": None,
            WORK_DIR / "rcpp_cache": None,
            Path(options.data_path): (WORK_DIR / "resources" / "datasets"),
            Path(options.inits_path): (WORK_DIR / "resources" / "inits"),
            Path(options.results_path, f"ct{ct}"): (WORK_DIR / "results"),
            reports_path / f"ct{ct}": (WORK_DIR / "reports"),
        }

        for directory, link in dirs_symlinks.items():
            directory.mkdir(parents=True, exist_ok=True)
            if link is not None and not link.is_symlink():
                link.symlink_to(directory.absolute())

        config_dict = {
            "num_inits": options.num_inits,
            "knn_filter": options.knn_filter,
            "filter_genes": options.filter_genes,
            "filter_samples": options.filter_samples,
            "scale_iterations": options.scale_iterations,
            "init_strategy": options.init_strategy,
            "cell_types": ct,
            "dt": DT_STAMP,
            "dataset": options.dataset,
            "analysis_name": f"{options.analysis_name}_ct{ct}",
            "blocks_pipeline": str((WORK_DIR / "config" / "blocks.csv").absolute()),
            "count": {
                "time": 6000,
                "mem_ram": 32,
                "threads": 8,
                "email": "aladyeva.e@wustl.edu",
                "nodes": -1,
                "docker": "dockerreg01.accounts.ad.wustl.edu/artyomov_lab/docker_linseed_snakemake:cpp",
            },
            "out_dir": str(WORK_DIR.absolute()),
            "rcpp_cache_dir": str((WORK_DIR / "rcpp_cache").absolute())
        }
        if options.top_mad is not None:
            config_dict["top_mad"] = options.top_mad
        if options.min_mad is not None:
            config_dict["min_mad"] = options.min_mad
        if options.min_median is not None:
            config_dict["min_median"] = options.min_median

        config_path = WORK_DIR / "config" / "config.yaml"
        with open(config_path, "w+") as ff:
            yaml.dump(config_dict, ff, allow_unicode=True, default_flow_style=False)

        shutil.copyfile(options.blocks_path, WORK_DIR / "config" / "blocks.csv")
        print(WORK_DIR)
        if options.local:
            cmd = [
                "snakemake",
                "-c4",
                "--rerun-incomplete",
                "--configfile",
                str(config_path.absolute()),
            ]
            print(' '.join(cmd))
            p = Popen(cmd, shell=False, stdin=None, stdout=None, stderr=None, close_fds=True)
        else:
            snk_cmd = """snakemake --profile lsf_demo --local-cores $L_CORES --jobs 50 \\
              -pr --conda-frontend conda --restart-times 3 --rerun-incomplete"""
            if options.apply_filters:
                snk_cmd += " -f apply_filters"
            cmd = f"""
            P_WD=`pwd`
            mkdir -p "$P_WD/tmp"
            echo "__LSF_JOB_CUSTOM_TMPDIR__=$P_WD/tmp" > lsf_docker_env_file.env
            chmod a+r lsf_docker_env_file.env
            LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env
            export SMK_DOCKER_IMG="dockerreg01.accounts.ad.wustl.edu/artyomov_lab/docker_linseed_snakemake:cpp"
            export P_LOG=$P_WD/logs/pipeline.log L_CORES=4 LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env 
            mkdir -p logs
            bsub -cwd $HOME -n $L_CORES -G compute-martyomov -q general -oo $P_LOG -R 'span[hosts=1]' \\
              -a "docker($SMK_DOCKER_IMG)" /usr/bin/script -fqe /dev/null \\
              -c "source /etc/bash.bashrc; cd $P_WD; export TMPDIR=$P_WD/tmp; {snk_cmd}"
            """
            print(cmd)
            p = Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True, cwd=WORK_DIR)
        processes.append(p)
        report.writelines(f"<a href='./ct{ct}/{DT_STAMP}.html'>{ct} cell types</a><br>")
    report.writelines("</body>\n</html>")
for p in processes:
    p.wait()
