from pathlib import Path
from datetime import datetime
import os

localrules: all, select_best, make_inits

INITS_NAMES = [f"{config['analysis_name']}_{i}" for i in range(1, config["num_inits"] + 1)]
if "dt" not in config.keys():
    DT_STAMP = datetime.now().strftime('%Y%m%d_%H%M%S')
else:
    DT_STAMP = config["dt"]
UNIQUE_PATH = os.path.join("results", DT_STAMP)
OUT_DIR = config["out_dir"]

rule all:
    input:
        expand(
            [
                OUT_DIR + "/results/{run_id}/best/trajectory.png",
                OUT_DIR + "/results/{run_id}/best/init_id.stamp",
                OUT_DIR + "/results/{run_id}/stats/merged_stats.tsv",
                OUT_DIR + "/results/{run_id}/plots/{sample}/trajectory.png",
                OUT_DIR + "/reports/{run_id}.html"
            ],
            run_id=DT_STAMP,
            sample=INITS_NAMES
        )


rule apply_filters:
    input:
        hk_cc_genes=OUT_DIR + "/resources/datasets/HK_CC_genes_list.rds"
    output:
        dataset=OUT_DIR + "/resources/preprocess/dataset.rds",
        top_mad=OUT_DIR + "/resources/preprocess/topMAD.png",
        top_mad_joined=OUT_DIR + "/resources/preprocess/topMAD_joined.png",
        mad_med_scatter=OUT_DIR + "/resources/preprocess/mad-med.png",
        distance_before=OUT_DIR + "/resources/preprocess/distancesBefore.png",
        distance_after=OUT_DIR + "/resources/preprocess/distancesAfter.png",
        svd_before=OUT_DIR + "/resources/preprocess/svdBefore.png",
        svd_before_plot=OUT_DIR + "/resources/preprocess/svdBeforeToPlot.rds",
        svd_after=OUT_DIR + "/resources/preprocess/svdAfter.png",
        distances_to_zero_before=OUT_DIR + "/resources/preprocess/distances_to_zero_before.png",
        distances_to_zero_after=OUT_DIR + "/resources/preprocess/distances_to_zero_after.png",
        metadata=OUT_DIR + "/resources/preprocess/dataset.txt"
    params:
        dataset=f"{OUT_DIR}/resources/datasets/{config['dataset']}.rds"
    threads: config['count']['threads']
    resources:
        mem_ram=config['count']['mem_ram'],
        time=config['count']['time'],
        email=config['count']['email'],
        nodes=config['count']['nodes'],
        docker=config['count']['docker']
    script:
        "scripts/PrepareDataset.R"

rule make_inits:
    input:
        OUT_DIR + "/resources/preprocess/dataset.rds"
    output:
        expand(
            OUT_DIR + "/resources/inits/{sample}.rds",
            sample=INITS_NAMES
        )
    threads: config['count']['threads']
    resources:
        mem_ram=config['count']['mem_ram'],
        time=config['count']['time'],
        email=config['count']['email'],
        nodes=config['count']['nodes'],
        docker=config['count']['docker']
    script:
        "scripts/MakeInits.R"

Path(UNIQUE_PATH).mkdir(parents=True,exist_ok=True)

rule run_optimization:
    input:
        init_file=OUT_DIR + "/resources/inits/{sample}.rds",
        dataset=OUT_DIR + "/resources/preprocess/dataset.rds"
    output:
        meta=OUT_DIR + "/results/{run_id}/meta/{sample}.meta",
        proportions=OUT_DIR + "/results/{run_id}/props/{sample}_proportions.tsv",
        basis_row=OUT_DIR + "/results/{run_id}/basis_row/{sample}_basis_fc.tsv",
        basis_column=OUT_DIR + "/results/{run_id}/basis_col/{sample}_basis_fc_clmn.tsv",
        stats=OUT_DIR + "/results/{run_id}/stats/{sample}_stats.tsv"
    threads: config['count']['threads']
    resources:
        mem_ram=config['count']['mem_ram'],
        time=config['count']['time'],
        email=config['count']['email'],
        nodes=config['count']['nodes'],
        docker=config['count']['docker']
    params:
        blocks_pipeline=config["blocks_pipeline"]
    log:
        OUT_DIR + "/logs/{run_id}/optimization/{sample}.log"
    script:
        "scripts/RunDeconvolution.R"

rule make_plots:
    input:
        OUT_DIR + "/results/{run_id}/meta/{sample}.meta"
    output:
        OUT_DIR + "/results/{run_id}/plots/{sample}/trajectory.png",
        OUT_DIR + "/results/{run_id}/plots/{sample}/negative_proportions.png",
        OUT_DIR + "/results/{run_id}/plots/{sample}/negative_basis.png",
        OUT_DIR + "/results/{run_id}/plots/{sample}/sum_to_one_const.png"
    threads: config['count']['threads']
    resources:
        mem_ram=config['count']['mem_ram'],
        time=config['count']['time'],
        email=config['count']['email'],
        nodes=config['count']['nodes'],
        docker=config['count']['docker']
    log:
        OUT_DIR + "/logs/{run_id}/plots/{sample}.log"
    script:
        "scripts/MakePlots.R"

rule select_best:
    input:
        plots=[
            OUT_DIR + "/results/{run_id}/plots/{sample}/trajectory.png".format(run_id=DT_STAMP, sample=sample)
            for sample in INITS_NAMES
        ],
        stats=[
            OUT_DIR + "/results/{run_id}/stats/{sample}_stats.tsv".format(run_id=DT_STAMP, sample=sample)
            for sample in INITS_NAMES
        ],
        basis_col=[
            OUT_DIR + "/results/{run_id}/basis_col/{sample}_basis_fc_clmn.tsv".format(run_id=DT_STAMP, sample=sample)
            for sample in INITS_NAMES
        ],
    params:
        run_id=DT_STAMP
    threads: config['count']['threads']
    resources:
        mem_ram=config['count']['mem_ram'],
        time=config['count']['time'],
        email=config['count']['email'],
        nodes=config['count']['nodes'],
        docker=config['count']['docker']
    output:
        OUT_DIR + "/results/{run_id}/stats/merged_stats.tsv",
        OUT_DIR + "/results/{run_id}/best/trajectory.png",
        OUT_DIR + "/results/{run_id}/best/init_id.stamp",
        OUT_DIR + "/results/{run_id}/best/negative_basis.png",
        OUT_DIR + "/results/{run_id}/best/negative_proportions.png",
        OUT_DIR + "/results/{run_id}/best/sum_to_one_const.png",
        OUT_DIR + "/results/{run_id}/best/init_points.png",
        OUT_DIR + "/results/{run_id}/best/final_points.png",
        OUT_DIR + "/results/{run_id}/best/metafile.meta",
        OUT_DIR + "/results/{run_id}/best/proportions.tsv",
        OUT_DIR + "/results/{run_id}/best/basis_row.tsv",
        OUT_DIR + "/results/{run_id}/best/basis_column.tsv",
        OUT_DIR + "/results/{run_id}/best/abundance.png",
        UMAP=OUT_DIR + "/results/{run_id}/plots/UMAP.png"

    script:
        "scripts/ProcessBestRun.R"


rule prepare_report:
    input:
        expand([OUT_DIR + "/results/{run_id}/stats/merged_stats.tsv",
                OUT_DIR + "/results/{run_id}/best/trajectory.png",
                OUT_DIR + "/results/{run_id}/best/negative_basis.png",
                OUT_DIR + "/results/{run_id}/best/negative_proportions.png",
                OUT_DIR + "/results/{run_id}/best/sum_to_one_const.png",
                OUT_DIR + "/results/{run_id}/best/metafile.meta",
                OUT_DIR + "/results/{run_id}/best/proportions.tsv",
                OUT_DIR + "/results/{run_id}/best/basis_column.tsv",
                OUT_DIR + "/results/{run_id}/plots/UMAP.png",
                OUT_DIR + "/results/{run_id}/best/abundance.png"],
            run_id=DT_STAMP)
    params:
        run_id=DT_STAMP
    threads: config['count']['threads']
    resources:
        mem_ram=config['count']['mem_ram'],
        time=config['count']['time'],
        email=config['count']['email'],
        nodes=config['count']['nodes'],
        docker=config['count']['docker']
    output:
        OUT_DIR + "/reports/{run_id}.html"

    script:
        "scripts/prepareReport.Rmd"
