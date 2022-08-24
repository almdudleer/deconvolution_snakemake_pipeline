export ANALYSIS_NAME="HNSC_HK_CC_MAD0_GENES_2000_RANDOM_INIT"
export OUT_DIR="out/$ANALYSIS_NAME"
python scripts/run_linseedv2.py -l --min_ct 3 --max_ct 9 --num_inits 10 \
  --snakemake_path $OUT_DIR/snakemake \
  --data_path datasets --dataset HNSC_HK_CC \
  --inits_path $OUT_DIR/inits \
  --results_path $OUT_DIR/results \
  --reports_path $OUT_DIR/reports \
  --blocks_path templates/blocks.csv \
  --min_mad 0 --filter_genes 2000 --knn_filter \
  --init_strategy SelectRandom \
  --analysis_name $ANALYSIS_NAME --dt 20220825_124300