MIN_MAD=2
INIT=SelectX
DATASET="HNSC"

ANALYSIS_NAME=$(echo "${DATASET}_MAD${MIN_MAD}_${INIT}_INIT" | tr '[:lower:]' '[:upper:]')
OUT_DIR="out/$ANALYSIS_NAME"
export ANALYSIS_NAME MIN_MAD INIT DATASET OUT_DIR

python scripts/run_linseedv2.py -l --min_ct 6 --max_ct 10 --num_inits 3 \
  --snakemake_path "$OUT_DIR"/snakemake \
  --data_path datasets --dataset $DATASET \
  --inits_path "$OUT_DIR"/inits \
  --results_path "$OUT_DIR"/results \
  --reports_path "$OUT_DIR"/reports \
  --blocks_path templates/blocks.csv \
  --min_mad $MIN_MAD --filter_genes 2000 --knn_filter \
  --init_strategy $INIT \
  --analysis_name "$ANALYSIS_NAME" --dt 20220825_124300