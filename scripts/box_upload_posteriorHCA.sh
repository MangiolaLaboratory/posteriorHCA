#!/usr/bin/env bash
# Upload posteriorHCA to UofA Box (posteriorHCA folder).
# Run inside screen: screen -r posteriorHCA_box
set -euo pipefail

LOG_DIR="/home/a1237163/lab/chen/posteriorHCA/logs"
mkdir -p "$LOG_DIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$LOG_DIR/box_upload_${STAMP}.log"
SRC="/home/a1237163/lab/chen/posteriorHCA"
REMOTE="UofA_Box:posteriorHCA"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rclone_env

EX=(
  --exclude "**/_cache/**"
  --exclude "**/*_files/**"
  --exclude "**/intermediate/GSE293189_pre_doublet_checkpoint.rds"
  --exclude "**/count_matrix_sparse.mtx"
  --exclude "**/formula_bench/_targets/**"
  --exclude "wget-log*"
  --exclude ".git/**"
)

RCLONE_FLAGS=(
  --progress
  --stats 1m
  --stats-one-line
  --transfers 4
  --checkers 8
  --retries 5
  --low-level-retries 10
  --log-file "$LOG"
  --log-level INFO
)

exec > >(tee -a "$LOG") 2>&1

echo "=== posteriorHCA Box upload started $(date -Is) ==="
echo "Log: $LOG"
echo "Remote: $REMOTE"
rclone version | head -1
echo

echo "--- [1/2] dev/ case studies ---"
rclone copy "$SRC/dev" "$REMOTE/dev" "${EX[@]}" "${RCLONE_FLAGS[@]}"
echo "dev/ done $(date -Is)"
echo

echo "--- [2/2] package source (no dev, no .git) ---"
rclone copy "$SRC" "$REMOTE/package" \
  --exclude "dev/**" \
  --exclude ".git/**" \
  --exclude "wget-log*" \
  --exclude "**/_cache/**" \
  --exclude "**/*_files/**" \
  "${RCLONE_FLAGS[@]}"
echo "package/ done $(date -Is)"
echo

echo "--- remote size ---"
rclone size "$REMOTE"
echo
echo "=== upload finished $(date -Is) ==="
