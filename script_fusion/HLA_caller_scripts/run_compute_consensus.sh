#!/bin/bash

module load python/3.6.1

python ${PIPELINE_DIR}/compute_consensus.py $RESULTS_DIR

sleep 1m