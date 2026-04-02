#!/bin/bash
# to be run each day, plotting sbd data
echo "---------- plot_simba.sh -----------"
echo " "
date

# load conda environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate /opt/anaconda3/envs/buoy-data

python3 /path/to/plot-simba/plot_simba_403539995324566/plot_simba.py
