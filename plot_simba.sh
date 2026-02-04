#!/bin/bash
# to be run each day, plotting sbd data
echo "---------- plot_simba.sh -----------"
echo " "
date

# load conda environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate /opt/anaconda3/envs/buoy-data

python3 /storage/common/buoy-data/plot-simba/plot_simba.py

cp /storage/common/buoy-data/plot-simba/300340657110080_top.png /storage2/website/public_html/images/
cp /storage/common/buoy-data/plot-simba/300340657110080_top_fr.png /storage2/website/public_html/images/
cp /storage/common/buoy-data/plot-simba/300340657110080_bottom.png /storage2/website/public_html/images/
cp /storage/common/buoy-data/plot-simba/300340657110080_bottom_fr.png /storage2/website/public_html/images/

cp /storage/common/buoy-data/plot-simba/300340657110080_edata.json /storage2/website/public_html/data/
cp /storage/common/buoy-data/plot-simba/300340657110080_colormap.json /storage2/website/public_html/data/
