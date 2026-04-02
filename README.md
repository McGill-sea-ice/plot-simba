# plot-simba
Code to detect ice interfaces and create plots as shown on https://web.meteo.mcgill.ca/tremblay/mistissini

# Documentation and usage

## Preparations

### conda environment
The python code in `plot_simba.py` requires certain modules included in `environment.yml`. We thus create a [conda](https://anaconda.org/channels/anaconda/packages/conda/overview) environment with those modules that we can later load to run the code. `conda env create -f environment.yml` will create the required environment.

### input data requirements
`plot-simba` requires input that has been created with [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) which in turn has been downloaded using [`get-buoy-data`](https://github.com/McGill-sea-ice/get-buoy-data).  If you create your data in another way, make sure to set up the directory structure and naming conventions as in [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) to make sure that `plot-simba` functions properly.

## [Usage](documentation.ipynb) and [documentation of the method to detect interfaces](documentation_detect.ipynb)
First, make sure you set the variable `imei` in `plot_simba.py` to the IMEI of your instrument. Then adjust `datapath` and `plotpath` in `plot_simba.py` as well as the path in `plot_simba.sh`.  
After loading the python environment, run `python plot_simba.py`. The script will verify that the required temperature date is present and perform the [detection of the interfaces](https://github.com/McGill-sea-ice/plot-simba/blob/08f210c09a43b525764115642018a4e448afdee5/plot_simba.py#L20) if ice is present. Once this is done data will be exported in `json` format for use in the interactive [Apache echarts](https://echarts.apache.org/examples/en/index.html) graphs on a [website like this one.](https://web.meteo.mcgill.ca/tremblay/mistissini). Then the non-interactive graphs will be plotted and saved.
For details on the usage see the [documentation notebook](documentation.ipynb). The method to detect the interfaces is described in further detail in another [notebook](documentation_detect.ipynb).

### Automation
The file `plot_simba.sh` contains bash code that handles loading the correct conda environment and executing `plot_simba.py`. Note that `source /opt/anaconda3/etc/profile.d/conda.sh` is necessary due to the way the conda environments are set up on the machine that this code was developped on but almost certainly needs to be adjusted or removed depending on your local machine.  
An example of how use [cron](https://en.wikipedia.org/wiki/Cron) to automatically run `plot_simba.sh` every day to plot the data created by [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) is shown in `to_crontab`. Including this in your crontab will create a log file `plot_simba.log`. Don't forget to adjust paths and directories.
