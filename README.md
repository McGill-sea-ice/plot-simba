# plot-simba
Code to detect ice interfaces and create plots as shown on https://web.meteo.mcgill.ca/tremblay/mistissini

# Documentation and usage

## Preparations

### conda environment
The python code in `plot_simba.py` requires certain modules included in `environment.yml`. We thus create a [conda](https://anaconda.org/channels/anaconda/packages/conda/overview) environment with those modules that we can later load to run the code. `conda env create -f environment.yml` will create the required environment.

### input data requirements
`plot-simba` requires input that has been created with [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) which in turn has been downloaded using [`get-buoy-data`](https://github.com/McGill-sea-ice/get-buoy-data).  If you create your data in another way, make sure to set up the directory structure and naming conventions as in [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) to make sure that `plot-simba` functions properly.

### Copernicus SentinelHub
To include satellite images in the produced plot, the user needs to have an account on the [Copernicus Data Space Ecosystem (CDSE)](https://documentation.dataspace.copernicus.eu/Registration.html) in order to use the [SentinelHub Process API](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Process.html) for the download of those images.

### `locdate.json`
For downloading of the satellite image and to plot the text and pointers on the downloaded image, we need to provide some coordinates and text. This is done in the file `IMEI_locdate.json`, where `IMEI` has to be replaced by your instrument's IMEI. An example is included in the repository as `locate.json`.

## Usage

### Automation
The file `plot_simba.sh` contains bash code that handles loading the correct conda environment and executing `plot_simba.py`. Note that `source /opt/anaconda3/etc/profile.d/conda.sh` is necessary due to the way the conda environments are set up on the machine that this code was developped on but almost certainly needs to be adjusted or removed depending on your local machine.  
An example of how use [cron](https://en.wikipedia.org/wiki/Cron) to automatically run `plot_simba.sh` every day to plot the data created by [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) is shown in `to_crontab`. Including this in your crontab will create a log file `plot_simba.log`. Don't forget to adjust paths and directories.
