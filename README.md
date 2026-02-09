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

## Usage and [documentation of the method to detect interfaces](documentation.ipynb)
After loading the python environment, run `python plot-simba.py`. In case it has not been done already, this will download a satellite image from [Copernicus](https://browser.dataspace.copernicus.eu/?zoom=7&lat=47.40579&lng=-71.80664&themeId=DEFAULT-THEME&visualizationUrl=U2FsdGVkX1%2BhknApBktOR3T9tLuf8eHf2vvMNVOxfbcd8ETVg1LHjRpAq1sZzr4IylHK7fdylcYbFK4cXCrbFJ1AevZuBe3RI3beTjwVu0tUN30WWFj1jisRYR8f2xLc&datasetId=S2_L2A_CDAS&fromTime=2026-02-03T00%3A00%3A00.000Z&toTime=2026-02-03T23%3A59%3A59.999Z&layerId=1_TRUE_COLOR&demSource3D=%22MAPZEN%22&cloudCoverage=80&dateMode=SINGLE) with the location and date specified as in [`locdate.json`](locdate.json). The script will then verify that the required temperature date is present and perform the [detection of the interfaces](https://github.com/McGill-sea-ice/plot-simba/blob/08f210c09a43b525764115642018a4e448afdee5/plot_simba.py#L20) if ice is present. Once this is done data will be exported in `json` format for use in the interactive [Apache echarts](https://echarts.apache.org/examples/en/index.html) graphs on a [website like this one.](https://web.meteo.mcgill.ca/tremblay/mistissini). Then the non-interactive graphs will be plotted and saved.
For details on the methods, see the [documentation notebook](documentation.ipynb).

### Automation
The file `plot_simba.sh` contains bash code that handles loading the correct conda environment and executing `plot_simba.py`. Note that `source /opt/anaconda3/etc/profile.d/conda.sh` is necessary due to the way the conda environments are set up on the machine that this code was developped on but almost certainly needs to be adjusted or removed depending on your local machine.  
An example of how use [cron](https://en.wikipedia.org/wiki/Cron) to automatically run `plot_simba.sh` every day to plot the data created by [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) is shown in `to_crontab`. Including this in your crontab will create a log file `plot_simba.log`. Don't forget to adjust paths and directories.
