# plot-simba
Code to detect ice interfaces and create plots as shown on https://web.meteo.mcgill.ca/tremblay/mistissini

# Documentation and usage

## Preparations

### conda environment
The python code in `plot_simba.py` requires certain modules included in `environment.yml`. We thus create a [conda](https://anaconda.org/channels/anaconda/packages/conda/overview) environment with those modules that we can later load to run the code. `conda env create -f environment.yml` will create the required environment.

### input data requirements
`plot-simba` requires input that has been created with [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) which in turn has been downloaded using `get-buoy-data](https://github.com/McGill-sea-ice/get-buoy-data).  If you create your data in another way, make sure to set up the directory structure and naming conventions as in [`convert-buoy-data`](https://github.com/McGill-sea-ice/convert-buoy-data) to make sure that `plot-simba` functions properly.

### Copernicus SentinelHub
To include satellite images in the produced plot, the user needs to have an account on the [Copernicus Data Space Ecosystem (CDSE)](https://documentation.dataspace.copernicus.eu/Registration.html) in order to use the [SentinelHub Process API](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Process.html) for the download of those images.

## Usage
