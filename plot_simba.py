import xarray as xr
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.patheffects as path_effects
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.patches import Rectangle
import cmocean.cm as cmo
from datetime import datetime
import numpy as np
import os
import sentinelhub as sh
import json

imei = "300340657110080"
datapath = "/storage/common/buoy-data/convert-buoy-data/decode_convert_" + imei
plotpath = "/storage/common/buoy-data/plot-simba/plot_simba_" + imei
config = sh.SHConfig("cdse")

def detect_interfaces(da, t_air, isurf, frozendate):
    '''Detect air-snow, snow-ice, and ice-water interfaces in temperature profile. Developed 
    for use with thermistor string data from a SAMS Enterprise Snow and Ice Mass Balance Apparatus (SIMBA)
    and NOT tested with any other data.

    Parameters
    ----------
    da : xarray Dataset
        2D xarray Dataset created by `convert-buoy-data` (https://github.com/McGill-sea-ice/convert-buoy-data)
        Dimensions are `time` and `pos` (position along thermistor string) and the temperature variable
        must be named `temp`. The last index along `pos` (air temperature) has to be removed.

    t_air : xarray DataArray
        1D xarray DataArray containing the air temperature associated with `da` (the last index along `pos`
        in the original dataset)

    isurf : int
        Integer describing at which index along `pos` the initial (non-frozen) water-air interface is located.

    forzendate : str
        Date on which the water first froze. Format is "YYYY-MM-DD"

    Returns
    -------
    idx : dict
        Dictionary containing the keys ("snowTop", "snowMid", "snowBot", "iceTop", "iceBot"), each with a
        time series of the indeces of the air-snow interface, the mid-point of the snow cover, the bottom
        of the snow cover, the top of the ice, and the ice-water interface. 

    dep : dict
        Same as `idx` but containing the "depths", i.e. the distance in cm along the thermistor string for
        each variable instead of the indeces.
    '''
    # Create boolean, True when air temperature is > 0
    t_air_gt_0 = (t_air >= 0).compute()
    # Get first and second vertical derivative
    dTdz = da.temp.differentiate("pos")
    d2Tdz2 = dTdz.differentiate("pos")
    # GENERAL: most interface detections apply some simple sanity checks like "snow temperature must be
    # below zero", "ice bottom must be below snow bottom", etc. In the following, those will not be 
    # commented.
    #
    # ice-water interface (bottom of the ice) is where temperature drops below zero for the first
    # time when starting from the bottom of the thermistor string
    iiceBot = xr.where(da.temp < 0, True, False)[:, ::-1].cumsum("pos")[:, ::-1].argmin("pos") 
    iceBot = da.pos[iiceBot].where(((iiceBot != 0) & (iiceBot >= isurf)), other=da.pos.isel(pos=isurf).values)
    # the maximum of the first derivative is located in the snow
    isnowMid = dTdz.where(((da.pos < iceBot) & (da.temp <= 0)), other=0).argmax(dim="pos")
    snowMid = da.pos[isnowMid].where(isnowMid != 0, other=da.pos.isel(pos=isurf).values)
    # the top of the snow is where the first derivative drops below 0.2 times its maximum value
    # for the first time when going upwards, starting for the middle of the snow pack
    isnowTop = xr.where(((da.pos < snowMid) & (dTdz < (0.2 * dTdz.max("pos"))) & (da.temp <= 0)), True, False)[::-1, :].cumsum("pos")[::-1, :].argmin("pos")
    snowTop = da.pos[isnowTop].where(isnowTop != 0, other=da.pos.isel(pos=isurf).values)
    # the bottom of the snow is where the second derivative drops below 0.5 times its minimum
    # value looking down
    isnowBot = xr.where(((da.pos >= snowTop) & (da.pos <= iceBot) & (da.temp <= 0) 
                         & (d2Tdz2 < (0.5 * d2Tdz2.min("pos")))), True, False).cumsum("pos").argmax("pos")
    snowBot = da.pos[isnowBot].where(isnowBot != 0, other=da.pos.isel(pos=isurf).values)
    # for the ice top, first find the maximum first derivative below the snow bottom
    iiceTop1st = dTdz.where(((da.pos >= snowBot) & (da.temp <= 0)), other=0).argmax(dim="pos")
    iceTop1st = da.pos[iiceTop1st].where(iiceTop1st != 0, other=da.pos.isel(pos=isurf).values)
    # the mimimum of the second derivative below the maximum of the first derivative
    # is the ice top (in most cases snow bottom == ice top)
    iiceTop = d2Tdz2.where(((da.pos >= iceTop1st) & (da.temp <= 0)), other=0).argmin(dim="pos")
    iceTop = da.pos[iiceTop].where(iiceTop != 0, other=da.pos.isel(pos=isurf).values)
    # redo for air temperatures >= 0 because most of the gradients reverse in that case
    # but the algorithm still works quite well unless temperature is close to zero
    isnowMid[t_air_gt_0] = dTdz.where(((da.pos < iceBot) & (da.temp <= 0)), other=0).argmin(dim="pos")[t_air_gt_0]
    snowMid = da.pos[isnowMid].where(isnowMid != 0, other=da.pos.isel(pos=isurf).values)
    isnowTop[t_air_gt_0] = xr.where(((da.pos < snowMid) & (dTdz > (0.2 * dTdz.min("pos"))) & (da.temp <= 0)), True, False)[::-1, :].cumsum("pos")[::-1, :].argmin("pos")[t_air_gt_0]
    snowTop = da.pos[isnowTop].where(isnowTop != 0, other=da.pos.isel(pos=isurf).values)
    isnowBot[t_air_gt_0] = d2Tdz2.where(((da.pos >= snowTop) & (da.pos <= iceBot) & (da.temp <= 0)
                                        & (d2Tdz2 > (0.5 * d2Tdz2.max("pos")))), True, False).cumsum("pos").argmax("pos")[t_air_gt_0]
    snowBot = da.pos[isnowBot].where(isnowBot != 0, other=da.pos.isel(pos=isurf).values)
    iiceTop1st[t_air_gt_0]  = dTdz.where(((da.pos >= snowBot) & (da.temp <= 0)), other=0).argmin(dim="pos")[t_air_gt_0] 
    iceTop1st = da.pos[iiceTop1st].where(iiceTop1st != 0, other=da.pos.isel(pos=isurf).values)
    iiceTop[t_air_gt_0]  = d2Tdz2.where(((da.pos >= iceTop1st) & (da.temp <= 0)), other=0).argmax(dim="pos")[t_air_gt_0] 
    iceTop = da.pos[iiceTop].where(iiceTop != 0, other=da.pos.isel(pos=isurf).values)
    # set everything to NaN before the freezing date
    iiceBot = iiceBot.where(iiceBot.time > np.datetime64(frozendate), other=np.nan)
    iceBot = iceBot.where(iceBot.time > np.datetime64(frozendate), other=np.nan)
    iiceTop = iiceTop.where(iiceTop.time > np.datetime64(frozendate), other=np.nan)
    iceTop = iceTop.where(iceTop.time > np.datetime64(frozendate), other=np.nan)
    isnowBot = isnowBot.where(isnowBot.time > np.datetime64(frozendate), other=np.nan)
    snowBot = snowBot.where(snowBot.time > np.datetime64(frozendate), other=np.nan)
    isnowMid = isnowMid.where(isnowMid.time > np.datetime64(frozendate), other=np.nan)
    snowMid = snowMid.where(snowMid.time > np.datetime64(frozendate), other=np.nan)
    isnowTop = isnowTop.where(isnowTop.time > np.datetime64(frozendate), other=np.nan)
    snowTop = snowTop.where(snowTop.time > np.datetime64(frozendate), other=np.nan)
    # get d/dt of ice top, if interface changes too much within one day, the algorithm likely failed and we take the values from the day before or after
    ddticeTop = np.zeros(len(iceTop))
    ddticeTop[1::] = iceTop.values[1::] - iceTop.values[0:-1]
    t = 1
    while t < len(iceTop):
        if np.abs(ddticeTop[t]) > 20:
            if t == len(iceTop) - 1:
                iceTop[t] = iceTop[t-1]
                t += 1
            else:
                tt = t
                while ((iceTop[tt] > iceTop[t-1] + 20) | (iceTop[tt] < iceTop[t-1] - 20)):
                    tt += 1
                iceTop[t:tt] = np.interp(np.arange(t, tt), [t-1, tt], [iceTop[t-1].values, iceTop[tt].values])
                t = tt + 1
        else:
            t += 1
    # redo snow bottom and top after corrections to ice top (they cannot be below ice top)
    snowTop[(snowTop > iceTop).compute()] = iceTop[(snowTop > iceTop).compute()]
    snowBot[(snowBot > iceTop).compute()] = iceTop[(snowBot > iceTop).compute()]
    # set everything to NaN where temperature is NaN near the surface
    if frozendate != 0:
        yaxmax = int(np.max(isurf + iiceBot - isnowTop).values)
        if np.min(isnowTop) < 5:
            yaxmin = 5
        else:
            yaxmin = int(np.min(isnowTop).values)
    else:
        yaxmax = isurf+10
        yaxmin = isurf-10
    isnowTop = isnowTop.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    isnowMid = isnowMid.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    isnowBot = isnowBot.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    iiceTop = iiceTop.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    iiceBot = iiceBot.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    snowTop = snowTop.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    snowMid = snowMid.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    snowBot = snowBot.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    iceTop = iceTop.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    iceBot = iceBot.where(~np.isnan(da.temp.isel(pos=slice(yaxmin-5, yaxmax+6)).mean("pos")), other=np.nan)
    # create dictionaries with the computed variables
    idx = {"snowTop": isnowTop, "snowMid": isnowMid, "snowBot": isnowBot, "iceTop": iiceTop, "iceBot": iiceBot}
    dep = {"snowTop": snowTop, "snowMid": snowMid, "snowBot": snowBot, "iceTop": iceTop, "iceBot": iceBot}
    return idx, dep

# download satellite image for upper plot is not present
for fr in [True, False]:
    if fr:
        if not os.path.isfile(plotpath + "/" + imei + "_location_fr.png"):
            with open(plotpath + "/" + imei + '_locdate.json', 'r') as ld:
                locdate = json.load(ld)
            box = (locdate["lon1"], locdate["lat1"], locdate["lon2"], locdate["lat2"])
            resolution = 20
            bbox = sh.BBox(bbox=box, crs=sh.CRS.WGS84)
            size = sh.bbox_to_dimensions(bbox, resolution=resolution)
            evalscript_true_color = """
                //VERSION=3
                
                function setup() {
                    return {
                        input: [{
                            bands: ["B02", "B03", "B04"]
                        }],
                        output: {
                            bands: 3
                        }
                    };
                }
            
                function evaluatePixel(sample) {
                    return [sample.B04, sample.B03, sample.B02];
                }
            """
            request_true_color = sh.SentinelHubRequest(
                evalscript=evalscript_true_color,
                input_data=[
                    sh.SentinelHubRequest.input_data(
                        data_collection=sh.DataCollection.SENTINEL2_L2A.define_from(
                            "s2l2a", service_url=config.sh_base_url
                        ),
                        time_interval=(locdate["date1"], locdate["date2"])
                    )
                ],
                responses=[
                    {
                        "identifier": "default",
                        "format": {"type": "image/png"},
                    },
                ],
                bbox=bbox,
                size=size,
                config=config,
                )
            response = request_true_color.get_data()
            true_color_imgs = response[0]
            dpi = 80
            height, width, nbands = true_color_imgs.shape
            figsize = width / float(dpi), height / float(dpi)
            fig = plt.figure(figsize=figsize)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.axis('off')
            ax.imshow(true_color_imgs/255*10, interpolation='nearest')
            ax.text(locdate["xname"], locdate["yname"], locdate["name"], fontsize=32, color="snow", bbox={"edgecolor": "gray", "facecolor": "gray", "alpha": 0.8})
            ax.plot(locdate["xmark"], locdate["ymark"], marker=">", markersize=36, color="orangered")
            ax.text(locdate["xmark"] - 30, locdate["ymark"] + 1, "La bouée", fontsize=26, fontweight="bold", color="dimgray", ha="right", va="center")
            r = Rectangle((locdate["xmark"] - 188, locdate["ymark"] - 18), 166, 36,  edgecolor="orangered", facecolor="snow", linewidth=5)
            ax.add_patch(r)
            plt.savefig(plotpath + "/" + imei + '_location_fr.png', dpi=dpi)
    else:
        if not os.path.isfile(plotpath + "/" + imei + "_location.png"):
            with open(plotpath + "/" + imei + '_locdate.json', 'r') as ld:
                locdate = json.load(ld)
            box = (locdate["lon1"], locdate["lat1"], locdate["lon2"], locdate["lat2"])
            resolution = 20
            bbox = sh.BBox(bbox=box, crs=sh.CRS.WGS84)
            size = sh.bbox_to_dimensions(bbox, resolution=resolution)
            evalscript_true_color = """
                //VERSION=3
                
                function setup() {
                    return {
                        input: [{
                            bands: ["B02", "B03", "B04"]
                        }],
                        output: {
                            bands: 3
                        }
                    };
                }
            
                function evaluatePixel(sample) {
                    return [sample.B04, sample.B03, sample.B02];
                }
            """
            request_true_color = sh.SentinelHubRequest(
                evalscript=evalscript_true_color,
                input_data=[
                    sh.SentinelHubRequest.input_data(
                        data_collection=sh.DataCollection.SENTINEL2_L2A.define_from(
                            "s2l2a", service_url=config.sh_base_url
                        ),
                        time_interval=(locdate["date1"], locdate["date2"])
                    )
                ],
                responses=[
                    {
                        "identifier": "default",
                        "format": {"type": "image/png"},
                    },
                ],
                bbox=bbox,
                size=size,
                config=config,
                )
            response = request_true_color.get_data()
            true_color_imgs = response[0]
            dpi = 80
            height, width, nbands = true_color_imgs.shape
            figsize = width / float(dpi), height / float(dpi)
            fig = plt.figure(figsize=figsize)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.axis('off')
            ax.imshow(true_color_imgs/255*10, interpolation='nearest')
            ax.text(locdate["xname"], locdate["yname"], locdate["name"], fontsize=32, color="snow", bbox={"edgecolor": "gray", "facecolor": "gray", "alpha": 0.8})
            ax.plot(locdate["xmark"], locdate["ymark"], marker=">", markersize=36, color="orangered")
            ax.text(locdate["xmark"] - 30, locdate["ymark"] + 1, "Buoy location", fontsize=26, fontweight="bold", color="dimgray", ha="right", va="center")
            r = Rectangle((locdate["xmark"] - 268, locdate["ymark"] - 18), 246, 36,  edgecolor="orangered", facecolor="snow", linewidth=5)
            ax.add_patch(r)
            plt.savefig(plotpath + "/" + imei + '_location.png', dpi=dpi)

now = datetime.now()
thisyear = now.year
thismonth = now.month

# we cut the first two time steps becaues for this buoy they were tests in a different location
try:
    data_orig = xr.open_dataset(datapath + "/" + imei + ".nc").isel(time=slice(2, None))
    no_data = False
    if thismonth > 6:
        data_in = data_orig.sel(time=slice(str(now.year) + "-07-01", None))
    else:
        data_in = data_orig.sel(time=slice(str(now.year - 1) + "-07-01", None))
    if len(data_in.time) == 0:
        no_data = True
except:
    no_data = True

dg = "#4f4f4f"
ddg = "#333333"

if no_data:
    for fr in [True, False]:
        plt.figure(figsize=(10, 4))
        if fr:
            plt.text(0.5, 0.5, "À venir", ha="center", va="center", fontweight="bold", fontsize=32, color=dg)
        else:
            plt.text(0.5, 0.5, "Coming soon", ha="center", va="center", fontweight="bold", fontsize=32, color=dg)
        plt.gca().axis("off")
        if fr:
            plt.savefig(plotpath + "/plot_" + imei + "_top_fr.png", dpi=300)
        else:
            plt.savefig(plotpath + "/plot_" + imei + "_top.png", dpi=300)
else:
    t_air_in = data_in.temp.isel(pos=-1)
    t_air_0 = (t_air_in.groupby("time.dayofyear") - t_air_in.groupby("time.dayofyear").min()).compute()
    t_air = t_air_in.where(t_air_0.drop_vars(["pos", "dayofyear"])==0, drop=True)
    data = data_in.where(t_air_0.drop_vars(["pos", "dayofyear"])==0, drop=True).isel(pos=slice(0, -1))
    current_t_air = data_in.temp.isel(pos=-1, time=-1)
    current_datetime = str(str(data_in.time.isel(time=-1).values)[0:10] + " "
                       + str(data_in.time.isel(time=-1).values)[11:16] + "H")
    #
    if os.path.isfile(plotpath + "/" + "frozen"):
        with open(plotpath + "/" + "frozen", "r") as f:
            froze = f.read()
        f.close()
    else:
        with open(plotpath + "/" + "frozen", "w") as f:
            f.write("False")
        f.close()
        froze = "False"
    #
    if froze == "True":
        frozen = True
        if thismonth > 6:
            if os.path.isfile(plotpath + "/" + "isurf" + str(thisyear) + "-" + str(thisyear + 1)):
                with open(plotpath + "/" + "isurf" + str(thisyear) + "-" + str(thisyear + 1), "r") as f:
                    isurf = int(f.read())
                f.close()
        else:
            if os.path.isfile(plotpath + "/" + "isurf" + str(thisyear - 1) + "-" + str(thisyear)):
                with open(plotpath + "/" + "isurf" + str(thisyear - 1) + "-" + str(thisyear), "r") as f:
                    isurf = int(f.read())
                f.close()
    else:
        frozen = False
    #
    if not frozen:
        isurf = int(data.temp.std("time").differentiate("pos").argmin("pos").values + 2)
        if thismonth > 6:
            with open(plotpath + "/" + "isurf" + str(thisyear) + "-" + str(thisyear + 1), "w") as f:
                f.write(str(isurf))
            f.close()
        else:
            with open(plotpath + "/" + "isurf" + str(thisyear - 1) + "-" + str(thisyear), "w") as f:
                f.write(str(isurf))
            f.close()
        surftemp = data.temp.isel(pos=isurf)
        if surftemp.isel(time=-1) <= 0:
            frozen = True
            with open(plotpath + "/" + "frozen", "w") as f:
                f.write(str(frozen))
            f.close()
            with open(plotpath + "/" + "frozendate", "w") as f:
                f.write(str(data.time.isel(time=-1).values)[0:10])
            f.close()
    else:
        if thismonth > 6:
            with open(plotpath + "/" + "isurf" + str(thisyear) + "-" + str(thisyear + 1), "r") as f:
                isurf = int(f.read())
            f.close()
        else:
            with open(plotpath + "/" + "isurf" + str(thisyear - 1) + "-" + str(thisyear), "r") as f:
                isurf = int(f.read())
            f.close()
    #
    if not frozen:
        da = data #.isel(time=slice(-14, None))
        t_a = t_air #.isel(time=slice(-14, None))
        frozendate = 0
    else:
        with open(plotpath + "/" + "frozendate", "r") as f:
            frozendate = f.read()
        f.close()
        da = data #.sel(time=slice(str((np.datetime64(frozendate) - np.timedelta64(2, "W"))), None))
        t_a = t_air #.sel(time=slice(str((np.datetime64(frozendate) - np.timedelta64(2, "W"))), None))
    #
    if frozen:
        idx, dep = detect_interfaces(da, t_a, isurf, frozendate)
    else:
        dummy = (np.ones(len(da.time)) * da.pos.isel(pos=isurf).values).astype(int)
        idx = {"snowTop": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "snowMid": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "snowBot": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "iceTop": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "iceBot": xr.DataArray(dummy, dims=["time"], coords={"time": da.time})}
        dep = {"snowTop": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "snowMid": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "snowBot": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "iceTop": xr.DataArray(dummy, dims=["time"], coords={"time": da.time}),
               "iceBot": xr.DataArray(dummy, dims=["time"], coords={"time": da.time})}
    # put data into json for use with the echarts on the website
    edata = {"twater": {}, "tair": {}, "icethick": {}, "snowicethick": {}, "totalicethick": {}, "snowthick": {}}
    colormap = {"cmap1": {}, "cmap2": {}, "cmap3": {}, "cmap4": {}}
    cmap1 = cmo.amp
    cmap2 = cmo.dense
    cmap3 = plt.cm.Blues
    cmap4 = plt.cm.Greys
    # use a random, non-leap year to create a time axis for the plot. only the
    # month and day will be used later, the year is unimportant
    climtime = pd.date_range(
        start=datetime.strptime("2001 182", "%Y %j"),
        end=datetime.strptime("2002 181", "%Y %j"),
        freq="1D"
        )
    edata['date'] = [str(climtime[i])[5:10] for i in np.arange(0, len(climtime))]
    firstyear = int(data_in.time.isel(time=0).dt.year.values)
    # depending on the month we are in, a different year will be needed to define
    # the colormap
    if thismonth > 6:
        tyear = thisyear + 1
    else:
        tyear = thisyear
    if firstyear == tyear:
        edata['years'] = str(firstyear) + "/" + str(firstyear + 1)
    else:
        edata['years'] = [str(i) + "/" + str(i + 1) for i in np.arange(firstyear, tyear)]
    # define range to read colors from colormap in the range (0.2, 0.8), avoiding
    # values too close to 0 or 1 that would result in very light or dark colors and
    # thus not be visible on the graph (depending on light- or dark-mode)
    cpos = np.linspace(0.4, 0.8, tyear - firstyear)
    # extract data for each season (Jul. 1 - Jun. 30)
    for y in np.arange(firstyear, tyear):
        year = str(y)
        if os.path.isfile(plotpath + "/" + "isurf" + str(y) + "-" + str(y + 1)):
            with open(plotpath + "/" + "isurf" + str(y) + "-" + str(y + 1), "r") as f:
                isurf = int(f.read())
            f.close()
        else:
            print("Error: " + "isurf" + str(y) + "-" + str(y + 1) + " not found")
            isurf = 20
        # start one day later in leap years, i.e. keeping the delta t to the last
        # day of the year the same
        if y % 4 == 0 and (y % 100 != 0 or y % 400 == 0):
            ctime = pd.date_range(
                start=datetime.strptime(year + " " + str(183), "%Y %j"),
                end=datetime.strptime(str(int(year) + 1) + " " + str(181),
                                               "%Y %j"),
                freq="1D"
                )
        else:
            ctime = pd.date_range(
                start=datetime.strptime(year + " " + str(182), "%Y %j"),
                end=datetime.strptime(str(int(year) + 1) + " " + str(181),
                                               "%Y %j"),
                freq="1D"
                )
        #
        icethick_tmp = xr.merge([(dep["iceBot"] - da.pos.isel(pos=isurf).values).resample(time="1D").mean().sel(time=slice(ctime[0], ctime[-1])).rename("pos"),
                                 xr.DataArray(np.zeros(365) + np.nan, dims=["time"], coords={"time": ctime}, name="pos")], compat="override", join="outer")
        #
        snowicethick_tmp = xr.merge([(da.pos.isel(pos=isurf).values - dep["snowBot"]).resample(time="1D").mean().sel(time=slice(ctime[0], ctime[-1])).rename("pos"),
                                 xr.DataArray(np.zeros(365) + np.nan, dims=["time"], coords={"time": ctime}, name="pos")], compat="override", join="outer")
        snowicethick_tmp = snowicethick_tmp.where(((snowicethick_tmp.pos >= 0) | (np.isnan(snowicethick_tmp.pos))), other=0)
        #
        totalicethick_tmp = list(icethick_tmp.pos.values + snowicethick_tmp.pos.values)
        #
        snow_tmp = (dep["snowBot"] - dep["snowTop"]).where(idx["snowBot"] < isurf, other=(da.pos.isel(pos=isurf).values - dep["snowTop"]))
        snowthick_tmp = xr.merge([snow_tmp.resample(time="1D").mean().sel(time=slice(ctime[0], ctime[-1])).rename("pos"),
                                  xr.DataArray(np.zeros(365) + np.nan, dims=["time"], coords={"time": ctime}, name="pos")], compat="override", join="outer")
        snowthick_tmp = snowthick_tmp.where(((snowthick_tmp.pos >= 0) | (np.isnan(snowthick_tmp.pos))), other=0)
        #
        edata["totalicethick"][str(y) + "/" + str(y + 1)] = ["null" if np.isnan(x) else x for x in totalicethick_tmp]
        edata["icethick"][str(y) + "/" + str(y + 1)] = ["null" if np.isnan(x) else x for x in list(icethick_tmp.pos.values)]
        edata["snowicethick"][str(y) + "/" + str(y + 1)] = ["null" if np.isnan(x) else x for x in list(snowicethick_tmp.pos.values)]
        edata["snowthick"][str(y) + "/" + str(y + 1)] = ["null" if np.isnan(x) else x for x in list(snowthick_tmp.pos.values)]
        #
        twater_tmp1 = xr.merge([data_in.isel(pos=isurf+5).resample(time="1D").mean().sel(time=slice(ctime[0], ctime[-1])),
                               xr.DataArray(np.zeros(365) + np.nan, dims=["time"], coords={"time": ctime}, name="temp")], compat="override", join="outer")
        twater_tmp = twater_tmp1.where(twater_tmp1.temp >= 0, other=0)
        twater_tmp = twater_tmp.where(~np.isnan(twater_tmp1), other=np.nan)
        edata["twater"][str(y) + "/" + str(y + 1)] = ["null" if np.isnan(x) else x for x in list(twater_tmp.temp.values)]
        #
        tair_tmp = xr.merge([t_air_in.resample(time="1D").mean().sel(time=slice(ctime[0], ctime[-1])),
                             xr.DataArray(np.zeros(365) + np.nan, dims=["time"], coords={"time": ctime}, name="temp")], compat="override", join="outer")
        edata["tair"][str(y) + "/" + str(y + 1)] = ["null" if np.isnan(x) else x for x in list(tair_tmp.temp.values)]
        # add the colors from cmocean colormap to the color data
        colormap["cmap1"][str(y) + "/" + str(y+1)] = '#%02x%02x%02x' % tuple([
            int(cmap1(cpos[y-firstyear])[i] * 255) for i in [0, 1, 2]
            ])
        colormap["cmap2"][str(y) + "/" + str(y+1)] = '#%02x%02x%02x' % tuple([
            int(cmap2(cpos[y-firstyear])[i] * 255) for i in [0, 1, 2]
            ])
        colormap["cmap3"][str(y) + "/" + str(y+1)] = '#%02x%02x%02x' % tuple([
            int(cmap3(cpos[y-firstyear])[i] * 255) for i in [0, 1, 2]
            ])
        colormap["cmap4"][str(y) + "/" + str(y+1)] = '#%02x%02x%02x' % tuple([
            int(cmap4(cpos[y-firstyear])[i] * 255) for i in [0, 1, 2]
            ])
    with open(plotpath + "/" + imei + '_edata.json', 'w') as fp:
        json.dump(edata, fp)
    with open(plotpath + "/" + imei + '_colormap.json', 'w') as fp:
        json.dump(colormap, fp)

    if frozen:
        yaxmax = int(np.max(isurf + idx["iceBot"] - idx["snowTop"]).values)
        if np.min(idx["snowTop"]) < 5:
            yaxmin = 5
        else:
            yaxmin = int(np.min(idx["snowTop"]).values)
    else:
        yaxmax = isurf + 10
        yaxmin = isurf - 10
        if yaxmin < 5:
            yaxmin = 5
    #
    yax = -(da.pos.isel(pos=slice(yaxmin-5, yaxmax+6)) - da.pos.isel(pos=isurf))
    lastIdx = int(len(idx["iceBot"]) - 1 - (idx["iceBot"] / idx["iceBot"]).where(~np.isnan(idx["iceBot"]), other=0)[::-1].argmax("time"))
    if frozen:
        current_t_water = data_in.temp.isel(time=-1).isel(pos=slice(int(idx["iceBot"][lastIdx].values), yaxmax+6)).mean("pos")
    else:
        current_t_water = data_in.temp.isel(time=-1).isel(pos=slice(isurf, yaxmax+6)).mean("pos")
    # fill in some NaNs when there are too large time gaps in the data so the plots look nicer
    dayinns = 3600 * 24 * 1e9
    dtdt = np.array([int(data_in.time.diff("time").values[i]) for i in range(0, len(data_in.time)-1)])
    if (dtdt >= 7 * dayinns).any():
        #create time axis
        newtime = []
        da_tmp = data_in.copy()
        for t in range(0, len(data_in.time)-1):
            newtime.append(data_in.time[t].values)
            if dtdt[t] >= 7 * dayinns:
                filltime = np.arange(data_in.time[t].values + np.timedelta64(6, "h"), data_in.time[t+1].values, np.timedelta64(6, "h"))
                newtime = newtime + [filltime[j] for j in range(0, len(filltime))]
        da_tmp = da_tmp.interp(time=newtime, method="nearest")
        da_tmp = da_tmp.where(da_tmp.time.isin(data_in.time), other=np.nan)
    else:
        da_tmp = data_in
    #
    water_evo = da_tmp.temp.isel(pos=slice(isurf+1, yaxmax+6)).values.T
    water_evo_yax = yax.sel(pos=slice(data_in.pos.isel(pos=isurf+1), None))
    water_max = np.floor(np.nanmax(water_evo))
    water_bounds = np.linspace(0, water_max, 25)
    water_norm = colors.BoundaryNorm(boundaries=water_bounds, ncolors=256)
    #
    air_evo = da_tmp.temp.isel(pos=slice(yaxmin-5, isurf)).values.T
    air_evo_yax = yax.sel(pos=slice(None, data_in.pos.isel(pos=isurf-1)))
    air_min = np.ceil(np.nanmin(air_evo))
    air_max = np.floor(np.nanmax(air_evo))
    air_lim = np.max(np.abs([air_min, air_max]))
    air_bounds = np.linspace(-air_lim, air_lim, 31)
    air_norm = colors.BoundaryNorm(boundaries=air_bounds, ncolors=256)
    #
    shadow_effect = [
        path_effects.Stroke(linewidth=2, foreground="linen"),
        path_effects.Normal()
    ]
    # plot the figures
    for fr in [True, False]:
        # top figure first
        fig = plt.figure(figsize=(9.2, 5))
        gs = fig.add_gridspec(3, 4, width_ratios=[5, 1.0, 0.1, 1], height_ratios=[6, 1, 6])
        ax1 = fig.add_subplot(gs[0:3, 0])
        caxa = fig.add_subplot(gs[0, 1])
        caxw = fig.add_subplot(gs[2, 1])
        ax2 = fig.add_subplot(gs[0:3, 3])
        if fr:
            locim = plt.imread(plotpath + "/" + imei + '_location_fr.png')
            air = "air"
            water = "eau"
            snow = "neige"
            ice = "glace"
            snowice = "neige gelée"
            snowice2 = "neige\ngelée"
            totalice = "glace totale"
        else:
            locim = plt.imread(plotpath + "/" + imei + '_location.png')
            air = "air"
            water = "water"
            snow = "snow"
            ice = "ice" 
            snowice= "snow-ice"
            snowice2 = "snow\n-ice"
            totalice = "total ice"
        ax1.imshow(locim, interpolation='nearest')
        if frozen:
            idep = -(dep["iceBot"] - da.pos.isel(pos=isurf))
            sidep = -(dep["snowBot"] - da.pos.isel(pos=isurf))
            sidep = sidep.where(sidep >= 0, 0)
            sidepnan = sidep.where(~np.isnan(idep), np.nan)
            tidep = idep - sidepnan
            sdep = -(dep["snowTop"] - da.pos.isel(pos=isurf)) - sidepnan
            ax2.fill_between([0, 1], np.min(yax), y2=-1, color=cmo.thermal(water_norm(0)))
            ax2.fill_between([0, 1], tidep[lastIdx], y2=-sidep[lastIdx], color="cornflowerblue", zorder=5)
            ax2.fill_between([0, 1], -sidep[lastIdx], y2=0, color="skyblue", zorder=4)
            ax2.fill_between([0, 1], 0, y2=sdep[lastIdx], color="ghostwhite", zorder=3)
            ax2.fill_between([0, 1], 1, y2=np.max(yax), color=cmo.balance(air_norm(current_t_air)))
            airthick = dep["snowTop"][lastIdx]
            if airthick < 5:
                ax2.annotate("", xytext=(1.55, -(0 - da.pos.isel(pos=isurf)) - (airthick / 2)),
                             xy=(1.0, -(0 - da.pos.isel(pos=isurf)) - (airthick / 2)),
                             arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                ax2.text(1.6, -(0 - da.pos.isel(pos=isurf)) - (airthick / 2), air,
                         ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
            else:
                ax2.text(0.5, np.max(yax) - ((np.max(yax) - np.min(yax)) / 95), air,
                         ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6, path_effects=shadow_effect)
            ax2.text(0.5, np.max(yax) + ((np.max(yax) - np.min(yax)) / 14), str(np.around(current_t_air.values, decimals=1)) + r"$^{\circ}$C",
                     ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6)
            #
            #
            if idx["snowBot"][lastIdx] < isurf:
                snowthick = (dep["snowBot"][lastIdx] - dep["snowTop"][lastIdx]).values
            else:
                snowthick = (da.pos.isel(pos=isurf) - dep["snowTop"][lastIdx]).values
            if snowthick > 0:
                if snowthick < 5:
                    ax2.annotate("", xytext=(1.55, snowthick - (snowthick / 2)), 
                                 xy=(1.0, snowthick - (snowthick / 2)), 
                                 arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                    ax2.text(1.6, -(dep["snowTop"] - da.pos.isel(pos=isurf))[-1] - (snowthick / 2), str(int(snowthick)) + r"$\,$cm" + " " + snow, 
                             ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
                else:
                    ax2.text(0.5, snowthick - (snowthick / 2), snow, 
                             ha="center", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
                    ax2.annotate("", xytext=(1.55, snowthick - (snowthick / 2)), 
                                 xy=(1.0, snowthick - (snowthick / 2)), 
                                 arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                    ax2.text(1.6, snowthick - (snowthick / 2), str(int(snowthick)) + r"$\,$cm", 
                             ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
            #
            snowicethick = (da.pos.isel(pos=isurf) - dep["snowBot"][lastIdx]).values
            if snowicethick > 0:
                if snowicethick < 10:
                    ax2.annotate("", xytext=(1.55, -snowicethick + (snowicethick / 2)),
                                 xy=(1.0, -snowicethick + (snowicethick / 2)),
                                arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                    ax2.text(1.6, -snowicethick + (snowicethick / 2), str(int(snowicethick)) + r"$\,$cm" + " " + snowice,
                             ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
                else:
                    ax2.text(0.5, -snowicethick + (snowicethick / 2), snowice2,
                             ha="center", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
                    ax2.annotate("", xytext=(1.55, -snowicethick + (snowicethick / 2)),
                                 xy=(1.0, -snowicethick + (snowicethick / 2)),
                                 arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                    ax2.text(1.6, -snowicethick + (snowicethick / 2), str(int(snowicethick)) + r"$\,$cm",
                             ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
            #
            #
            icethick = (dep["iceBot"][lastIdx] - da.pos.isel(pos=isurf)).values
            if icethick < 5:
                ax2.annotate("", xytext=(1.55, -(snowicethick + icethick) + (icethick / 2)),
                             xy=(1.0, -(snowicethick + icethick) + (icethick / 2)),
                             arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                ax2.text(1.6, -(snowicethick + icethick) + (icethick / 2), str(int(icethick)) + r"$\,$cm" + " " + ice,
                         ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
            else:
                ax2.text(0.5, -(snowicethick + icethick) + (icethick / 2), ice,
                         ha="center", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
                ax2.annotate("", xytext=(1.55, -(snowicethick + icethick) + (icethick / 2)),
                             xy=(1.0, -(snowicethick + icethick) + (icethick / 2)),
                             arrowprops=dict(arrowstyle="->", color=dg), zorder=6)
                ax2.text(1.6, -(snowicethick + icethick) + (icethick / 2), str(int(icethick)) + r"$\,$cm",
                         ha="left", va="center", fontsize=14, fontweight="bold", color=dg, zorder=6)
            totalicethick = icethick + snowicethick
            ax2.annotate("", xytext=(-0.5, -(snowicethick + icethick) + ((snowicethick + icethick) / 2)),
                         xy=(0.0, 0),
                         arrowprops=dict(arrowstyle="-", color=dg), zorder=6)
            ax2.annotate("", xytext=(-0.5, -(snowicethick + icethick) + ((snowicethick + icethick) / 2)),
                         xy=(0.0, -(snowicethick + icethick)),
                         arrowprops=dict(arrowstyle="-", color=dg), zorder=6)
            ax2.text(-0.55, -(snowicethick + icethick) + ((snowicethick + icethick) / 2.3), totalice + ":\n" + r"$\mathdefault{\bf{" + str(snowicethick + icethick) + "\,cm" + "}}$",
                     ha="right", va="bottom", fontsize=14, color=dg)
            ax2.text(0.5, np.min(yax) + ((np.max(yax) - np.min(yax)) / 14), water,
                     ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6, path_effects=shadow_effect)
            ax2.text(0.5, np.min(yax) - ((np.max(yax) - np.min(yax)) / 45), str(np.around(0.0, decimals=1)) + r"$^{\circ}$C",
                     ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6)
        else:
            ax2.fill_between([0, 1], np.min(yax), y2=-1, color=cmo.thermal(water_norm(current_t_water)))
            ax2.fill_between([0, 1], 1, y2=np.max(yax), color=cmo.balance(air_norm(current_t_air)))
            ax2.text(0.5, np.max(yax) - ((np.max(yax) - np.min(yax)) / 45), air, ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6, path_effects=shadow_effect)
            ax2.text(0.5, np.max(yax) + ((np.max(yax) - np.min(yax)) / 14), str(np.around(current_t_air.values, decimals=1)) + r"$^{\circ}$C",
                     ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6)
            ax2.text(0.5, np.min(yax) + ((np.max(yax) - np.min(yax)) / 14), water, ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6, path_effects=shadow_effect)
            ax2.text(0.5, np.min(yax) - ((np.max(yax) - np.min(yax)) / 45), str(np.around(current_t_water.values, decimals=1)) + r"$^{\circ}$C",
                     ha="center", va="top", fontsize=14, fontweight="bold", color=dg, zorder=6)
        ax2.text(0.5, np.min(yax) - ((np.max(yax) - np.min(yax))/ 7), current_datetime[0:10] + "\n" + current_datetime[11:17],
                 ha="center", va="top", fontsize=12, fontweight="bold", color=dg, zorder=6)
        ax2.set_ylim(np.min(yax), np.max(yax))
        ax2.axis('off')
        ax1.axis('off')
        caxw.axis('off')
        caxa.axis('off')
        if fr:
            fig.text(0.05, 0.9, "Épaisseur de glace actuelle estimée: " + r"$\mathdefault{\bf{" + str(snowicethick + icethick) + "\,cm" + "}}$", fontsize=24, ha="left", color="#525252",
                     bbox={"facecolor": "w", "edgecolor": "#787878", "linewidth": 3, "pad": 10})
        else:
            fig.text(0.05, 0.9, "Current estimated ice thickness: " + r"$\mathdefault{\bf{" + str(snowicethick + icethick) + "\,cm" + "}}$", fontsize=24, ha="left", color="#525252",
                     bbox={"facecolor": "w", "edgecolor": "#787878", "linewidth": 3, "pad": 10})
        plt.subplots_adjust(wspace=0.1, left=0.05, right=0.75, bottom=0.2, top=0.77)
        if fr:
            plt.savefig(plotpath + "/" + imei + "_top_fr.png", dpi=300)
        else:
            plt.savefig(plotpath + "/" + imei + "_top.png", dpi=300)
        # bottom figure
        fig = plt.figure(figsize=(10, 5))
        gs = fig.add_gridspec(3, 2, width_ratios=[6, 0.2], height_ratios=[6, 1.8, 6])
        ax1 = fig.add_subplot(gs[0:3, 0])
        caxa = fig.add_subplot(gs[0, 1])
        caxw = fig.add_subplot(gs[2, 1])
        if frozen:
            idep = -(dep["iceBot"] - da.pos.isel(pos=isurf))
            sidep = -(dep["snowBot"] - da.pos.isel(pos=isurf))
            sidep = sidep.where(sidep >= 0, 0)
            sidepnan = sidep.where(~np.isnan(idep), np.nan)
            tidep = idep - sidepnan
            sdep = -(dep["snowTop"] - da.pos.isel(pos=isurf)) - sidepnan
            airzero = (sidep.copy().interp(time=da_tmp.time) * 0).bfill("time")
            airzero = airzero.where(airzero.time <= np.datetime64(frozendate), other=sidepnan.interp(time=da_tmp.time).ffill("time").bfill("time"))
            a1w = ax1.pcolormesh(da_tmp.time, water_evo_yax, water_evo, cmap=cmo.thermal, norm=water_norm)
            plt.colorbar(a1w, ax=ax1, cax=caxw)
            f1 = ax1.fill_between(da.time, 0, y2=sdep.ffill("time"), color="ghostwhite", zorder=3, label=snow)
            f2 = ax1.fill_between(da.time, -sidepnan.ffill("time"), y2=0, color="skyblue", zorder=5, label=snowice)
            f3 = ax1.fill_between(da.time, tidep.ffill("time"), y2=0, color="cornflowerblue", zorder=4, label=ice)
            a1a = ax1.pcolormesh(da_tmp.time, air_evo_yax - airzero, air_evo, cmap=cmo.balance, norm=air_norm)
            plt.colorbar(a1a, ax=ax1, cax=caxa)
            ax1.legend(ncols=1, facecolor="grey", loc="lower left", bbox_to_anchor=(0.997, 0.385))
        else:
            snowicethick = 0
            icethick = 0
            a1w = ax1.pcolormesh(da_tmp.time, water_evo_yax, water_evo, cmap=cmo.thermal, norm=water_norm)
            plt.colorbar(a1w, ax=ax1, cax=caxw)
            a1a = ax1.pcolormesh(da_tmp.time, air_evo_yax, air_evo, cmap=cmo.balance, norm=air_norm)
            plt.colorbar(a1a, ax=ax1, cax=caxa)
        ax1.set_ylim(np.min(yax), np.max(yax))
        if fr:
            ax1.set_ylabel("cm du surface de la glace")
        else:
            ax1.set_ylabel("cm from ice surface")
        t = ax1.get_xticks()
        l = ax1.get_xticklabels()
        ax1.set_xticks(t)
        ax1.set_xticklabels(l, rotation=45, ha="right", rotation_mode="anchor")
        ax1.set_xlim(da.time[0], da.time[-1])
        ax1.xaxis.label.set_color(dg)
        ax1.yaxis.label.set_color(dg)
        ax1.tick_params(axis='x', colors=ddg)
        ax1.tick_params(axis='y', colors=ddg)
        ax1.spines['bottom'].set_color(ddg)
        ax1.spines['left'].set_color(ddg)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        air_tick_diff = np.floor(np.mean(np.diff(np.unique(np.floor(caxa.get_yticks())))))
        air_ticks = np.arange(-((air_lim // air_tick_diff) * air_tick_diff), air_lim, air_tick_diff)
        caxa.set_yticks(air_ticks)
        caxa.set_yticklabels([int(air_ticks[i]) for i in range(0, len(air_ticks))])
        caxa.yaxis.set_minor_locator(MultipleLocator(1))
        if fr:
            caxa.set_ylabel(r"temp. d'air $^{\circ}$C")
        else:
            caxa.set_ylabel(r"air temp. $^{\circ}$C")
        water_ticks = np.arange(0, water_max, np.floor(np.mean(np.diff(np.unique(np.floor(caxw.get_yticks()))))))
        caxw.set_yticks(water_ticks)
        caxw.set_yticklabels([int(water_ticks[i]) for i in range(0, len(water_ticks))])
        caxw.yaxis.set_minor_locator(MultipleLocator(0.5))
        if fr:
            caxw.set_ylabel(r"temp. d'eau $^{\circ}$C")
        else:
            caxw.set_ylabel(r"water temp. $^{\circ}$C")
        plt.subplots_adjust(wspace=0.02, left=0.1, right=0.85, bottom=0.2, top=0.9)
        if fr:
            plt.savefig(plotpath + "/" + imei + "_bottom_fr.png", dpi=300)
        else:
            plt.savefig(plotpath + "/" + imei + "_bottom.png", dpi=300)
