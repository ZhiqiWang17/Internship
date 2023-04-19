# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:27:37 2021

@author: d17827
"""

from osgeo import gdal
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio as rio
import matplotlib.pyplot as plt
import os
import fiona
import glob
from pyproj import Transformer
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.crs import CRS
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.interpolate import interp1d
from tqdm import tqdm
import dask
#%% define directories
cwd = os.getcwd()

#%%
pilots = gpd.read_file("C:Users/d17827/Desktop/cpns_clip_VS_kisters/pilots/pilootgebieden.shp")
dijle = pilots.loc[pilots["gebied"]=="Dijle"]

# transformer from lambert to wgs
dijle = dijle.to_crs("epsg:4326")

flanders = gpd.read_file("Q:/admin/admin_grenzen_vla/admin_grenzen_19/RefgewG10.shp")
flanders = flanders.to_crs("epsg:4326")

dijle.plot()
flanders.plot()

              
#combine flanders and dijle shp in one image
base = flanders.plot(color='white', edgecolor='black')
dijle.plot(ax=base, color='red', markersize=5)

#%%
#clip copernicus geotiff to shp
in_raster_folder = "C:/Users/d17827/Desktop/cpns_clip_VS_kisters/cpns/"
out_clip_folder = "C:/Users/d17827/Desktop/cpns_clip_VS_kisters/cpns_clip/"

rasterlist=glob.glob(in_raster_folder+"/*.tiff")

for raster in rasterlist:
    date= os.path.basename(raster)[10:18]
    with rio.open(raster) as src:
        out_image, out_transform = rio.mask.mask(src, flanders.geometry,crop=True)
        out_meta = src.meta
                     
    out_meta.update({"driver": "GTiff","height": out_image.shape[1], "width": out_image.shape[2],
                     "transform": out_transform})
   
    with rio.open(out_clip_folder + f"cpns_{date}.tif", "w", **out_meta) as dest:
        dest.write(out_image)

#%%
#copernicus
cpns = rio.open("C:/Users/d17827/Desktop/cpns_clip_VS_kisters/cpns_clip/cpns_20180101.tif")
cpns=cpns.read()
print(cpns.shape)
plt.hist(cpns_data.ravel(), bins=256, range=(0, 250), fc='k', ec='k')

# imgplot = plt.imshow(cpns_data,clim=(100, 180),extent=[0,350,10,120])
cpns_data = cpns[cpns <= 200]
ecdf_cpns = ECDF(cpns_data)
# plot the cdf
plt.plot(ecdf_cpns.x, ecdf_cpns.y)
#rescaling-remove nodata
cpns_rescale = cpns.copy().astype(float).squeeze()
cpns_rescale[cpns_rescale > 200] = np.nan
plt.imshow(cpns_rescale)
#rescaling with ecdf axis 
cpns_rescale[cpns_rescale <= 200] = ecdf_cpns(cpns_rescale[cpns_rescale <= 200])
plt.imshow(cpns_rescale)
plt.colorbar(orientation='horizontal')

#kister
kister = rio.open("C:/Users/d17827/Desktop/cpns_clip_VS_kisters/kisters/image20180101T00.geotiff")
kister=kister.read()
print(kister.shape)
plt.hist(kister.ravel(), bins=256, range=(0.0, 0.4), fc='k', ec='k')

kister_data = kister[kister <= 65534]
ecdf_kister = ECDF(kister_data)
# plot the cdf
plt.plot(ecdf_kister.x, ecdf_kister.y)
#rescaling-remove nodata
kister_rescale = kister.copy().astype(float).squeeze()
kister_rescale[kister_rescale > 65534] = np.nan
plt.imshow(kister_rescale)
#rescaling with ecdf axis 
kister_rescale[kister_rescale <= 65534] = ecdf_kister(kister_rescale[kister_rescale <= 65534])
plt.imshow(kister_rescale)
plt.colorbar(orientation='horizontal')

#compare in one plot
fig, ax = plt.subplots(1, 2,figsize = (10, 5))
ax[0].imshow(cpns_rescale)
ax[0].set_title('Copernicus')
ax[0].axis('off')
ax[1].imshow(kister_rescale)
ax[1].set_title('Kister')
ax[1].axis('off')
#plt.colorbar(ax=ax[1],location='bottom')
#fig.colorbar(fig,ax=ax,location='bottom')
#plt.colorbar(plt, ax=ax,location='bottom',shrink=0.75)
#plt.colorbar(plt,orientation='horizontal', ax = ax)

# fig = plt.figure(constrained_layout=True)
# (subfig_l, subfig_r) = fig.subfigures(nrows=1, ncols=2)
# axes_l = subfig_l.subplots(nrows=1, ncols=2, sharey=True)
# for ax in axes_l:
#     im_l = ax.imshow(cpns_rescale)
#     im_r = ax.imshow(kister_rescale)
# # shared colorbar for left subfigure
# subfig_l.colorbar(im_l, ax=axes_l, location='bottom')
plt.savefig("compare_2018.01.01.png")

#histogram
fig, ax = plt.subplots(figsize = (10, 10))
dt = np.mean(np.diff(np.sort(np.unique(kister[kister < 65535]))))
dnoise = np.random.rand(np.sum(kister < 65535)) * dt - dt/2
y_kisters, x_kisters = np.histogram(kister[kister < 65535] + dnoise, bins = 60)
ax.plot((x_kisters[1:] + x_kisters[:-1])/2, y_kisters/np.sum(y_kisters), label = "Kisters",color='r')

ax1 = ax.twiny()
dt = np.mean(np.diff(np.sort(np.unique(cpns[cpns <= 200]))))
dnoise = np.random.rand(np.sum(cpns <= 200)) * dt - dt/2
y_cpns, x_cpns = np.histogram(cpns[cpns <= 200] + dnoise , bins = 60)
ax1.plot((x_cpns[1:] + x_cpns[:-1])/2, y_cpns/np.sum(y_cpns), label = "copernicus")
fig.legend()
plt.title('Histogram')
vals = np.arange(120, 170)
cnts = np.zeros(len(vals))

for i in range(len(vals)):
    cnts[i] = np.sum(cpns == vals[i])
    
plt.savefig("His_2018.01.01.png")    
#%%
#read nc file of copernicus data
pilots = gpd.read_file(cwd + "/pilots/pilootgebieden.shp")
dijle = pilots.loc[pilots["gebied"]=="Dijle"]

#%% read pydov download and only keep measurements in Dijle valley, in 2018
df1 = pd.read_csv(cwd + "/pydov/gpd_query_large_1.csv", index_col=0)
df2 = pd.read_csv(cwd + "/pydov/gpd_query_large_2.csv", index_col=0)

# stitch dataframes and put in desired format
df = pd.concat([df1, df2])
del(df1, df2)

# drop duplicates and empy values and reset the index
df = df.drop_duplicates()
df = df.dropna()
df = df.reset_index(drop=True)

# put date column in datetime format and only keep data from 2018
df["datum"] = pd.to_datetime(df["datum"])
df["year"] = df["datum"].apply(lambda x: x.year)
df = df[df['year'] == 2018]

# %% only keep measurements/piezometers in dijle valley
# transform to geodataframe and pass the coordinate system of mask layer
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y))
gdf = gdf.set_crs("epsg:31370")

# clip geodataframe to the dijle shapefile
gdf = gpd.clip(gdf, dijle.geometry)

#%% transform coordinates from lambert to wgs
transformer = Transformer.from_crs('epsg:31370', 'epsg:4326')
gdf["lat"], gdf["lon"] = transformer.transform(gdf["x"].tolist(),
                                               gdf["y"].tolist())

#read copernicus nc file 
SWI = xr.open_dataset(cwd + "/copernicus/SWI.nc")
dates = np.unique(gdf["datum"].tolist())
depth = 0.60
#date = dates[0]
for date in tqdm(dates):
    tmpSWI = SWI.SWI.sel(depth = depth, 
                  lon = xr.DataArray(gdf[gdf['datum'] == date].lon, dims = "points"),
                  lat = xr.DataArray(gdf[gdf['datum'] == date].lat, dims = "points"),
                  time = date, method = 'ffill')
    gdf.loc[gdf['datum'] == date, "SWI"] = tmpSWI.values
#some_variable = tmpSWI.values

xvalues = np.array(x_cpns[1:]);
yvalues = np.array(y_cpns);
xx, yy = np.meshgrid(xvalues, yvalues, sparse=True)
plt.plot(xx,yy,marker='.',linestyle="none")
SWI_f = interp1d(xvalues,yvalues)
plt.plot(x, y, 'o', xx, yy, '-')
plt.show()

#read Kisters nc file 
Kisters=xr.open_mfdataset("T:/R466000/WAT/466674/8/4666748002/1_Ruwe_data/Kisters/Kisters.nc",concat_dim="time", combine="nested",
                  data_vars='minimal', coords='minimal', compat='override')

dates = np.unique(gdf["datum"].tolist())
date = dates[0]
tmpKisters = Kisters.Kisters.sel(
                  lon = xr.DataArray(gdf[gdf['datum'] == date].lon, dims = "points"),
                  lat = xr.DataArray(gdf[gdf['datum'] == date].lat, dims = "points"),
                  time = date, method = 'ffill')
some_variable = tmpKisters.values

#plot scatter
plt.scatter(SWI-f, Kisters)
plt.show()





# imgplot = plt.imshow(cpns)


