# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 16:01:27 2022

@author: Austin Abreu
"""
#This isn't designed to run all at once, and you probably shouldn't. Run in 
# Sections to prevent memory overload. Also, the helper functions automatically
# save plots.
#%% Libraries
from obspy import UTCDateTime
from obspy.core.event.catalog import read_events
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt
from matplotlib.transforms import offset_copy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
#%% DATA IMPORT
###############################################################################
#Provider Catalogs
ncedc = pd.read_csv("Data\\geysers_NCEDC_catalog2006_2016.csv",header = 0, sep = ",")
lbnl = pd.read_csv("Data\\geysers_LBL_catalog2006_2016.csv",header = 0, sep = ",")
ncedc["DateTime"] = pd.to_datetime(ncedc["DateTime"]) #The times are held as strings so we must convert them
lbnl["DateTime"] = pd.to_datetime(lbnl["DateTime"])

#Well Locations
wells = pd.read_csv("Data\\GEYSERS_WELLS.csv",header = 0, sep = "\t")
wells = wells[wells["Data\\Longitude"] < -122]

#Injector API ID, then parse out development wells.
injectors = pd.read_csv('inj_wells.csv',header = 0,sep='\t')
wells = wells[wells['APINumber'].isin(injectors['API'])] #Filtering Argument

#Station locations
stations = pd.read_csv('Data\\station_list.csv',header=0,sep = "|")

#EQT Picks & catalog
EQTcat = read_events('Data\\associations.xml')

#Misc
ncedcColor = 'xkcd:medium blue'
lbnlColor = 'xkcd:brick red'

keys = set(['DEB','ACR','AL6','EPR','GAXB','LCK','PFR','RGP','SQK','TCH'])
# %% PLOTTING FUNCTIONS
##############################################################################
def scatterplot (x,y,title,xlabel,ylabel,size,pointColor='k',label=None,
                 marker = 'o',colorMap = 'None', normal = 'None', mS = 7):
    fig=plt.figure(title,figsize=size)
    ax=plt.subplot()
    ax.scatter(x,y,c=pointColor,s=mS,label = label,marker = marker,
               cmap = colorMap, norm = normal)
    ax.set_title(name,fontsize=30)
    ax.set_ylabel(ylabel,fontsize=20)
    ax.set_xlabel(xlabel,fontsize=20)
    ax.tick_params('both',labelsize=15,width=3)
    title = title.replace('\n','')
    plt.savefig(title.replace(".",","))

def histplot (x,bins,name,xlabel,ylabel,size,label = 'None',fC = 'r'):
    fig=plt.figure(name,figsize=size)
    ax=plt.subplot()

    ax.hist(x,bins,facecolor=fC,edgecolor='k', label = label)
    ax.set_title(name,fontsize=20,color='k')
    ax.set_ylabel(ylabel,fontsize=15,color='k')
    ax.set_xlabel(xlabel,fontsize=15,color='k')
    ax.tick_params('both',labelsize=10,color='k',width=3)
    plt.savefig(name.replace("\n",""))
    
def b_val_prep(quantity,bins,correction = 0):
    #######################################################################
    # correction = the desired magnitude of completeness correction factor
    #######################################################################
    fig=plt.figure('zzz')
    ax=plt.subplot()    
    magHist = ax.hist(quantity,bins) #I need to create a miscellaneous plot for this method
    plt.close('zzz')

    Mc_correction=correction
    Mc_max_curv=magHist[1][np.argmax(magHist[0])] #Pulls maximum of histogram bins
    
    Mc=Mc_max_curv+Mc_correction
    print('Magnitude of completeness = '+str(Mc))

    N,M = np.histogram(quantity,bins)
    Nmax=np.array([len(quantity[quantity==max(quantity)])])
    N=np.append(N,Nmax)
    return N, M, Mc

def b_value(magnitudes,numberEvents,magCompleteness,title, lC):
    fig = plt.figure(title,figsize = (19,10))
    ax=plt.subplot()
    numberEvents=np.flip(np.cumsum(np.flip(numberEvents)))
    ax.scatter(magnitudes,numberEvents,c='k',s=5)
    ax.set_yscale('log')
    ax.set_xlabel('Magnitude of Events',fontsize=15)
    ax.set_ylabel('Number of Events, Log Scale',fontsize=15)
    ax.set_title(title,fontsize=20)
    #Creates a line along the magnitude of completeness
    ax.axvline(x=magCompleteness,c='m',
               label = 'Magnitude of Completeness=' + str(magCompleteness))
    
    Mcopy = magnitudes[0:-1]
    indexAboveComplete = [i for i,k in enumerate(Mcopy) if k >= magCompleteness] #returns INDICES on VALUE condition
    linearMags = [Mcopy[i] for i in indexAboveComplete]
    
    linearN = [numberEvents[i] for i in indexAboveComplete]
    
    for i,k in enumerate(linearN): #Sometimes we have a magnitude with 0 events that isn't culled
        if k == 0:
            del linearN[i]
            del linearMags[i]
    
    #Highlights events above the magnitude of completeness
    ax.scatter(linearMags,linearN,c='blue',s=5,label='Values Above Completeness')
    
    #We want to make a line model and extract the b-value from the dataset:
    b_val,a_val = np.polyfit(linearMags,np.log10(linearN),1)
    
    #Now we can solve a line with these coefficients. Here's two ways to do it:
        #the list "linearMags" isn't a np.array, and will return an error if it isn't transformed or indexed directly.
    modelN = 10**(b_val*np.array(linearMags) + a_val)
    
    #Finally, we can put it all together.
    ax.plot(linearMags,modelN,linewidth=1.5,color = lC,label=f'b = {str(np.round(b_val,4))}')
    
    #If we were curious how this deviates from a b-value of 1:
    model_1 = 10**(-1*np.array(linearMags) + a_val)
    ax.plot(linearMags,model_1,linewidth=1.5,color = 'green',linestyle='dotted',label = 'b = -1')
    
    
    #keep this at the bottom of the cell.
    ax.legend()
    plt.savefig(name.replace(".",","))

# I'll be using a few misc plotting functions repeatedly
def steamLine():
    plt.axhline(y = -0.7, label = 'Avg Depth of Steam Field', ls = ':', 
                color = 'b', lw = 3)

def felsiteLine():
    plt.axhline(y = -2, label = 'Avg Depth of Pluton', ls = ':', color = 'r', lw = 3)
    
def stationLabels(dataFrame):
    for i, text in enumerate(dataFrame.Station):
        plt.annotate(text, (dataFrame['Longitude'].iloc[i],
                            dataFrame['Latitude'].iloc[i]), fontsize = 'medium',
                            fontweight = 'bold')
        
#%% Catalog Disassembly Function
def reorganizer(cat): #requires an Obspy catalog object
    #Dictionaries allow for faster processing and easy mapping of station
        #names to data
    stationsCat={ #used to map station names to pick indices
        }
    channels={  #map station names to pick-relative channels
        }
    picks={ #map station names to pick time data
        }
    origins = {}
    p_time = {}
    s_time = {}
    for j,event in enumerate(cat):
        origins.setdefault(j,event.origins[0].time)
        for i,p in enumerate(event.picks):
            stationsCat.setdefault(p.waveform_id.station_code,[]).append((j,i))
            picks.setdefault(p.waveform_id.station_code,[]).append((j,i,p.time))
            channels.setdefault(p.waveform_id.station_code,[]).append((j,i,
                p.waveform_id.channel_code))
            if p.phase_hint == 'P':
                p_time.setdefault(p.waveform_id.station_code,[]).append((j,i,p.time))
            if p.phase_hint == 'S':
                s_time.setdefault(p.waveform_id.station_code,[]).append((j,i,p.time))
            
    key_dict = ['catTime','sta_Evnt','cha_pick','sta_times','P_times','S_times']
    dict_list = [origins,stationsCat,channels,picks,p_time,s_time]
    catDict = dict(zip(key_dict,dict_list))    
    
    return catDict
#%% MAG HIST
###############################################################################
name = 'Catalog Magnitude Distributions'
bins=np.unique(np.round(ncedc["Magnitude"],decimals = 1))
ylab = 'Number of Events'
xlab = 'Magnitude'
sz = (19,10)
histplot(ncedc["Magnitude"],bins,name,xlab,ylab,sz, fC = ncedcColor)
#%% B-Value Generation for NCEDC CAT
###############################################################################
bins=np.unique(np.round(ncedc["Magnitude"],decimals = 1))
N,M,Mc = b_val_prep(ncedc['Magnitude'],bins,correction = 0.2)

name = 'Cumulative Gutenburg-Richter Relationship Model'
b_value(M,N,Mc,name,lC = ncedcColor)

ncedcCat = ncedc[ncedc['Magnitude'] >= Mc] #Notice me! We don't analyze values below Mc.
#%% LBNL MAG HIST
###############################################################################
name = 'Catalog Magnitude Distributions'

bins=np.unique(np.round(lbnl["Magnitude"],decimals = 1))
ylab = 'Number of Events'
xlab = 'Magnitude'
sz = (19,10)

histplot(lbnl["Magnitude"],bins,name,xlab,ylab, size = sz,fC = lbnlColor)
#%% B-Value Generation for LBNL CAT
###############################################################################
bins=np.unique(np.round(lbnl["Magnitude"],decimals = 1))
N,M,Mc = b_val_prep(lbnl['Magnitude'],bins,correction = 0.2)

name = 'Cumulative Gutenburg-Richter Relationship Model'
b_value(M,N,Mc,name, lC = lbnlColor)

lbnlCat = lbnl[lbnl['Magnitude'] >= Mc] #Notice me! We don't analyze values below Mc.
#%% MAG OVER TIME
###############################################################################
name = 'Temporal Magnitude Distribution'
xlab = 'Time'
ylab = 'Magnitude'
sz = (19,10)
scatterplot(lbnlCat['DateTime'],lbnlCat['Magnitude'],name,xlab,ylab,sz,
            pointColor = lbnlColor,label = 'LBNL')
scatterplot(ncedcCat['DateTime'],ncedcCat['Magnitude'],name,xlab,ylab,sz,
            pointColor = ncedcColor,label = 'NCEDC')
plt.legend()
plt.savefig(name)
#%%
name = 'Mag - Long'
xlab = 'Latitude'
ylab = 'Magnitude'
sz = (19,10)
scatterplot(lbnlCat['Latitude'],lbnlCat['Magnitude'],name,xlab,ylab,sz,
            pointColor = lbnlColor,label = 'LBNL')
scatterplot(ncedcCat['Latitude'],ncedcCat['Magnitude'],name,xlab,ylab,sz,
            pointColor = ncedcColor,label = 'NCEDC')
plt.legend()
plt.savefig(name)
#%% LAT-LON
###############################################################################
name = 'Map-View Distribution of Events, LBNL'
xlab = 'Longitude'
ylab = 'Latitude'
sz = (10,10)

#normalize1 = mpl.colors.LogNorm(vmin = lbnlCat.Magnitude.min(), vmax = lbnlCat.Magnitude.max())
#normalize2 = mpl.colors.LogNorm(vmin = ncedcCat.Magnitude.min(), vmax = ncedcCat.Magnitude.max())

scatterplot(lbnlCat['Longitude'],lbnlCat['Latitude'],
            name,xlab,ylab,sz,pointColor= 'xkcd:light red',label = 'LBNL Events')
# scatterplot(ncedcCat['Longitude'],ncedcCat['Latitude'],
#             name,xlab,ylab,sz,pointColor= 'xkcd:light red', label = 'NCEDC Events')
scatterplot(wells['Longitude'],wells['Latitude'],name,xlab,ylab,sz,
            pointColor = 'b', label = 'Well', marker = '*')
scatterplot(stations['Longitude'],stations['Latitude'],name,xlab,ylab,sz,
        pointColor = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
        mS = 24)

stationLabels(stations)

plt.axis('equal')
plt.legend()

plt.savefig(name.replace(".",","))
#%% Cartopy utils
'''First plot, global position of area of study using Cartopy'''
fig = plt.figure(figsize=(19, 10))

ax=plt.subplot2grid((2,2),(0,0),projection=ccrs.Robinson(central_longitude=np.mean(ncedcCat['Longitude'])),colspan=2)
ax.set_title('Global position of area of study',fontsize=30)
ax.set_global()
ax.stock_img()
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle='dotted')
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.add_feature(cfeature.RIVERS)
ax.plot(np.mean(ncedcCat['Longitude']),np.mean(ncedcCat['Latitude']), 
        marker='*', color='r', markersize=15,alpha=0.8, transform=ccrs.Geodetic())

#%% Cartopy Map-View LBNL
# Create a Stamen terrain background instance.
stamen_terrain = cimgt.Stamen('terrain-background')
fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=ccrs.Miller(central_longitude=np.mean(ncedcCat['Longitude'])))
name = 'Regional Map of Study Area'
ax.set_title(name,fontsize=30)


# Limit the extent of the map to a small longitude/latitude range.
# ax2.set_extent([np.min(lons)-0.5,np.max(lons)+0.5,np.min(lats)-0.5,np.max(lats)+0.5], crs=ccrs.Geodetic())
ax.set_extent([np.min(ncedcCat['Longitude'])-0.03,
                np.max(ncedcCat['Longitude'])+0.03,
                np.min(ncedcCat['Latitude'])-0.03,
                np.max(ncedcCat['Latitude'])+0.03], crs=ccrs.Geodetic())


ax.gridlines(draw_labels=True,color='k')

# Add the Stamen data at zoom level 8.
ax.add_image(stamen_terrain, 10)

ax.set_aspect('equal')

#Add earthquakes with a scatter plot
ax.scatter(lbnlCat['Longitude'],lbnlCat['Latitude'],marker='o', color=lbnlColor,
            s=0.05,transform=ccrs.Geodetic(), label = 'LBNL Event')
ax.scatter(ncedcCat['Longitude'],ncedcCat['Latitude'],marker='o', color=ncedcColor,
            s=0.05,transform=ccrs.Geodetic(), label = 'NCEDC Event')
ax.scatter(wells['Longitude'],wells['Latitude'],s = 10,
            c = 'b', label = 'Well', marker = '*',transform=ccrs.Geodetic())
ax.scatter(stations['Longitude'],stations['Latitude'], 
        c = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
        s= 24,transform=ccrs.Geodetic())

# Use the cartopy interface to create a matplotlib transform object
# for the Geodetic coordinate system. We will use this along with
# matplotlib's offset_copy function to define a coordinate system which
# translates the text by 25 pixels to the left.
plt.legend()
geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
text_transform = offset_copy(geodetic_transform, units='dots', x=-50)

#%% Cartopy Map-View NCEDC
fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=ccrs.Miller(central_longitude=np.mean(lbnlCat['Longitude'])))
name = 'Regional Map of Study Area'
ax.set_title(name,fontsize=30)

ax.set_extent([np.min(lbnlCat['Longitude'])-0.03,
                np.max(lbnlCat['Longitude'])+0.03,
                np.min(lbnlCat['Latitude'])-0.03,
                np.max(lbnlCat['Latitude'])+0.03], crs=ccrs.Geodetic())


ax.gridlines(draw_labels=True,color='k')

ax.add_image(stamen_terrain, 10)

ax.set_aspect('equal')

ax.scatter(lbnlCat['Longitude'],lbnlCat['Latitude'],marker='o', color=lbnlColor,
            s=0.1,transform=ccrs.Geodetic(), label = 'LBNL Event')
ax.scatter(wells['Longitude'],wells['Latitude'],s = 10,
            c = 'b', label = 'Well', marker = '*',transform=ccrs.Geodetic())
ax.scatter(stations['Longitude'],stations['Latitude'], 
        c = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
        s= 24,transform=ccrs.Geodetic())

plt.legend()
geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
text_transform = offset_copy(geodetic_transform, units='dots', x=-50)

#%% DEPTH-LON XSXN
###############################################################################
name = 'Longitudinal Cross-Section, LBNL'
xlab = 'Longitude'
ylab = 'Depth (km)'
sz = (19,10)
scatterplot(lbnlCat['Longitude'][lbnlCat['Depth'] < 10],
            -lbnlCat['Depth'][lbnlCat['Depth'] < 10],name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

name = 'Longitudinal Cross-Section, NCEDC'
scatterplot(ncedcCat['Longitude'][ncedcCat['Depth'] < 10],
            -ncedcCat['Depth'][ncedcCat['Depth'] < 10],name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

#%% DEPTH-LAT XSXN
###############################################################################
name = 'Latitudinal Cross-Section, LBNL'
xlab = 'Latitude'
ylab = 'Depth (km)'
sz = (19,10)
scatterplot(lbnlCat['Latitude'][lbnlCat['Depth'] < 10],
            -lbnlCat['Depth'][lbnlCat['Depth'] < 10],name,xlab,ylab,sz)
scatterplot(wells['Latitude'],np.zeros(len(wells)),name,xlab,ylab,sz, pointColor = 'b')
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

name = 'Latitudinal Cross-Section, NCEDC'
scatterplot(ncedcCat['Latitude'][ncedcCat['Depth'] < 10],
            -(ncedcCat['Depth'][ncedcCat['Depth'] < 10] + 0.779),name,xlab,ylab,sz)
scatterplot(wells['Latitude'],np.zeros(len(wells)),name,xlab,ylab,sz, pointColor = 'b')
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))
#%% DEPTH / TIME
###############################################################################
name = 'Temporal Distribution of Depths, LBNL'
xlab = 'Time'
ylab = 'Depth (km)'
sz = (19,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime
windowLBNL = lbnlCat[(lbnlCat['DateTime'] > sTime) & (lbnlCat['DateTime'] < eTime)]
windowNCEDC = ncedcCat[(ncedcCat['DateTime'] > sTime) & (ncedcCat['DateTime'] < eTime)]
scatterplot(windowLBNL['DateTime'][windowLBNL['Depth'] < 10],
            -windowLBNL['Depth'][windowLBNL['Depth'] < 10],name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

name = 'Temporal Distribution of Depths, NCEDC'
scatterplot(windowNCEDC['DateTime'][windowNCEDC['Depth'] < 10],
            -(windowNCEDC['Depth'][windowNCEDC['Depth'] < 10] + 0.779),name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()

plt.savefig(name.replace('.',','))
#%% DEPTH / MAG
###############################################################################
xlab = 'Magnitude, M$_L$'
ylab = 'Depth (km)'
sz = (19,10)

# LBNL
name = 'Distribution of Magnitudes at Depth, LBNL'

scatterplot(lbnlCat['Magnitude'][lbnlCat['Depth'] < 10],
            -lbnlCat['Depth'][lbnlCat['Depth'] < 10],name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

# NCEDC
name = 'Distribution of Magnitudes at Depth, NCEDC'
xlab = 'Magnitude, M$_d$'
scatterplot(ncedcCat['Magnitude'][ncedcCat['Depth'] < 10],
            -(ncedcCat['Depth'][ncedcCat['Depth'] < 10] + 0.779),name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

#%% DEPTH HIST
###############################################################################
name = 'Event Depths, Whole Field'
fig=plt.figure(name,figsize=(19,10))
ax=plt.subplot()

ax.hist(-lbnlCat['Depth'][lbnlCat['Depth'] < 10],bins=50,facecolor=lbnlColor,
        edgecolor='k',label='LBNL Distribution',
        orientation='horizontal')
ax.hist(-ncedcCat['Depth'][ncedcCat['Depth'] < 10],bins=50,facecolor=ncedcColor,
        edgecolor='k',label='NCEDC Distribution',
        orientation='horizontal')
ax.set_title(name,fontsize=40,color='k')
ax.set_xlabel('Number of events',fontsize=30,color='k')
ax.set_ylabel('Depth (km)',fontsize=30,color='k')
ax.tick_params('both',labelsize=30,color='black',width=3)

steamLine()
felsiteLine()
ax.legend(fontsize=15)
plt.savefig(name)
#%% LONGITUDE OVER TIME
###############################################################################
name = 'Time Change in Longitude, LBNL'
xlab = 'Time'
ylab = 'Longitude'
sz = (19,10)
sTime = UTCDateTime(2006,1,1).datetime
eTime = UTCDateTime(2016,1,1).datetime
windowLBNL = lbnlCat[(lbnlCat['DateTime'] > sTime) & (lbnlCat['DateTime'] < eTime)]
windowNCEDC = ncedcCat[(ncedcCat['DateTime'] > sTime) & (ncedcCat['DateTime'] < eTime)]
scatterplot(windowLBNL['DateTime'],windowLBNL['Longitude'],name,xlab,ylab,sz)
name = 'Time Change in Longitude, NCEDC'
scatterplot(windowNCEDC['DateTime'],windowNCEDC['Longitude'],name,xlab,ylab,sz)
#%% LATITUDE OVER TIME
###############################################################################
name = 'Time Change in Latitude, LBNL'
xlab = 'Time'
ylab = 'Latitude'
sz = (19,10)
sTime = UTCDateTime(2006,1,1).datetime
eTime = UTCDateTime(2016,1,1).datetime
windowLBNL = lbnlCat[(lbnlCat['DateTime'] > sTime) & (lbnlCat['DateTime'] < eTime)]
windowNCEDC = ncedcCat[(ncedcCat['DateTime'] > sTime) & (ncedcCat['DateTime'] < eTime)]
scatterplot(windowLBNL['DateTime'],windowLBNL['Latitude'],name,xlab,ylab,sz)
name = 'Time Change in Latitude, NCEDC'
scatterplot(windowNCEDC['DateTime'],windowNCEDC['Latitude'],name,xlab,ylab,sz)
#%% Time Targeting
nnn = ncedcCat.resample('Q',on='DateTime').count()
lll = lbnlCat.resample('Q',on='DateTime').count()
#%% Quarterly Event Frequency Comparison Plot
name = 'Comparison of Events Recorded in Provider Catalogs'
xlab = 'Time, Weeks'
ylab = 'Number of Recorded Events'
fig=plt.figure(name,figsize=(19,10))
ax=plt.subplot()
ax.plot(list(nnn.index),nnn.Magnitude, color = ncedcColor,
        label = 'NCEDC')
ax.plot(list(lll.index),lll.Magnitude, color = lbnlColor,
        label = 'LBNL')
ax.set_title(name,fontsize=30)
ax.set_ylabel(ylab,fontsize=20)
ax.set_xlabel(xlab,fontsize=20)
ax.tick_params('both',labelsize=15,width=3)
plt.grid('on')
plt.legend()
plt.savefig(name)
#%% Geographic Targeting
###############################################################################
northBounds = [38.79, 38.90, -122.90, -122.75]
southBounds = [38.72, 38.80, -122.80, -122.66]
lbnl_NZ = lbnlCat[(lbnlCat['Latitude'] > northBounds[0]) &
                    (lbnlCat['Latitude'] < northBounds[1]) &
                    (lbnlCat['Longitude'] > northBounds[2]) &
                    (lbnlCat['Longitude'] < northBounds[3])]
lbnl_SZ = lbnlCat[(lbnlCat['Latitude'] > southBounds[0]) &
                (lbnlCat['Latitude'] < southBounds[1]) &
                (lbnlCat['Longitude'] > southBounds[2]) &
                (lbnlCat['Longitude'] < southBounds[3])]
ncedc_NZ = ncedcCat[(ncedcCat['Latitude'] > northBounds[0]) &
                    (ncedcCat['Latitude'] < northBounds[1]) &
                    (ncedcCat['Longitude'] > northBounds[2]) &
                    (ncedcCat['Longitude'] < northBounds[3])]
ncedc_SZ = ncedcCat[(ncedcCat['Latitude'] > southBounds[0]) &
                (ncedcCat['Latitude'] < southBounds[1]) &
                (ncedcCat['Longitude'] > southBounds[2]) &
                (ncedcCat['Longitude'] < southBounds[3])]
wells_NZ = wells[(wells['Latitude'] > northBounds[0]) &
                (wells['Latitude'] < northBounds[1]) &
                (wells['Longitude'] > northBounds[2]) &
                (wells['Longitude'] < northBounds[3])]
wells_SZ = wells[(wells['Latitude'] > southBounds[0]) &
                (wells['Latitude'] < southBounds[1]) &
                (wells['Longitude'] > southBounds[2]) &
                (wells['Longitude'] < southBounds[3])]
stations_NZ = stations[(stations['Latitude'] > northBounds[0]) &
                (stations['Latitude'] < northBounds[1]) &
                (stations['Longitude'] > northBounds[2]) &
                (stations['Longitude'] < northBounds[3])]
stations_SZ = stations[(stations['Latitude'] > southBounds[0]) &
                (stations['Latitude'] < southBounds[1]) &
                (stations['Longitude'] > southBounds[2]) &
                (stations['Longitude'] < southBounds[3])]

#%% NORTH ZONE
###############################################################################
###############################################################################
## MAP
stamen_terrain = cimgt.Stamen('terrain-background')
fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=ccrs.Miller(central_longitude=np.mean(ncedc_NZ['Longitude'])))
name = 'North-West Region Overview'
ax.set_title(name,fontsize=30)


# Limit the extent of the map to a small longitude/latitude range.
# ax2.set_extent([np.min(lons)-0.5,np.max(lons)+0.5,np.min(lats)-0.5,np.max(lats)+0.5], crs=ccrs.Geodetic())
ax.set_extent([np.min(ncedc_NZ['Longitude']),
                np.max(ncedc_NZ['Longitude']),
                np.min(ncedc_NZ['Latitude']),
                np.max(ncedc_NZ['Latitude'])], crs=ccrs.Geodetic())


ax.gridlines(draw_labels=True,color='k')

# Add the Stamen data at zoom level 8.
ax.add_image(stamen_terrain, 12)

ax.set_aspect('equal')

#Add earthquakes with a scatter plot
ax.scatter(lbnl_NZ['Longitude'],lbnl_NZ['Latitude'],marker='o', color=lbnlColor,
            s=0.05,transform=ccrs.Geodetic(), label = 'LBNL Event')
ax.scatter(ncedc_NZ['Longitude'],ncedc_NZ['Latitude'],marker='o', color=ncedcColor,
            s=0.07,transform=ccrs.Geodetic(), label = 'NCEDC Event')
ax.scatter(wells['Longitude'],wells['Latitude'],s = 10,
            c = 'b', label = 'Well', marker = '*',transform=ccrs.Geodetic())
ax.scatter(stations['Longitude'],stations['Latitude'], 
        c = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
        s= 24,transform=ccrs.Geodetic())

# Use the cartopy interface to create a matplotlib transform object
# for the Geodetic coordinate system. We will use this along with
# matplotlib's offset_copy function to define a coordinate system which
# translates the text by 25 pixels to the left.
plt.legend()
geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
text_transform = offset_copy(geodetic_transform, units='dots', x=-50)

## MAG HIST
###############################################################################
name = 'NW Region Magnitude Distribution'
bins1 = np.unique(np.round(ncedc_NZ["Magnitude"],decimals = 1))
bins2 = np.unique(np.round(lbnl_NZ["Magnitude"],decimals = 1))
ylab = 'Number of Events'
xlab = 'Magnitude'
sz = (19,10)
histplot(lbnl_NZ["Magnitude"],bins2,name,xlab,ylab,sz,label = 'LBNL',
         fC = lbnlColor)
histplot(ncedc_NZ["Magnitude"],bins1,name,xlab,ylab,sz, 
         label = 'NCEDC', fC = ncedcColor)
plt.legend()
plt.savefig(name)

## B-VALUES
###############################################################################
### NCEDC ###
name = 'Cumulative Gutenburg-Richter Model, NW Region, NCEDC'
N,M,Mc = b_val_prep(ncedc_NZ['Magnitude'],bins1,correction = 0)
b_value(M,N,Mc,name)

ncedc_NZ = ncedc_NZ[ncedc_NZ['Magnitude'] > Mc]

### LBNL ###
name = 'Cumulative Gutenburg-Richter Model, NW Region, LBNL'
N,M,Mc = b_val_prep(lbnl_NZ['Magnitude'],bins2,correction = 0)
b_value(M,N,Mc,name)


## MAG / TIME ##
###############################################################################
name = 'Temporal Magnitude Distribution in NW Region'
xlab = 'Time'
ylab = 'Magnitude'
sz = (19,10)
scatterplot(lbnl_NZ["DateTime"],lbnl_NZ["Depth"],name,xlab,ylab,
            sz, pointColor = 'xkcd:brick red', label = 'LBNL')
scatterplot(ncedc_NZ["DateTime"],ncedc_NZ["Depth"],name,xlab,
            ylab,sz, pointColor = 'xkcd:medium blue',label = 'NCEDC')
plt.legend()
plt.savefig(name.replace('.',','))

## DEPTH / TIME ##
###############################################################################
name = 'Temporal Distribution of Depths in NW Region'
xlab = 'Time'
ylab = 'Depth (km)'
sz = (19,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime

# I need to connect injection rates to increases in seismicity...
windowLBNL = lbnl_NZ[(lbnl_NZ['DateTime'] > sTime) & (lbnl_NZ['DateTime'] < eTime)]
windowNCEDC = ncedc_NZ[(ncedc_NZ['DateTime'] > sTime) & (ncedc_NZ['DateTime'] < eTime)]

### LBNL ###

scatterplot(windowLBNL['DateTime'][windowLBNL['Depth'] < 10],
            -windowLBNL['Depth'][windowLBNL['Depth'] < 10],name,xlab,ylab,sz,
            pointColor = lbnlColor, label = 'LBNL')

### NCEDC ###
scatterplot(windowNCEDC['DateTime'][windowNCEDC['Depth'] < 10],
            -(windowNCEDC['Depth'][windowNCEDC['Depth'] < 10] + 0.779),name,
            xlab,ylab,sz, pointColor = ncedcColor,label = 'NCEDC')
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

## DEPTH / MAG ##
###############################################################################

ylab = 'Depth (km)'
sz = (19,10)

### LBNL ###
name = 'Distribution of Magnitudes at Depth, LBNL'
xlab = 'Magnitude, M$_L$'
scatterplot(lbnl_NZ['Magnitude'][lbnl_NZ['Depth'] < 10],
            -lbnl_NZ['Depth'][lbnl_NZ['Depth'] < 10],name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

### NCEDC ###
name = 'Distribution of Magnitudes at Depth, NCEDC'
xlab = 'Magnitude, M$_d$'
scatterplot(ncedc_NZ['Magnitude'][ncedc_NZ['Depth'] < 10],
            -(ncedc_NZ['Depth'][ncedc_NZ['Depth'] < 10] + 0.779),name,xlab,
            ylab,sz)

steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

## DEPTH HIST ##
###############################################################################
name = 'Event Depths in NW Region'
fig=plt.figure(name,figsize=(11,11))
ax=plt.subplot()

ax.hist(-lbnl_NZ['Depth'][lbnl_NZ['Depth'] < 10],bins=50,facecolor=lbnlColor,
        edgecolor='k',label='LBNL Distribution',
        orientation='horizontal')
ax.hist(-(ncedc_NZ['Depth'][ncedc_NZ['Depth'] < 10]+ 1),bins=50,facecolor=ncedcColor,
        edgecolor='k',label='NCEDC Distribution',
        orientation='horizontal')
ax.set_title(name,fontsize=40,color='k')
ax.set_xlabel('Number of events',fontsize=20,color='k')
ax.set_ylabel('Depth (km)',fontsize=20,color='k')
ax.tick_params('both',labelsize=30,color='black',width=3)

steamLine()
felsiteLine()
ax.legend(fontsize=15)
plt.savefig(name)

## LONG/TIME ##
###############################################################################
name = 'Event Longitude Through Time in NW Region'
xlab = 'Time'
ylab = 'Longitude'
sz = (19,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime

windowLBNL = lbnl_NZ[(lbnl_NZ['DateTime'] > sTime) & (lbnl_NZ['DateTime'] < eTime)]
windowNCEDC = ncedc_NZ[(ncedc_NZ['DateTime'] > sTime) & (ncedc_NZ['DateTime'] < eTime)]

scatterplot(windowLBNL['DateTime'], windowLBNL['Longitude'], name, xlab, ylab, sz,
            pointColor = lbnlColor, label = 'LBNL')
scatterplot(windowNCEDC['DateTime'],windowNCEDC['Longitude'],name,xlab,ylab,sz,
            pointColor = ncedcColor, label = 'NCEDC')
plt.legend()
plt.savefig(name)

## LAT/TIME ##
###############################################################################
name = 'Event Latitude Through Time in NW Region'
xlab = 'Time'
ylab = 'Latitude'
sz = (10,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime

windowLBNL = lbnl_NZ[(lbnl_NZ['DateTime'] > sTime) & (lbnl_NZ['DateTime'] < eTime)]
windowNCEDC = ncedc_NZ[(ncedc_NZ['DateTime'] > sTime) & (ncedc_NZ['DateTime'] < eTime)]

scatterplot(windowLBNL['DateTime'],windowLBNL['Latitude'],name,xlab,ylab,sz,
            pointColor = lbnlColor,label = 'LBNL')
scatterplot(windowNCEDC['DateTime'],windowNCEDC['Latitude'],name,xlab,ylab,sz,
            pointColor = ncedcColor, label = 'NCEDC')
plt.legend()
plt.savefig(name)

## Temporally Constrained Map View
###############################################################################
name = 'Correlation Between Events in NW Region,\n 22-02-2015 to 27-02-2015'
xlab = 'Longitude'
ylab = 'Latitude'
sz = (10,10)
 
scatterplot(windowLBNL["Longitude"],windowLBNL["Latitude"],name,xlab,
            ylab,sz,pointColor = lbnlColor, label = 'LBNL')
scatterplot(windowNCEDC["Longitude"],windowNCEDC["Latitude"],name,xlab,
            ylab,sz,pointColor = ncedcColor, label = 'NCEDC')
scatterplot(wells_NZ['Longitude'],wells_NZ['Latitude'],name,xlab,ylab,sz,
            pointColor = 'xkcd:vibrant green', label = 'Injector',
            marker = '1', mS = 20)
# scatterplot(stations_NZ['Longitude'],stations_NZ['Latitude'],name,xlab,ylab,sz,
#         pointColor = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
#         mS = 20)

#stationLabels(stations_NZ)
plt.legend()
name = name.replace('\n','')
plt.savefig(name)
#%% SOUTH ZONE
###############################################################################
###############################################################################
### MAP VIEW
###############################################################################
stamen_terrain = cimgt.Stamen('terrain-background')
fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=ccrs.Miller(central_longitude=np.mean(ncedc_SZ['Longitude'])))
name = 'South-East Region Overview'
ax.set_title(name,fontsize=30)


# Limit the extent of the map to a small longitude/latitude range.
# ax2.set_extent([np.min(lons)-0.5,np.max(lons)+0.5,np.min(lats)-0.5,np.max(lats)+0.5], crs=ccrs.Geodetic())
ax.set_extent([np.min(ncedc_SZ['Longitude']),
                np.max(ncedc_SZ['Longitude']),
                np.min(ncedc_SZ['Latitude']),
                np.max(ncedc_SZ['Latitude'])], crs=ccrs.Geodetic())


ax.gridlines(draw_labels=True,color='k')

# Add the Stamen data at zoom level 8.
ax.add_image(stamen_terrain, 12)

ax.set_aspect('equal')

#Add earthquakes with a scatter plot
ax.scatter(lbnl_SZ['Longitude'],lbnl_SZ['Latitude'],marker='o', color=lbnlColor,
            s=0.05,transform=ccrs.Geodetic(), label = 'LBNL Event')
ax.scatter(ncedc_SZ['Longitude'],ncedc_SZ['Latitude'],marker='o', color=ncedcColor,
            s=0.1,transform=ccrs.Geodetic(), label = 'NCEDC Event')
ax.scatter(wells['Longitude'],wells['Latitude'],s = 10,
            c = 'b', label = 'Well', marker = '*',transform=ccrs.Geodetic())
ax.scatter(stations['Longitude'],stations['Latitude'], 
        c = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
        s= 24,transform=ccrs.Geodetic())

# Use the cartopy interface to create a matplotlib transform object
# for the Geodetic coordinate system. We will use this along with
# matplotlib's offset_copy function to define a coordinate system which
# translates the text by 25 pixels to the left.
plt.legend()
geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
text_transform = offset_copy(geodetic_transform, units='dots', x=-50)

## MAG HIST
###############################################################################
name = 'South-East Zone Magnitude Distribution'
bins1 = np.unique(np.round(ncedc_SZ["Magnitude"],decimals = 1))
bins2 = np.unique(np.round(lbnl_SZ["Magnitude"],decimals = 1))
ylab = 'Number of Events'
xlab = 'Magnitude'
sz = (19,10)
histplot(lbnl_SZ["Magnitude"],bins2,name,xlab,ylab,sz,label = 'LBNL',
         fC = lbnlColor)
histplot(ncedc_SZ["Magnitude"],bins1,name,xlab,ylab,sz, 
         label = 'NCEDC', fC = ncedcColor)
plt.legend()
plt.savefig(name)

## B-VALUES
###############################################################################
### NCEDC ###
name = 'Cumulative Gutenburg-Richter Model, SE Region, NCEDC'
N,M,Mc = b_val_prep(ncedc_SZ['Magnitude'],bins1,correction = 0)
b_value(M,N,Mc,name)

ncedc_SZ = ncedc_SZ[ncedc_SZ['Magnitude'] > Mc]

### LBNL ###
name = 'Cumulative Gutenburg-Richter Model, SE Region, LBNL'
N,M,Mc = b_val_prep(lbnl_SZ['Magnitude'],bins2,correction = 0)
b_value(M,N,Mc,name)


## MAG / TIME ##
###############################################################################
name = 'Temporal Magnitude Distribution in SE Zone'
xlab = 'Time'
ylab = 'Magnitude'
sz = (19,10)
scatterplot(lbnl_SZ["DateTime"],lbnl_SZ["Depth"],name,xlab,ylab,
            sz, pointColor = 'xkcd:brick red', label = 'LBNL')
scatterplot(ncedc_SZ["DateTime"],ncedc_SZ["Depth"],name,xlab,
            ylab,sz, pointColor = 'xkcd:medium blue',label = 'NCEDC')
plt.legend()
plt.savefig(name.replace('.',','))

## DEPTH / TIME ##
###############################################################################
name = 'Temporal Distribution of Depths in SE Region'
xlab = 'Time'
ylab = 'Depth (km)'
sz = (19,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime

# I need to connect injection rates to increases in seismicity...
windowLBNL = lbnl_SZ[(lbnl_SZ['DateTime'] > sTime) & (lbnl_SZ['DateTime'] < eTime)]
windowNCEDC = ncedc_SZ[(ncedc_SZ['DateTime'] > sTime) & (ncedc_SZ['DateTime'] < eTime)]

### LBNL ###

scatterplot(windowLBNL['DateTime'][windowLBNL['Depth'] < 10],
            -windowLBNL['Depth'][windowLBNL['Depth'] < 10],name,xlab,ylab,sz,
            pointColor = lbnlColor, label = 'LBNL')

### NCEDC ###
scatterplot(windowNCEDC['DateTime'][windowNCEDC['Depth'] < 10],
            -(windowNCEDC['Depth'][windowNCEDC['Depth'] < 10] + 0.779),name,
            xlab,ylab,sz, pointColor = ncedcColor,label = 'NCEDC')
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

## DEPTH / MAG ##
###############################################################################

ylab = 'Depth (km)'
sz = (19,10)

### LBNL ###
name = 'Distribution of Magnitudes at Depth in SE Region, LBNL'
xlab = 'Magnitude, M$_L$'
scatterplot(lbnl_SZ['Magnitude'][lbnl_SZ['Depth'] < 10],
            -lbnl_SZ['Depth'][lbnl_SZ['Depth'] < 10],name,xlab,ylab,sz)
steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

### NCEDC ###
name = 'Distribution of Magnitudes at Depth in SE Region, NCEDC'
xlab = 'Magnitude, M$_d$'
scatterplot(ncedc_SZ['Magnitude'][ncedc_SZ['Depth'] < 10],
            -(ncedc_SZ['Depth'][ncedc_SZ['Depth'] < 10] + 0.779),name,xlab,
            ylab,sz)

steamLine()
felsiteLine()
plt.legend()
plt.savefig(name.replace('.',','))

## DEPTH HIST ##
###############################################################################
name = 'Event Depths in SE Region'
fig=plt.figure(name,figsize=(19,10))
ax=plt.subplot()

ax.hist(-lbnl_SZ['Depth'][lbnl_SZ['Depth'] < 10],bins=50,facecolor=lbnlColor,
        edgecolor='k',label='LBNL Distribution',
        orientation='horizontal')
ax.hist(-(ncedc_SZ['Depth'][ncedc_SZ['Depth'] < 10]+ 0.779),bins=50,facecolor=ncedcColor,
        edgecolor='k',label='NCEDC Distribution',
        orientation='horizontal')
ax.set_title(name,fontsize=40,color='k')
ax.set_xlabel('Number of events',fontsize=30,color='k')
ax.set_ylabel('Depth (km)',fontsize=30,color='k')
ax.tick_params('both',labelsize=30,color='black',width=3)

steamLine()
felsiteLine()
ax.legend(fontsize=15)
plt.savefig(name)

## LONG/TIME ##
###############################################################################
name = 'Event Longitude Through Time in SE Region'
xlab = 'Time'
ylab = 'Longitude'
sz = (19,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime

windowLBNL = lbnl_SZ[(lbnl_SZ['DateTime'] > sTime) & (lbnl_SZ['DateTime'] < eTime)]
windowNCEDC = ncedc_SZ[(ncedc_SZ['DateTime'] > sTime) & (ncedc_SZ['DateTime'] < eTime)]

scatterplot(windowLBNL['DateTime'],windowLBNL['Longitude'],name,xlab,ylab,sz,
            pointColor = lbnlColor,label = 'LBNL')
scatterplot(windowNCEDC['DateTime'],windowNCEDC['Longitude'],name,xlab,ylab,sz,
            pointColor = ncedcColor, label = 'NCEDC')
plt.legend()
plt.savefig(name)

## LAT/TIME ##
###############################################################################
name = 'Event Latitude Through Time in SE Region'
xlab = 'Time'
ylab = 'Latitude'
sz = (19,10)
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime

windowLBNL = lbnl_SZ[(lbnl_SZ['DateTime'] > sTime) & (lbnl_SZ['DateTime'] < eTime)]
windowNCEDC = ncedc_SZ[(ncedc_SZ['DateTime'] > sTime) & (ncedc_SZ['DateTime'] < eTime)]

scatterplot(windowLBNL['DateTime'],windowLBNL['Latitude'],name,xlab,ylab,sz,
            pointColor = lbnlColor,label = 'LBNL')
scatterplot(windowNCEDC['DateTime'],windowNCEDC['Latitude'],name,xlab,ylab,sz,
            pointColor = ncedcColor, label = 'NCEDC')
plt.legend()
plt.savefig(name)

## Temporally Constrained Map View
###############################################################################
name = 'Correlation Between Events in SE Region,\n 22-02-2015 to 27-02-2015'
xlab = 'Longitude'
ylab = 'Latitude'
sz = (10,10)
 
scatterplot(windowLBNL["Longitude"],windowLBNL["Latitude"],name,xlab,
            ylab,sz,pointColor = lbnlColor, label = 'LBNL')
scatterplot(windowNCEDC["Longitude"],windowNCEDC["Latitude"],name,xlab,
            ylab,sz,pointColor = ncedcColor, label = 'NCEDC')
scatterplot(wells_SZ['Longitude'],wells_SZ['Latitude'],name,xlab,ylab,sz,
            pointColor = 'xkcd:vibrant green', label = 'Injector',
            marker = '1', mS = 20)
# scatterplot(stations_NZ['Longitude'],stations_NZ['Latitude'],name,xlab,ylab,sz,
#         pointColor = 'xkcd:very dark blue', label = 'Seismogram', marker = '^',
#         mS = 20)

#stationLabels(stations_NZ)
plt.legend()
name = name.replace('\n','')
plt.savefig(name)



#%% EQT Histogram
EQTDict = reorganizer(EQTcat)

#Resamples origin times so I can accumulate them and compare
originsDT = {key:value.datetime for key, value in EQTDict['catTime'].items()}
origins = pd.Series(originsDT.keys(),index = originsDT.values())
eventCounts = origins.resample('D').count()

    
sTime = UTCDateTime(2015,2,22).datetime
eTime = UTCDateTime(2015,2,27).datetime
# This is pretty redundant but it's not hurting anyone
lbnlPicks = lbnlCat['DateTime'][(lbnlCat['DateTime'] > sTime) & (lbnlCat['DateTime'] < eTime)]
ncedcPicks = ncedcCat['DateTime'][(ncedcCat['DateTime'] > sTime) & (ncedcCat['DateTime'] < eTime)]

windowLBNL = lbnlCat[(lbnlCat['DateTime'] > sTime) & (lbnlCat['DateTime'] < eTime)]
windowNCEDC = ncedcCat[(ncedcCat['DateTime'] > sTime) & (ncedcCat['DateTime'] < eTime)]

lbnlCounts = windowLBNL.resample('D', on = 'DateTime').count()
ncedcCounts = windowNCEDC.resample('D', on = 'DateTime').count()

# Self explanatory
name = 'Comparison of Events Recorded in Catalogs'
xlab = 'Time, Days'
ylabl = 'Number of Recorded Events'
fig=plt.figure(name,figsize=(10,10))
ax=plt.subplot()
ax.plot(list(eventCounts.index),eventCounts, color = 'g', label = 'EQT')
ax.plot(list(ncedcCounts.index),ncedcCounts.Magnitude, color = ncedcColor,
        label = 'NCEDC')
ax.plot(list(lbnlCounts.index),lbnlCounts.Magnitude, color = lbnlColor,
        label = 'LBNL')
ax.set_title(name,fontsize=30)
ax.set_ylabel(ylab,fontsize=20)
ax.set_xlabel(xlab,fontsize=20)
ax.tick_params('both',labelsize=15,width=3)
plt.legend()
plt.savefig(name)

#%% Event Correlation Algorithm
Corr = {'lbnl':[],'ncedc':[]}

#Pulls EQT and lbnl times and computes a time difference between them.
    # If the difference falls within a threshold, it is appended to a dictionary
    # Containing correlated times per station.
    
diffsL = []
diffsN = []
for i,EQTtime in EQTDict['catTime'].items():
    for j,lbnltime in enumerate(lbnlPicks):
        dateDiff = EQTtime.datetime - lbnltime
        #some weird thing happens where it loops around the start of the day,
            # hence the != -1 at the end.
        if dateDiff <= dt.timedelta(0,3) and dateDiff >= -dt.timedelta(0,3) and dateDiff._d == 0:
            Corr['lbnl'].append((i,j))
            diffsL += [dateDiff.total_seconds()]
    for k,ncedctime in enumerate(ncedcPicks):
        dateDiff = EQTtime.datetime - ncedctime
        if dateDiff <= dt.timedelta(0,3) and dateDiff >= -dt.timedelta(0,3) and dateDiff._d == 0:
            Corr['ncedc'].append((i,k))
            diffsN += [dateDiff.total_seconds()]
#A function of similar form can be used to correlate all three catalogs.
del dateDiff

#%% Quantify Differences in Event Matching

title = 'Absolute Differences of Correlated Events \n EQT vs NCEDC vs LBNL'
xlab = 'Difference in Time'
ylab = 'Number of Occurrences'
bins = 20
histplot(diffsL,bins,title,xlab,ylab,size = (10,10),label = 'LBNL', fC = lbnlColor)
histplot(diffsN,bins,title,xlab,ylab,size = (10,10),label = 'NCEDC', fC = ncedcColor)
plt.legend()
plt.savefig(name.replace('\n',''))
#%% Compute triangulations
# Utilizes Time difference beteen P and S wave arrivals to calculate the 
#distance between the source and the station. Aggregates distances into a 
#4-size tuple with the following indexing: 0 = Source Catalog Event index
#1 = station longitude | 2 = station latitude | 3 = calculated distance in km
##############################################################################
s_p = []
dist = []
triangulations = {}
for key in keys:
    for P_group in EQTDict['P_times'][key]:
        for S_group in EQTDict['S_times'][key]:
            if P_group[0] == S_group[0]:
                s_p = S_group[2] - P_group[2]
                dist = s_p/((1/4.2)-(1/8))
                for pair in EQTDict['sta_Evnt'][key]:
                        if (P_group[0],P_group[1]) == pair:
                            triangulations.setdefault(key,[]).append((P_group[0],
                                  list(stations['Longitude'][stations['Station'] == str(key)])[0],
                                  list(stations['Latitude'][stations['Station'] == str(key)])[0], dist))
#%% Cross-Correlate Events to Build Location Database
# Uses the correlation data from Corr dict to identify events by matching in-
# ices from trangulation tuple[0] to the appropriate index in Corr. The function
# organizes two lists and a dictionary for triangulation data that replaces the
# catalog index in triangulation with the catalog index for the corresponding event

# I'd like to use pandas to do this, but I wasn't able to figure out how.

lbnlLoc = [] #initialize variables
ncedcLoc = []
triMatch = {}
for key in keys: #match key set of stations
    for catKey, relation in Corr.items(): 
        for item in relation: #unfortunately items returns the bulk list
            for group in triangulations[key]:
                if group[0] == item[0]: #here the catalog index for correlated events is matched
                    if catKey == 'lbnl':
                        lbnlLoc.append((lbnlCat['Longitude'].iloc[item[1]],lbnlCat['Latitude'].iloc[item[1]]))
                        triMatch.setdefault(catKey,[]).append((item[1],group[1:4]))
                    if catKey == 'ncedc':
                        ncedcLoc.append((ncedcCat['Longitude'].iloc[item[1]],ncedcCat['Latitude'].iloc[item[1]]))
                        triMatch.setdefault(catKey,[]).append((item[1],group[1:4]))


#%% Plotting Triangulations
# This function, and the one below it, plots a schematic diagram in plan view
# with which to judge the accuracy of an triangulation effort. The catalog event
# is highlighted and label. Be careful with setting i to a large quantity of
# indices, as it may crash python.

#LBNL Triangulation Match
for j,group in enumerate(triMatch['lbnl']):
    for i in range(5):
        if group[0] == i:
            name = 'Triangulation Method for LBNL Catalog Event ' + str(i)
            fig=plt.figure(name,figsize=(10,10))
            ax=plt.subplot()
            ax.scatter(group[1][0], group[1][1], color = 'xkcd:very dark blue')
            ax.scatter(lbnlLoc[i][0],lbnlLoc[i][1], color = 'g')
            plt.annotate('Catalog Event',(lbnlLoc[i][0],lbnlLoc[i][1]))
            a = plt.Circle((group[1][0], group[1][1]),group[1][2]/111, color='r',fill=False)
            ax.set_aspect('equal', adjustable='datalim')
            ax.add_patch(a)
            ax.set_ylabel('Latitude',fontsize=10)
            ax.set_xlabel('Longitude',fontsize=10)
            ax.set_title(name.replace('\n',''), fontsize = 20)
        plt.show()
        plt.savefig(name)
            
#%% NCEDC Triangulation Match
for j,group in enumerate(triMatch['ncedc']):
    for i in range(5):
        if group[0] == i:
            name = 'Triangulation Method for NCEDC Catalog Event ' + str(i)
            fig=plt.figure(name,figsize=(10,10))
            ax=plt.subplot()
            ax.scatter(group[1][0], group[1][1], color = 'xkcd:very dark blue')
            ax.scatter(ncedcLoc[i][0],ncedcLoc[i][1], color = 'g')
            plt.annotate('Catalog Event',(ncedcLoc[i][0],ncedcLoc[i][1]))
            a = plt.Circle((group[1][0], group[1][1]),group[1][2]/111, color='r',fill=False)
            ax.set_aspect('equal', adjustable='datalim')
            ax.add_patch(a)
            ax.set_ylabel('Latitude',fontsize=10)
            ax.set_xlabel('Longitude',fontsize=10)
            ax.set_title(name.replace('\n',''), fontsize = 20)
        plt.show()
        plt.savefig(name)
            