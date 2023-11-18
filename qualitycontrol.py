# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:36:44 2022

@author: Austin Abreu
"""
###############################################################################
# This function analyzes a set of waveforms relative to a station, then uses
# simple fuzzy matching methods to correlate pick times between catalogs.

#  BE CAREFUL WHEN USING THIS FUNCTION, IT LIKES TO ALLOCATE LOTS OF MEMORY
#   This occurs primarily in figure generation and obspy stream handling.
#%% Libraries
from obspy import UTCDateTime, Stream, read, Catalog
from obspy.clients.fdsn import Client
from obspy.core.event.catalog import read_events
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt



def reorganizer(cat, keys): #requires an Obspy catalog object
    # keys is defined as the list of station names that are relevent to the task.
        #All stations not in keys will be parsed out to improve processing time.
    #Dictionaries allow for faster processing and easy mapping of station
        #names to data
    stations={ #used to map station names to pick indices
        }
    channels={  #map station names to pick-relative channels
        }
    picks={ #map station names to pick time data
        }
    origins = {
        }
    for j,event in enumerate(cat):
        origins.setdefault(j,event.origins[0].time)
        for i,p in enumerate(event.picks):
            stations.setdefault(p.waveform_id.station_code,[]).append((j,i))
            picks.setdefault(p.waveform_id.station_code,[]).append(p.time)
            channels.setdefault(p.waveform_id.station_code,[]).append(
                p.waveform_id.channel_code)
            
    bestStations = stations.copy() #copying the dictionary is probably the best way
                                        # to handle this
    bestPicks = picks.copy()
    bestChannels = channels.copy()
    for key in stations:
        if key not in keys:
            bestStations.pop(key) #since we copied the dictionaries we can use keys 
                                    #interchangeably
            bestPicks.pop(key)
            bestChannels.pop(key)
            
    del stations, channels, picks #discarding unused vars
    return bestStations, bestChannels, bestPicks, origins


#%%Catalog Imports
ncedc = pd.read_csv("geysers_NCEDC_catalog2006_2016.csv",header = 0, sep = ",")
lbnl = pd.read_csv("geysers_LBL_catalog2006_2016.csv",header = 0, sep = ",")
ncedc["DateTime"] = pd.to_datetime(ncedc["DateTime"]) #The times are held as strings so we must convert them
lbnl["DateTime"] = pd.to_datetime(lbnl["DateTime"])

#EQT Picks & catalog
EQTcat = read_events('associations.xml')

#Misc
ncedcColor = 'xkcd:medium blue'
lbnlColor = 'xkcd:brick red'


#%% Waveform Import
st = Stream()
st = read('downloads_mseeds\\DEB\\*.mseed') #Set explicity or use wildcard.
keys = set(['DEB']) #adjust appropriate to your use-case
#%% Waveform pre-processing
sTime = UTCDateTime(2015,2,23)
eTime = UTCDateTime(2015,2,24)

st = st.trim(sTime,eTime)

stW = st.copy()
stW = stW.detrend('linear') #detrend
for tr in stW: #demean
    tr.data = tr.data-np.mean(tr.data)
del tr
stW = stW.taper(max_percentage=0.05) #taper
stW = stW.filter("bandpass",freqmin=0.5,freqmax=20,zerophase=True)

stW = stW.merge()

del st
#%% APPEND EQT CAT TO WAVEFORM
stations, channels, picks, origins = reorganizer(EQTcat, keys)

relPicks = {}
relStations = {}
relChannels = {}
# for key, time in picks.items()):
#     if sTime < time < eTime:
#         relPicks.setdefault(key,[]).append(time)
#         relStations.setdefault(key,[]).append(time)

for key,i in picks.items():
    for k,val in enumerate(i):    
        if sTime < val < eTime:
            relPicks.setdefault(key,[]).append(val) 
            relChannels.setdefault(key,[]).append(channels[key][k])
            relStations.setdefault(key,[]).append(stations[key][k])
            
for key, val in origins.items():
    if val < sTime and val > eTime:
        origins.pop(key)
            
lbnlPicks = lbnl['DateTime'][(lbnl['DateTime'] > sTime.datetime) & (lbnl['DateTime'] < eTime.datetime)]
ncedcPicks = ncedc['DateTime'][(ncedc['DateTime'] > sTime.datetime) & (ncedc['DateTime'] < eTime.datetime)]

for tr in stW:
    tr.stats.setdefault('EQTPicks',[])
    tr.stats.EQTPicks.clear()
    tr.stats.setdefault('EQTPhase',[])
    tr.stats.EQTPhase.clear()
    #These catalogs are station agnostic, so we'll just slap em in
    tr.stats.setdefault('LBNLPicks',list(lbnlPicks)) 
    tr.stats.setdefault('NCEDCPicks',list(ncedcPicks))
del tr
    
#Now we can append EQT picks to the traces
for key in list(relStations.keys()):
        #matching each key to its channels from the catalog
        for i,cha in enumerate(relChannels[key]):
            #Now I match the downloaded stream to the station (key) AND 
                #respective channels.
            try:
                for tr in stW.select(station=key,channel=cha):
                #Nothing happens if it doesn't match; if it does:
                    tr.stats.setdefault('EQTPicks').append(relPicks[key][i])
                    tr.stats.setdefault('EQTPhase').append(
                        EQTcat[relStations[key][i][0]].picks[
                            relStations[key][i][0]].phase_hint)
                #I attach a dictionary to tr.stats containing a list of the 
                    #picks for the proper channel!
            except Exception:
                for tr in stW:
                    tr.stats.setdefault('EQTPicks').append(relPicks[key][i])
                    tr.stats.setdefault('EQTPhase').append(
                        EQTcat[relStations[key][i][0]].picks[
                            relStations[key][i][1]].phase_hint)
#This function is flexible to any number of station names, channels
    # it can be modified to use any catalog, as long as that catalog is an Obspy object
del stations, channels, EQTcat,lbnl,ncedc, key, i, cha, tr
#%%
fig, ax = plt.subplots(1,figsize=(19,10))

chan = stW[0].stats.channel
sta = stW[0].stats.station
net = stW[0].stats.network
stationID = net+'.'+sta+'$_{'+chan+'}$'

#Waveform Plotting
ax.plot(stW[0].times("matplotlib"),stW[0].data,label=stationID + ' Signal',color = 'k')

#Pick Plotting
for j in range(len(stW[0].stats.EQTPicks)):
    ax.axvline(x=stW[0].stats.EQTPicks[j].datetime, color='g',
                   linewidth = 1)
    
for k in range(len(stW[0].stats.LBNLPicks)):
    ax.axvline(x=stW[0].stats.LBNLPicks[k], color=lbnlColor,
                   linewidth = 0.5)
    
for a in range(len(stW[0].stats.NCEDCPicks)):
    ax.axvline(x=stW[0].stats.NCEDCPicks[a], color=ncedcColor,
                   linewidth = 0.5)
    
# Plot Traits
ax.set_ylabel("Counts",fontsize=10,color = 'k')
ax.set_title(stationID,
                 fontsize=20,color='k');
ax.tick_params('both',labelsize=10,color='k',labelcolor='k',width=3)
    
#%% Event Correlation Algorithm
Corr = {'lbnl':[],'ncedc':[]}

#Pulls EQT and lbnl times and computes a time difference between them.
    # If the difference falls within a threshold, it is appended to a dictionary
    # Containing correlated times per station.
dateDiff = []
for i,EQTtime in origins.items():
    for j,lbnltime in enumerate(lbnlPicks):
        dateDiff += [EQTtime.datetime - lbnltime]
        #some weird thing happens where it loops around the start of the day,
            # hence the != -1 at the end.
        if dateDiff <= dt.timedelta(0,3) and dateDiff >= -dt.timedelta(0,3) and dateDiff._d == 0:
            Corr['lbnl'].append((i,j))
    for k,ncedctime in enumerate(ncedcPicks):
        dateDiff = EQTtime.datetime - ncedctime
        if dateDiff <= dt.timedelta(0,3) and dateDiff >= -dt.timedelta(0,3) and dateDiff._d == 0:
            Corr['ncedc'].append((i,k))
#A function of similar form can be used to correlate all three catalogs.

del dateDiff, lbnltime, ncedctime
#%%%    
fig, ax = plt.subplots(1,figsize=(19,10))

#Metadata
chan = stW[0].stats.channel
sta = stW[0].stats.station
net = stW[0].stats.network
stationID = net + '.' + sta + '.' + chan

#Waveform Plotting
ax.plot(stW[0].times("matplotlib"),stW[0].data,label=stationID + ' Signal',color = 'k')

#CORRELATED Pick Plotting
for pair in Corr['lbnl']:
    ax.axvline(x=origins[pair[0]].datetime, color=lbnlColor,
                   linewidth = .75)
    ax.axvline(x=stW[0].stats.LBNLPicks[pair[1]], color=lbnlColor,
                   linewidth = .75)
        
for pair in Corr['ncedc']:
    ax.axvline(x=origins[pair[0]].datetime, color=ncedcColor,
                   linewidth = .75)
    ax.axvline(x=stW[0].stats.NCEDCPicks[pair[1]], color=ncedcColor,
                   linewidth = .75)
    
# Plot Traits
ax.set_ylabel("Counts",fontsize=10,color = 'k')
ax.set_title(stationID,
                 fontsize=20,color='k');
ax.tick_params('both',labelsize=10,color='k',labelcolor='k',width=3)
