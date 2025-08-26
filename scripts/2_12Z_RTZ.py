from obspy import read_inventory, read, UTCDateTime
from obspy.signal import rotate
from obspy.geodetics.base import gps2dist_azimuth
from obspy import Stream
import sys
import os
import glob
import numpy as np

# ------------------- Auxiliar functions ------------------- # 
def isset(v):
    try:
        type (eval(v))
    except:
        return 0
    else:
        return 1
    
#Set up inputs
infile_orientations="./YO_sta_orientations.txt" #Format: stanm, bh1_azimuth, bh1_azimuth_stdev
event_id='2014-06-15_20.14.51'

# Prepare output waveform and station directories
wav_dir = event_id+'/waveforms'
sac_rot_dir = event_id+'/sac_rotated'
station_dir = event_id+'/stations'
station_rot_dir =event_id+'/stations_rotated'
sac_dir = event_id+'/sac'

#Options: "DISP", "VEL" or "ACC"
final_record_type = "VEL" 


# # Sampling rate (in Hz)
# fs = 40
eventfile="evlist.txt"
with open('evlist.txt', "r") as infile:
    for line in infile:
        items = line.split()
        if items[10] == event_id:
            elat=float(items[7])
            elon=float(items[8])
            edep=float(items[9])
        else:
            print('need to change the event id')

# #Phase data
# model1d="prem"
# phasename="SKKS"

# if len(glob.glob(station_dir+'/*.xml'))==0:
# 	print ("ERROR: no xml files in "+station_dir+". Quitting")
# 	sys.exit()

# # Create folder for rotated station files
# if os.path.exists(station_rot_dir):
#     os.system('rm -r '+station_rot_dir)
# os.mkdir( station_rot_dir )


# Create folder for rotated sac files
if os.path.exists(sac_rot_dir):
    os.system('rm -r '+sac_rot_dir)
os.mkdir( sac_rot_dir )


# fname_stations = station_dir+'/*.xml'
# fname_sac= sac_dir+'/*.sac'

# #read in station data
# allstations=read_inventory(fname_stations)

# #update station xml
# with open (infile_orientations) as azfile:
#     for line in azfile:
#         sta,azi=line.split(',',2)
#         azi=float(azi)
#         for n in allstations:
#             for s in n.stations:
#                 for c in s.channels:
#                     if c.code.startswith("BH") or c.code.startswith('HH'):
#                         if c.code.endswith('1'):
#                             c.azimuth=azi
#                             c.dip = 0
#                         elif c.code.endswith('2'):
#                             c.azimuth=np.round(((azi+90) % 360),1)
#                             c.dip=0
#                 allstations.write(station_rot_dir+"/"+n.code+"."+s.code+".xml", format="STATIONXML")

#rotate sac files now with update xml
fname_H1sac= sac_dir+'/*.sac'
st=read(fname_H1sac)
allxml=read_inventory(station_rot_dir+'/*.xml')
NErotStream=st.rotate('->ZNE',inventory=allxml)

RTrotStream=Stream()
stations=set(tr.stats.station for tr in NErotStream)

#get backazi and rotate to RT
for sta in stations:
    trace=NErotStream.select(station=sta)
    code=trace[0].id
    srttime=trace[0].stats.starttime
    loc=allxml.get_coordinates(code,srttime)
    stalat=loc['latitude']
    stalon=loc['longitude']
    stael=loc['elevation']
    [res,az,backazi]=gps2dist_azimuth( elat, elon,stalat,stalon)
    #print(res,az,baz)
    rotated=trace.rotate('NE->RT',inventory=allxml,back_azimuth=backazi)
    RTrotStream+=rotated

for trace in RTrotStream:
    sac_file_name = sac_rot_dir+'/'+trace.stats['network']+'.'+trace.stats['station']+'.'
    trace.write(sac_file_name+trace.stats['channel']+'.sac', format='sac')