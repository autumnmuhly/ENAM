from obspy import read_inventory, read, UTCDateTime
from obspy.signal import rotate
from obspy.geodetics.base import gps2dist_azimuth
import sys
import os
import glob
import numpy as np
import subprocess

event_id='2014-06-25_11.52.02'
station_dir = event_id+'/stations'
station_rot_dir =event_id+'/stations_rotated'
sac_dir = event_id+'/sac'
sac_rot_dir = event_id+'/sac_rotated'

infile_orientations="./YO_sta_orientations.txt"


# ---- Rotate to ENZ data -----------
# Write data to (sac) files
sac_file_name = sac_rot_dir+'/'+tn.stats['network']+'.'+\
                tn.stats['station']+'.'

#Rotate to ENZ
if isset('t1') and isset('t2'):
    print(inv[0][0][0].azimuth,inv[0][0][0].dip)
    print(inv[0][0][1].azimuth,inv[0][0][1].dip)
    print(inv[0][0][2].azimuth,inv[0][0][2].dip)
    [tz.data,tn.data,te.data] = rotate.rotate2zne(tz.data,inv[0][0][2].azimuth,inv[0][0][2].dip,
                                        t1.data,inv[0][0][0].azimuth,inv[0][0][0].dip,
                                        t2.data,inv[0][0][1].azimuth,inv[0][0][1].dip)
    te.stats['channel'] = te.stats['channel'][:-1]+'E'
    tn.stats['channel'] = tn.stats['channel'][:-1]+'N'


#Write out ENZ to SAC
tn.stats.sac=sachdr
te.stats.sac=sachdr
tz.stats.sac=sachdr

tn.write(sac_file_name+tn.stats['channel']+'.sac', format='sac')
te.write(sac_file_name+te.stats['channel']+'.sac', format='sac')
tz.write(sac_file_name+tz.stats['channel']+'.sac', format='sac')



# ---- Rotate to RTZ data -----------
# Write data to (sac) files
sac_file_name = sac_rot_dir+'/'+tn.stats['network']+'.'+\
                tn.stats['station']+'.'

tr = tn.copy()
tt = te.copy()

# Rotate N-E to R-T
[res,az,baz]=gps2dist_azimuth( elat, elon,\
                                inv[0][0].latitude,inv[0][0].longitude)

[tr.data,tt.data] = rotate.rotate_ne_rt(tr.data,tt.data,baz)
tr.stats['channel'] = tr.stats['channel'][:-1]+'R'
tt.stats['channel'] = tt.stats['channel'][:-1]+'T'


tr.stats.sac=sachdr
tt.stats.sac=sachdr
tz.stats.sac=sachdr
tr.write(sac_file_name+tr.stats['channel']+'.sac', format='sac')
tt.write(sac_file_name+tt.stats['channel']+'.sac', format='sac')
tz.write(sac_file_name+tz.stats['channel']+'.sac', format='sac')
