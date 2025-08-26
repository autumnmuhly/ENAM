#Convert from miniseed to sac 
from obspy import UTCDateTime, read, read_inventory
from obspy.core.event import read_events
import os 

print(" >> Converting MSEED to SAC files")

event_id='2014-06-15_20.14.51'


# Create folder for sac files


station_dir = event_id+'/stations'
station_rot_dir =event_id+'/stations_rotated'
sac_dir = event_id+'/sac'
sac_rot_dir = event_id+'/sac_rotated'
if os.path.exists(sac_dir):
    os.system('rm -r '+sac_dir)
os.mkdir( sac_dir )

with open('evlist.txt', "r") as infile:
    for line in infile:
        items = line.split()
        if items[10] == event_id:
            elat=items[7]
            elon=items[8]
            edep=items[9]
        else:
            print('need to change the event id')

# If using quakeml format to get evetn data
# cat = read_events('2014-06-25_11.52.02/2014-06-25_11.52.02.qml')
# for evt in cat:
#     event=evt.origins[0]
#     elat=event.latitude
#     elon=event.longitude
#     edep=event.depth
    
def isset(v):
    try:
        type (eval(v))
    except:
        return 0
    else:
        return 1

final_record_type = "VEL" #Options: "DISP", "VEL" or "ACC"
# Sampling rate (in Hz)
fs = 40
# record length
rlength =3600


wav_dir = event_id+'/waveforms'
station_dir = event_id+'/stations'
sac_dir = event_id+'/sac'

fname_stations = station_dir+'/*.xml'
fname_mseed    = wav_dir+'/*.mseed'

# Read traces and station information
trs = read(fname_mseed)
sts = read_inventory(fname_stations)

    

# Get station names
receivers = sts.get_contents()['stations']
for rcv in receivers:
    #print((  "%s station" %(str(rcv.split('.')[1].split(' ')[0])) ))
    # Read the corresponding traces to the receiver
    network = str(rcv.split('.')[0])
    rcvname = str(rcv.split('.')[1].split(' ')[0])

    ts=trs.select(network=network,station=rcvname)
    #print(ts)
    inv=sts.select(network=network,station=rcvname)

        
    # Check that all 3 components were downloaded
    if len(ts) != 3:
        print(( "Not all the 3 components are read for %s %s" %(network,rcvname)))
        print( "Skipping")
    # Remove instrument response
    for ti in ts:
        #print(ti)
        pre_filt = [0.008, 0.016, 16, 32]
        ti.remove_response(inventory=inv, output=final_record_type,
                        pre_filt=pre_filt,
                        water_level=60, plot=False)
#            fstop = 1./2.*fs
#            fpass = 1./2.*fs - 0.0001
#            ti.remove_response(inventory=inv, output=final_record_type,
#                               pre_filt=(0.,1./3600,fpass,fstop),
#                               water_level=60, plot=False)
        tname=(ti.get_id())
        chan=(ti.stats.channel)
        slat=(inv[0][0].latitude)
        slon=(inv[0][0].longitude)
        selv=(inv[0][0].elevation)
                        
        # Resample data
        ts = ts.resample( sampling_rate = fs )

        # Trim data
        stime0 = UTCDateTime(ts[0].stats.starttime)
        stime1 = UTCDateTime(ts[1].stats.starttime)
        stime2 = UTCDateTime(ts[2].stats.starttime)
        stmax  = max(stime0,stime1,stime2)

        etime0 = UTCDateTime(ts[0].stats.endtime)
        etime1 = UTCDateTime(ts[1].stats.endtime)
        etime2 = UTCDateTime(ts[2].stats.endtime)
        etmin  = min(etime0,etime1,etime2)

        ref_time = stmax        
        etmin = ref_time + rlength
        ts = ts.trim(stmax,etmin,pad=True,fill_value=0)

        # No rotation for OBS data (to rotate use script OBS_rotate_12ZtoNEZ.py)
    if isset('t1') and isset('t2'):
        del t1,t2
    t1=5
    for ti in ts:
        print(f' this is ti {ti.id[-1]}')
        if str(ti.id[-1]) == 'E':
            te = ti
        elif str(ti.id[-1]) == 'N':
            print('equal to N')
            tn = ti
            print(f't1 doesnt exist for {sac_dir}+{tn.stats['network']}')
        elif str(ti.id[-1]) == 'Z':
            tz = ti
        elif str(ti.id[-1]) == '1':
            t1 = ti
            tn = ti
        elif str(ti.id[-1]) == '2':  
            t2 = ti
            te = ti
        else:
            print( 'Unknown component')
            skip = 1
    #SAC header
    sachdr = {
        'evla': elat,
        'evlo': elon,
        'evdp': edep,
        'stla': slat,
        'stlo': slon,
        'stel': selv
    }
    # ---- 12Z data -----------
    # Write data to (sac) files
    if t1==5:
        sac_file_name = sac_dir+'/'+tn.stats['network']+'.'+\
                        tn.stats['station']+'.'
        #Write out NEZ to SAC
        tn.stats.sac=sachdr
        te.stats.sac=sachdr
        tz.stats.sac=sachdr
        tn.write(sac_file_name+tn.stats['channel']+'.sac', format='sac')
        te.write(sac_file_name+te.stats['channel']+'.sac', format='sac')
        tz.write(sac_file_name+tz.stats['channel']+'.sac', format='sac')
        # want to copy those that are aleady in coord to the sac directory
        #os.system(f"cp /usc/data/ADEPT/{evt_name}/{evt_name}.{sacfile.stanm}.{sacfile.network}*R.D.sac .")
    if t1 !=5 :
        sac_file_name = sac_dir+'/'+tn.stats['network']+'.'+\
                        t1.stats['station']+'.'

        #Write out 12Z to SAC
        t1.stats.sac=sachdr
        t2.stats.sac=sachdr
        tz.stats.sac=sachdr
        print(f'writing files to .sac to {sac_dir} {t1}')
        t1.write(sac_file_name+t1.stats['channel']+'.sac', format='sac')
        t2.write(sac_file_name+t2.stats['channel']+'.sac', format='sac')
        tz.write(sac_file_name+tz.stats['channel']+'.sac', format='sac')

