#
# Copyright (C) 2008
# ASTRON (Netherlands Foundation for Research in Astronomy)
# and The MeqTree Foundation
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, seg@astron.nl
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# CASA script that uses the CASA simulator to create
# UV tracks for SKA simulations

import random
import numpy
import os

# basic parameters - things which are most likely to be changed
observatory  = 'ASKAP'         # where the telescope is situated
RA0          = "0h0m1.0"     # field centre RA
#DEC0         = "-28d00m00s"   # field centre DEC
#Freq         = "1400MHz"     # frequency at lower edge of band
DEC0         = "-45d00m00s"   # field centre DEC
#DEC0         = "-80d00m00s"   # field centre DEC
Freq         = "1400MHz"     # frequency at lower edge of band
num_channels = 1             # number of channels
channel_inc  = '50.0MHz'      # channel increment
Diameter     = 12.0          # dish diameter in metres
num_dishes   = 12            # number of dishes in the array
Stokes       = 'XX XY YX YY' # Stokes parameters - other option is 'RR RL LR LL'
int_time     = '30.0s'         # integration period
#int_time     = '6s'         # integration period
ref_time = me.epoch('UT','18jun2014/20:00:00')  # a reference time
a = qa.quantity('21dec2013/20:00:00')
print 'a',a
#ref_time = me.epoch('UTC',a)  # a reference time
print 'ref_time', ref_time
starttime    = -21600.0        # start of obs, in seconds relative to ref_time
scanlength   =  2*21600.0      # length of observation, in sec
noise        =  '0.0Jy'      # add some noise. Use '0.0Jy' if no noise wanted
FOV          =  30            # field of view in arcmin
num_sources  = 5             # number of sources to observe, randomized over FOV

clname = 'mymodel.cl'        # file for component list
msname = 'ASKAP_12_test.MS'  # name of measurement set to create
# first delete old MS, test images, etc
print '*** deleting previous stuff ***'
system_call = "/bin/rm -rf " + msname + " " + clname
os.system(system_call)

# create sources with randomized offsets from field centre
dec_offset=[]
ra_offset=[]
for i in range(num_sources):
  offset = random.uniform(-FOV, FOV)
  dec_offset.append(str(offset)+'arcmin')
  offset = random.uniform(-FOV, FOV)
  ra_offset.append(str(offset)+'arcmin')

# now add above offsets to field centre
refRA = qa.unit(RA0);
refDec = qa.unit(DEC0);
for i in range(num_sources):
  ra = qa.add(refRA,ra_offset[i]);
  dec = qa.add(refDec,dec_offset[i]);
  direction = me.direction(rf='J2000',v0=ra,v1=dec)
  print 'source dir ', direction
  cl.addcomponent(dir=direction,flux=1.0, freq=Freq)

# the setspectrum function seems to be broken
# cl.setspectrum(0, 'spectral index', [-0.5, 0, 0, 0])

cl.rename(filename=clname)
cl.done()

#  Define ASKAP array, antenna locations from Tim Cornwell, local coordinates
    
xx =  [-2556112.856429219, -2556090.5197019787, -2556031.7010004967, -2555975.3071098486, -2556062.352267588, -2556500.22800737, -2555961.8252890203, -2555326.623093444, -2557352.698163088, -2556409.20407202, -2555597.447614275, -2556559.13952461]

yy =  [5097394.421311565, 5097429.818865128, 5097457.645149432, 5097238.697250006, 5097566.995432242, 5097341.658547914, 5096984.5323678795, 5098269.789240807, 5097178.811844573, 5097068.70966771, 5097844.084160582, 5097779.361367356]

zz =  [-2848443.383871815, -2848400.3722200543, -2848403.3368191277, -2848842.763482064, -2848181.6233792487, -2848191.8455700525, -2849306.44495641, -2847587.817717512, -2847721.1231166106, -2848758.1500504687, -2848103.3976325844, -2847361.063454385]


diam = numpy.zeros((num_dishes,), numpy.float64) + Diameter;
ant_names = [ str(i) for i in range(1,num_dishes+1) ]

#Following simulates the ASKAP configuration
#askapwgs=me.position(rf='wgs84', v0='116d39m32.0', v1='-26d42m15', v2='125.0m')
#askapitrf=me.measure(v=askapwgs, rf='ITRF')

sm.open(msname)
print 'opened MS'
sm.setspwindow(spwname='SKA', freq=Freq,
		      deltafreq=channel_inc, freqresolution=channel_inc, 
		      nchannels=num_channels, stokes=Stokes)
pos_obs = me.observatory(observatory)
sm.setconfig(telescopename=observatory, 
       x=xx, y=yy,z=zz, dishdiameter=diam.tolist(), mount='ALT-AZ', antname = ant_names, padname=ant_names,
       coordsystem='global', referencelocation=pos_obs)

print 'setting up simulator specifications'
dir0 = me.direction(rf="J2000",  v0=RA0, v1=DEC0)
print 'field centre', dir0 
sm.setfield(sourcename='ASKAP_12', sourcedirection=dir0)
print 'passed sm.setfield'
sm.settimes(integrationtime=int_time, usehourangle=True, referencetime=ref_time)
print 'passed sm.settimes'
#sm.settimes(integrationtime=int_time, referencetime=ref_time)
sm.setlimits(shadowlimit=0.001, elevationlimit='8.0deg')
print 'passed sm.setlimits'
sm.setauto(autocorrwt=0.0)
print 'setting up for w projection'
#sm.setoptions(ftmachine='wproject',wprojplanes=512)
print 'specifying model image'
#model_image = 'wproj_random_test_image'

scan=0
endtime = starttime + scanlength
midtime = starttime + scanlength * 0.5
while(starttime<endtime):
  print ' **** observing'
  sm.observe('ASKAP_12', 'SKA', starttime=str(starttime)+'s', 
                                  stoptime=str(starttime+scanlength)+'s')
  print 'called sm.observe' 
  sm.setdata(msselect='SCAN_NUMBER=='+ str(scan))
  print 'calculated sm.setdata'
  print '**** predicting'
  sm.predict(complist=clname)
# sm.predict(imagename=model_image)
  starttime=starttime+scanlength
  scan=scan+1

if noise != '0.0Jy':
  try:
    sm.setnoise(mode='simplenoise', simplenoise=noise)
    print ' **** corrupting'
    sm.corrupt();
  except:
    print ' **** failure when setting noise for corruption'
    pass

sm.done()
