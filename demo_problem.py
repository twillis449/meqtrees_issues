# standard preamble
#
#% $Id$ 
#
#
# Copyright (C) 2002-2007
# The MeqTree Foundation & 
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

# This script uses a standard pointing model for AzEl mounted telescopes
# to generate an observation which has been corrupted by pointing errors.
# Adapted from Oleg's workshop script example10-tracking.py
# The script includes the possibility of a parallactic angle track tilt - 
# i.e we simulate a 'sky rotator' error.


from Timba.TDL import *
from Timba.Meq import meq
from Timba import pynode
import math
import random
import numpy

import Meow

from Meow import Bookmarks
import Meow.StdTrees
import sky_models

# defining MS at compile time doesn't seem to help ...
TDLCompileOptionSeparator("MS selection");
# MS options first
mssel = Meow.Context.mssel = Meow.MSUtils.MSSelector(has_input=True,tile_sizes=None,
                  read_flags=True,write_flags=True)
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
TDLCompileMenu('L and M position of phased-up beam',
  TDLOption('l_beam','L position of beam (in units of FWHM)',[0,1,2,3],more=float),
  TDLOption('m_beam','M position of beam (in units of FWHM)',[0,1,2,3],more=float),
);
# get Azimuth / Elevation telescope sinusoidal random tracking errors
TDLCompileMenu('Additional Telescope Tracking Errors',
TDLOption('min_random_error','Minimum Random Error (arcsec)',[0,5,10,20],more=float),
TDLOption('max_random_error','Maximum Random Error (arcsec)',[0,5,10,20],more=float),
# Hold RMS jitter constant over tile?
TDLOption('hold_rms_constant','hold pointing jitter constant over tile?',[True, False]),

TDLOption('max_thermal_deflection','Maximum Thermal Deflection (arcsec)',[0,5,10,20],more=float),
);

# if not None, a per-ifr noise term with the given stddev will be added
TDLCompileOption('noise_stddev',"Add background noise (Jy)",[0,1e-9,1e-8,1e-6,5e-6,1e-3,0.05,0.1,0.33,0.5],more=float);
TDLCompileOption('subtract_perfect_obs','Subtract perfect observation from error observation?',[False,True])


# MS run-time options
TDLRuntimeMenu("Data selection & flag handling",*mssel.runtime_options());
# some GUI options
Meow.Utils.include_ms_options(has_input=False,tile_sizes=[16,32,48,96]);

TDLRuntimeMenu("Imaging options",
    *Meow.Utils.imaging_options(npix=256,arcmin=sky_models.imagesize(),channels=[[32,1,1]]));

# useful constant: 1 deg in radians
DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = ARCMIN/60;

# define desired half-intensity width of power pattern (HPBW)
# as we are fitting total intensity I pattern (here .021 rad = 74.8 arcmin)
fwhm  = 0.021747 # beam FWHM
fwhm1  = fwhm * 1.01
fwhm_sq = fwhm * fwhm
fwhm1_sq = fwhm1 * fwhm1


# for offset FPA
#fwhm_offset = 2.06 * fwhm
fwhm_offset = 0.0

ln_16  = -2.7725887
def ASKAP_elliptical_voltage_response(E, lm):
  """This makes the nodes to compute the beam gain, E, given an lm position.
  'lm' is an input node, giving sqr of position of source
  'E' is an output node to which the gain will be assigned""";
  # gaussian to which we want to optimize beams
  lmsq = E('lmsq') << Meq.Sqr(lm);
  lsq = E('lsq') << Meq.Selector(lmsq,index=0);
  msq = E('msq') << Meq.Selector(lmsq,index=1);
  lm_sq  = lsq/fwhm_sq + msq/fwhm1_sq
  E << Meq.Sqrt(Meq.Exp(lm_sq * ln_16))
  return E

def ASKAP_voltage_response(E, lm,):
  """This makes the nodes to compute the beam gain, E, given an lm position.
  'lm' is an input node, giving position of source 
  'E' is an output node to which the gain will be assigned""";
  # gaussian to which we want to optimize beams
  lmsq = E('lmsq') << Meq.Sqr(lm);
  lsq = E('lsq') << Meq.Selector(lmsq,index=0);
  msq = E('msq') << Meq.Selector(lmsq,index=1);
  lm_sq  = (lsq + msq) / fwhm_sq
  E << Meq.Sqrt(Meq.Exp(lm_sq * ln_16))
  return E

def noise_matrix (stddev=0.1):
  """helper function to create a 2x2 complex gaussian noise matrix""";
  noise = Meq.GaussNoise(stddev=stddev);
  return Meq.Matrix22(
    Meq.ToComplex(noise,noise),Meq.ToComplex(noise,noise),
    Meq.ToComplex(noise,noise),Meq.ToComplex(noise,noise)
  );

def _define_forest (ns):

  array,observation = mssel.setup_observation_context(ns);

# Get RA, Dec of Sun
  ns.sun <<Meq.ObjectRADec(obj_name="SUN")

  if noise_stddev > 0.0:
    # create per-antenna noise term
    for sta in array.stations():
      ns.noise(sta) << noise_matrix(noise_stddev/2);

  if abs(min_random_error) > 0.0:
    ns.min_random << Meq.Constant(min_random_error)
    ns.max_random << Meq.Constant(max_random_error)
    for p in array.stations():
      ns.random_Az(p)  << Meq.RandomNoise(ns.min_random, ns.max_random)
      ns.random_El(p)  << Meq.RandomNoise(ns.min_random, ns.max_random)
      if hold_rms_constant:
        ns.const_random_Az(p) << Meq.PyNode(ns.random_Az(p),ns.min_random, ns.max_random,class_name="PyRandomConst",module_name="PYRandomConst");
        ns.const_random_El(p) << Meq.PyNode(ns.random_El(p),ns.min_random, ns.max_random,class_name="PyRandomConst",module_name="PYRandomConst");
        ns.ConstRandomPoint(p) << Meq.Composer(ns.const_random_Az(p), ns.const_random_El(p)) * ARCSEC
      else:
        ns.RandomPoint(p) << Meq.Composer(ns.random_Az(p), ns.random_El(p)) * ARCSEC
  else:
    for p in array.stations():
      ns.RandomPoint(p) << Meq.Composer(0.0, 0.0)

# get its elevation
  for p in array.stations():
    ns.AzEl_stn(p) << Meq.AzEl(radec=observation.phase_centre.radec(), xyz=array.xyz(p))
    ns.AzEl_sun(p) << Meq.AzEl(radec=ns.sun, xyz=array.xyz(p))
    ns.Az_diff(p) << Meq.PyNode(ns.AzEl_sun(p),ns.AzEl_stn(p),class_name="PyDiffAzimuth",module_name="PYDiffAzimuth");
    ns.pos_elev_sun(p) << Meq.PyNode(ns.AzEl_sun(p),class_name="PyPosElevation",module_name="PYPosElevation");
#   ns.pos_elev_sun(p) << Meq.PyNode(ns.AzEl_sun(p),class_name="PyPosElevation",module_name=__file__);
    

  # create a source model and make list of corrupted sources
  allsky = Meow.Patch(ns,'all',observation.phase_centre);
  if subtract_perfect_obs:
    perfectsky = Meow.Patch(ns,'perfect',observation.phase_centre);
  sources = sky_models.make_model(ns,"S",l_beam*fwhm,m_beam*fwhm);
  for src in sources:
    lm = src.direction.lm();
    E = ns.E(src.name);
    if subtract_perfect_obs:
      EP = ns.EP(src.name);

    for p in array.stations():
      # thermal deviation is ~ 5 arcsec when Sun is at maximum elevation
      ns.sine_stn_elev(p) << Meq.Sin(ns.pos_elev_sun(p))
      ns.sine_az_diff(p) << Meq.Sin(ns.Az_diff(p))
      ns.thermal(p) << Meq.Abs(ns.sine_az_diff(p)) * max_thermal_deflection * ARCSEC * ns.sine_stn_elev(p)

      if hold_rms_constant:
        lm1 = ns.lm1(src.name,p) << lm + ns.ConstRandomPoint(p) + ns.thermal(p)
      else:
        lm1 = ns.lm1(src.name,p) << lm + ns.RandomPoint(p) + ns.thermal(p)

      # compute E for apparent position
      ASKAP_voltage_response(E(p),lm1);

      if subtract_perfect_obs:
        ASKAP_voltage_response(EP(p),lm);
    allsky.add(src.corrupt(E));
    if subtract_perfect_obs:
      perfectsky.add(src.corrupt(EP));

  observed = allsky.visibilities();

 # throw in a bit of noise?
  if noise_stddev > 0:
    for p,q in array.ifrs():
      ns.noisy_predict(p,q) << observed(p,q) + ns.noise(p) + ns.noise(q);
    observed = ns.noisy_predict;
  if subtract_perfect_obs:
    perfect = perfectsky.visibilities()
    for p,q in array.ifrs():
      ns.diff(p,q) << observed(p,q) - perfect(p,q);
     # make sinks and vdm. Note that we don't want to make any spigots...
    Meow.StdTrees.make_sinks(ns,ns.diff,spigots=False);


  # make some useful inspectors. Collect them into a list, since we need
  # to give a list of 'post' nodes to make_sinks() below
  pg = Bookmarks.Page("Inspectors",1,2);
  inspectors = [];
  inspectors.append(
    Meow.StdTrees.vis_inspector(ns.inspect_observed,observed) );
  pg.add(ns.inspect_observed,viewer="Collections Plotter");

 # create spigots, condeqs, residuals
  if not subtract_perfect_obs:
    Meow.StdTrees.make_sinks(ns,observed,spigots=False,post=inspectors);


def _test_forest(mqs,parent,wait=False):

  req = Meow.Utils.create_io_request();
  # execute    
  mqs.execute('VisDataMux',req,wait=wait);
  
# this is a useful thing to have at the bottom of the script,  
# it allows us to check the tree for consistency simply by 
# running 'python script.tdl'

if __name__ == '__main__':
# Timba.TDL._dbg.set_verbose(5);
  ns=NodeScope()
  _define_forest(ns)
  ns.Resolve()
  print "Added %d nodes" % len(ns.AllNodes())
