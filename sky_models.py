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

# This script is just Oleg's sky_models.py script from the January 2007
# workshop 

from Timba.TDL import *
import Meow
import math
import random

# define desired half-intensity width of power pattern (HPBW)
# as we are fitting total intensity I pattern (here .021 rad = 74.8 arcmin)
FWHM  = 0.021747 # beam FWHM
#FWHM  = 0.02 # cortes beam 1 FWHM
#FWHM  = 0.023 #  Bruce central ant pat
#FWHM  = 0.02425  # from visibility obs
#FWHM  = 0.0219597 # cortes beam 0 FWHM
#FWHM = 0.021816595

DEG = math.pi/180.0
ARCMIN = DEG/60;
ARCSEC = ARCMIN/60;

ARCSEC_8 = 8.0 * ARCSEC

def imagesize ():
  """Returns current image size, based on grid size and step""";
  return grid_size*grid_step;

#def point_source (ns,name,l,m,I=1):
def point_source (ns,name,l,m,I=2):
  """shortcut for making a pointsource with a direction""";
  l = int(l / ARCSEC_8) * ARCSEC_8
  m = int(m / ARCSEC_8) * ARCSEC_8
  srcdir = Meow.LMDirection(ns,name,l,m);
  return Meow.PointSource(ns,name,srcdir,I=I);

def point_source_poln (ns,name,l,m,I,Q,U,V):
  """shortcut for making a pointsource with a direction""";
  srcdir = Meow.LMDirection(ns,name,l,m);
  return Meow.PointSource(ns,name,srcdir,I=I,Q=Q,U=U,V=V);


def cross_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a cross""";
  model = [point_source(ns,basename+"+0+0",l0,m0,I)];
  dy = 0;
  for dx in range(-nsrc,nsrc+1):
    if dx:
      name = "%s%+d%+d" % (basename,dx,dy);
      model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  dx = 0;
  for dy in range(-nsrc,nsrc+1):
    if dy:
      name = "%s%+d%+d" % (basename,dx,dy);
      model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  return model;

def mbar_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a line along the m axis""";
  model = [];
  model.append(point_source(ns,basename+"+0",l0,m0,I));
  for dy in range(-nsrc,nsrc+1):
    if dy:
      name = "%s%+d" % (basename,dy);
      model.append(point_source(ns,name,l0,m0+dm*dy,I));
  return model;

def lbar_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a line along the m axis""";
  model = [];
  model.append(point_source(ns,basename+"+0",l0,m0,I));
  for dx in range(-nsrc,nsrc+1):
    if dx:
      name = "%s%+d" % (basename,dx);
      model.append(point_source(ns,name,l0+dl*dx,m0,I));
  return model;
  
def star8_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in an 8-armed star""";
  model = [point_source(ns,basename+"+0+0",l0,m0,I)];
  for n in range(1,nsrc+1):
    for dx in (-n,0,n):
      for dy in (-n,0,n):
        if dx or dy:
          name = "%s%+d%+d" % (basename,dx,dy);
          model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  return model;

def grid_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a grid""";
# model = [point_source(ns,basename+"+0+0",l0,m0,I)];
  model = [point_source(ns,basename+"+0+0",l0,m0,2)];
  for dx in range(-nsrc,nsrc+1):
    for dy in range(-nsrc,nsrc+1):
      if dx or dy:
        name = "%s%+d%+d" % (basename,dx,dy);
        model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,2));
  return model;

def random_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged randomly in the field""";
  model = []
  num_sources = 2 + nsrc*2
  print 'constructing random field with sources = ', num_sources
  for i in range(num_sources):
#   x_offset = random.uniform(-4*dl, 4*dl)
#   y_offset = random.uniform(-4*dm, 4*dm)
    x_offset = random.uniform(-1*dl, 1*dl)
    y_offset = random.uniform(-1*dm, 1*dm)
#   I = random.uniform(0.02,0.2)
    I = 2.0
    print 'iteration, signal ', i,I
    name = "%s%+d" % (basename,i)
    model.append(point_source(ns,name,l0+x_offset,m0+y_offset,I));
#   model.append(point_source(ns,name,l0+x_offset,m0+y_offset,0.0));
  return model;

def ASKAP_power_response(l, m):
  """This makes the nodes to compute the beam gain, E, given an lm position.
  'lm' is an input node, giving position of source
  'E' is an output node to which the gain will be assigned""";
  fwhm  = FWHM

  fwhm_sq = fwhm * fwhm
  ln_16  = -2.7725887
  lsq = l * l
  msq = m * m
  lm_sq  = (lsq + msq) / fwhm_sq
  E = math.exp(lm_sq * ln_16)
  return E

def source_count(ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources from Stil model""";
  from string import split, strip
  filename = 'Tony.Willis_1.dat'
  # note we 'reset' ref 0 RA to -1 deg so we have extra 'room' in + L direction
  text = open(filename, 'r').readlines()
  source_num = -1
  model = []
  ref_x = l0 / DEG
  ref_y = m0 / DEG
  print 'ref_x ref_y ', ref_x,ref_y
  for i in range(len(text)):
    info = split(strip(text[i]))
    x_offset = ARCSEC * 8.0 * int((float(info[6]) +1.0 - ref_x) * 3600.0 / 8.0)
    y_offset = ARCSEC * 8.0 * int((float(info[7]) - ref_y) * 3600.0 / 8.0)
    beam = ASKAP_power_response(x_offset,y_offset)
    flux = float(info[15]) 
    I = flux * beam
#   if I >= 0.0015 and I < 0.02:
#   if I > 0.001 and I < 0.01:
#   if I >= 0.005 and I < 0.05:
#   if I >= 0.2 and I < 0.4:
# gives 49 sources at l0 = 0, m0 = 0 offset
# gives 32 sources at l0 = 1, m0 = 0 offset
# gives 40 sources at l0 = 2, m0 = 0 offset
#   if I >= 0.0075:
    if I >= 0.01:
      source_num = source_num + 1
      x_offset = ARCSEC * 8.0 * int((float(info[6])+1.0) * 3600.0 / 8.0)
      y_offset = ARCSEC * 8.0 * int((float(info[7])) * 3600.0 / 8.0)
      print 'iteration, signal ', source_num,x_offset, y_offset,I
      name = "%s%+d" % (basename,source_num)
      # adjust for CASA flux densities
      I = flux * 2
      model.append(point_source(ns,name,x_offset,y_offset,I));

  return model;

def source_count_pointing(ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources from Stil model""";
  from string import split, strip
  filename = 'Tony.Willis_1.dat'
  # note we 'reset' ref 0 RA to -1 deg so we have extra 'room' in + L direction
  text = open(filename, 'r').readlines()
  source_num = -1
  model = []
  ref_x = l0 / DEG
  ref_y = m0 / DEG
  print 'ref_x ref_y ', ref_x,ref_y
  for i in range(len(text)):
    info = split(strip(text[i]))
    x_offset = ARCSEC * 8.0 * int((float(info[6]) +1.0 - ref_x) * 3600.0 / 8.0)
    y_offset = ARCSEC * 8.0 * int((float(info[7]) - ref_y) * 3600.0 / 8.0)
    beam = ASKAP_power_response(x_offset,y_offset)
    flux = float(info[15]) 
    I = flux * beam
#   if I >= 0.0015 and I < 0.02:
#   if I > 0.001 and I < 0.01:
#   if I >= 0.005 and I < 0.05:
#   if I >= 0.2 and I < 0.4:
# gives 49 sources at l0 = 0, m0 = 0 offset
# gives 32 sources at l0 = 1, m0 = 0 offset
# gives 40 sources at l0 = 2, m0 = 0 offset
#   if I >= 0.0075:
    if I >= 0.01:
      source_num = source_num + 1
      print 'iteration, atten_signal orig_signal ', source_num,x_offset, y_offset,I,flux
      name = "%s%+d" % (basename,source_num)
      # adjust for CASA flux densities
      I = flux * 2
      model.append(point_source(ns,name,x_offset,y_offset,I));

  return model;

def source_count_poln(ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources from Stil model""";
  from string import split, strip
  filename = 'Tony.Willis_1.dat'
  # note we 'reset' ref 0 RA to -1 deg so we have extra 'room' in + L direction
  text = open(filename, 'r').readlines()
  source_num = -1
  model = []
  ref_x = l0 / DEG
  ref_y = m0 / DEG
  print 'ref_x ref_y ', ref_x,ref_y
  for i in range(len(text)):
    info = split(strip(text[i]))
    x_offset = ARCSEC * 8.0 * int((float(info[6]) +1.0 - ref_x) * 3600.0 / 8.0)
    y_offset = ARCSEC * 8.0 * int((float(info[7]) - ref_y) * 3600.0 / 8.0)
    beam = ASKAP_power_response(x_offset,y_offset)
    flux = float(info[15]) 
    I = flux * beam
#   if I >= 0.0015 and I < 0.02:
#   if I > 0.001 and I < 0.01:
#   if I >= 0.005 and I < 0.05:
#   if I >= 0.2 and I < 0.4:
    if I >= 0.020:
      source_num = source_num + 1
      x_offset = ARCSEC * 8.0 * int((float(info[6])+1.0) * 3600.0 / 8.0)
      y_offset = ARCSEC * 8.0 * int((float(info[7])) * 3600.0 / 8.0)
      print 'iteration, signal ', source_num,x_offset, y_offset,I
      name = "%s%+d" % (basename,source_num)
      # adjust for CASA flux densities
      I = flux * 2
      Q= float(info[16]) * 2
      U= float(info[17]) * 2
      V = 0.0
      model.append(point_source_poln(ns,name,x_offset,y_offset,I,Q,U,V));
  return model;

def point_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns single offset point source""";
  dl = 0.0 
  dm = 0.0
  I = 2.0
# I = 0.5 * 2
  model = [point_source(ns,basename+"+0+0",l0+dl,m0+dm,I)];
  return model;

def make_model (ns,basename="S",l0=5.0e-10,m0=0,I=1):
  """Creates and returns selected model""";
  return model_func(ns,basename,l0,m0,grid_step*FWHM,grid_step*FWHM,(grid_size-1)/2,I);


# model options
TDLCompileOption("model_func","Model type",[point_model,cross_model,grid_model,star8_model,lbar_model,mbar_model,random_model,source_count,source_count_pointing,source_count_poln]);
TDLCompileOption("grid_size","Number of sources in each direction",[1,3,5,7],more=int);
TDLCompileOption("grid_step","Grid step, in units of FWHM",[.1,0.25,.5,1],more=float);


