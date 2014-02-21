from Timba.Meq import meq
from Timba import pynode
import numpy

class PyDiffAzimuth (pynode.PyNode):
  """calculated difference in azimuth between direction of sun and
  direction of antenna pointing / phase centre """

  def get_result (self,request,*children):
    try:
      if len(children) != 2:
        raise TypeError,"two child expected";
      res0 = children[0];
      res1 = children[1];
      vellsets = [];
      cells = res0.cells;
      solar_azimuth = res0.vellsets[0].value;
      shape = solar_azimuth.shape
# there's probably a numpy function for this
      for i in range(shape[0]):
        if solar_azimuth[i,0] < 0.0:
          solar_azimuth[i,0] = solar_azimuth[i,0] + (math.pi * 2.0)

      obs_azimuth = res1.vellsets[0].value;
      shape = obs_azimuth.shape
      for i in range(shape[0]):
        if obs_azimuth[i,0] < 0.0:
          obs_azimuth[i,0] = obs_azimuth[i,0] + (math.pi * 2.0)
        obs_azimuth[i,0] = obs_azimuth[i,0] - solar_azimuth[i,0]
      vellsets.append(meq.vellset(obs_azimuth))
      res = meq.result(cells = cells)
      res.vellsets = vellsets
      return res
    except:
      traceback.print_exc();
      raise;

