from Timba.Meq import meq
from Timba import pynode
import numpy

class PyPosElevation (pynode.PyNode):
  """sets elevation to zero if it actually goes negative """
  def get_result (self,request,*children):
    try:
      if len(children) != 1:
        raise TypeError,"one child expected";
      res0 = children[0];
      vellsets = [];
      cells = res0.cells;
      elevation = res0.vellsets[1].value;
      shape = elevation.shape
# there's probably a numpy function for this
      for i in range(shape[0]):
        if elevation[i,0] < 0.0:
          elevation[i,0] = 0.0
      vellsets.append(meq.vellset(elevation))
      res = meq.result(cells = cells)
      res.vellsets = vellsets
      return res
    except:
      traceback.print_exc();
      raise;


