
from Timba.Meq import meq
from Timba import pynode

class PyRandomConst (pynode.PyNode):
  """generate a random number that holds for the life of a tile"""

  def get_result (self,request,*children):
    try:
      if len(children) != 3:
        raise TypeError,"two children expected";
      res0 = children[0];
      vellsets = [];
      cells = res0.cells;
      result_vector = res0.vellsets[0].value;
      shape = result_vector.shape

      res1 = children[1];
      min_random = res1.vellsets[0].value
      res2 = children[2];
      max_random = res2.vellsets[0].value
      noise = random.uniform(min_random, max_random)
      for i in range(shape[0]):
          result_vector[i,0] = noise
      vellsets.append(meq.vellset(result_vector))
      res = meq.result(cells = cells)
      res.vellsets = vellsets
      return res
    except:
      traceback.print_exc();
      raise;

