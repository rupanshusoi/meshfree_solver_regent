import "regent"
require "config"
require "point"

local pow = regentlib.pow(double)
local fabs = regentlib.fabs(double)

__demand(__inline)
task venkat_limiter(qtilde : double[4], pgp : region(ispace(int1d), Point), idx : int, config : Config)
where
  reads(pgp.{q, min_dist, minq, maxq})
do
  var VL_CONST = config.vl_const
  var phi : double[4]
  var count : int = 0
  for i = 0, 4 do
    var q = pgp[idx].q[i]
    var del_neg = qtilde[i] - q
    var del_pos : double
    if fabs(del_neg) <= 1e-5 then
      phi[count] = 1;
      count += 1
    else
      if del_neg > 0 then
        var max_q = pgp[idx].maxq[i]
        del_pos = max_q - q
      elseif del_neg < 0 then
        var min_q = pgp[idx].minq[i]
        del_pos = min_q - q
      end      
    
      var ds = pgp[idx].min_dist
      var epsi = VL_CONST * ds
      epsi = pow(epsi, 3)
    
      var num = (del_pos * del_pos) + (epsi * epsi)
      num = num * del_neg + 2.0 * del_neg * del_neg * del_pos

      var den = del_pos*del_pos + 2.0*del_neg*del_neg
      den = den + del_neg*del_pos + epsi*epsi
      den = den*del_neg

      var temp = num / den

      if temp < 1 then
        phi[count] = temp
        count += 1
      else
        phi[count] = 1
        count += 1
      end    
    end
  end
  return phi
end
