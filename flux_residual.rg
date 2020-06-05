import "regent"
require "config"
require "point"
require "outer_fluxes"
require "wall_fluxes"
require "interior_fluxes"

local C = regentlib.c

terra pprint(a : double[4])
  C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])
end

__demand(__cuda)
task cal_flux_residual(compact : region(ispace(int1d), Point), pmap : region(ispace(int1d), int), config : Config)
where 
  reads(compact.{localID, flag_1, x, y, nx, ny, q, dq0, dq1, xpos_conn, xneg_conn,                      ypos_conn, yneg_conn, min_dist, minq, maxq}, pmap),
  writes(compact.flux_res)
do
  for point in compact do
    if point.flag_1 == 0 then
      var Gxp = wall_dGx_pos(compact, point.localID, pmap, config)
      var Gxn = wall_dGx_neg(compact, point.localID, pmap, config)
      var Gyn = wall_dGy_neg(compact, point.localID, pmap, config)
      var GTemp = array(2.0, 2.0, 2.0, 2.0) * (Gxp + Gxn + Gyn)
      point.flux_res = GTemp
    end
    if point.flag_1 == 1 then
      var Gxp = interior_dGx_pos(compact, point.localID, pmap, config)
      var Gxn = interior_dGx_neg(compact, point.localID, pmap, config)
      var Gyp = interior_dGy_pos(compact, point.localID, pmap, config)
      var Gyn = interior_dGy_neg(compact, point.localID, pmap, config)
      var GTemp = Gxp + Gxn + Gyp + Gyn
      point.flux_res = GTemp
    end
    if point.flag_1 == 2 then
      var Gxp = outer_dGx_pos(compact, point.localID, pmap, config)
      var Gxn = outer_dGx_neg(compact, point.localID, pmap, config)
      var Gyp = outer_dGy_pos(compact, point.localID, pmap, config)
      var GTemp = Gxp + Gxn + Gyp  
      point.flux_res = GTemp
    end
  end
end
