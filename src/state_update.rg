import "regent"
require "point"
require "config"

local C = regentlib.c
local sqrt = regentlib.sqrt(double)
local exp = regentlib.exp(double)
local erf = regentlib.erf(double)
local sin = regentlib.sin(double)
local cos = regentlib.cos(double)
local PI = 3.1415926535897932

terra pprint(a : double[4])
  C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])        
end

__demand(__cuda)
task func_delta(pe : region(ispace(int1d), Point), 
    pn : region(ispace(int1d), Point), config : Config)
where
  reads(pn.{x, y, prim, conn}), writes(pe.{delta, prim_old})
do
  __demand(__openmp)
  for point in pe do 
    point.prim_old = pn[point].prim

    var conn = pn[point].conn
    var x = pn[point].x
    var y = pn[point].y

    var min_delt : double = 1.0

    for i = 0, 20 do
      var itm : int = conn[i]
      if itm == 0 then break
      else
        var rho = pn[itm].prim[0]
        var u1 = pn[itm].prim[1]
        var u2 = pn[itm].prim[2]
        var pr = pn[itm].prim[3]
      
        var x_k = pn[itm].x
        var y_k = pn[itm].y

        var dist = sqrt((x_k - x) * (x_k - x) + (y_k - y) * (y_k - y))

        var mod_u : double = sqrt(u1 * u1 + u2 * u2)
    
        var delta_t = config.cfl * (dist / (mod_u + 3 * sqrt(pr / rho)))

        if (min_delt > delta_t) then min_delt = delta_t end
      end
    end
    point.delta = min_delt
  end
end

__demand(__inline)
task primitive_to_conserved(nx : double, ny : double, prim : double[4])
  var U : double[4]

  var rho = prim[0]
  U[0] = rho

  var temp1 = rho * prim[1]
  var temp2 = rho * prim[2]

  U[1] = temp1 * ny - temp2 * nx

  U[2] = temp1 * nx + temp2 * ny

  U[3] = 2.5 * prim[3] + 0.5 * (temp1 * temp1 + temp2 * temp2) / rho

  return U
end

__demand(__inline)
task conserved_vector_Ubar(nx : double, ny : double, prim : double[4], config : Config)
  var mach = config.mach
  var gamma = config.gamma
  var pr_inf = config.pr_inf
  var rho_inf = config.rho_inf

  var Ubar : double[4]

  var theta = config.aoa * PI / 180.0

  var u1_inf = mach * cos(theta)
  var u2_inf = mach * sin(theta)

  var tx = ny
  var ty = -nx

  var u1_inf_rot = u1_inf * tx + u2_inf * ty
  var u2_inf_rot = u1_inf * nx + u2_inf * ny

  var temp1 = u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot
  var e_inf = pr_inf / (rho_inf * (gamma - 1)) + 0.5 * temp1

  var beta = (0.5 * rho_inf) / pr_inf
  var S2 = u2_inf_rot * sqrt(beta)
  var B2_inf = exp(-S2*S2) / (2 * sqrt(PI * beta))
  var A2n_inf = 0.5 * (1 - erf(S2))

  var rho = prim[0]
  var u1 = prim[1]
  var u2 = prim[2]
  var pr = prim[3]

  var u1_rot = u1*tx + u2*ty
  var u2_rot = u1*nx + u2*ny

  temp1 = u1_rot*u1_rot + u2_rot*u2_rot
  var e = pr / (rho * (gamma - 1)) + 0.5 * temp1

  beta = rho / (2 * pr)
  S2 = u2_rot * sqrt(beta)
  var B2 = exp(-S2*S2)/(2 * sqrt(PI * beta))
  var A2p = 0.5 * (1 + erf(S2))  

  Ubar[0] = (rho_inf * A2n_inf) + (rho * A2p)

  Ubar[1] = (rho_inf*u1_inf_rot*A2n_inf) + (rho*u1_rot*A2p)

  temp1 = rho_inf * (u2_inf_rot * A2n_inf - B2_inf)
  var temp2 = rho*(u2_rot*A2p + B2)

  Ubar[2] = temp1 + temp2

  temp1 = (rho_inf * A2n_inf * e_inf - 0.5 * rho_inf * u2_inf_rot * B2_inf)
  temp2 = rho * A2p * e + 0.5 * rho * u2_rot * B2

  Ubar[3] = temp1 + temp2

  return Ubar
end

__demand(__cuda)
task state_update(pe : region(ispace(int1d), Point), iter : int, rk : int, 
      eu : int, res_old : double, config: Config)
where
  reads(pe.{localID, flag_1, nx, ny, prim, prim_old, delta, flux_res}), writes(pe.prim)
do
  var sum_res_sqr : double = 0.0

  var obt : double = 1.0 / 3.0
  var tbt : double = 2.0 / 3.0
  
  var itm : int

  __demand(__openmp)
  for point in pe do
    if point.localID > 0 then
      if point.flag_1 == 0 then
        var nx = point.nx
        var ny = point.ny

        var U : double[4] = primitive_to_conserved(nx, ny, point.prim)
        var U_old : double[4] = primitive_to_conserved(nx, ny, point.prim_old)
        var temp = U[0]

        if rk ~= 3 then
          for j = 0, 4 do
            U[j] = U[j] - (0.5 * eu * point.flux_res[j])
          end
        else
          for j = 0, 4 do
            U[j] = (tbt * U_old[j]) + obt * (U[j] - (0.5 * point.flux_res[j]))
          end  
        end

        U[2] = 0  
        var U2_rot = U[1]
        var U3_rot = U[2]
        
        U[1] = U2_rot * ny + U3_rot * nx
        U[2] = U3_rot * ny - U2_rot * nx
        
        var res_sqr = (U[0] - temp) * (U[0] - temp)
        sum_res_sqr += res_sqr
        
        var tempU : double[4]
        tempU[0] = U[0]
        temp = 1 / U[0]
        tempU[1] = U[1] * temp
        tempU[2] = U[2] * temp
    
        tempU[3] = 0.4 * U[3] - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]))
        point.prim = tempU
      end
      if point.flag_1 == 1 then
        var nx = point.nx
        var ny = point.ny

        var U : double[4] = primitive_to_conserved(nx, ny, point.prim)
        var U_old : double[4] = primitive_to_conserved(nx, ny, point.prim_old)
        var temp = U[0]

        if rk ~= 3 then
          for j = 0, 4 do
            U[j] = U[j] - (0.5 * eu * point.flux_res[j])
          end
        else
          for j = 0, 4 do
            U[j] = (tbt * U_old[j]) + obt * (U[j] - (0.5 * point.flux_res[j]))
          end  
        end
        
        var U2_rot = U[1]
        var U3_rot = U[2]
        U[1] = U2_rot * ny + U3_rot * nx
        U[2] = U3_rot * ny - U2_rot * nx

        var res_sqr = (U[0] - temp) * (U[0] - temp)

        sum_res_sqr += res_sqr
        
        var tempU : double[4]
        tempU[0] = U[0]
        temp = 1 / U[0]
        tempU[1] = U[1] * temp
        tempU[2] = U[2] * temp
        tempU[3] = 0.4 * U[3] - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]))
        point.prim = tempU
      end
      if point.flag_1 == 2 then
        var nx = point.nx
        var ny = point.ny

        var U : double[4] = conserved_vector_Ubar(nx, ny, point.prim, config)
        var U_old : double[4] = conserved_vector_Ubar(nx, ny, point.prim_old, config)
        var temp = U[0]

        if rk ~= 3 then
          for j = 0, 4 do
            U[j] = U[j] - (0.5 * eu * point.flux_res[j])
          end
        else
          for j = 0, 4 do
            U[j] = (tbt * U_old[j]) + obt * (U[j] - (0.5 * point.flux_res[j]))
          end  
        end
        
        var U2_rot = U[1]
        var U3_rot = U[2]
        U[1] = U2_rot * ny + U3_rot * nx
        U[2] = U3_rot * ny - U2_rot * nx

        var tempU : double[4]
        tempU[0] = U[0]
        temp = 1 / U[0]
        tempU[1] = U[1] * temp
        tempU[2] = U[2] * temp

        tempU[3] = 0.4 * U[3] - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]))
        point.prim = tempU
      end
    end
  end

  return sum_res_sqr
end
