import "regent"
require "config"
require "point"
require "state_update"
require "flux_residual"

local C = regentlib.c
local cos = regentlib.cos(double)
local sin = regentlib.sin(double)
local sqrt = regentlib.sqrt(double)
local log = regentlib.log(double)
local pow = regentlib.pow(double)
local PI = 3.1415926535898

terra pprint(a : double[4])
  C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])
end

__demand(__inline)
task getInitialPrimitive()
  var rho_inf : double = 1.0
  var mach : double = 0.85
  var machcos = mach * cos(PI/180) -- removed call to calculateTheta
  var machsin = mach * sin(PI/180) -- removed call to calculateTheta
  var pr_inf = 0.7142857142857143
  var primal : double[4]
  primal[0] = rho_inf
  primal[1] = machcos
  primal[2] = machsin
  primal[3] = pr_inf
  return primal
end

task calculateNormals(left : double[2], right : double[2], mx : double, my : double)
  var lx = left[0]
  var ly = left[1]

  var rx = right[0]
  var ry = right[1]

  var nx1 = my - ly
  var nx2 = ry - my

  var ny1 = mx - lx
  var ny2 = rx - mx
  
  var nx = 0.5 * (nx1 + nx2)
  var ny = 0.5 * (ny1 + ny2)

  var det = sqrt(nx*nx + ny*ny)
  
  nx = -nx / det
  ny = ny / det

  var arr : double[2]
  arr[0] = nx
  arr[1] = ny

  return arr
end

__demand(__inline)
task calculateConnectivity(globaldata : region(ispace(int1d), Point), idx : int)
where
  reads writes(globaldata)
do
  var curr = globaldata[idx]
  var currx = curr.x
  var curry = curr.y
  var nx = curr.nx
  var ny = curr.ny
  var flag = curr.flag_1

  var tx = ny
  var ty = -nx
    
  var bigarr : int[80]  -- this array will contain all four of xpos_conn,
  for i = 0, 80 do  -- xneg_conn, ypos_conn and yneg_conn back to back
    bigarr[i] = 0  -- with each subarray of size 20
  end

  var xpos_count = 0
  var xneg_count = 0
  var ypos_count = 0
  var yneg_count = 0
  
  for i = 0, 20 do
    if curr.conn[i] == 0 then
      break
    else
      var itm = curr.conn[i]
      var itmx = globaldata[itm].x
      var itmy = globaldata[itm].y
      
      var delx = itmx - currx
      var dely = itmy - curry
      
      var dels = delx*tx + dely*ty
      var deln = delx*nx + dely*ny
      
      if dels <= 0 then
        bigarr[xpos_count] = itm
        xpos_count += 1
      end
      if dels >= 0 then
        bigarr[20 + xneg_count] = itm  
        xneg_count += 1
      end
      
      if flag == 1 then
        if deln <= 0 then
          bigarr[40 + ypos_count] = itm
          ypos_count += 1
        end
        if deln >= 0 then
          bigarr[60 + yneg_count] = itm
          yneg_count += 1
        end
      elseif flag == 0 then
        bigarr[60 + yneg_count] = itm
        yneg_count += 1
      elseif flag == 2 then
        bigarr[40 + ypos_count] = itm
        ypos_count += 1
      end
    end
  end   
  return bigarr
end

task setq(pe: region(ispace(int1d), Point))
where
  reads(pe.{localID, prim}), writes(pe.q)
do
  for itm in pe do
    if itm.localID > 0 then
      var rho = itm.prim[0]
      var u1 = itm.prim[1]
      var u2 = itm.prim[2]
      var pr = itm.prim[3]
      
      var beta = 0.5 * (rho / pr)

      var tempq : double[4]

      tempq[0] = (log(rho) + (log(beta) * 2.5) - (beta * ((u1*u1) + (u2*u2))))

      tempq[1] = 2 * beta * u1
      tempq[2] = 2 * beta * u2
      tempq[3] = -2 * beta

      itm.q = tempq  
    end
  end
end

task setdq(ppr : region(ispace(int1d), Point), 
     pgp : region(ispace(int1d), Point),
     config : Config)
where
  reads(ppr.{x, y, q, localID, conn}, pgp.{x, y, q}), 
  writes(ppr.{dq0, dq1, minq, maxq})
do
  var power = config.power

  for itm in ppr do
    if itm.localID > 0 then
      var x_i = itm.x
      var y_i = itm.y
      
      var sum_delx_sqr : double = 0
      var sum_dely_sqr : double = 0
      var sum_delx_dely : double = 0
    
      var sum_delx_delq : double[4]
      var sum_dely_delq : double[4]

      var minq : double[4]
      var maxq : double[4]

      for i = 0, 4 do
        sum_delx_delq[i] = 0
        sum_dely_delq[i] = 0
        minq[i] = 0
        maxq[i] = 0
      end

      for i = 0, 20 do
        if itm.conn[i] == 0 then
          break
        else
          var conn = itm.conn[i]
          for j = 0, 4 do
            if minq[j] > pgp[conn].q[j] then
              minq[j] = pgp[conn].q[j]
            end
            if maxq[j] < pgp[conn].q[j] then
              maxq[j] = pgp[conn].q[j]
            end
          end
          
          var x_k = pgp[conn].x
          var y_k = pgp[conn].y

          var delx = x_k - x_i
          var dely = y_k - y_i

          var dist : double = sqrt(delx*delx + dely*dely)
          var weights : double = pow(dist, power)

          sum_delx_sqr = sum_delx_sqr + ((delx * delx) * weights)
          sum_dely_sqr = sum_dely_sqr + ((dely * dely) * weights)
          sum_delx_dely = sum_delx_dely + ((delx * dely) * weights)
          
          for j = 0, 4 do
            sum_delx_delq[j] = sum_delx_delq[j] + (weights * delx * (pgp[conn].q[j] - pgp[itm].q[j]))
            sum_dely_delq[j] = sum_dely_delq[j] + (weights * dely * (pgp[conn].q[j] - pgp[itm].q[j]))
          end
        end  
      end

      var det : double = (sum_delx_sqr * sum_dely_sqr) - (sum_delx_dely * sum_delx_dely)
      var tempdq0 : double[4]
      var tempdq1 : double[4]

      var sum_delx_delq1 : double[4]
      var sum_dely_delq1 : double[4]

      for i = 0, 4 do
        sum_delx_delq1[i] = sum_delx_delq[i] * sum_dely_sqr
        sum_dely_delq1[i] = sum_dely_delq[i] * sum_delx_dely
      end

      var tempsumx : double[4] 
      for i = 0, 4 do
        tempsumx[i] = (1 / det) * (sum_delx_delq1[i] - sum_dely_delq1[i])
      end

      var sum_dely_delq2 : double[4]
      for i = 0, 4 do
        sum_dely_delq2[i] = sum_dely_delq[i] * sum_delx_sqr
      end
    
      var sum_delx_delq2 : double[4]
      for i = 0, 4 do
        sum_delx_delq2[i] = sum_delx_delq[i] * sum_delx_dely
      end

      var tempsumy : double[4]
      for i = 0, 4 do
        tempsumy[i] = (1 / det) * (sum_dely_delq2[i] - sum_delx_delq2[i])
      end

      tempdq0 = tempsumx
      tempdq1 = tempsumy

      itm.dq0 = tempdq0
      itm.dq1 = tempdq1
      itm.minq = minq
      itm.maxq = maxq        
    end
  end  

end

task setqinner(ppr : region(ispace(int1d), Point), 
     pgp : region(ispace(int1d), Point),
     config : Config)
  
where
  reads(ppr.{localID, x, y, q, dq0, dq1, conn}, pgp.{x, y, q, dq0, dq1}), 
  writes(ppr.{inner0, inner1})
do
  var power = config.power

  for itm in ppr do
    if itm.localID > 0 then
      var x_i = itm.x
      var y_i = itm.y
      
      var sum_delx_sqr : double = 0.0
      var sum_dely_sqr : double = 0.0
      var sum_delx_dely : double = 0.0

      var tmp1 : double = 0.0
      var tmp2 : double = 0.0

      var sum_delx_delq1 : double = 0.0
      var sum_delx_delq2 : double = 0.0
      var sum_delx_delq3 : double = 0.0
      var sum_delx_delq4 : double = 0.0

      var sum_dely_delq1 : double = 0.0
      var sum_dely_delq2 : double = 0.0
      var sum_dely_delq3 : double = 0.0
      var sum_dely_delq4 : double = 0.0

      var q1 : double = itm.q[0]
      var q2 : double = itm.q[1]
      var q3 : double = itm.q[2]
      var q4 : double = itm.q[3]

      for i = 0, 20 do
        if itm.conn[i] == 0 then
          break
        else
          var conn = itm.conn[i]
          var x_k = pgp[conn].x
          var y_k = pgp[conn].y

          var delx = x_k - x_i
          var dely = y_k - y_i

          var dist : double = sqrt(delx*delx + dely*dely)
          var weights : double = pow(dist, power)

          sum_delx_sqr = sum_delx_sqr + ((delx * delx) * weights)
          sum_dely_sqr = sum_dely_sqr + ((dely * dely) * weights)
          sum_delx_dely = sum_delx_dely + ((delx * dely) * weights)
          
          tmp1 = q1 - 0.5 * (delx * itm.dq0[0] + dely * itm.dq1[0])
          tmp2 = q1 - 0.5 * (delx * pgp[conn].dq0[0] + dely * pgp[conn].dq1[0])

          sum_delx_delq1 += weights * delx * (tmp2 - tmp1)
          sum_dely_delq1 += weights * dely * (tmp2 - tmp1)

          tmp1 = q2 - 0.5 * (delx * itm.dq0[1] + dely * itm.dq1[1])
          tmp2 = q2 - 0.5 * (delx * pgp[conn].dq0[1] + dely * pgp[conn].dq1[1])

          sum_delx_delq2 += weights * delx * (tmp2 - tmp1)
          sum_dely_delq2 += weights * dely * (tmp2 - tmp1)

          tmp1 = q3 - 0.5 * (delx * itm.dq0[2] + dely * itm.dq1[2])
          tmp2 = q3 - 0.5 * (delx * pgp[conn].dq0[2] + dely * pgp[conn].dq1[2])

          sum_delx_delq3 += weights * delx * (tmp2 - tmp1)
          sum_dely_delq3 += weights * dely * (tmp2 - tmp1)

          tmp1 = q4 - 0.5 * (delx * itm.dq0[3] + dely * itm.dq1[3])
          tmp2 = q4 - 0.5 * (delx * pgp[conn].dq0[3] + dely * pgp[conn].dq1[3])

          sum_delx_delq4 += weights * delx * (tmp2 - tmp1)
          sum_dely_delq4 += weights * dely * (tmp2 - tmp1)

        end
      end

      var det : double = (sum_delx_sqr * sum_dely_sqr) - (sum_delx_dely * sum_delx_dely)
      
      itm.inner0[0] = (1/det) * (sum_delx_delq1 * sum_dely_sqr - sum_dely_delq1 * sum_delx_dely)
      itm.inner0[1] = (1/det) * (sum_delx_delq2 * sum_dely_sqr - sum_dely_delq2 * sum_delx_dely)
      itm.inner0[2] = (1/det) * (sum_delx_delq3 * sum_dely_sqr - sum_dely_delq3 * sum_delx_dely)
      itm.inner0[3] = (1/det) * (sum_delx_delq4 * sum_dely_sqr - sum_dely_delq4 * sum_delx_dely)

      itm.inner1[0] = (1/det) * (sum_dely_delq1 * sum_delx_sqr - sum_delx_delq1 * sum_delx_dely)  
      itm.inner1[1] = (1/det) * (sum_dely_delq2 * sum_delx_sqr - sum_delx_delq2 * sum_delx_dely)  
      itm.inner1[2] = (1/det) * (sum_dely_delq3 * sum_delx_sqr - sum_delx_delq3 * sum_delx_dely)  
      itm.inner1[3] = (1/det) * (sum_dely_delq4 * sum_delx_sqr - sum_delx_delq4 * sum_delx_dely)  
    end
  end
end

task updateqinner(ppr : region(ispace(int1d), Point))
where
  reads(ppr.{inner0, inner1}), writes(ppr.{dq0, dq1})
do
  for point in ppr do
    point.dq0 = point.inner0
    point.dq1 = point.inner1
  end  
end
