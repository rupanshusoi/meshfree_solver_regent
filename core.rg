import "regent"
require "point"
require "state_update"
require "flux_residual"

local C = regentlib.c
local Cmath = terralib.includec("math.h")


terra getInitialPrimitive()
	var rho_inf : double = 1
	var mach : double = 0.85
	var machcos = mach * Cmath.cos(Cmath.M_PI/180) -- removed call to calculateTheta
	var machsin = mach * Cmath.sin(Cmath.M_PI/180) -- removed call to calculateTheta
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

	var det = Cmath.sqrt(nx*nx + ny*ny)
	
	nx = -nx / det
	ny = ny / det

	var arr : double[2]
	arr[0] = nx
	arr[1] = ny

	return arr
end

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
		
	var bigarr : int[80]	-- this array will contain all four of xpos_conn,
	for i = 0, 80 do	-- xneg_conn, ypos_conn and yneg_conn back to back
		bigarr[i] = 0	-- with each subarray of size 20
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

task q_var_derivatives(globaldata : region(ispace(int1d), Point), size : int)
where 
	reads(globaldata), writes(globaldata.{q, dq, minq, maxq})
do
	var power : double = 0
	for idx = 1, size + 1 do
		var itm : Point = globaldata[idx]
		var rho = itm.prim[0]
		var u1 = itm.prim[1]
		var u2 = itm.prim[2]
		var pr = itm.prim[3]
		
		var beta = 0.5 * (rho / pr)

		var tempq : double[4]

		tempq[0] = (Cmath.log(rho) + (Cmath.log(beta) * 2.5) - (beta * ((u1*u1) + (u2*u2))))

		tempq[1] = 2 * beta * u1
		tempq[2] = 2 * beta * u2
		tempq[3] = -2 * beta

		globaldata[idx].q = tempq	
	end
		
	for idx = 1, size + 1 do
		var itm = globaldata[idx]
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
					if minq[j] > globaldata[conn].q[j] then
						minq[j] = globaldata[conn].q[j]
					end
					if maxq[j] < globaldata[conn].q[j] then
						maxq[j] = globaldata[conn].q[j]
					end
				end
				
				var x_k = globaldata[conn].x
				var y_k = globaldata[conn].y

				var delx = x_k - x_i
				var dely = y_k - y_i

				var dist = Cmath.sqrt(delx*delx + dely*dely)
				var weights = Cmath.pow(dist, power)

				sum_delx_sqr = sum_delx_sqr + ((delx * delx) * weights)
				sum_dely_sqr = sum_dely_sqr + ((dely * dely) * weights)
				sum_delx_dely = sum_delx_dely + ((delx * dely) * weights)
				
				for j = 0, 4 do -- CHECK THIS
					sum_delx_delq[j] = sum_delx_delq[j] + (weights * delx * (globaldata[conn].q[j] - globaldata[idx].q[j]))
					sum_dely_delq[j] = sum_dely_delq[j] + (weights * dely * (globaldata[conn].q[j] - globaldata[idx].q[j]))
				end
			end	
		end
		
		var det : double = (sum_delx_sqr * sum_dely_sqr) - (sum_delx_dely * sum_delx_dely)
		var tempdq : double[2][4]

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
			sum_delx_delq2[i] = sum_delx_delq[i] * sum_delx_sqr
		end

		var tempsumy : double[4]
		for i = 0, 4 do
			tempsumy[i] = (1 / det) * (sum_dely_delq2[i] - sum_delx_delq2[i])
		end

		for i = 0, 4 do
			tempdq[0][i] = tempsumx[i]
		end
		for i = 0, 4 do
			tempdq[1][i]  = tempsumy[i]
		end

		globaldata[idx].dq = tempdq
			
		globaldata[idx].minq = minq
		globaldata[idx].maxq = maxq				

	end
end

task qtilde_to_primitive(qtilde : double[4]) -- MOVED TO OUTER FLUXES DUE TO
	var gamma : double = 1.4	     -- DUE TO CIRCULAR IMPORT ERROR	
	
	var q1 = qtilde[0]
    	var q2 = qtilde[1]
    	var q3 = qtilde[2]
	var q4 = qtilde[3]
	
	var beta = -q4 * 0.5
	var temp = 0.5 / beta
	
	var u1 = q2 * temp
	var u2 = q3 * temp

	var temp1 = q1 + beta * (u1*u1 + u2*u2)
	var temp2 = temp1 - (Cmath.log(beta)/(gamma - 1))
	var rho = Cmath.exp(temp2)
	var pr = rho * temp

	var arr : double[4]
	arr[0] = u1
	arr[1] = u2
	arr[2] = rho
	arr[3] = pr

	return arr
end

task fpi_solver(iter : int, globaldata : region(ispace(int1d), Point), wallindices : region(ispace(int1d), int), outerindices : region(ispace(int1d), int), interiorindices : region(ispace(int1d), int), res_old : int)
where
	reads (outerindices, interiorindices, wallindices, globaldata), writes(globaldata)
do
	var rks : int = 5
	var eu : int = 1

	for i = 1, iter do
		func_delta(globaldata, 48738)
		for rk = 1, rks do
			q_var_derivatives(globaldata, 48738)
			cal_flux_residual(globaldata, wallindices, outerindices, interiorindices)
			var res_old = state_update(globaldata, wallindices, outerindices, interiorindices, i, rk, eu, res_old)
			--todo: file writing here
			return res_old
		end
	end	

end
