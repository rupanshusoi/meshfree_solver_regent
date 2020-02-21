import "regent"
require "point"

local C = regentlib.c
local Cmath = terralib.includec("math.h")

terra pprint(a : double[4])
	C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])        
end

task func_delta(ppr : region(ispace(int1d), Point), 
		pgp : region(ispace(int1d), Point))
where
	reads(pgp.{x, y, prim, conn}, ppr.{x, y, prim, conn}), 
	writes(ppr.{delta, prim_old})
do
	var cfl : double = 0.2
	for point in ppr do
		var min_delt : double = 1
		point.prim_old = point.prim

		for i = 0, 20 do
			var itm : int = point.conn[i]
			if (itm == 0) then
				break
			else
				var rho = pgp[itm].prim[0]
				var u1 = pgp[itm].prim[1]
                		var u2 = pgp[itm].prim[2]
                		var pr = pgp[itm].prim[3]
			
				var x_i = point.x
                		var y_i = point.y

               			var x_k = pgp[itm].x
               			var y_k = pgp[itm].y

				var dist = (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i)
                		dist = Cmath.sqrt(dist)

				var mod_u : double = Cmath.sqrt(u1*u1 + u2*u2)
		
				var delta_t = dist/(mod_u + 3 * Cmath.sqrt(pr/rho))
				delta_t = cfl*delta_t

				if (min_delt > delta_t) then
					min_delt = delta_t
				end
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
task conserved_vector_Ubar(nx : double, ny : double, prim : double[4])
	var Mach : double = 0.85
	var gamma : double = 1.4
	var pr_inf : double = 0.7142857142857143
	var rho_inf : double = 1.0

	var Ubar : double[4]

	var theta = Cmath.M_PI / 180

	var u1_inf = Mach * Cmath.cos(theta)
	var u2_inf = Mach * Cmath.sin(theta)

	var tx = ny
	var ty = -nx

	var u1_inf_rot = u1_inf * tx + u2_inf * ty
	var u2_inf_rot = u1_inf * nx + u2_inf * ny

	var temp1 = u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot
	var e_inf = pr_inf / (rho_inf * (gamma - 1)) + 0.5 * temp1

	var beta = (0.5 * rho_inf) / pr_inf
	var S2 = u2_inf_rot * Cmath.sqrt(beta)
	var B2_inf = Cmath.exp(-S2*S2) / (2 * Cmath.sqrt(Cmath.M_PI * beta))
	var A2n_inf = 0.5 * (1 - Cmath.erf(S2))

	var rho = prim[0]
	var u1 = prim[1]
	var u2 = prim[2]
	var pr = prim[3]

	var u1_rot = u1*tx + u2*ty
	var u2_rot = u1*nx + u2*ny

	temp1 = u1_rot*u1_rot + u2_rot*u2_rot
	var e = pr / (rho * (gamma - 1)) + 0.5 * temp1

	beta = rho / (2 * pr)
	S2 = u2_rot * Cmath.sqrt(beta)
	var B2 = Cmath.exp(-S2*S2)/(2 * Cmath.sqrt(Cmath.M_PI * beta))
	var A2p = 0.5 * (1 + Cmath.erf(S2))	

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

task state_update(pe : region(ispace(int1d), Point), iter : int, rk : int, eu : int, res_old : double)
where
	reads(pe.{localID, flag_1, nx, ny, prim, prim_old, delta, flux_res}), writes(pe.prim)
do
	var max_res : double = 0.0
	var sum_res_sqr : double = 0.0

	var obt : double = 1.0 / 3.0
	var tbt : double = 2.0 / 3.0
	
	var itm : int

	C.printf("First Sum %lf\n", sum_res_sqr)

	for point in pe do
		if point.localID > 0 then
			if point.flag_1 == 0 then

				if point.localID == 1 then
					C.printf("Sum %lf\n", sum_res_sqr)
				end

				var nx = point.nx
				var ny = point.ny
				var Utemp : double[4] = primitive_to_conserved(nx, ny, point.prim)
					
				-- making some ugly changes as a workaround to an unfixed
				-- Regent bug

				var U : double[4]
				for i = 0, 4 do
					U[i] = Utemp[i]
				end

				var U_old : double[4] = primitive_to_conserved(nx, ny, point.prim_old)
				var temp = U[0]

				if rk == 1 or rk == 2 or rk == 4 then
					for j = 0, 4 do
						U[j] = U[j] - (0.5 * eu * point.delta * point.flux_res[j])
					end
				else
					for j = 0, 4 do
						U[j] = (tbt * U_old[j]) + obt * (U[j] - (0.5 * point.delta * point.flux_res[j]))
					end	
				end

				U[2] = 0	
				var U2_rot = U[1]
				var U3_rot = U[2]
				
				U[1] = U2_rot * ny + U3_rot * nx
				U[2] = U3_rot * ny - U2_rot * nx
				
				var res_sqr = (U[0] - temp) * (U[0] - temp)
	--[[
				if res_sqr > max_res then

					max_res = res_sqr
					var max_res_point = itm
				end
	--]]
				sum_res_sqr = sum_res_sqr + res_sqr
				
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
				var Utemp : double[4] = primitive_to_conserved(nx, ny, point.prim)

				-- making some ugly changes as a workaround to an unfixed
				-- Regent bug

				var U : double[4]
				for i = 0, 4 do
					U[i] = Utemp[i]
				end

				var U_old : double[4] = primitive_to_conserved(nx, ny, point.prim_old)
				var temp = U[0]

				if rk == 1 or rk == 2 or rk == 4 then
					for j = 0, 4 do
						U[j] = U[j] - (0.5 * eu * point.delta * point.flux_res[j])
					end
				else
					for j = 0, 4 do
						U[j] = (tbt * U_old[j]) + obt * (U[j] - (0.5 * point.delta * point.flux_res[j]))
					end	
				end
				
				var U2_rot = U[1]
				var U3_rot = U[2]
				U[1] = U2_rot * ny + U3_rot * nx
				U[2] = U3_rot * ny - U2_rot * nx

				var res_sqr = (U[0] - temp) * (U[0] - temp)

				if res_sqr > max_res then

					max_res = res_sqr
					var max_res_point = itm
				end

				sum_res_sqr = sum_res_sqr + res_sqr
				
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
				var Utemp : double[4] = conserved_vector_Ubar(nx, ny, point.prim)
				
				-- making some ugly changes as a workaround to an unfixed
				-- Regent bug

				var U: double[4]
				for i = 0, 4 do
					U[i] = Utemp[i]
				end

				var U_old : double[4] = conserved_vector_Ubar(nx, ny, point.prim_old)
				var temp = U[0]

				if rk == 1 or rk == 2 or rk == 4 then
					for j = 0, 4 do
						U[j] = U[j] - (0.5 * eu * point.delta * point.flux_res[j])
					end
				else
					for j = 0, 4 do
						U[j] = (tbt * U_old[j]) + obt * (U[j] - (0.5 * point.delta * point.flux_res[j]))
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

	var res_new : double = Cmath.sqrt(sum_res_sqr) / 48738.0
	var residue : double	
	if iter <= 2 then
		res_old = res_new
		residue = 0
	else
		residue = Cmath.log10(res_new / res_old)
	end
	
	-- todo : put file writing here
	C.printf("\x1b[32m\nIteration number: %d, %d\n", iter, rk)
	C.printf("Residue: %0.15lf\n \x1b[0m", residue)

	return res_old
end
