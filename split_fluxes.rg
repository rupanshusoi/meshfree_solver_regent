import "regent"

local C = regentlib.c
local Cmath = terralib.includec("math.h")

task flux_Gxp(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)
	var Gxp : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1 * tx + u2 * ty
	var un = u1 * nx + u2*ny

	var beta = 0.5 * rho / pr
	var S1 = ut * Cmath.sqrt(beta)
	var B1 = 0.5 * Cmath.exp(-S1 * S1) / Cmath.sqrt(Cmath.M_PI * beta)
	var A1pos = 0.5 * (1 + Cmath.erf(S1))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un

	var count : int = 0
	Gxp[count] = rho * (ut * A1pos + B1)
	count += 1	

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1pos + ut * B1

	Gxp[count] = rho * temp2
	count += 1

	temp1 = ut * un * A1pos + un * B1
	Gxp[count]= rho * temp1
	count += 1

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * ut * temp1 * A1pos
	temp1 = 6 * pr_by_rho + u_sqr
	Gxp[count] = rho * (temp2 + 0.5 * temp1 * B1)
	count += 1

	return Gxp
end

task flux_Gxn(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)
	var Gxn : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1 * tx + u2 * ty
	var un = u1 * nx + u2*ny

	var beta = 0.5 * rho / pr
	var S1 = ut * Cmath.sqrt(beta)
	var B1 = 0.5 * Cmath.exp(-S1 * S1) / Cmath.sqrt(Cmath.M_PI * beta)
	var A1neg = 0.5 * (1 - Cmath.erf(S1))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un

	var count : int = 0
	Gxn[count] = rho * (ut * A1neg - B1)
	count += 1	

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1neg - ut * B1

	Gxn[count] = rho * temp2
	count += 1

	temp1 = ut * un * A1neg - un * B1
	Gxn[count]= rho * temp1
	count += 1

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * ut * temp1 * A1neg
	temp1 = 6 * pr_by_rho + u_sqr
	Gxn[count] = rho * (temp2 - 0.5 * temp1 * B1)
	count += 1

	return Gxn
end

task flux_Gyp(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)
	var Gyp : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1 * tx + u2 * ty
	var un = u1 * nx + u2 * ny

	var beta = 0.5 * rho / pr
	var S2 = un * Cmath.sqrt(beta)
	var B2 = 0.5 * Cmath.exp(-S2 * S2) / Cmath.sqrt(Cmath.M_PI * beta)
	var A2pos = 0.5 * (1 + Cmath.erf(S2))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un

	var count : int = 0
	Gyp[count] = rho * (un * A2pos + B2)
	count += 1	

	var temp1 = pr_by_rho + un * un
	var temp2 = temp1 * A2pos + un * B2

	Gyp[count] = rho * temp2
	count += 1

	temp1 = ut * un * A2pos + ut * B2
	Gyp[count]= rho * temp1
	count += 1

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * un * temp1 * A2pos
	temp1 = 6 * pr_by_rho + u_sqr
	Gyp[count] = rho * (temp2 + 0.5 * temp1 * B2)
	count += 1

	return Gyp
end

task flux_Gyn(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)
	var Gyn : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1 * tx + u2 * ty
	var un = u1 * nx + u2 * ny

	var beta = 0.5 * rho / pr
	var S2 = un * Cmath.sqrt(beta)
	var B2 = 0.5 * Cmath.exp(-S2 * S2) / Cmath.sqrt(Cmath.M_PI * beta)
	var A2neg = 0.5 * (1 - Cmath.erf(S2))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un

	var count : int = 0
	Gyn[count] = rho * (un * A2neg - B2)
	count += 1	

	var temp1 = pr_by_rho + un * un
	var temp2 = temp1 * A2neg - un * B2

	temp1 = ut * un * A2neg - ut * B2
	Gyn[count] = rho * temp1
	count += 1

	Gyn[count]= rho * temp2
	count += 1

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * un * temp1 * A2neg
	temp1 = 6 * pr_by_rho + u_sqr
	Gyn[count] = rho * (temp2 - 0.5 * temp1 * B2)
	count += 1

	return Gyn
end
