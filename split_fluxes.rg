import "regent"

local C = regentlib.c
local Cmath = terralib.includec("math.h")

__demand(__inline)
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

	Gxp[0] = rho * (ut * A1pos + B1)

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1pos + ut * B1

	Gxp[1] = rho * temp2

	temp1 = ut * un * A1pos + un * B1
	Gxp[2]= rho * temp1

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * ut * (7 * pr_by_rho + u_sqr)  * A1pos
	temp1 = 6 * pr_by_rho + u_sqr
	Gxp[3] = rho * (temp2 + 0.5 * temp1 * B1)

	return Gxp
end

__demand(__inline)
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

	Gxn[0] = rho * (ut * A1neg - B1)

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1neg - ut * B1

	Gxn[1] = rho * temp2

	temp1 = ut * un * A1neg - un * B1
	Gxn[2]= rho * temp1

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * ut * (7 * pr_by_rho + u_sqr)  * A1neg
	temp1 = 6 * pr_by_rho + u_sqr
	Gxn[3] = rho * (temp2 - 0.5 * temp1 * B1)

	return Gxn
end

__demand(__inline)
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

	Gyp[0] = rho * (un * A2pos + B2)

	var temp1 = pr_by_rho + un * un
	var temp2 = temp1 * A2pos + un * B2

	temp1 = ut * un * A2pos + ut * B2
	Gyp[1] = rho * temp1

	Gyp[2] = rho * temp2

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * un * (7 * pr_by_rho + u_sqr)  * A2pos
	temp1 = 6 * pr_by_rho + u_sqr
	Gyp[3] = rho * (temp2 + 0.5 * temp1 * B2)

	return Gyp
end

__demand(__inline)
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

	Gyn[0] = rho * (un * A2neg - B2)

	var temp1 = pr_by_rho + un * un
	var temp2 = temp1 * A2neg - un * B2

	temp1 = ut * un * A2neg - ut * B2
	Gyn[1] = rho * temp1

	Gyn[2]= rho * temp2

	temp1 = 7 * pr_by_rho + u_sqr
	temp2 = 0.5 * un * (7 * pr_by_rho + u_sqr)  * A2neg
	temp1 = 6 * pr_by_rho + u_sqr
	Gyn[3] = rho * (temp2 - 0.5 * temp1 * B2)

	return Gyn
end
