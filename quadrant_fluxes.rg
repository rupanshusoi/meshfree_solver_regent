import "regent"

local C = regentlib.c
local Cmath = terralib.includec("math.h")

terra printArr(a : double[4])
        C.printf("[%0.15lf, %0.15lf, %0.15lf, %0.15lf]\n", a[0], a[1], a[2], a[3])
end

task flux_quad_GxI(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)

	var G : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1*tx + u2*ty
	var un = u1*nx + u2*ny

	var beta = 0.5 * rho / pr
	var S1 = ut * Cmath.sqrt(beta)
	var S2 = un * Cmath.sqrt(beta)
	var B1 = 0.5 * Cmath.exp(-S1 * S1)/Cmath.sqrt(Cmath.M_PI * beta)
	var B2 = 0.5 * Cmath.exp(-S2 * S2)/Cmath.sqrt(Cmath.M_PI * beta)
	var A1neg = 0.5*(1 - Cmath.erf(S1))
	var A2neg = 0.5*(1 - Cmath.erf(S2))
	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un


	G[0] = rho * A2neg * (ut * A1neg - B1)

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1neg - ut * B1
	G[1] = rho * A2neg * temp2
	temp1 = ut * A1neg - B1
	temp2 = un * A2neg - B2
	G[2] = rho * temp1 * temp2

	temp1 = (7 * pr_by_rho) + u_sqr
	temp2 = 0.5 * ut * ((7 * pr_by_rho) + u_sqr) * A1neg
	temp1 = (6 * pr_by_rho) + u_sqr
	var temp3 = 0.5 * B1 * temp1

	temp1 = ut * A1neg - B1
	var temp4 = 0.5 * rho * un * B2 * temp1
	G[3] = rho * A2neg * (temp2 - temp3) - temp4

	return G

end

task flux_quad_GxII(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)

	var G : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1*tx + u2*ty
	var un = u1*nx + u2*ny

	var beta = 0.5 * rho / pr
	var S1 = ut * Cmath.sqrt(beta)
	var S2 = un * Cmath.sqrt(beta)
	var B1 = 0.5 * Cmath.exp(-S1 * S1)/Cmath.sqrt(Cmath.M_PI * beta)
	var B2 = 0.5 * Cmath.exp(-S2 * S2)/Cmath.sqrt(Cmath.M_PI * beta)
	var A1pos = 0.5*(1 + Cmath.erf(S1))
	var A2neg = 0.5*(1 - Cmath.erf(S2))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un


	G[0] = rho * A2neg * (ut * A1pos + B1)

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1pos + ut * B1
	G[1] = rho * A2neg * temp2

	temp1 = ut * A1pos + B1
	temp2 = un * A2neg - B2
	G[2] = rho * temp1 * temp2

	temp1 = (7 * pr_by_rho) + u_sqr
	temp2 = 0.5 * ut * ((7 * pr_by_rho) + u_sqr) * A1pos

	temp1 = (6 * pr_by_rho) + u_sqr
	var temp3 = 0.5 * B1 * temp1

	temp1 = ut * A1pos + B1
	var temp4 = 0.5 * rho * un * B2 * temp1

	G[3] = rho * A2neg * (temp2 + temp3) - temp4

	return G

end
task flux_quad_GxIII(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)

	var G : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1*tx + u2*ty
	var un = u1*nx + u2*ny

	var beta = 0.5 * rho / pr
	var S1 = ut * Cmath.sqrt(beta)
	var S2 = un * Cmath.sqrt(beta)
	var B1 = 0.5 * Cmath.exp(-S1 * S1)/Cmath.sqrt(Cmath.M_PI * beta)
	var B2 = 0.5 * Cmath.exp(-S2 * S2)/Cmath.sqrt(Cmath.M_PI * beta)
	var A1pos = 0.5*(1 + Cmath.erf(S1))
	var A2pos = 0.5*(1 + Cmath.erf(S2))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un


	G[0] = rho * A2pos * (ut * A1pos + B1)

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1pos + ut * B1
	G[1] = rho * A2pos * temp2

	temp1 = ut * A1pos + B1
	temp2 = un * A2pos + B2
	G[2] = rho * temp1 * temp2

	temp1 = (7 * pr_by_rho) + u_sqr
	temp2 = 0.5 * ut * ((7 * pr_by_rho) + u_sqr) * A1pos

	temp1 = (6 * pr_by_rho) + u_sqr
	var temp3 = 0.5 * B1 * temp1

	temp1 = ut * A1pos + B1
	var temp4 = 0.5 * rho * un * B2 * temp1

	G[3] = rho * A2pos * (temp2 + temp3) + temp4

	return G

end
	
task flux_quad_GxIV(nx : double, ny : double, u1 : double, u2 : double, rho : double, pr : double)

	var G : double[4]

	var tx = ny
	var ty = -nx

	var ut = u1*tx + u2*ty
	var un = u1*nx + u2*ny

	var beta = 0.5 * rho / pr
	var S1 = ut * Cmath.sqrt(beta)
	var S2 = un * Cmath.sqrt(beta)
	var B1 = 0.5 * Cmath.exp(-S1 * S1)/Cmath.sqrt(Cmath.M_PI * beta)
	var B2 = 0.5 * Cmath.exp(-S2 * S2)/Cmath.sqrt(Cmath.M_PI * beta)
	var A1neg = 0.5*(1 - Cmath.erf(S1))
	var A2pos = 0.5*(1 + Cmath.erf(S2))

	var pr_by_rho = pr / rho
	var u_sqr = ut * ut + un * un


	G[0] = rho * A2pos * (ut * A1neg - B1)

	var temp1 = pr_by_rho + ut * ut
	var temp2 = temp1 * A1neg - ut * B1
	G[1] = rho * A2pos * temp2

	temp1 = ut * A1neg - B1
	temp2 = un * A2pos + B2
	G[2] = rho * temp1 * temp2

	temp1 = ((7 * pr_by_rho) + u_sqr)
	temp2 = 0.5 * ut * ((7 * pr_by_rho) + u_sqr) * A1neg
	temp1 = (6 * pr_by_rho) + u_sqr
	var temp3 = 0.5 * B1 * temp1
	temp1 = ut * A1neg - B1
	var temp4 = 0.5 * rho * un * B2 * temp1
	G[3] = rho * A2pos * (temp2 - temp3) + temp4

	return G

end
