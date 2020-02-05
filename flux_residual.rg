import "regent"
require "point"
require "outer_fluxes"
require "wall_fluxes"
require "interior_fluxes"

local C = regentlib.c

terra printArr(a : double[4])
	C.printf("[%lf, %lf, %lf, %lf]\n", a[0], a[1], a[2], a[3])
end

task cal_flux_residual(globaldata : region(ispace(int1d), Point), wallindices : region(ispace(int1d), int), outerindices : region(ispace(int1d), int), interiorindices : region(ispace(int1d), int))
where
	reads(wallindices, outerindices, interiorindices, globaldata), writes(globaldata)
do
	var itm : int
	for i = 0, 48738 do
		itm = wallindices[i]
		if itm == 0 then
			break
		else
			var Gxp = wall_dGx_pos(globaldata, itm)
			var Gxn = wall_dGx_neg(globaldata, itm)
			var Gyn = wall_dGy_neg(globaldata, itm)
			var GTemp : double[4]
			for j = 0, 4 do
				GTemp[j] = 2 * (Gxp[j] + Gxn[j] + Gyn[j])
			end
			globaldata[itm].flux_res = GTemp
		end
	end	

	for i = 0, 48738 do
		itm = outerindices[i]
		if itm == 0 then
			break
		else
			var Gxp = outer_dGx_pos(globaldata, itm)
			var Gxn = outer_dGx_neg(globaldata, itm)
			var Gyp = outer_dGy_pos(globaldata, itm)
			var GTemp : double[4]
			for j = 0, 4 do
				GTemp[j] = Gxp[j] + Gxn[j] + Gyp[j]	
			end
			globaldata[itm].flux_res = GTemp
		end
	end	

	for i = 0, 48738 do
		itm = interiorindices[i]
		if itm == 0 then
			break
		else
			var Gxp = interior_dGx_pos(globaldata, itm)
			var Gxn = interior_dGx_neg(globaldata, itm)
			var Gyp = interior_dGy_pos(globaldata, itm)
			var Gyn = interior_dGy_neg(globaldata, itm)
			var GTemp : double[4]
			if itm == 46053 then
				C.printf("printing Gxp etc for 46053\n")
				printArr(Gxp)
				printArr(Gxn)
				printArr(Gyp)
				printArr(Gyn)
			end
			for j = 0, 4 do
				GTemp[j] = Gxp[j] + Gxn[j] + Gyp[j] + Gyn[j]	
			end
			globaldata[itm].flux_res = GTemp
		end
	end
end
