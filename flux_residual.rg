import "regent"
require "point"
require "outer_fluxes"

local C = regentlib.c

task cal_flux_residual(globaldata : region(ispace(int1d), Point), wallindices : region(ispace(int1d), int), outerindices : region(ispace(int1d), int), interiorindices : region(ispace(int1d), int))
where
	reads writes(globaldata, outerindices) -- MINIMISE THESE PRIVILEGES
do
	var itm : int
	C.printf("***cal_flux_residual called***\n")
	for i = 0, 48738 do
		itm = outerindices[i]
		if itm == 0 then
			break
		else
			var Gxp = outer_dGx_pos(globaldata, itm)
			var Gxn = outer_dGx_neg(globaldata, itm)
			var Gyp = outer_dGy_pos(globaldata, itm)
			var GTemp = double[4]
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
			var GTemp = double[4]
			for j = 0, 4 do
				GTemp[j] = Gxp[j] + Gxn[j] + Gyp[j] + Gyn[j]	
			end
			globaldata[itm].flux_res = GTemp
		end
	end
end