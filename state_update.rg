import "regent"
require "point"

local C = regentlib.c
local Cmath = terralib.includec("math.h")

task func_delta(globaldata : region(ispace(int1d), Point), size : int)
where
	reads(globaldata.{x, y, prim, conn}), writes(globaldata.{delta, prim_old})
do
	var cfl : double = 0.2
	for idx = 1, size + 1 do
		var min_delt : double = 1
		globaldata[idx].prim_old = globaldata[idx].prim
		
		for i = 0, 20 do
			var itm : int = globaldata[idx].conn[i]
			if (itm == 0) then
				break
			else
				var rho = globaldata[itm].prim[0]
				var u1 = globaldata[itm].prim[1]
                		var u2 = globaldata[itm].prim[2]
                		var pr = globaldata[itm].prim[3]
			
				var x_i = globaldata[idx].x
                		var y_i = globaldata[idx].y

               			var x_k = globaldata[itm].x
               			var y_k = globaldata[itm].y

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
		globaldata[idx].delta = min_delt
		--C.printf("min_delt = %0.16lf\n", min_delt)
	end
end
