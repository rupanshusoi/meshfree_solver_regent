import "regent"
require "point"

local C = regentlib.c
local Cmath = terralib.includec("math.h")

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
