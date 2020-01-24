import "regent"

local C = regentlib.c

struct Point
{
	localID : int;
	x : double;
	y : double;
	left : int;
	right : int;
	flag_1 : int;
	flag_2 : int;
	nbhs : int;
	conn : int[20];
	nx : double;
	ny : double;
	prim : double[4];
	prim_old : double[4];
	flux_res : double[4];
	q : double[4];
	dq : double[2][4];
	entropy : double;
	xpos_nbhs : int;
	xneg_nbhs : int;
	ypos_nbhs : int;
	yneg_nbhs : int;
	xpos_conn : int[20];
	xneg_conn : int[20];
	ypos_conn : int[20];
	yneg_conn : int[20];
	delta : double;
	min_dist : double;
	minq : double[4];
	maxq : double[4];
}

terra Point : getxy()
	var a : double[2]
	a[0] = self.x
	a[1] = self.y
	return a
end

task setNormals(globaldata : region(ispace(int1d), Point), 
			idx : int, arr : double[2])
where 
	writes(globaldata.{nx, ny})
do
	globaldata[idx].nx = arr[0]
	globaldata[idx].ny = arr[1]
end

task setConnectivity(globaldata : region(ispace(int1d), Point), idx : int, bigarr : int[80])
where
	reads(globaldata.{localID, xpos_conn}), writes(globaldata.{xpos_conn, xneg_conn, ypos_conn, yneg_conn, xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs})
do

	var xpos : int[20]
	var xneg : int[20]
	var ypos : int[20]
	var yneg : int[20]
	
	var xpos_count = 0
	var xneg_count = 0
	var ypos_count = 0
	var yneg_count = 0

	for i = 0, 20 do
		xpos[i] = 0
		xneg[i] = 0
		ypos[i] = 0
		yneg[i] = 0
	end
	for i = 0, 20 do
		if bigarr[i] == 0 then
			break
		else
			xpos[xpos_count] = bigarr[i]
			xpos_count = xpos_count + 1
		end
	end	
	for i = 0, 20 do
		if bigarr[20 + i] == 0 then
			break
		else
			xneg[xneg_count] = bigarr[20 + i]
			xneg_count = xneg_count + 1 
		end
	end
	for i = 0, 20 do
		if bigarr[40 + i] == 0 then
			break
		else
			ypos[ypos_count] = bigarr[40 + i]
			ypos_count = ypos_count + 1
		end
	end
	for i = 0, 20 do
		if bigarr[60 + i] == 0 then
			break
		else
			yneg[yneg_count] = bigarr[60 + i]
			yneg_count = yneg_count + 1
		end
	end

	globaldata[idx].xpos_conn = xpos
	globaldata[idx].xneg_conn = xneg
	globaldata[idx].ypos_conn = ypos
	globaldata[idx].yneg_conn = yneg

	globaldata[idx].xpos_nbhs = xpos_count
	globaldata[idx].xneg_nbhs = xneg_count
	globaldata[idx].ypos_nbhs = ypos_count
	globaldata[idx].yneg_nbhs = yneg_count

end
