import "regent"
require "arr"

local C = terralib.includec("stdio.h")

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
	conn : region(ispace(int1d, 20), int_arr); 
	nx : double;
	ny : double;
	prim : region(ispace(int1d, 4), double_arr);
	flux_res : region(ispace(int1d, 4), double_arr);
	q : region(ispace(int1d, 4), double_arr);
	dq : region(ispace(int2d, {2,4}), double_arr);
	entropy : double;
	xpos_nbhs : int;
	xneg_nbhs : int;
	ypos_nbhs : int;
	yneg_nbhs : int;
	xpos_conn : region(ispace(int1d, 20), int_arr);
	xneg_conn : region(ispace(int1d, 20), int_arr);
	ypos_conn : region(ispace(int1d, 20), int_arr);
	yneg_conn : region(ispace(int1d, 20), int_arr);
	delta : double;
	min_dist : double;
	minq : region(ispace(int1d, 4), double_arr);
	maxq : region(ispace(int1d, 4), double_arr);
}

terra Point : setNormals(x : double, y : double)
	self.nx = x
	self.ny = y
end

terra Point : getxy()
	return self.x
end

terra Point : gety()
	return self.y
end

terra Point : setConnectivity()
	
end
