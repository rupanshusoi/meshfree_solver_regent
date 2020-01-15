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
