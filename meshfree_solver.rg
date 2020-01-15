import "regent"
require "point"

local C = regentlib.c

task main()
	var file = C.fopen("partGrid40K", "r")
	
	var size : int
	C.fscanf(file, "%d", &size)
	
	var globaldata = region(ispace(int1d, size), Point)	
	
	var x : double
	var y : double
	var left : int
	var right : int
	var flag1 : int
	var flag2 : int
	var nx : double
	var ny : double
	var qt_depth : int
	var min_dist : double
	var nbhs : int
	var dummy_int : int[20]
	var dummy_double : double[4]
	var dummy_double2 : double[2][4]
	
	for count = 0, size do
		C.fscanf(file, "%lf %lf %d %d %d %d %lf %lf %d %lf %d", 
			&x, &y, &left, &right, &flag1, &flag2, &nx, &ny, 
			&qt_depth, &min_dist, &nbhs)
		
		var nbhs_arr : int[20]
		for i = 0, 20 do
			nbhs_arr[i] = 0
		end

		for j = 0, nbhs do
			var temp : int
			C.fscanf(file, "%d", &temp)
			nbhs_arr[j] = temp
		end

		var p = Point {0, x, y, left, right, flag1, flag2, 
			 nbhs, nbhs_arr, nx, ny, dummy_double, dummy_double,
			 dummy_double, dummy_double2, 0, 0, 0, 0, 0,
			 dummy_int, dummy_int, dummy_int, dummy_int, 0, min_dist, 
			 dummy_double, dummy_double}
		
		globaldata[count] = p
	end
	
	var p : Point
	
	var wallpts = 0
	var outerpts = 0
	var interiorpts = 0
	
	var wallptsidx : int[48738]
	var outerptsidx : int[48738]
	var interiorptsidx : int[48738]
	
	for i = 0, 48738 do
		wallptsidx[i] = 0
		outerptsidx[i] = 0
		interiorptsidx[i] = 0
	end	

	for i = [int](globaldata.bounds.lo), [int](globaldata.bounds.hi) do
		p = globaldata[i]
		if p.flag_1 == 0 then
			wallptsidx[wallpts] = i
			wallpts += 1
		elseif p.flag_1 == 1 then
			interiorptsidx[interiorpts] = i
			interiorpts += 1
		else
			outerptsidx[outerpts] = i
			outerpts += 1
		end
	end
end
regentlib.start(main)
