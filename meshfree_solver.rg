import "regent"
require "arr"
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
	var dummy_int = region(ispace(int1d, 0), int_arr)
	var dummy_double = region(ispace(int1d, 0), double_arr)
	
	size = 5 -- REMOVE BEFORE RUNNING

	for count = 0, size do
		C.fscanf(file, "%lf %lf %d %d %d %d %lf %lf %d %lf %d", &x, &y, &left, &right, &flag1, &flag2, &nx, &ny, &qt_depth, &min_dist, &nbhs)
		var nbhs_arr = region(ispace(int1d, nbhs), int_arr)
		for j = 0, nbhs do
			var temp : int
			C.fscanf(file, "%d", &temp)
			nbhs_arr[j].n = temp
		end
		globaldata[count] = Point {0, x, y, left, right, flag1, flag2, nbhs, nbhs_arr, nx, ny, dummy_double, dummy_double, dummy_double, dummy_double, 0, 0, 0, 0, 0, dummy_int, dummy_int, dummy_int, dummy_int, 0, min_dist, dummy_double, dummy_double}
		C.printf("%lf %lf %lf\n", globaldata[count].x, globaldata[count].y, globaldata[count].min_dist)
	end
	
end
regentlib.start(main)
