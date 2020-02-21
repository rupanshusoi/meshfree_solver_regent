import "regent"
require "point"
require "core"
local C = regentlib.c

terra pprint(a : double[4])
	C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])
end

task blocker(r : region(ispace(int1d), Point))
where reads writes(r) do
end

task main()
	var file = C.fopen("partGrid40K", "r")
	var size : int
	C.fscanf(file, "%d", &size)
	
	var globaldata = region(ispace(int1d, size + 1), Point)	
	var edges = region(ispace(int1d, 393993 + 1), Edge)
	
	var defprimal = getInitialPrimitive()
	
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
	
	var wallpts = 0
	var outerpts = 0
	var interiorpts = 0
	
	C.printf("Populating globaldata\n")	

	globaldata[0].localID = 0

	var edgecount = 0
	for count = 0, size do
		C.fscanf(file, "%lf %lf %d %d %d %d %lf %lf %d %lf %d", 
			&x, &y, &left, &right, &flag1, &flag2, &nx, &ny, 
			&qt_depth, &min_dist, &nbhs)
		
		var nbhs_arr : int[20]
		for i = 0, 20 do
			nbhs_arr[i] = 0
		end

		var temp : int
		for j = 0, nbhs do
			C.fscanf(file, "%d", &temp)
			nbhs_arr[j] = temp
			edges[edgecount].in_ptr = count + 1
			edges[edgecount].out_ptr = temp
			edgecount += 1
		end

		var p = Point {count + 1, x, y, left, right, flag1, flag2, nbhs, 
				nbhs_arr, nx, ny, defprimal, dummy_double, 
				dummy_double, dummy_double, dummy_double, 
				dummy_double, 0, 0, 0, 0, 0, dummy_int, 
				dummy_int, dummy_int, dummy_int, 0, min_dist, 
				dummy_double, dummy_double}

		globaldata[count + 1] = p
	end
	
	C.fclose(file)
	
	-- making partitions

	var points_equal = partition(equal, globaldata, ispace(int1d, 15))
	var edges_out = preimage(edges, points_equal, edges.in_ptr)
	var points_out = image(globaldata, edges_out, edges.out_ptr)
	var points_ghost = points_out - points_equal
	var points_allnbhs = points_equal | points_ghost
	
	var idx : int
	var curr : double[2]
	var leftpt : double[2]
	var rightpt : double[2]
	var normals : double[2]
	
	C.printf("Setting normals\n")
	for point in globaldata do
		if point.flag_1 == 0 or point.flag_1 == 2 then
			var curr : double[2]
			curr[0] = point.x
			curr[1] = point.y
			var leftpt : double[2]
			var rightpt : double[2]
			leftpt[0] = globaldata[point.left].x
			leftpt[1] = globaldata[point.left].y
			rightpt[0] = globaldata[point.right].x
			rightpt[1] = globaldata[point.right].y
			normals = calculateNormals(leftpt, rightpt, curr[0], curr[1])
			point.nx = normals[0]
			point.ny = normals[1]
		end
	end

	C.printf("Calculating connectivity\n")
	var connectivity : int[80]
	for count = 0, size do
		idx = count + 1
		connectivity = calculateConnectivity(globaldata, idx)
		setConnectivity(globaldata, idx, connectivity)
	end

	var res_old : double = 0.0
	var eu : int = 1
	var rks : int = 5
	var iter : int = 200

	C.printf("Starting FPI solver\n")
	
	-- refactoring to make FPI solver run from here

	for i = 1, iter do	
		__demand(__index_launch)
		for color in points_equal.colors do
			func_delta(points_equal[color], points_allnbhs[color])
		end
		for rk = 1, rks do
			__demand(__index_launch)
			for color in points_equal.colors do
				setq(points_equal[color])
			end

			__demand(__index_launch)
			for color in points_equal.colors do
				setdq(points_equal[color], points_allnbhs[color])
			end
			
			__demand(__index_launch)
			for color in points_equal.colors do
				cal_flux_residual(points_equal[color], points_allnbhs[color])
			end
			res_old = state_update(globaldata, i, rk, eu, res_old)
		end
	end
end
regentlib.start(main)
