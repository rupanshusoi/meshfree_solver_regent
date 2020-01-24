import "regent"
require "point"
require "core"
local C = regentlib.c

task main()
	var file = C.fopen("partGrid40K", "r")
	var size : int
	C.fscanf(file, "%d", &size)
	
	var globaldata = region(ispace(int1d, size + 1), Point)	
	
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
	var dummy_double2 : double[2][4]
	
	var wallpts = 0
	var outerpts = 0
	var interiorpts = 0
	
	var wallptsidx = region(ispace(int1d, 48738), int)
	var outerptsidx = region(ispace(int1d, 48738), int)
	var interiorptsidx = region(ispace(int1d, 48738), int)

	fill(wallptsidx, 0)
	fill(outerptsidx, 0)
	fill(interiorptsidx, 0)

	C.printf("Populating globaldata\n")	
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
		var p = Point {count + 1, x, y, left, right, flag1, flag2, nbhs, 
				nbhs_arr, nx, ny, defprimal, dummy_double, 
				dummy_double, dummy_double, dummy_double2, 0, 0, 
				0, 0, 0, dummy_int, dummy_int, dummy_int, dummy_int, 
				0, min_dist, dummy_double, dummy_double}

		globaldata[count + 1] = p

		if p.flag_1 == 0 then
			wallptsidx[wallpts] = count + 1
			wallpts += 1
		elseif p.flag_1 == 1 then
			interiorptsidx[interiorpts] = count + 1
			interiorpts += 1
		else
			outerptsidx[outerpts] = count + 1
			outerpts += 1
		end

	end
	
	var idx : int
	var curr : double[2]
	var leftpt : double[2]
	var rightpt : double[2]
	var normals : double[2]

	C.printf("Setting normals\n")
	for count = 0, 48738 do
		idx = count 
		idx = wallptsidx[idx]
		if idx == 0 then
			break
		else
			curr = globaldata[idx]:getxy()
			leftpt = globaldata[globaldata[idx].left]:getxy()
			rightpt = globaldata[globaldata[idx].right]:getxy()
			normals = calculateNormals(leftpt, rightpt, curr[0], curr[1])
			setNormals(globaldata, idx, normals)
		end
	end
	for count = 0, 48738 do
		idx = count
		idx = outerptsidx[idx]
		if idx == 0 then
			break
		else
			curr = globaldata[idx]:getxy()
			leftpt = globaldata[globaldata[idx].left]:getxy()
			rightpt = globaldata[globaldata[idx].right]:getxy()
			normals = calculateNormals(leftpt, rightpt, curr[0], curr[1])
			setNormals(globaldata, idx, normals)
		end
	end

	C.printf("Calculating connectivity\n")
	var connectivity : int[80]
	for count = 0, size do
		idx = count + 1
		connectivity = calculateConnectivity(globaldata, idx)
		setConnectivity(globaldata, idx, connectivity)
	end

	var res_old = 0

	C.printf("Starting FPI solver\n")
	fpi_solver(1 + 1, globaldata, wallptsidx, outerptsidx, interiorptsidx, 0)

end
regentlib.start(main)
