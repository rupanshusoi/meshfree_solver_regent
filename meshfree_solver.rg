import "regent"
require "config"
require "point"
require "core"

--[[
local MAPPER
do
  local root_dir = arg[0]:match(".*/") or "./"

  local include_path = ""
  local include_dirs = terralib.newlist()
  include_dirs:insert("-I")
  include_dirs:insert(root_dir)
  for path in string.gmatch(os.getenv("INCLUDE_PATH"), "[^;]+") do
    include_path = include_path .. " -I " .. path
    include_dirs:insert("-I")
    include_dirs:insert(path)
  end

  local mapper_cc = root_dir .. "meshfree_mapper.cc"
  local mapper_so
  if os.getenv('OBJNAME') then
    local out_dir = os.getenv('OBJNAME'):match('.*/') or './'
    mapper_so = out_dir .. "meshfree_mapper.so"
  elseif os.getenv('SAVEOBJ') == '1' then
    mapper_so = root_dir .. "meshfree_mapper.so"
  else
    mapper_so = os.tmpname() .. ".so" -- root_dir .. "circuit_mapper.so"
  end
  local cxx = os.getenv('CXX') or 'c++'

  local cxx_flags = os.getenv('CC_FLAGS') or ''
  cxx_flags = cxx_flags .. " -O2 -Wall -Werror"
  if os.execute('test "$(uname)" = Darwin') == 0 then
    cxx_flags =
      (cxx_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cxx_flags = cxx_flags .. " -shared -fPIC"
  end

  cxx_flags = cxx_flags .. " " .. "-std=c++11"
  local cmd = (cxx .. " " .. cxx_flags .. " " .. include_path .. " " ..
                 mapper_cc .. " -o " .. mapper_so)
  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. mapper_cc)
    assert(false)
  end
  terralib.linklibrary(mapper_so)
  MAPPER = terralib.includec("meshfree_mapper.h", include_dirs)
end
--]]

local C = regentlib.c
local sqrt = regentlib.sqrt(double)
local log10 = regentlib.log10(double)

terra pprint(a : double[4])
  C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])
end  

task run_setq(p : region(ispace(int1d), Point), part : partition(disjoint, p, ispace(int1d)))
where
  reads(p.{localID, prim}), writes(p.q)
do
  __demand(__index_launch)
  for color in part.colors do
    setq(part[color])
  end
end

task main()
  var file = C.fopen("grids/partGrid40K", "r")

  var size : int
  C.fscanf(file, "%d", &size)
  
  var config : Config
  config : initConfig(size)

  var globaldata = region(ispace(int1d, size + 1), Point)  
  var edges = region(ispace(int1d, config.totalnbhs + 1), Edge)
  
  var defprimal = getInitialPrimitive()
  
  var localID : int
  var part_number : int
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
  globaldata[0].part_number = 0

  var edgecount = 0
  for count = 0, size do
    if not config.isMETIS then
      C.fscanf(file, "%lf %lf %d %d %d %d %lf %lf %d %lf %d", &x, &y, &left, &right, &flag1, &flag2, &nx, &ny, &qt_depth, &min_dist, &nbhs)
    else
      C.fscanf(file, "%d %d %lf %lf %d %d %d %d %lf %lf %d %lf %d", &part_number, &localID, &x, &y, &left, &right, &flag1, &flag2, &nx, &ny, &qt_depth, &min_dist, &nbhs)
    end
    
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

    var p : Point
    if not config.isMETIS then
      p = Point {0, count + 1, x, y, left, right, flag1, flag2, nbhs, nbhs_arr, nx, ny, defprimal, dummy_double, dummy_double, dummy_double, dummy_double, dummy_double, dummy_double, dummy_double, 0, 0, 0, 0, 0, dummy_int, dummy_int, dummy_int, dummy_int, 0, min_dist, dummy_double, dummy_double}
    else
      p = Point {part_number, localID, x, y, left, right, flag1, flag2, nbhs, nbhs_arr, nx, ny, defprimal, dummy_double, dummy_double, dummy_double, dummy_double, dummy_double, dummy_double, dummy_double, 0, 0, 0, 0, 0, dummy_int, dummy_int, dummy_int, dummy_int, 0, min_dist, dummy_double, dummy_double}
    end

    globaldata[count + 1] = p
  end
  
  C.fclose(file)
  
  -- making partitions

  var points_equal = partition(equal, globaldata, ispace(int1d, config.partitions))
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
  var rks : int = config.rks
  var iter : int = config.iter

  var itime = C.legion_get_current_time_in_micros()
  C.printf("Starting FPI solver\n")
  
  __demand(__trace)
  for i = 1, iter + 1 do  
    __demand(__index_launch, __trace)
    for color in points_equal.colors do
      func_delta(points_equal[color], points_allnbhs[color],
           config)
    end
    for rk = 1, rks do

      run_setq(globaldata, points_equal)

      __demand(__index_launch, __trace)
      for color in points_equal.colors do
        setdq(points_equal[color], points_allnbhs[color],
              config)
      end

      for j = 0, config.inner_iter do
        __demand(__index_launch, __trace)
        for color in points_equal.colors do
          setqinner(points_equal[color], points_allnbhs[color], config)
        end

        __demand(__index_launch, __trace)
        for color in points_equal.colors do
          updateqinner(points_equal[color])
        end
      end

      __demand(__index_launch, __trace)
      for color in points_equal.colors do
        cal_flux_residual(points_equal[color], points_allnbhs[color], config)
      end

      var sum_res_sqr : double = 0.0

      __demand(__index_launch, __trace)
      for color in points_equal.colors do
        sum_res_sqr += state_update(points_equal[color], i, rk, eu, res_old)
      end

      var res_new : double = sqrt(sum_res_sqr) / config.size
      var residue : double
      if i <= 2 then
        res_old = res_new
        residue = 0
      else residue = log10(res_new / res_old) end

      if i % 100 == 0 then
        C.printf("Residue = %0.13lf for iteration %d, %d\n", residue, i, rk)
      end

    end
  end
  var ftime = C.legion_get_current_time_in_micros()
  regentlib.c.printf("***Time = %lld***\n", ftime - itime)
end
regentlib.start(main)
