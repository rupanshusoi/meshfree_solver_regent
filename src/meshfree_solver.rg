import "regent"
require "config"
require "point"
require "core"
require "fpi_solver"

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

local C = regentlib.c
local Cstring = terralib.includec("string.h")

terra pprint(a : double[4])
  C.printf("[\x1b[33m %0.15lf, %0.15lf, %0.15lf, %0.15lf]\n \x1b[0m", a[0], a[1], a[2], a[3])
end  

task read_grid(pt_distr : region(ispace(int1d), Point), edges : region(ispace(int1d), Edge), config : Config)
where writes(pt_distr, edges) do
  var defprimal = get_initial_primitive(config)
  
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
  
  var file = C.fopen([rawstring](config.filename), "r")
  var tmp : int
  C.fscanf(file, "%d", &tmp)

  --C.printf("Populating pt_distr\n")  
  pt_distr[0].localID = 0
  pt_distr[0].part_number = 0

  var edgecount = 0
  for count = 0, config.size do
    if not config.isMETIS then
      C.fscanf(file, "%lf %lf %d %d %d %d %lf %lf %d %lf %d", &x, &y, &left, &right, &flag1, &flag2, &nx, &ny, &qt_depth, &min_dist, &nbhs)
    else
      C.fscanf(file, "%d %d %lf %lf %d %d %d %d %lf %lf %d %lf %d", &part_number, &localID, &x, &y, &left, &right, &flag1, &flag2, &nx, &ny, &qt_depth, &min_dist, &nbhs)
    end
    
    -- Runge-Kutta point
    -- if count == 0 then flag1 = 1 end

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

    pt_distr[count + 1] = p
  end
  
  C.fclose(file)
end

task get_args()
  var args = C.legion_runtime_get_input_args()
  
  var iter, inner_iter = 200, 3

  var i = 1
  while i < args.argc do
    if Cstring.strcmp(args.argv[i], "--iter") == 0 then
      iter = C.atoi(args.argv[i+1])
    elseif Cstring.strcmp(args.argv[i], "--inner-iter") == 0 then
      inner_iter = C.atoi(args.argv[i+1])
    end
    i += 1
  end

  var config = initConfig(iter, inner_iter)
  return config
end

task output_file(pt_distr : region(ispace(int1d), Point), config : Config)
where reads(pt_distr) do

  var file = C.fopen("output.dat", "w")
  C.fprintf(file, "Size = %d Iter = %d Inner iter = %d\n", config.size, config.iter, config.inner_iter)
  for pt in pt_distr do
    C.fprintf(file, "%0.13lf %0.13lf %0.13lf %0.13lf %0.13lf %0.13lf\n", pt.x, pt.y, pt.prim[0], pt.prim[1], pt.prim[2], pt.prim[3])
    --C.fprintf(file, "%0.13lf %0.13lf %0.13lf %0.13lf %0.13lf %0.13lf\n", pt.x, pt.y, pt.flux_res[0], pt.flux_res[1], pt.flux_res[2], pt.flux_res[3])
  end

end

__demand(__replicable)
task main()
  var config = get_args()

  var pt_distr = region(ispace(int1d, config.size + 1), Point)  
  var edges = region(ispace(int1d, config.totalnbhs + 1), Edge)
  
  read_grid(pt_distr, edges, config)
  
  var idx : int
  var curr : double[2]
  var leftpt : double[2]
  var rightpt : double[2]
  var normals : double[2]

  --C.printf("Setting normals\n")
  for point in pt_distr do
    if point.flag_1 == 0 or point.flag_1 == 2 then
      var curr : double[2]
      curr[0] = point.x
      curr[1] = point.y
      var leftpt : double[2]
      var rightpt : double[2]
      leftpt[0] = pt_distr[point.left].x
      leftpt[1] = pt_distr[point.left].y
      rightpt[0] = pt_distr[point.right].x
      rightpt[1] = pt_distr[point.right].y
      normals = cal_normals(leftpt, rightpt, curr[0], curr[1])
      point.nx = normals[0]
      point.ny = normals[1]
    end
  end

  --C.printf("Calculating connectivity\n")
  var connectivity : int[80]
  for count = 0, config.size do
    idx = count + 1
    connectivity = cal_connectivity(pt_distr, idx)
    set_connectivity(pt_distr, idx, connectivity)
  end

  fpi_solver(pt_distr, edges, config)
  --output_file(pt_distr, config)
end
regentlib.start(main) -- , MAPPER.register_mappers)
