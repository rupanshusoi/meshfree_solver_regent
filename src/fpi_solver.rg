import "regent"
require "point"
require "core"
require "config"

local C = regentlib.c
local sqrt = regentlib.sqrt(double)
local log10 = regentlib.log10(double)

task print_residue(residue : double, i : int, rk : int)
  C.printf("Residue = %0.13lf for iteration %d, %d\n", residue, i, rk)
end

task print_details(itime : int64, ftime : int64, config : Config)
  C.printf("*** Time = %lld ***\n", ftime - itime)
  C.printf("*** Filename = %s, Partitions = %d, Grid size = %lld points, Total Neighbours = %lld ***\n", [rawstring](config.filename), config.partitions, config.size, config.totalnbhs)
end

__demand(__inline, __openmp)
task fpi_solver(pt_distr : region(ispace(int1d), Point), edges : region(ispace(int1d), Edge), config : Config)
where reads writes(pt_distr, edges) do

  --var points_equal = partition(equal, pt_distr, ispace(int1d, config.partitions))
  var points_equal = partition(complete, pt_distr.part_number, ispace(int1d, config.partitions))
  var edges_out = preimage(edges, points_equal, edges.in_ptr)
  var points_out = image(pt_distr, edges_out, edges.out_ptr)
  var points_ghost = points_out - points_equal
  var points_allnbhs = points_equal | points_ghost

  var res_old : double = 0.0
  var eu : int = 1

  var itime = C.legion_get_current_time_in_micros()

  for i = 1, config.iter + 1 do
    __demand(__index_launch)
    for color in points_equal.colors do
      func_delta(points_equal[color], points_allnbhs[color],
           config)
    end
    for rk = 1, config.rks do
      __demand(__index_launch)
      for color in points_equal.colors do
        set_q(points_equal[color])
      end

      __demand(__index_launch)
      for color in points_equal.colors do
        set_dq(points_equal[color], points_allnbhs[color],
              config)
      end

      for j = 0, config.inner_iter do
        __demand(__index_launch)
        for color in points_equal.colors do
          set_q_inner(points_equal[color], points_allnbhs[color], config)
        end

        __demand(__index_launch)
        for color in points_equal.colors do
          update_q_inner(points_equal[color])
        end
      end

      __demand(__index_launch)
      for color in points_equal.colors do
        cal_flux_residual(points_equal[color], points_allnbhs[color], config)
      end

      var sum_res_sqr : double = 0.0

      __demand(__index_launch)
      for color in points_equal.colors do
        sum_res_sqr += state_update(points_equal[color], i, rk, eu, res_old, config)
      end

      var res_new : double = sqrt(sum_res_sqr) / config.size
      var residue : double
      if i <= 2 then
        res_old = res_new
        residue = 0
      else residue = log10(res_new / res_old) end

      if true then
        print_residue(residue, i, rk)
      end
    end
  end

  var ftime = C.legion_get_current_time_in_micros()
  print_details(itime, ftime, config)
end
