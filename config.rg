import "regent"

fspace Config
{
  filename : regentlib.string,
  size : int,
  totalnbhs : int,
  partitions : int,
  cfl : double,
  mach : double,
  aoa : int,
  power : int,
  vl_const : int,
  rho_inf : int,
  pr_inf : double,
  gamma : double,
  rks : int,
  iter : int,
  inner_iter : int,
  isMETIS : bool
}

--[[
800 K : 804824, 6460047
2.5 M : 2642264, 21172800
10 M  : 9992000, 79997009
25 M  : 25330172, 202730842
40 M  : 39381464, 315166328
--]]

task initConfig()
  var c = Config {
    filename = "grids/partGrid2.5M",
    size = 2642264,
    totalnbhs = 21172800,
    partitions = 16,
    cfl = 0.01,
    mach = 0.85,
    aoa = 1,
    power = 0,
    vl_const = 20,
    rho_inf = 1,
    pr_inf = 0.7142857142857143,
    gamma = 1.4,
    rks = 5,
    iter = 10,
    inner_iter = 0,
    isMETIS = false
  }
  return c
end
