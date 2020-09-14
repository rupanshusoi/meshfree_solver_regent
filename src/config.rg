import "regent"

fspace Config
{
  filename : regentlib.string,
  size : int,
  totalnbhs : int,
  partitions : int,
  cfl : double,
  mach : double,
  aoa : double,
  power : double,
  vl_const : double,
  rho_inf : double,
  pr_inf : double,
  gamma : double,
  rks : int,
  iter : int,
  inner_iter : int,
  isMETIS : bool
}

task initConfig(iter : int, inner_iter : int)
  var grid_size = 40
  var size : int, totalnbhs : int

  if grid_size == 0.04 then
    size = 48738
    totalnbhs = 393993
  elseif grid_size == 0.8 then
    size = 804824
    totalnbhs = 6460047
  elseif grid_size == 2.5 then
    size = 2642264
    totalnbhs = 21172800
  elseif grid_size == 10 then
    size = 9992000
    totalnbhs = 79997009
  elseif grid_size == 25 then
    size = 25330172
    totalnbhs = 202730842
  elseif grid_size == 40 then
    size = 39381464
    totalnbhs = 315166328
  end

  var c = Config {
    filename = "../grids/partGrid40M_16_b",
    size = size,
    totalnbhs = totalnbhs,
    partitions = 16,
    cfl = 0.2,
    mach = 0.63,
    aoa = 2.0,
    power = 0.0,
    vl_const = 50.0,
    rho_inf = 1.0,
    pr_inf = 1.0 / 1.4,
    gamma = 1.4,
    rks = 5,
    iter = iter,
    inner_iter = inner_iter,
    isMETIS = true
  }
  return c
end
