import "regent"

struct Config
{
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
  gamma : double
  rks : int,
  iter : int,
  inner_iter : int,
  isMETIS : bool
}

terra Config : initConfig(size : int)
  self.size = size
  self.totalnbhs = 393993
  self.partitions = 2
  self.cfl = 0.01
  self.mach = 0.85
  self.aoa = 1
  self.power = 0
  self.vl_const = 20
  self.rho_inf = 1
  self.pr_inf = 0.7142857142857143
  self.gamma = 1.4
  self.rks = 5
  self.iter = 1
  self.inner_iter = 0
  self.isMETIS = false
end
