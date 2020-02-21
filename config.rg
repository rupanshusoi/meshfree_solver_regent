import "regent"

struct Config
{
	size : int,
	totalnbhs : int,
	partitions : int,
	cfl : double,
	max_iters : double,
	mach : double,
	aoa : int,
	power : int,
	vl_const : int,
	rho_inf : int,
	pr_inf : double,
	gamma : double
	rks : int,
	iter : int,
}

terra Config : initConfig(size : int)
	self.size = size
	self.totalnbhs = 393993
	self.partitions = 16
	self.cfl = 0.2
	self.max_iters = 20000
	self.mach = 0.85
	self.aoa = 1
	self.power = 0
	self.vl_const = 150
	self.rho_inf = 1
	self.pr_inf = 0.7142857142857143
	self.gamma = 1.4
	self.rks = 5
	self.iter = 200
end
