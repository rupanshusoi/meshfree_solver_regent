import "regent"

__demand(__inline)
task foo()

	var sum_delx_sqr : double = 0
    	var sum_dely_sqr : double = 0
    	var sum_delx_dely : double = 0

	var nx = 1.0
	var ny = 1.0
	var tx = 1.0
	var ty = 1.0

	var delx = 1.0
	var dely = 1.0

	var dels = delx*tx + dely*ty
	var deln = delx*nx + dely*ny
	var dels_weights = dels
	var deln_weights = deln

	sum_delx_sqr = sum_delx_sqr + dels*dels_weights
	sum_dely_sqr = sum_dely_sqr + deln*deln_weights

	regentlib.c.printf("Sums %0.13lf %0.13lf\n", sum_delx_sqr, sum_dely_sqr)
end

task main()
	foo()
end
regentlib.start(main)
