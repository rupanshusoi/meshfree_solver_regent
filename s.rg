import "regent"

require "config"

task main()
	var conf : Config
	conf : initConfig()
	regentlib.c.printf("%d\n", conf.vl_const)
end
regentlib.start(main)
