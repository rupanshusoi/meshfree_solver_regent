import "regent"

local C = regentlib.c

task t()
	var c : int[1]
	c[0] = 1
	return c
end

task main()
	var a = t()
	a[0] = 22
	C.printf("%d\n", a[0])
end
regentlib.start(main)
