import "regent"
local C = regentlib.c
terra printArr(a : double[4])
	C.printf("from child task %lf\n", a[0])
end
task main()
	var a = 0
	if a == 1 or a == 2 or a == 3 then
		C.printf("yes\n")
	else
		C.printf("NO\n")
	end
end
regentlib.start(main)
