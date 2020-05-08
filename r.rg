import "regent"

local C = regentlib.c

struct Point
{
  color : int1d,
  num : int
}

task main()
  var data = region(ispace(int1d, 3), Point)
  var i = 0
  for pt in data do
     pt.color = 2 - i
     pt.num = 117 + i
     i += 1
  end
  var p = partition(data.color, ispace(int1d, 3))
  for color in p.colors do
    C.printf("color = %d ", color)
    for pt in p[color] do
      C.printf("%d ", pt.num)
      C.printf("%d ", pt)
    end
    C.printf("; ")
  end
  C.printf("%d\n", p[0][2].num)
end
regentlib.start(main)
