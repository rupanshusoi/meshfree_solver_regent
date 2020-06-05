import "regent"
require "point"

local C = regentlib.c

__demand(__inline) __forbid(__cuda)
task update_compact_instance(pe : region(ispace(int1d), Point), pg : region(ispace(int1d), Point), compact : region(ispace(int1d), Point)) 
where reads(pe, pg), writes(compact) do

  var num = compact.ispace.volume
  var pmap = region(ispace(int1d, num), int)
  for i = 0, num do pmap[i] = i end
  
  var count = 0
  for pt in pe do
    compact[count] = pe[pt]
    pmap[pt] = count
    count += 1
  end 
  for pt in pg do
    compact[count] = pg[pt]
    pmap[pt] = count
    count += 1
  end 

  regentlib.assert(num == count, "assertion failed!")
  return pmap
end

__demand(__inline) __forbid(__cuda)
task update_sparse_instance(pe : region(ispace(int1d), Point), pg : region(ispace(int1d), Point), compact : region(ispace(int1d), Point), pmap : region(ispace(int1d), int)) 
where reads(compact, pmap), writes(pe, pg) do

  for pt in pe do
    pe[pt] = compact[pmap[pt]]
  end
  for pt in pg do
    pg[pt] = compact[pmap[pt]]
  end
end
