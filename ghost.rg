import "regent"
require "point"
require "config"

local C = regentlib.c

__demand(__inline) --__forbid(__cuda)
task update_compact_instance(pe : region(ispace(int1d), Point), pg : region(ispace(int1d), Point), compact : region(ispace(int1d), Point), config : Config) 
where reads(pe, pg), writes(compact) do

  var pmap = region(ispace(int1d, config.size + 2), int)
  for i = 0, config.size + 1 do pmap[i] = i end
  
  var count = 1
  for pt in pe do
    if pt.localID ~= 0 then
      compact[count] = pe[pt]
      pmap[pt] = count
      count += 1
    end
  end 
  for pt in pg do
    if pt.localID ~= 0 then
      compact[count] = pg[pt]
      pmap[pt] = count
      count += 1
    end
  end 

  return pmap
end

__demand(__inline) --__forbid(__cuda)
task update_sparse_instance(pe : region(ispace(int1d), Point), pg : region(ispace(int1d), Point), compact : region(ispace(int1d), Point), pmap : region(ispace(int1d), int)) 
where reads(compact, pmap), writes(pe, pg) do

  for pt in pe do
    pe[pt] = compact[pmap[pt]]
  end
  for pt in pg do
    pg[pt] = compact[pmap[pt]]
  end
end
