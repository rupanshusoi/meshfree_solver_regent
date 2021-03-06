import "regent"

local C = regentlib.c

fspace Edge
{
  in_ptr : int1d;
  out_ptr : int1d;
}

struct Point
{
  part_number : int1d;
  localID : int;
  x : double;
  y : double;
  left : int;
  right : int;
  flag_1 : int;
  flag_2 : int;
  nbhs : int;
  conn : int[20];
  nx : double;
  ny : double;
  prim : double[4];
  prim_old : double[4];
  flux_res : double[4];
  q : double[4];
  dq0 : double[4];
  dq1 : double[4];
  inner0 : double[4];
  inner1 : double[4];
  entropy : double;
  xpos_nbhs : int;
  xneg_nbhs : int;
  ypos_nbhs : int;
  yneg_nbhs : int;
  xpos_conn : int[20];
  xneg_conn : int[20];
  ypos_conn : int[20];
  yneg_conn : int[20];
  delta : double;
  min_dist : double;
  minq : double[4];
  maxq : double[4];
}

task set_normals(pt_distr : region(ispace(int1d), Point), 
      idx : int, arr : double[2])
where 
  writes(pt_distr.{nx, ny})
do
  pt_distr[idx].nx = arr[0]
  pt_distr[idx].ny = arr[1]
end

__demand(__inline)
task set_connectivity(pt_distr : region(ispace(int1d), Point), idx : int, bigarr : int[80])
where
  writes(pt_distr.{xpos_conn, xneg_conn, ypos_conn, yneg_conn, 
         xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs})
do

  var xpos : int[20]
  var xneg : int[20]
  var ypos : int[20]
  var yneg : int[20]
  
  var xpos_count = 0
  var xneg_count = 0
  var ypos_count = 0
  var yneg_count = 0

  for i = 0, 20 do
    xpos[i] = 0
    xneg[i] = 0
    ypos[i] = 0
    yneg[i] = 0
  end
  for i = 0, 20 do
    if bigarr[i] == 0 then
      break
    else
      xpos[xpos_count] = bigarr[i]
      xpos_count = xpos_count + 1
    end
  end  
  for i = 0, 20 do
    if bigarr[20 + i] == 0 then
      break
    else
      xneg[xneg_count] = bigarr[20 + i]
      xneg_count = xneg_count + 1 
    end
  end
  for i = 0, 20 do
    if bigarr[40 + i] == 0 then
      break
    else
      ypos[ypos_count] = bigarr[40 + i]
      ypos_count = ypos_count + 1
    end
  end
  for i = 0, 20 do
    if bigarr[60 + i] == 0 then
      break
    else
      yneg[yneg_count] = bigarr[60 + i]
      yneg_count = yneg_count + 1
    end
  end

  pt_distr[idx].xpos_conn = xpos
  pt_distr[idx].xneg_conn = xneg
  pt_distr[idx].ypos_conn = ypos
  pt_distr[idx].yneg_conn = yneg

  pt_distr[idx].xpos_nbhs = xpos_count
  pt_distr[idx].xneg_nbhs = xneg_count
  pt_distr[idx].ypos_nbhs = ypos_count
  pt_distr[idx].yneg_nbhs = yneg_count

end
