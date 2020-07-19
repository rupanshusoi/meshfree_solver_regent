import argparse
from functools import cmp_to_key

# This takes a Regent-METIS format file and sorts it by bitmask

def compare(a, b):
  a_bits = a[0].count(1)
  b_bits = b[0].count(1)
  
  if a_bits > b_bits:
    return 1
  elif a_bits < b_bits:
    return -1
  else:
    a_ = int(''.join(map(str, a[0])), 2)    
    b_ = int(''.join(map(str, b[0])), 2)    
    if a_ > b_:
      return 1
    elif a_ < b_:
      return -1
    else:
      a_part = int(a[1].split()[0])
      b_part = int(b[1].split()[0])
      if a_part > b_part:
        return 1
      elif a_part < b_part:
        return -1
      else:
        return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Regent-METIS format file", type=str)
    parser.add_argument("-b", "--bitmask", help="Bitmask file", type=str)
    args = parser.parse_args()
    grid_file_path = args.file
    bitmask_file_path = args.bitmask

    with open(grid_file_path, "r") as grid_file:
      data = grid_file.read().split("\n")
      data.pop()
      size = int(data.pop(0))

      arr = [i for i in range(size)]
      with open(bitmask_file_path, "r") as bitmask_file:
        bitmask = bitmask_file.read().split("\n")
        for line in bitmask[:-1]: 
          line = line.split()
          idx = int(line[0])
          bits = list(map(int, line[1:]))
          arr[idx - 1] = (bits, data[idx - 1]) 

        arr = sorted(arr, key=cmp_to_key(compare))
        order = [0 for i in range(size + 1)]
        for idx, t in enumerate(arr):
          order[int(t[1].split()[1])] = idx + 1

      fin = [i[1].split() for i in arr]    
     
      with open(grid_file_path + "_sb", "w+") as fin_file_path:
        fin_file_path.write(str(size) + "\n")
        for idx, line in enumerate(fin):
          line[1] = str(idx + 1)
          line[4] = str(order[int(line[4])])
          line[5] = str(order[int(line[5])])

          for index, i in enumerate(line[13:]):
            line[index + 13] = str(order[int(i)])

          fin_file_path.write(" ".join(line))
          fin_file_path.write("\n")
          

    print("Done")

if __name__ == "__main__":
    main()
