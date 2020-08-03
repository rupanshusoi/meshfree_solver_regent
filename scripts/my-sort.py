import os
import argparse

# This takes a Regent-METIS format file and sorts it by partition number

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Regent-METIS format file", type=str)
    args = parser.parse_args()
    grid_file = args.file

    with open(grid_file, "r") as grid:
      data = grid.read().split("\n")
      data.pop()
      size = int(data.pop(0))

      fin = sorted(data, key=lambda x:int(x.split()[0]))
      fin = [line.split() for line in fin]

      order = [i for i in range(size + 1)]
      for idx, line in enumerate(fin):
        order[int(line[1])] = idx + 1

      with open(grid_file + "_s", "w+") as fin_file:
        for idx, line in enumerate(fin):
          line[1] = str(idx + 1)
          line[4] = str(order[int(line[4])])
          line[5] = str(order[int(line[5])])

          for index, i in enumerate(line[13:]):
            line[index + 13] = str(order[int(i)])

          fin_file.write(" ".join(line))
          fin_file.write("\n")
          

    print("Done")

if __name__ == "__main__":
    main()
