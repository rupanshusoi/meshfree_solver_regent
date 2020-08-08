# Gets total runtime for each task in q-LSKUM from Legion Prof generated tsv files

import os

def main():
  data = []
  num_files = 0
  for file in os.listdir("."):
    # OpenMP procs end with 3
    if file.endswith("3.tsv"):
      with open(file, "r") as the_file:
        num_files += 1
        data += the_file.read().split("\n")[1:-1]
        break

  times = {"func_delta" : 0, "setq" : 0, "setdq" : 0, "setqinner" : 0, "updateqinner" : 0, "cal_flux_residual" : 0, "state_update" : 0}
  
  exec_start = 0
  exec_end = 0

  for idx, line in enumerate(data):
    line = line.split("\t")

    start = float(line[3])
    end = float(line[4])
    task = line[7].split()[0]

    if idx == 0:
      exec_start = start
      exec_end = end
    else:
      exec_start = min(start, exec_start)
      exec_end = max(end, exec_end)

    times[task] += end - start

  print("Total execution time = {} s\n".format((exec_end - exec_start) / 1e6))

  for key in times:
    times[key] /= (exec_end - exec_start) * num_files
    print("{} : {:.2f} %".format(key, 100 * times[key]))


if __name__ == "__main__":
  main()
