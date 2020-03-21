import os
import argparse
import glob
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Grid File Folder", type=str, default="grids")
    args = parser.parse_args()
    gridPath = args.folder
    if not os.path.isdir(gridPath):
        print("Grid Folder does not exist")
        exit()

    gridPoints = {}

    gridFiles = glob.glob(os.path.join(gridPath, "*"))
    totalpoints = 0
    for gridFile in gridFiles:
        partValue = int(os.path.basename(gridFile).replace("partGrid",""))
        data = open(gridFile).read().split("\n")
        readLimit = 0
        for idx, itm in enumerate(data[:-1]):
            if idx == 0:
                readLimit = int(itm.split(" ")[2])
                totalpoints += readLimit
            if idx > 0 and idx <= readLimit:
                pointIdx = int(itm.split(" ")[0])
                gridPoints[pointIdx] = "{} {}".format(str(partValue),itm)

    assert len(gridPoints) == totalpoints, "Something went wrong"
    print("Writing output file")
    with open("partGrid{}".format(str(totalpoints)), "w+") as the_file:
        the_file.write("{}\n".format(str(totalpoints)))
        for i in range(1, totalpoints + 1):
                the_file.write("{}\n".format(gridPoints[i]))
    print("Done")

if __name__ == "__main__":
    main()
