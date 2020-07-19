import argparse
import os
import glob
from natsort import natsorted

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input Grid Folder", type=str, default="dummy")
    args = parser.parse_args()
    inputFile = args.input
    if not os.path.isdir(inputFile):
        print("Input file not found or invalid Input file")
        exit()
    files = natsorted(glob.glob(os.path.join(inputFile, "*")))
    globaldata = {}
    bitmask = {}
    for idx, itm in enumerate(files):
        data = open(itm).read().split("\n")
        data.pop(-1)
        smalldata = {"local":[], "ghost":[]}
        totalpts = int(data[0].split(" ")[1])
        localpts = int(data[0].split(" ")[2])
        ghostpts = int(data[0].split(" ")[3])
        for itm2 in data[1:localpts+1]:
            ptidx = int(itm2.split(" ")[0])
            smalldata["local"].append(ptidx)
            bitmask[ptidx] = [0] * len(files)
        for itm2 in data[localpts+1:]:
            ptidx = int(itm2.split(" ")[0])
            smalldata["ghost"].append(ptidx)
        assert localpts == len(smalldata["local"])
        assert ghostpts == len(smalldata["ghost"])
        globaldata[idx] = smalldata
    
    for itm in globaldata.keys():
        # for itm2 in globaldata[itm]["local"]:
        #     if itm2 not in bitmask.keys():
        #         bitmask[itm2] = [0] * len(globaldata)
        #         bitmask[itm2][itm] = 0
        
        for itm2 in globaldata[itm]["ghost"]:
            # if itm2 not in bitmask.keys():
            #     bitmask[itm2] = [0] * len(globaldata)
            #     bitmask[itm2][itm] = 1
            #     # print(itm2)
            # else:
            bitmask[itm2][itm] = 1

    with open("bitmask", "w+") as the_file:
        # for itm in range(1, len(bitmask.keys()) + 1):
        for itm in bitmask.keys():
            the_file.write("{} {}\n".format(itm, " ".join(map(str,bitmask[itm]))))
    print("Done")
        
        

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    main()
