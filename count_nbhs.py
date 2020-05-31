import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Grid File Location", type=str, default="partGridNew")
    args = parser.parse_args()

    data = open(args.file).read()
    data = data.split("\n")[1:]
    data.pop()
    totalnbh = 0
    for itm in data:
        itm = " ".join(itm.split(" ")).split(" ")
        totalnbh += int(itm[10])
    
    print(totalnbh)

if __name__ == "__main__":
    main()