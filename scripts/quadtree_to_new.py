import math
def calculateNormals(left, right, mx, my, flag):
    lx = left[0]
    ly = left[1]

    rx = right[0]
    ry = right[1]

    nx1 = my - ly
    nx2 = ry - my

    ny1 = mx - lx
    ny2 = rx - mx

    nx = 0.5*(nx1 + nx2)
    ny = 0.5*(ny1 + ny2)

    det = math.sqrt(nx*nx + ny*ny)

    nx = -nx/det
    ny = ny/det

    if flag == 1:
    	return (1, 0)
    return (nx,ny)

file = open("partGrid9600", "r")
data = file.read()
data = data.split('\n')
data.pop(-1)
file.close()

with open("partGrid10K", "w") as file:
	for line in data[1:]:
		itm = line.split()
		lx = float(data[int(itm[2])].split()[0])
		ly = float(data[int(itm[2])].split()[1])
		rx = float(data[int(itm[3])].split()[0])
		ry = float(data[int(itm[3])].split()[1])
		flag = int(itm[4])
		nx, ny = calculateNormals((lx, ly), (rx, ry), float(itm[0]), float(itm[1]), flag)
		for i in range(0, 6):
			file.write(str(itm[i]) + ' ')
		file.write("{0:.20f}".format(nx) + ' ')
		file.write("{0:.20f}".format(ny) + ' ')
		file.write('0' + ' ') # quadtree depth
		file.write(str(itm[6]) + ' ')
		for i in range(7, len(itm) - 1):
			file.write(str(itm[i]) + ' ')
		file.write(itm[len(itm) - 1])	
		file.write('\n')
