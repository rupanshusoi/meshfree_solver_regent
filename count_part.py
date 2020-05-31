arr = [0, 0, 0, 0, 0, 0, 0, 0]
with open("grids/partGrid48738", "r") as file:
    data = file.read()
    data = data.split('\n')
    data.pop()
    data.pop(0)
    for line in data:
        for i in range(8):
            if int(line[0]) == i:
                arr[i] += 1
    print(arr)
    print(sum(arr))
