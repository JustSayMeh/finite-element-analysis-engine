import math
x = []
y = []

def F(r, z):
    return -10*math.log(r) + 10 * math.log(100)

with open('axisX.txt', 'r') as f :
    for line in f.readlines():
        x.append(float(line))

with open('axisY.txt', 'r') as f :
    for line in f.readlines():
        y.append(float(line))

with open('nodes.txt', 'w') as f :
    f.write("{} 2 1\n".format(len(x)*len(y)))
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            f.write("{} {} 1.0\n".format(x[j], y[i]))

with open('elems.txt', 'w') as f :
    elems = []
    lenx = len(x)
    leny = len(y)
    print(lenx)
    for i in range(0, leny - 1):
        for j in range(0, lenx - 1):
            jbase = j + i * lenx
            elems.append([jbase, jbase + 1, jbase + lenx, jbase + lenx + 1])
    f.write("{} 4\n".format(len(elems)))
    for i in range(0, len(elems)):
        f.write("{} {} {} {}\n".format(elems[i][0], elems[i][1], elems[i][2], elems[i][3]))

secondL = []
with open('secondL.txt', 'r') as f :
    for line in f.readlines():
        secondL.append(float(line))


with open('firstB.txt', 'w') as f :
    elems = []
    lenx = len(x)
    leny = len(y)
    f.write("{}\n".format(2*len(x) + 2*(len(y) - 2) - len(secondL)))
    for i in range(0, lenx):
        f.write("{} {}\n".format(i, F(x[i], y[0])))
        f.write("{} {}\n".format(lenx*(leny - 1) + i, F(x[i], y[-1])))

    for i in range(1, leny - 1):
        if i not in secondL:
            f.write("{} {}\n".format(i*lenx, F(x[0], y[i])))
    for i in range(1, leny - 1):
        f.write("{} {}\n".format(i*lenx -1, F(x[-1], y[i])))


secondL.append(secondL[-1] + 1)
secondL.append(secondL[0] - 1)
for i in secondL:
    print("{}".format(i*lenx))