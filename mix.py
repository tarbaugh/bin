from typing import get_type_hints
import numpy as np
import math

# ge = "1.93 2.181 1.8 31 1.2 109.5 7.049556277 0.6022245584 4 0 0"
# te = "1.03 2.51 1.80 25.0 1.20 90.0 8.1415 0.6671 4.0 0.0 0.0"

# ge = "1.918	2.181 1.8 21 1.2 109.5 7.049556277 0.6022245584 4 0 0"
# te = "1.03 2.51 1.80 25.0 1.20 90.0 8.1415 0.6671 4.0 0.0 0.0"

# ge = "1.918	2.181 1.8 21 1.2 90.0 7.049556277 0.6022245584 4 0 0"
# te = "1.03 2.51 1.80 25.0 1.20 90.0 8.1415 0.6671 4.0 0.0 0.0"

# ge = "1.918	2.181 1.8 21 1.2 109.5 7.049556277 0.6022245584 4 0 0"
# te = "1.849775   2.905254   1.594353   32.50000   1.200000  109.5   7.917000   .7307283   4.000000   0.000000   0.000000"

# ge = "2.1683  2.0951  1.80  21.0  1.20  109.5 7.049556277  0.6022245584  4.0  0.0 0.0" #Si
# te = "1.849775   2.905254   1.594353   32.50000   1.200000  109.5   8.1415   .7307283   4.000000   0.000000   0.000000"

ge = "1.2000000 2.100 1.6 32.5 1.2 109.5 7.917 0.72 4.0 0.0 0.0" # Ga
te = "1.03 2.51 1.80 25.0 1.20 109.5 8.1415 0.6671 4.0 0.0 0.0"

ge_params = ge.split()
te_params = te.split()

params = np.zeros((2,2,2,11)).astype(np.float64)

for i, p in enumerate(ge_params):
    params[0][0][0][i] = float(p)

for i, p in enumerate(te_params):
    params[1][1][1][i] = float(p)

types = ["Ge", "Te"]

def rule0(i):
    return math.sqrt(params[0][0][0][i]*params[1][1][1][i])

def rule1(i, j, k):
    return math.sqrt(math.sqrt(i*j*j*k))

def rule2(i):
    two = math.sqrt(params[0][0][0][i]*params[1][1][1][i])
    return 

for i in range(2):
    for j in range(2):
        for k in range(2):
            if i == j and j == k:
                s = types[i]+" "+types[j]+" "+types[k]+" "
                for p in range(11):
                    if p == 5:
                        s += str(-0.3333333)+" "
                    else:
                        s += str(params[i][j][k][p])+" "
                print(s)
            else:
                # params[i][j][k][0] = rule0(0)
                # params[i][j][k][3] = rule0(1)

                params[i][j][k][0] = rule1(params[i][i][i][0], params[j][j][j][0], params[k][k][k][0])
                params[i][j][k][3] = rule1(params[i][i][i][3], params[j][j][j][3], params[k][k][k][3])
                params[i][j][k][5] = rule1(params[i][i][i][5], params[j][j][j][5], params[k][k][k][5])
                params[i][j][k][6] = rule1(params[i][i][i][6], params[j][j][j][6], params[k][k][k][6])
                params[i][j][k][7] = rule1(params[i][i][i][7], params[j][j][j][7], params[k][k][k][7])
                
                params[i][j][k][1] = (params[i][i][i][1] + params[j][j][j][1])/2
                params[i][j][k][2] = params[0][0][0][2]
                params[i][j][k][4] = params[0][0][0][4]
                params[i][j][k][8] = params[0][0][0][8]
                params[i][j][k][9] = params[0][0][0][9]
                params[i][j][k][10] = params[0][0][0][10]
                s = types[i]+" "+types[j]+" "+types[k]+" "
                for p in range(11):
                    if p == 5:
                        s += str(round(math.cos(math.pi*params[i][j][k][p]/180), 5))+" "
                    else:
                        s += str(params[i][j][k][p])+" "
                print(s)

