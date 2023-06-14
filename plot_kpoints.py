import numpy as np
import matplotlib.pyplot as plt

with open("E.txt", "r") as f:
    arr = []
    for i, l in enumerate(f):
        arr.append(list(map(float, l.split())))

arr = np.array(arr)
plt.plot(arr[:,0], arr[:,1]/200)

with open("train.cfg", "r") as f:
    arr = []
    for i, l in enumerate(f):
        if l.

for r in rs:
    plt.axvline(x=r)
plt.show()