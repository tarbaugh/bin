import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 6))

exp = np.loadtxt("te_exper_gr.txt")
ax.plot(exp[:,0], exp[:,1], label="Experimental")

mtp = np.loadtxt("te_gr.txt")
ax.plot(mtp[:,0], mtp[:,1], label="MTP")

ax.legend()

plt.show()