import numpy as np
from scipy import interpolate

exp = np.loadtxt("newest")
mintemp = np.amin(exp[:,0])
mindens = np.amin(exp[:,1])
print("minimum temp: {}".format(mintemp))
print("minimum dens: {}".format(mindens))


tck = interpolate.interp1d(exp[:,0], exp[:,1])

t = np.linspace(900, 600, num=31, endpoint=True)
d = t.copy()
d[d < mintemp] = mindens

for i in range(len(d)):
    if d[i] != mindens:
        d[i] = tck(d[i])

print("initial dens: {}".format(d[0]))

# import matplotlib.pyplot as plt
# plt.plot(t, d)
# plt.show()

with open("lammps.train", "w") as f:
    f.write("""units			metal
boundary		p p p
atom_style		atomic

read_data		otf_input
replicate       2 2 2

mass			1 72.63 
mass			2 127.603

pair_style		mlip mlip.ini
pair_coeff		* * mlip

thermo_style		custom step temp pe ke etotal press density

thermo			5000
timestep		0.001
velocity		all create 900 900 mom yes rot yes dist gaussian

fix			    997 all nvt temp 900 900 0.1
run			    1000000
dump			1 all custom 1000 md.xyz id type x y z fx fy fz
unfix           997\n\n""")

    for i in range(len(d)-1):
        f.write("fix			{} all nvt/sllod temp {} {} 0.1\n".format(2*i + 1, t[i], t[i+1]))
        f.write("fix			{} all deform 1 x scale {} y scale {} z scale {} remap v\n".format(2*i + 2, np.cbrt(d[i]/d[i+1]), np.cbrt(d[i]/d[i+1]), np.cbrt(d[i]/d[i+1])))
        f.write("run			200000\n")
        f.write("unfix			{}\n".format(2*i + 1))
        f.write("unfix			{}\n".format(2*i + 2))
        f.write("write_restart	{}.restart\n\n".format(int(t[i+1])))

    f.write("fix			998 all nvt/sllod temp {} {} 0.1\n".format(t[-1], t[-1]))
    f.write("fix			999 all deform 1 x scale {} y scale {} z scale {} remap v\n".format(np.cbrt(d[i]/d[i+1]), np.cbrt(d[i]/d[i+1]), np.cbrt(d[i]/d[i+1])))
    f.write("run			1000000\n")