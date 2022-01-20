import os
import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file, export_file
from ovito.modifiers import CreateBondsModifier, WrapPeriodicImagesModifier, TimeAveragingModifier, BondAnalysisModifier

filename = '1100KMSD.xyz'

def change_pbc(frame, data):
    data.cell_.pbc = (True,True,True)

def getBAD(filename, dir):
    # Load input data.
  pipeline = import_file(filename)

  pipeline.modifiers.append(change_pbc)
  pipeline.modifiers.append(WrapPeriodicImagesModifier())

  # Calculate partial RDFs:
  mod = CreateBondsModifier(mode = CreateBondsModifier.Mode.Pairwise)
  mod.set_pairwise_cutoff(1,1,3.0)
  mod.set_pairwise_cutoff(1,2,3.22)
  mod.set_pairwise_cutoff(2,2,3.0)
  pipeline.modifiers.append(mod)

  pipeline.modifiers.append(BondAnalysisModifier(partition=BondAnalysisModifier.Partition.ByParticleType))

  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:bond-angle-distr'))
  # Access the output DataTable:
  bond_table = pipeline.compute().tables['bond-angle-distr[average]']
  for component, name in enumerate(bond_table.y.component_names):
    print(name, component)
  bond_table = bond_table.xy()

  return bond_table[:,0], bond_table[:,1], bond_table[:,2], bond_table[:,3], bond_table[:,4], bond_table[:,5], bond_table[:,6]

def getFiles():
  for _,dirs,_ in os.walk(os.getcwd()):
    for dir in dirs:
      for root,_,files in os.walk(dir):
        for file in files:
            if file.endswith(".xyz"):
              print('Computing '+os.path.join(root, file)+'...')
              
rtx, rty1, rty2, rty3, rty4, rty5, rty6 = getBAD(filename, 'root')

sum1 = []
sum2 = []
sumTot = []

for i, _ in enumerate(rty1):
  sum1.append(rty1[i] + rty2[i] + rty3[i])
  sum2.append(rty4[i] + rty5[i] + rty6[i])
  sumTot.append(rty1[i] + rty2[i] + rty3[i] + rty4[i] + rty5[i] + rty6[i])

sum1 = np.array(sum1)
sum2 = np.array(sum2)
sumTot = np.array(sumTot)

sum1 /= 94
sum2 /= 80
sumTot /= 88.5

fig = plt.figure()
fig, axs = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(hspace=0)


# axs[0].plot(rtx, sum1)
axs[0].plot(rtx, sum1)
# axs[0].set_yticks(np.arange(0.0, 1.5, 0.5))
# axs[0].set_ylim(0, 2.5)
axs[0].annotate('X-Ge-Y', xy=(160, 0.8), xycoords="data")



axs[1].plot(rtx, sum2)
# axs[1].set_yticks(np.arange(0.0, 0.7, 0.3))
axs[1].annotate('X-Te-Y', xy=(160, 0.71), xycoords="data")


axs[2].plot(rtx, sumTot)
axs[2].set_yticks(np.arange(0.0, 1.5, 0.5))
# axs[2].set_ylim(0, 2)
axs[2].annotate('Total', xy=(160, 1.2), xycoords="data")

axs[2].set_xlabel(r"Angle ${\Theta}$ [deg.]")
axs[1].set_ylabel(r"P(${\Theta}$) [arb. units]")
plt.show()




