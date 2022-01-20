import os
import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file, export_file
from ovito.modifiers import CreateBondsModifier, WrapPeriodicImagesModifier, TimeAveragingModifier, BondAnalysisModifier

def change_pbc(frame, data):
    data.cell_.pbc = (True,True,True)

def getBAD(filename, dir):
    # Load input data.
  pipeline = import_file(filename)

  pipeline.modifiers.append(change_pbc)
  pipeline.modifiers.append(WrapPeriodicImagesModifier())

  # Calculate partial RDFs:
  mod = CreateBondsModifier(mode = CreateBondsModifier.Mode.Uniform, cutoff=3.2)
  pipeline.modifiers.append(mod)

  pipeline.modifiers.append(BondAnalysisModifier(bins=80))

  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:bond-angle-distr'))
  # Access the output DataTable:
  bond_table = pipeline.compute().tables['bond-angle-distr[average]']
  for component, name in enumerate(bond_table.y.component_names):
    print(name, component)
  bond_table = bond_table.xy()

  return bond_table[:,0], bond_table[:,1]

def getFiles():
  for _,dirs,_ in os.walk(os.getcwd()):
    for dir in dirs:
      for root,_,files in os.walk(dir):
        for file in files:
            if file.endswith(".xyz"):
              print('Computing '+os.path.join(root, file)+'...')
              
liquid_x, liquid_y = getBAD('liquid.xyz', 'root')
amor_x, amor_y = getBAD('amor.xyz', 'root')
liq_max = liquid_y.max()
amor_max = amor_y.max()
amor_y *= 1.32

GAP_l = np.loadtxt("gap_liquid_ba.txt", dtype=float)
GAP_a = np.loadtxt("gap_amor_ba.txt", dtype=float)
GAP_l = GAP_l[GAP_l[:, 0].argsort()]
GAP_a = GAP_a[GAP_l[:, 0].argsort()]

_, gap_liq_max = GAP_l.max(axis=0)
_, gap_amor_max = GAP_a.max(axis=0)

liquid_y /= liq_max
amor_y /= amor_max
GAP_l[:, 1] /= gap_liq_max
GAP_a[:, 1] /= gap_amor_max

fig = plt.figure()
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)

axs[0].plot(liquid_x, liquid_y, label='MTP')
axs[0].plot(GAP_l[:,0],GAP_l[:,1], dashes=[6, 2], color='black', label='DFT')
# axs[0].plot(GeGeDFT[:,0],GeGeDFT[:,1])
axs[0].set_yticks(np.arange(0.0, 1.5, 0.5))
# axs[0].set_ylim(0, 2.5)
axs[0].annotate('Liquid 1200K', xy=(20, 0.5), xycoords="data")
axs[0].legend()



axs[1].plot(amor_x, amor_y, label='MTP')
axs[1].plot(GAP_a[:,0],GAP_a[:,1], dashes=[6, 2], color='black', label='DFT')
# axs[1].plot(GeTeDFT[:,0],GeTeDFT[:,1])
axs[1].set_yticks(np.arange(0.0, 0.7, 0.3))
# axs[1].set_ylim(0, 1)
axs[1].annotate('Amorphous 300K', xy=(20, 0.5), xycoords="data")

axs[1].set_xlabel(r"Angle ${\Theta}$ [deg.]")
axs[0].set_ylabel(r"P(${\Theta}$) [arb. units]")
axs[0].yaxis.set_label_coords(-.1, 0.0)

plt.show()




