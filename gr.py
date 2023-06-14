import os
import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, WrapPeriodicImagesModifier, TimeAveragingModifier

def change_pbc(frame, data):
    data.cell_.pbc = (True,True,True)

def getRDF(filename, dir):
    # Load input data.
  pipeline = import_file(filename)

  pipeline.modifiers.append(change_pbc)
  pipeline.modifiers.append(WrapPeriodicImagesModifier())

  # Calculate partial RDFs:
  pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=9.0, number_of_bins=100))

  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
  # Access the output DataTable:
  rdf_table = pipeline.compute().tables['coordination-rdf[average]'].xy()
  return rdf_table[:,0], rdf_table[:,1]

def getFiles():
  for _,dirs,_ in os.walk(os.getcwd()):
    for dir in dirs:
      for root,_,files in os.walk(dir):
        for file in files:
            if file.endswith(".xyz"):
              print('Computing '+os.path.join(root, file)+'...')
              
liquid_x, liquid_y = getRDF('liquid.xyz', 'root')
amor_x, amor_y = getRDF('amor.xyz', 'root')

# fig, axs = plt.subplots()
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)

GAP_l = np.loadtxt("gap_liquid_gr.txt", dtype=float)
GAP_a = np.loadtxt("gap_amor_gr.txt", dtype=float)
GAP_l = GAP_l[GAP_l[:, 0].argsort()]
GAP_a = GAP_a[GAP_l[:, 0].argsort()]

axs[0].plot(liquid_x, liquid_y, label='MTP')
axs[0].plot(GAP_l[:,0],GAP_l[:,1], dashes=[6, 2], color='black', label='DFT')
axs[0].set_yticks(np.arange(0.0, 2.1, 0.5))
axs[0].set_ylim(0, 2.5)
axs[0].set_xlim(1,6.5)
axs[0].annotate('Liquid 1200K',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
axs[0].legend()
# axs.axvline(x=3.0, ymin=0, ymax=1, linestyle='--')


axs[1].plot(amor_x, amor_y, label='MTP')
axs[1].plot(GAP_a[:,0],GAP_a[:,1], dashes=[6, 2], color='black', label='DFT')
axs[1].set_yticks(np.arange(0.0, 2.6, 0.5))
axs[1].set_ylim(0, 3.0)
axs[1].set_xlim(1,6.5)
axs[1].annotate('Amorphous 300K',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
# axs[1].axvline(x=3.22, ymin=0, ymax=1, linestyle='--')

axs[1].set_xlabel(r"r [$\AA$]")
axs[0].set_ylabel(r"g(r)")
axs[0].yaxis.set_label_coords(-.1, 0.0)

plt.show()