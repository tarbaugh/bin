import os
import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, WrapPeriodicImagesModifier, TimeAveragingModifier

filename = 'zipoli900.xyz'

def change_pbc(frame, data):
    data.cell_.pbc = (True,True,True)

def getRDF(filename, dir):
    # Load input data.
  pipeline = import_file(filename)

  pipeline.modifiers.append(change_pbc)
  pipeline.modifiers.append(WrapPeriodicImagesModifier())

  # Calculate partial RDFs:
  pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=9.0, number_of_bins=100, partial=True))

  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
  # Access the output DataTable:
  rdf_table = pipeline.compute().tables['coordination-rdf[average]'].xy()
  return rdf_table[:,0], rdf_table[:,1], rdf_table[:,2], rdf_table[:,3]

def getFiles():
  for _,dirs,_ in os.walk(os.getcwd()):
    for dir in dirs:
      for root,_,files in os.walk(dir):
        for file in files:
            if file.endswith(".xyz"):
              print('Computing '+os.path.join(root, file)+'...')
              
rtx, rty1, rty2, rty3 = getRDF(filename, 'root')

fig = plt.figure()
fig, axs = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(hspace=0)

# GeGeZ = np.loadtxt("ZipoliGeGe.txt", dtype=float)
# GeTeZ = np.loadtxt("ZipoliGeTe.txt", dtype=float)
# TeTeZ = np.loadtxt("ZipoliTeTe.txt", dtype=float)
# GeGeDFT = np.loadtxt("ZipoliDFTGeGe.txt", dtype=float)
# GeTeDFT = np.loadtxt("ZipoliDFTGeTe.txt", dtype=float)
# TeTeDFT = np.loadtxt("ZipoliDFTTeTe.txt", dtype=float)


# GeGeZ = np.loadtxt("ZipoliGeGe.txt", dtype=float)
# GeTeZ = np.loadtxt("ZipoliGeTe.txt", dtype=float)
# TeTeZ = np.loadtxt("ZipoliTeTe.txt", dtype=float)
GeGeDFT = np.loadtxt("gege.txt", dtype=float)
GeTeDFT = np.loadtxt("gete.txt", dtype=float)
TeTeDFT = np.loadtxt("tete.txt", dtype=float)

GeTeDFT[:,1] -= 2
TeTeDFT[:,1] -= 1

axs[0].plot(rtx, rty1)
# axs[0].plot(GeGeZ[:,0],GeGeZ[:,1])
axs[0].plot(GeGeDFT[:,0],GeGeDFT[:,1], dashes=[6, 2], color='black')
axs[0].set_yticks(np.arange(0.0, 2.0, 0.5))
axs[0].set_ylim(0, 2.5)
axs[0].set_xlim(0,9)
axs[0].annotate('Ge-Ge',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
axs[0].axvline(x=3.0, ymin=0, ymax=1, linestyle='--')


axs[1].plot(rtx, rty2)
# axs[1].plot(GeTeZ[:,0],GeTeZ[:,1])
axs[1].plot(GeTeDFT[:,0],GeTeDFT[:,1], dashes=[6, 2], color='black')
axs[1].set_yticks(np.arange(0.0, 3.5, 0.5))
axs[1].set_ylim(0, 3.5)
axs[1].set_xlim(0,9)
axs[1].annotate('Ge-Te',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
axs[1].axvline(x=3.22, ymin=0, ymax=1, linestyle='--')


axs[2].plot(rtx, rty3)
# axs[2].plot(TeTeZ[:,0],TeTeZ[:,1])
axs[2].plot(TeTeDFT[:,0],TeTeDFT[:,1], dashes=[6, 2], color='black')
axs[2].set_yticks(np.arange(0.0, 2.0, 0.5))
axs[2].set_ylim(0, 2)
axs[2].set_xlim(0,9)
axs[2].annotate('Te-Te',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
axs[2].axvline(x=3.0, ymin=0, ymax=1, linestyle='--')


axs[2].set_xlabel(r"r [$\AA$]")
axs[1].set_ylabel(r"g(r)")

plt.show()




