import os
import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, WrapPeriodicImagesModifier, TimeAveragingModifier

Reduced = False

def change_pbc(frame, data):
    data.cell_.pbc = (True,True,True)

def getRDF(filename, start=0, skip=1, end=-1):
    # Load input data.
  pipeline = import_file(filename)

  pipeline.modifiers.append(change_pbc)
  pipeline.modifiers.append(WrapPeriodicImagesModifier())

  # Calculate partial RDFs:
  pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=12.0, number_of_bins=125))

  # print(pipeline.source.num_frames)
  if end == -1:
    end = pipeline.source.num_frames-1
  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf', interval=(start, end), sampling_frequency=skip))
  # Access the output DataTable:
  rdf_table = pipeline.compute().tables['coordination-rdf[average]'].xy()
  x = rdf_table[:,0]
  y = rdf_table[:,1]

  if Reduced:
    ref = pipeline.source.compute(0)
    cell = ref.cell
    a = cell[:,0][0]
    b = cell[:,1][1]
    c = cell[:,2][2]
    # o = cell[:,3]

    num_atoms = len(ref.particles.positions)
    num_dens = num_atoms/(a*b*c)
    # print(num_dens)
    for i, r in enumerate(x):
      y[i] = (y[i]-1)*4*np.pi*num_dens*r
    return x, y, num_dens
  else:
    return x, y
  # return x, y, num_dens

def getFiles():
  ctot = []
  tot = []
  for _,dirs,_ in os.walk(os.getcwd()):
    for dir in dirs:
      for root,_,files in os.walk(dir):
        for file in files:
            if file.endswith(".xyz"):
              print('Computing '+os.path.join(root, file)+'...')
              x, y = getRDF(os.path.join(root, file))
              ctot.append([x, y])
        print(np.average(np.array(ctot), axis=1))
        break

  for c in tot:
    fig, axs = plt.subplots(1, 1)
    axs.plot(c[0], c[1])
  plt.show()

getFiles()

