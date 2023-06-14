# from ovito.io import import_file, export_file
# from ovito.modifiers import CoordinationAnalysisModifier, WrapPeriodicImagesModifier, TimeAveragingModifier

# def change_pbc(frame, data):
#     data.cell_.pbc = (True,True,True)

# def getRDF(filename, start=0, skip=1, end=-1):
#     # Load input data.
#   pipeline = import_file(filename)

#   pipeline.modifiers.append(change_pbc)
#   pipeline.modifiers.append(WrapPeriodicImagesModifier())

#   # Calculate partial RDFs:
#   pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=12.0, number_of_bins=175))

#   # print(pipeline.source.num_frames)
#   if end == -1:
#     end = pipeline.source.num_frames-1
#   pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf', interval=(start, end), sampling_frequency=skip))
#   # Access the output DataTable:
#   rdf_table = pipeline.compute().tables['coordination-rdf[average]'].xy()
#   x = rdf_table[:,0]
#   y = rdf_table[:,1]

#   print(list(y))

# getRDF('novdw/vdw_1150.xyz')

from ase.io import read
from asap3.analysis.rdf import RadialDistributionFunction

traj = read('vdw_1150.xyz', ':')
RDFobj = None
for atoms in traj:
    if RDFobj is None:
        RDFobj = RadialDistributionFunction(atoms, 10.0, 100)
    else:
        RDFobj.atoms = atoms  # Fool RDFobj to use the new atoms
    RDFobj.update()           # Collect data
rdf = RDFobj.get_rdf()

import matplotlib.pyplot as plt
plt.plot(rdf)
plt.show()
print(rdf)