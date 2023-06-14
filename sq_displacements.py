import numpy
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

from ovito.io import import_file, export_file
from ovito.modifiers import CalculateDisplacementsModifier
import numpy

# Load input data and create a data pipeline.
pipeline = import_file(file)
length = pipeline.source.num_frames-1

def calculate_msd(frame, data):

    # Access the per-particle displacement magnitudes computed by the 
    # CalculateDisplacementsModifier that precedes this user-defined modifier in the 
    # data pipeline:
    types = -1
    if types != -1:
        displacement_magnitudes = data.particles['Displacement Magnitude']
        types = data.particles['Particle Type']
        bool = types == 2
        disp = displacement_magnitudes[bool]
    else:
        disp = data.particles['Displacement Magnitude']

    # Compute MSD:
    #msd = numpy.sum(displacement_magnitudes ** 2) / len(displacement_magnitudes)
    msd = numpy.sum(disp ** 2) / len(disp)

    # Output MSD value as a global attribute: 
    data.attributes["Disp"] = msd 

for i in range(1, length, 100):
    pipeline = import_file(file)

    pipeline.modifiers.append(CalculateDisplacementsModifier(use_frame_offset=True, frame_offset=-i))
    # pipeline.modifiers.append(CalculateDisplacementsModifier())
    # Insert user-defined modifier function into the data pipeline.
    pipeline.modifiers.append(calculate_msd)

    # Export calculated MSD value to a text file and let OVITO's data pipeline do the rest:
    export_file(pipeline, "sqdisp{}.txt".format(i), 
        format = "txt/attr",
        columns = ["Timestep", "Disp"],
        multiple_frames = True,
        start_frame = i)

avgs = []
for i in range(1, length, 100):
    with open("sqdisp{}.txt".format(i), "r") as f:
        tot = 0.0
        count = 0
        for iter, l in enumerate(f):
            if iter != 0:
                tot += float(l.split()[-1])
                count += 1
        avgs.append((i, tot/count))

with open("msd.txt", "w") as f:
    for (iter, x) in avgs:
        f.write("{}\t{}\n".format(iter, x))
