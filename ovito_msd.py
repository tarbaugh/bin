import numpy
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

def modify(frame, data):
    # Access the per-particle displacement magnitudes computed by the 
    # 'Displacement Vectors' modifier preceding this user-defined modifier in the 
    # data pipeline:
    displacement_magnitudes = data.particles['Displacement Magnitude']

    # Compute MSD:
    msd = numpy.sum(displacement_magnitudes ** 2) / len(displacement_magnitudes)

    # Output value as a global attribute to make it available within OVITO:
    data.attributes["MSD"] = msd 

from ovito.io import import_file, export_file
from ovito.modifiers import CalculateDisplacementsModifier
import numpy

# Load input data and create a data pipeline.
pipeline = import_file(file)

# Calculate per-particle displacements with respect to initial simulation frame:
pipeline.modifiers.append(CalculateDisplacementsModifier())

# Define the custom modifier function:
def calculate_msd(frame, data):

    # Access the per-particle displacement magnitudes computed by the 
    # CalculateDisplacementsModifier that precedes this user-defined modifier in the 
    # data pipeline:
    displacement_magnitudes = data.particles['Displacement Magnitude']

    # Compute MSD:
    msd = numpy.sum(displacement_magnitudes ** 2) / len(displacement_magnitudes)

    # Output MSD value as a global attribute: 
    data.attributes["MSD"] = msd 

# Insert user-defined modifier function into the data pipeline.
pipeline.modifiers.append(calculate_msd)

# Export calculated MSD value to a text file and let OVITO's data pipeline do the rest:
export_file(pipeline, "output/msd_data.txt", 
    format = "txt/attr",
    columns = ["Timestep", "MSD"],
    multiple_frames = True)