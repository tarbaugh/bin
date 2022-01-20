import random

# experimental amorphous density of 5.88 g cm^–3
# -> 3.54 amu angstrom^-3

# Atomic mass of Germanium is 72.64 u. -> 108 * 72.64 = 7845.12 u
# Atomic mass of Tellurium is 127.6 u. -> 108 * 127.6 = 13780.8 u

# Total u = 21625.92

ByDensity = True
ByBoxSize = False
ByADensity = False

Ge = 108
Te = 108

N = Ge + Te

filename = 'data_gete5050'

if ByDensity:
    Density = 5.5
    # Density = 5.323 # Germanium
    # Density = 6.5 # Te

    GeM = Ge * 72.64 * 1.6605402e-24
    TeM = Te * 127.6 * 1.6605402e-24

    BoxSize = (100000000*((GeM+TeM)*(1/Density))**(1/3))
    BoxSizeCubed = BoxSize*BoxSize*BoxSize

if ByADensity:
    Amorphous = True
    Crystal = False
    if Crystal:
        Density = 0.0357
    if Amorphous:
        Density = 0.0334
        # Density = 0.0282 # Experimentail Ge1Te6
    
    BoxSizeCubed = N*(1/Density)
    BoxSize = BoxSizeCubed**(1/3)

if ByBoxSize:
    BoxSize = 18.6
    BoxSize = 37.047181672448836
    BoxSizeCubed = BoxSize*BoxSize*BoxSize


NperL = 0
while ((NperL*NperL*NperL)<N):
    NperL += 1
Nl = BoxSize/(NperL)

types = []
for _ in range(Ge):
    types.append(1)
for _ in range(Te):
    types.append(2)

random.shuffle(types)

with open(filename, 'w') as f:
    f.write('GeTe Input Data\n\n')

    f.write('\t{}\tatoms\n'.format(N))
    f.write('\t2\tatom types\n\n')
    f.write('0.00\t{}\txlo xhi\n'.format(BoxSize))
    f.write('0.00\t{}\tylo yhi\n'.format(BoxSize))
    f.write('0.00\t{}\tzlo zhi\n\n'.format(BoxSize))

    f.write('Masses\n\n')

    f.write('1 72.63\n')
    f.write('2 127.603\n\n')

    f.write('Atoms\n\n')


    x=0
    y=0
    z=0
    atom_num=0
    for k in range(NperL):
        y=0
        for j in range(NperL):
            x = 0
            for i in range(NperL):
                atom_num += 1
                if types:
                    t = types.pop()
                else:
                    break
                f.write('{}\t{}\t{}\t{}\t{}\n'.format(atom_num, t, x, y, z))
                x += Nl
            y += Nl
        z+= Nl

print("File {} created with ....\n".format(filename))
print("\t{} Ge Atoms\n".format(Ge))
print("\t{} Te Atoms\n".format(Te))
print("A Box Size of {} and Volume of {}\n".format(BoxSize, BoxSizeCubed))
