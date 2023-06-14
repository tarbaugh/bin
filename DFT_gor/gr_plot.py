import os
import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, WrapPeriodicImagesModifier, TimeAveragingModifier

comp796 = False
Temp = False
comp20 = False
DFT = False
Reduced = False
Reduced_Amor = True
Neutron = False
Reduced_Amor_Cutoff = False

def change_pbc(frame, data):
    data.cell_.pbc = (True,True,True)

def getRDF(filename, dir, start=0, skip=1, end=-1):
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

  if Reduced or Reduced_Amor or Reduced_Amor_Cutoff:
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
  for _,dirs,_ in os.walk(os.getcwd()):
    for dir in dirs:
      for root,_,files in os.walk(dir):
        for file in files:
            if file.endswith(".xyz"):
              print('Computing '+os.path.join(root, file)+'...')

if DFT:
  dft_x, dft_y = getRDF('md.nb', 'root', start=5000, skip=2)
  vdw_1150_x, vdw_1150_y = getRDF('vdw/vdw_1150.xyz', 'root', start=0, skip=1)
  novdw_1150_x, novdw_1150_y = getRDF('novdw/vdw_1150.xyz', 'root', start=0, skip=1)
  vdw_20_x, vdw_20_y = getRDF('20/vdw20.xyz', 'root', start=0, skip=1)
  novdw_20_x, novdw_20_y = getRDF('20/novdw20.xyz', 'root', start=0, skip=1)

  fig, axs = plt.subplots(2, 1)
  fig.subplots_adjust(hspace=0)
  axs[0].plot(dft_x, dft_y, label='DFT', color='black')
  axs[1].plot(dft_x, dft_y, label='DFT', color='black')
  axs[0].plot(novdw_1150_x, novdw_1150_y, label='MTP', color='blue')
  axs[1].plot(vdw_1150_x, vdw_1150_y, label='MTP + VdW', color='blue')
  axs[0].plot(novdw_20_x, novdw_20_y, label='MTP20', color='red')
  axs[1].plot(vdw_20_x, vdw_20_y, label='MTP20 + VdW', color='red')

  for i in range(2):
    axs[i].annotate('BLYP Functional (1150K)',
              xy=(1, 0), xycoords='axes fraction',
              xytext=(-20, 20), textcoords='offset pixels',
              horizontalalignment='right',
              verticalalignment='bottom')               
  axs[0].legend()
  axs[1].legend()
  plt.show()

if Neutron:
  vdw_653_x, vdw_653_y = getRDF('650.xyz', 'root')
  vdw_685_x, vdw_685_y = getRDF('690.xyz', 'root')
  vdw_796_x, vdw_796_y = getRDF('800.xyz', 'root')

  novdw_653_x, novdw_653_y = getRDF('pbe.xyz', 'root', start=0, end=2, skip=1)
  novdw_685_x, novdw_685_y = getRDF('pbe.xyz', 'root', start=0, end=2, skip=1)
  novdw_796_x, novdw_796_y = getRDF('pbe.xyz', 'root', start=0, end=2, skip=1)
  
  neut_796 = np.loadtxt("neutron/ge15te85_neutron_796K", dtype=float)
  neut_796 = neut_796[neut_796[:, 0].argsort()]

  neut_685 = np.loadtxt("neutron/ge15te85_neutron_685K", dtype=float)
  neut_685 = neut_685[neut_685[:, 0].argsort()]

  neut_653 = np.loadtxt("neutron/ge15te85_neutron_653K", dtype=float)
  neut_653 = neut_653[neut_653[:, 0].argsort()]

  fig, axs = plt.subplots(2, 2)
  # fig.subplots_adjust(hspace=0)
  
  for i in range(2):
    axs[i,0].plot(neut_796[:,0],neut_796[:,1], alpha=0.6, color='red', label='Neutron Scattering (796 K)')
    axs[i,0].plot(neut_685[:,0],neut_685[:,1], alpha=0.6, color='orange', label='Neutron Scattering (685 K)')
    axs[i,0].plot(neut_653[:,0],neut_653[:,1], alpha=0.6, color='green', label='Neutron Scattering (653 K)')

  axs[0,1].plot(vdw_653_x, vdw_653_y, alpha=0.6, label='MTP BLYP (653 K)', color='green')
  axs[0,1].plot(vdw_685_x, vdw_685_y, alpha=0.6, label='MTP BLYP (685 K)', color='orange')
  axs[0,1].plot(vdw_796_x, vdw_796_y, alpha=0.6, label='MTP BLYP (796 K)', color='red')

  axs[1,1].plot(novdw_653_x, novdw_653_y, alpha=0.6, label='MTP PBE (653 K)', color='green')
  axs[1,1].plot(novdw_685_x, novdw_685_y, alpha=0.6, label='MTP PBE (685 K)', color='orange')
  axs[1,1].plot(novdw_796_x, novdw_796_y, alpha=0.6, label='MTP PBE (796 K)', color='red')

  for i in range(2):
    for j in range(2):
      # axs[i,j].annotate('BLYP Functional',
      #           xy=(1, 0), xycoords='axes fraction',
      #           xytext=(-20, 20), textcoords='offset pixels',
      #           horizontalalignment='right',
      #           verticalalignment='bottom')      
      axs[i,j].legend()
      axs[i,j].set_ylim(-0.1, 2.5)
      axs[i,j].set_yticks(np.arange(0.0, 2.6, 0.5))
      axs[i,j].set_xlim(1,8)
  
  plt.show()

if Reduced:
  vdw_592_x, vdw_592_y, num_dens = getRDF('690.xyz', 'root')
  # for i, r in enumerate(vdw_592_x):
  #   vdw_592_y[i] = (vdw_592_y[i]-1)*4*np.pi*num_dens*r
  vdw_359_x, vdw_359_y, num_dens = getRDF('650.xyz', 'root')
  # for i, r in enumerate(vdw_359_x):
  #   vdw_359_y[i] = (vdw_359_y[i]-1)*4*np.pi*num_dens*r
  vdw_833_x, vdw_833_y, num_dens = getRDF('800.xyz', 'root')
  # for i, r in enumerate(vdw_833_x):
  #   vdw_833_y[i] = (vdw_833_y[i]-1)*4*np.pi*num_dens*r
  # novdw_592_x, novdw_592_y = getRDF('novdw/vdw_592.xyz', 'root', start=0, skip=1)
  # novdw_833_x, novdw_833_y = getRDF('novdw/vdw_833.xyz', 'root', start=0, skip=1)
  
  # neut_592 = np.loadtxt("Gr/exper_592.txt", dtype=float)
  # neut_592 = neut_592[neut_592[:, 0].argsort()]

  neut_833 = np.loadtxt("Gr/exper_833.txt", dtype=float)
  neut_833 = neut_833[neut_833[:, 0].argsort()]

  fig, axs = plt.subplots(1, 1)
  fig.subplots_adjust(hspace=0)
  
  axs.plot(neut_833[:,0],neut_833[:,1], alpha=0.6, color='blue', label='Neutron Scattering (833 K)')

  axs.plot(vdw_359_x, vdw_359_y, alpha=0.6, label='(650 K)', color='yellow')
  axs.plot(vdw_592_x, vdw_592_y, alpha=0.6, label='(690 K)', color='orange')
  axs.plot(vdw_833_x, vdw_833_y, alpha=0.6, label='(800 K)', color='red')

  axs.annotate('BLYP Functional',
    xy=(1, 0), xycoords='axes fraction',
    xytext=(-20, 20), textcoords='offset pixels',
    horizontalalignment='right',
    verticalalignment='bottom')      
  axs.legend()
  axs.set_ylim(-1, 1.2)
  # axs[i,j].set_yticks(np.arange(0.0, 2.6, 0.5))
  axs.set_xlim(1,12)

  plt.show()

if Reduced_Amor:
  blyp_359_x, blyp_359_y, num_dens = getRDF('660.xyz', 'root', skip=1)
  blyp_359_x2, blyp_359_y2, num_dens = getRDF('800.xyz', 'root', skip=1)
  # pbe_359_x, pbe_359_y, num_dens = getRDF('pbe.xyz', 'root', skip=1)

  amor_Gr = np.loadtxt("Gr/exper_660.txt", dtype=float)
  amor_Gr = amor_Gr[amor_Gr[:, 0].argsort()]

  amor_Gr_high = np.loadtxt("Gr/exper_833.txt", dtype=float)
  amor_Gr_high = amor_Gr_high[amor_Gr_high[:, 0].argsort()]
  
  fig, axs = plt.subplots(2, 1)
  fig.subplots_adjust(hspace=0.4)

  # neut_653 = np.loadtxt("neutron/ge15te85_neutron_653K", dtype=float)
  # neut_653 = neut_653[neut_653[:, 0].argsort()]
  # for i, r in enumerate(neut_653[:,0]):
  #   neut_653[i,1] = (neut_653[i,1]-1)*4*np.pi*num_dens*r

  axs[1].plot(blyp_359_x, blyp_359_y, label='MTP (660 K)', color='blue')
  axs[1].plot(amor_Gr[:,0], amor_Gr[:,1], label='Exp. (~660 K)', color='black')
  # axs[1].plot(neut_653[:,0], neut_653[:,1], label='Exp. (592 K)', color='black')
  axs[0].plot(amor_Gr_high[:,0], amor_Gr_high[:,1], label='Exp. (833 K)', color='black')
  axs[0].plot(blyp_359_x2, blyp_359_y2, label='MTP (830 K)', color='red')

  for ax in axs.flat:
    ax.set(xlabel='r (Angstroms)', ylabel='')
  for ax in axs.flat:
    ax.label_outer()

  axs[1].set_title('Low Temperature G(r)')
  axs[1].set_ylim(-1.1, 1.8)
  axs[1].set_yticks(np.arange(-1, 1.8, 0.4))
  axs[1].set_xlim(1,12)
  axs[1].set_xticks(np.arange(0.0, 12, 1.0))
  axs[1].legend()

  axs[0].set_title('High Temperature G(r)')
  axs[0].set_ylim(-1.1, 1.8)
  axs[0].set_yticks(np.arange(-1, 1.8, 0.4))
  axs[0].set_xlim(1,12)
  axs[0].set_xticks(np.arange(0.0, 12, 1.0))
  axs[0].legend()

  plt.show()
  
if Reduced_Amor_Cutoff:
  # blyp_359_x1, blyp_359_y1, num_dens = getRDF('d3blyp.xyz', 'root', start=200, end=400, skip=1)
  blyp_359_x1, blyp_359_y1, num_dens = getRDF('800.xyz', 'root', start=50, end=290, skip=1)
  blyp_359_x, blyp_359_y, num_dens = getRDF('d3blyp.xyz', 'root', start=6800, end=7000, skip=1)
  pbe_359_x1, pbe_359_y1, num_dens = getRDF('187.xyz', 'root', start=200, end=400, skip=1)
  pbe_359_x, pbe_359_y, num_dens = getRDF('187.xyz', 'root', start=6800, end=7000, skip=1)

  amor_Gr = np.loadtxt("Gr/exper_592.txt", dtype=float)
  amor_Gr = amor_Gr[amor_Gr[:, 0].argsort()]

  amor_Gr_high = np.loadtxt("Gr/exper_833.txt", dtype=float)
  amor_Gr_high = amor_Gr_high[amor_Gr_high[:, 0].argsort()]
  
  fig, axs = plt.subplots(2, 1)
  fig.subplots_adjust(hspace=0)
  axs[0].plot(pbe_359_x1, pbe_359_y1, label='~900K', color='red')
  axs[0].plot(pbe_359_x, pbe_359_y, label='~600K', color='blue')
  axs[1].plot(blyp_359_x1, blyp_359_y1, label='~900K', color='red')
  axs[1].plot(blyp_359_x, blyp_359_y, label='~600K', color='blue')

  axs[0].annotate('7 A cutoff',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')  
  axs[1].annotate('8 A cutoff',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')  
  for i in range(2):
    # axs[i].plot(amor_Gr[:,0], amor_Gr[:,1], label='Experiment (592 K)', color='blue')
    # axs[i].plot(amor_Gr_high[:,0], amor_Gr_high[:,1], label='Experiment (833 K)', color='red')
    axs[i].fill(np.append(amor_Gr[:,0], amor_Gr_high[:,0][::-1]), np.append(amor_Gr[:,1], amor_Gr_high[:,1][::-1]), color='lightblue')
    axs[i].set_ylim(-1.1, 1.8)
    axs[i].set_yticks(np.arange(-1, 1.8, 0.4))
    axs[i].set_xlim(1,12)
    axs[i].set_xticks(np.arange(0.0, 12, 1.0))
    axs[i].legend()

  plt.show()

if comp20:
  vdw_20_x, vdw_20_y = getRDF('20/vdw/vdw_685.xyz', 'root', start=0, skip=1)
  novdw_20_x, novdw_20_y = getRDF('20/novdw/vdw_685.xyz', 'root', start=0, skip=1)
  vdw_x, vdw_y = getRDF('vdw/vdw_685.xyz', 'root', start=0, skip=1)
  novdw_x, novdw_y = getRDF('novdw/vdw_685.xyz', 'root', start=0, skip=1)

  fig, axs = plt.subplots(2, 1)
  fig.subplots_adjust(hspace=0)
  axs[0].plot(novdw_x, novdw_y, label='MTP', color='blue')
  axs[1].plot(vdw_x, vdw_y, label='MTP + VdW', color='blue')
  axs[0].plot(novdw_20_x, novdw_20_y, label='MTP20', color='red')
  axs[1].plot(vdw_20_x, vdw_20_y, label='MTP20 + VdW', color='red')

  for i in range(2):
    axs[i].annotate('BLYP Functional (1150K)',
              xy=(1, 0), xycoords='axes fraction',
              xytext=(-20, 20), textcoords='offset pixels',
              horizontalalignment='right',
              verticalalignment='bottom')               
  axs[0].legend()
  axs[1].legend()
  plt.show()

if Temp:
  vdw_685_x, vdw_685_y = getRDF('vdw/vdw_685.xyz', 'root', start=0, skip=1)
  vdw_653_x, vdw_653_y = getRDF('vdw/vdw_653.xyz', 'root', start=0, skip=1)
  vdw_592_x, vdw_592_y = getRDF('vdw/vdw_592.xyz', 'root', start=0, skip=1)
  vdw_359_x, vdw_359_y = getRDF('vdw/vdw_359.xyz', 'root', start=0, skip=1)
  novdw_685_x, novdw_685_y = getRDF('novdw/vdw_685.xyz', 'root', start=0, skip=1)
  novdw_653_x, novdw_653_y = getRDF('novdw/vdw_653.xyz', 'root', start=0, skip=1)
  novdw_592_x, novdw_592_y = getRDF('novdw/vdw_592.xyz', 'root', start=0, skip=1)
  novdw_359_x, novdw_359_y = getRDF('novdw/vdw_359.xyz', 'root', start=0, skip=1)

  fig, axs = plt.subplots(2, 1)
  fig.subplots_adjust(hspace=0)
  axs[0].plot(novdw_685_x, novdw_685_y, label='MTP', color='red')
  axs[0].plot(novdw_653_x, novdw_653_y, label='MTP', color='orange')
  axs[0].plot(novdw_592_x, novdw_592_y, label='MTP', color='lightblue')
  axs[0].plot(novdw_359_x, novdw_359_y, label='MTP', color='blue')
  axs[1].plot(vdw_685_x, vdw_685_y, label='MTP + Vdw', color='red')
  axs[1].plot(vdw_653_x, vdw_653_y, label='MTP + Vdw', color='orange')
  axs[1].plot(vdw_592_x, vdw_592_y, label='MTP + Vdw', color='lightblue')
  axs[1].plot(vdw_359_x, vdw_359_y, label='MTP + Vdw', color='blue')

  for i in range(2):
    axs[i].annotate('BLYP Functional (1150K)',
              xy=(1, 0), xycoords='axes fraction',
              xytext=(-20, 20), textcoords='offset pixels',
              horizontalalignment='right',
              verticalalignment='bottom')               
  axs[0].legend()
  axs[1].legend()
  plt.show()

if comp796:
  # mod_x, mod_y = getRDF('blyp.xyz', 'root', start=0, end=1)
  mod2_x, mod2_y = getRDF('blyp650.xyz', 'root')
  
  neut_685 = np.loadtxt("NewPaper/653.txt", dtype=float)
  neut_685 = neut_685[neut_685[:, 0].argsort()]
  neut_796 = np.loadtxt("NewPaper/653.txt", dtype=float)
  neut_796 = neut_796[neut_796[:, 0].argsort()]

  fig, axs = plt.subplots(1, 1)
  fig.subplots_adjust(hspace=0)
  
  # axs.plot(neut_653[:,0],neut_653[:,1], color='blue', label='Neutron Scattering (653 K)')
  # axs.plot(neut_685[:,0],neut_685[:,1], color='orange', label='Neutron Scattering (685 K)')
  axs.plot(neut_796[:,0],neut_796[:,1], color='red', label='Neutron Scattering (796 K)')

  # axs.plot(mod_x, mod_y, label='blyp', color='green')
  axs.plot(mod2_x, mod2_y, label='~800K', color='blue')

  # axs.annotate('BLYP Functional -D3',
  #           xy=(1, 0), xycoords='axes fraction',
  #           xytext=(-20, 20), textcoords='offset pixels',
  #           horizontalalignment='right',
  #           verticalalignment='bottom')      
  axs.legend()
  axs.set_ylim(-0.1, 2.5)
  axs.set_yticks(np.arange(0.0, 2.6, 0.5))
  axs.set_xlim(1,12)
  
  plt.show()
  # mod_x, mod_y = getRDF('mod.xyz', 'root', start=0, skip=1)
  # vdw_796_x, vdw_796_y = getRDF('vdw/vdw_796.xyz', 'root', start=0, skip=1)
  # novdw_796_x, novdw_796_y = getRDF('novdw/vdw_796.xyz', 'root', start=0, skip=1)
  # xyz36_x, xyz36_y = getRDF('md.nb', 'root', start=3500, skip=2)
  
  # neut_796 = np.loadtxt("neutron/ge15te85_neutron_796K", dtype=float)
  # neut_796 = neut_796[neut_796[:, 0].argsort()]

  # fig, axs = plt.subplots(1, 1)
  # fig.subplots_adjust(hspace=0)
  
  # axs.plot(neut_796[:,0],neut_796[:,1], color='black', label='Neutron Scattering (796 K)')

  # axs.plot(mod_x, mod_y, label='MTP + Corrections (796 K)', color='purple')
  # axs.plot(vdw_796_x, vdw_796_y, label='MTP + VdW (796 K)', color='blue')
  # axs.plot(xyz36_x, xyz36_y, label='DFT (796 K)', color='red')
  # axs.plot(novdw_796_x, novdw_796_y, label='MTP (796 K)', color='green')

  # axs.annotate('BLYP Functional',
  #           xy=(1, 0), xycoords='axes fraction',
  #           xytext=(-20, 20), textcoords='offset pixels',
  #           horizontalalignment='right',
  #           verticalalignment='bottom')      
  # axs.legend()
  # axs.set_ylim(-0.1, 2.5)
  # axs.set_yticks(np.arange(0.0, 2.6, 0.5))
  # axs.set_xlim(1,12)
  
  # plt.show()
  
