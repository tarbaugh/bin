import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

def readCFGEnergies(file, mtp=False):
    with open(file) as f:
        mark = -2
        energies = []
        for i, l in enumerate(f):
            if l.startswith(" Energy"):
                mark = int(i)
            if i == mark+1:
                energies.append(float(l))
                mark = -2

    return energies

def readCFGForces(file, mtp=False):
    with open(file) as f:
        mark = 99999999999
        forces = []
        for i, l in enumerate(f):
            if l.startswith(" AtomData:"):
                mark = int(i)
            if i == mark+201:
                mark = 99999999999
            if i >= mark+1:
                tmp_forces = l.split()[-3:]
                forces.append([float(j) for j in tmp_forces])

    return forces

def r2_plot(train_cfg, mtp_train_cfg, valid_cfg, mtp_valid_cfg, t):
    fig, ax = plt.subplots()
    ax.scatter(train_cfg, mtp_train_cfg, color="blue", label="Train")
    ax.scatter(valid_cfg, mtp_valid_cfg, color="orange", label="Test")
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("DFT "+t)
    ax.set_ylabel("MTP Predicted "+t)

    r2_train = r2_score(train_cfg, mtp_train_cfg)
    r2_valid = r2_score(valid_cfg, mtp_valid_cfg)
    ax.annotate('Train R2: {}'.format(r2_train),
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 35), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
    ax.annotate('Valid R2: {}'.format(r2_valid),
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
    ax.legend()

    plt.savefig("r2_"+t+".png")
    plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', action='store_true', default=False)
    parser.add_argument('-f', action='store_true', default=False)
    args = parser.parse_args()

    if not args.e:
        valid_cfg = readCFGForces("valid.cfg")
        train_cfg = readCFGForces("train.cfg")
        mtp_valid_cfg = readCFGForces("mtp_valid.cfg", mtp=True)
        mtp_train_cfg = readCFGForces("mtp_train.cfg", mtp=True)
        r2_plot(train_cfg, mtp_train_cfg, valid_cfg, mtp_valid_cfg, "Forces")

    if not args.f:
        valid_cfg = readCFGEnergies("valid.cfg")
        train_cfg = readCFGEnergies("train.cfg")
        mtp_valid_cfg = readCFGEnergies("mtp_valid.cfg", mtp=True)
        mtp_train_cfg = readCFGEnergies("mtp_train.cfg", mtp=True)
        r2_plot(train_cfg, mtp_train_cfg, valid_cfg, mtp_valid_cfg, "Energies")