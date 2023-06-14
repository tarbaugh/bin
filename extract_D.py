import numpy as np
from sklearn.linear_model import LinearRegression
import os

all_dat = []
for _, dirs, _ in os.walk(os.getcwd()):
  for dir in dirs:
    dat = np.loadtxt(os.path.join(dir, "msd.dat"))

    threshold = 1.2

    mask = dat[:,0] <= threshold
    new_matrix = np.delete(dat, np.where(mask), axis=0)

    x = new_matrix[:,0].reshape(-1, 1)
    y = new_matrix[:,1]

    model = LinearRegression()
    model.fit(x, y)
    all_dat.append([int(dir), model.coef_[0]])

    # r_sq = model.score(x, y)
    # print("slope: {} \nintercept: {}".format(model.coef_, model.intercept_))
    # print("r2: {}".format(r_sq))
all_dat.sort(key=lambda x: x[0])
all_dat = np.array(all_dat)
all_dat[:,1] *= 0.0001
np.savetxt("D.txt", all_dat)