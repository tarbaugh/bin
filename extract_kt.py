import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import os

all_dat = []
for _, dirs, _ in os.walk(os.getcwd()):
  for dir in dirs:
    if int(dir) > 640 and int(dir) != 690:
        with open(os.path.join(dir, "{}.01.xyz".format(dir)), "r") as f:
            for i, l in enumerate(f):
                if i == 6:
                    lo, hi = l.split()
                    vol = ((float(hi) - float(lo))**3) * 1e-30
                    break
        dat = np.loadtxt(os.path.join(dir, "sq.dat"))

        threshold = 1.0

        mask = dat[:,0] >= threshold
        new_matrix = np.delete(dat, np.where(mask), axis=0)

        x = new_matrix[:,0].reshape(-1, 1)
        y = new_matrix[:,1]

        model = LinearRegression()
        poly = PolynomialFeatures(degree=2)
        poly_features = poly.fit_transform(x)
        model.fit(poly_features, y)

        y_predicted = model.predict(poly.fit_transform(np.array([0]).reshape(-1,1)))
        print(dir,"\t",y_predicted[0]*vol/(1600*int(dir)*1.380649e-23)) # NEED BOLTZMANS CONST and CORRECT VOL UNITS