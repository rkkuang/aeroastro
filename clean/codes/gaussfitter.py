# *************** https://scipy-cookbook.readthedocs.io/items/FittingData.html#Fitting-a-2D-gaussian ******
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fitbeam(coresize, data):
    # coresize = 512
    Xin, Yin = np.mgrid[0:coresize, 0:coresize]
    # data = tele1.dirty_beam
    # plt.matshow(data, cmap=plt.cm.gist_earth_r)

    params = fitgaussian(data)
    fit = gaussian(*params)
    fit_res = fit(*np.indices(data.shape))
    # plt.contour(fit_res, cmap=plt.cm.copper)
    # ax = plt.gca()
    # (height, x, y, width_x, width_y) = params

    # plt.text(0.95, 0.05, """
    # x : %.1f
    # y : %.1f
    # width_x : %.1f
    # width_y : %.1f""" %(x, y, width_x, width_y),
    #         fontsize=16, horizontalalignment='right',
    #         verticalalignment='bottom', transform=ax.transAxes)

    # plt.matshow(fit_res, cmap=plt.cm.gist_earth_r)
    return fit_res

if __name__ == '__main__':
    # Create the gaussian data
    Xin, Yin = np.mgrid[0:201, 0:201]
    data = gaussian(3, 100, 100, 20, 40)(Xin, Yin) + np.random.random(Xin.shape)

    plt.matshow(data, cmap=plt.cm.gist_earth_r)

    params = fitgaussian(data)
    fit = gaussian(*params)

    plt.contour(fit(*np.indices(data.shape)), cmap=plt.cm.copper)
    ax = plt.gca()
    (height, x, y, width_x, width_y) = params

    plt.text(0.95, 0.05, """
    x : %.1f
    y : %.1f
    width_x : %.1f
    width_y : %.1f""" %(x, y, width_x, width_y),
            fontsize=16, horizontalalignment='right',
            verticalalignment='bottom', transform=ax.transAxes)
    plt.show()