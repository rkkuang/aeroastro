import matplotlib.pyplot as plt
# from skimage import exposure
# import cv2
import numpy as np
import scipy.optimize as opt
from scipy.optimize import curve_fit

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
if __name__ == '__main__':
    # # 2d gaussian fitting 
    # # Create x and y indices
    # x = np.linspace(0, 200, 201)
    # y = np.linspace(0, 200, 201)
    # x, y = np.meshgrid(x, y)

    # #create data
    # data = twoD_Gaussian((x, y), 3, 100, 100, 20, 40, 0, 10)

    # # plot twoD_Gaussian data generated above
    # plt.figure()
    # plt.imshow(data.reshape(201, 201))
    # plt.colorbar()

    # # add some noise to the data and try to fit the data generated beforehand
    # initial_guess = (3,100,100,20,40,0,10)

    # data_noisy = data + 0.2*np.random.normal(size=data.shape)
    # # print("size: ",data_noisy.shape)
    # # input(">>>>>>>>>>>>>>>>>")
    # popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_noisy, p0=initial_guess)

    # data_fitted = twoD_Gaussian((x, y), *popt)

    # fig, ax = plt.subplots(1, 1)
    # # ax.hold(True)
    # ax.imshow(data_noisy.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
    #     extent=(x.min(), x.max(), y.min(), y.max()))
    # ax.contour(x, y, data_fitted.reshape(201, 201), 8, colors='w')
    # plt.figure()
    # plt.imshow(data_fitted.reshape(201, 201))
    # plt.show()


    # one dimension:
    # Create data:
    # x0, sigma = 0, 0.1
    # y, xe  = np.histogram(np.random.normal(x0, sigma, 1000))
    # x = .5 * (xe[:-1] + xe[1:])

    # # Function to be fitted
    # def gauss(x, x0, y0, sigma):
    #     p = [x0, y0, sigma]
    #     return p[1]* np.exp(-((x-p[0])/p[2])**2)

    # # Initialization parameters
    # p0 = [1., 1., 1.]
    # # Fit the data with the function
    # fit, tmp = curve_fit(gauss, x, y, p0=p0)

    # # Plot the results
    # plt.title('Fit parameters:\n x0=%.2e y0=%.2e sigma=%.2e' % (fit[0], fit[1], fit[2]))
    # # Data
    # plt.plot(x, y, 'r--')
    # # Fitted function
    # x_fine = np.linspace(xe[0], xe[-1], 100)
    # plt.plot(x_fine, gauss(x_fine, fit[0], fit[1], fit[2]), 'b-')
    # plt.savefig('Gaussian_fit.png')
    # plt.show()

    # http://www.pianshen.com/article/5879285638/


    # NumPy: Generate a generic 2D Gaussian-like array
    # https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-79.php
    # import numpy as np
    # x, y = np.meshgrid(np.linspace(-1,1,10), np.linspace(-1,1,10))
    # d = np.sqrt(x*x+y*y)
    # sigma, mu = 1.0, 0.0
    # g = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
    # print("2D Gaussian-like array:")
    # print(g)

    # *************** https://scipy-cookbook.readthedocs.io/items/FittingData.html#Fitting-a-2D-gaussian ******
    # 

    pass



