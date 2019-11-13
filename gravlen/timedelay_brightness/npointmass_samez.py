from scipy.optimize import fsolve
import numpy as np

arcsec2rad = 1/3600/180*np.pi
rad2arcsec = 180/np.pi*3600

# def func(xy, a, px, py, betaxy):
def func(xy):
    #https://blog.csdn.net/jayloncheng/article/details/80003182
    x,y=xy[0],xy[1]
    r2 = (x-px)**2 + (y-py)**2
    return [ x - np.sum(a*(x-px)/r2) - betax,
            y - np.sum(a*(y-py)/r2) - betay]
if __name__ == '__main__':

    # position of the N point mass
    px = np.array([1,-1,2]).astype("float32") # in arcsec
    py = np.array([0,0,0]).astype("float32")
    px *= arcsec2rad
    py *= arcsec2rad
    # Einstein radius of the N point mass
    thetaE = np.array([1,1,1]).astype("float32") # in arcsec
    thetaE *= arcsec2rad
    # a = thetaE**arcsec2rad # bug !!!
    a = thetaE**2


    betaxy = np.array([1,0]).astype("float32") # in arcsec
    betaxy *= arcsec2rad

    betax, betay = betaxy[0], betaxy[1]

    s=fsolve(func,[1.6*arcsec2rad,0*arcsec2rad])# 1.6, -0.6
    print(s[0]*rad2arcsec,s[1]*rad2arcsec)

    s=fsolve(func,[-0.6*arcsec2rad,0*arcsec2rad])# 1.6, -0.6
    print(s[0]*rad2arcsec,s[1]*rad2arcsec)