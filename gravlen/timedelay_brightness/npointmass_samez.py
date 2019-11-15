from scipy.optimize import fsolve,leastsq
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

def commagnitude(a, px, py, Solux, Soluy):
    Mag = []
    for thetax,thetay in zip(Solux,Soluy):
        # print(thetax,thetay)
        r2 = (thetax-px)**2 + (thetay-py)**2
        A11 = 1 - np.sum( a*(1/r2 - 2*(thetax-px)**2/r2**2) )
        A12 = 2*np.sum( a*((thetax-px)*(thetay-py)/r2**2) )
        A22 = 1 - np.sum( a*(1/r2 - 2*(thetay-py)**2/r2**2) )
        mag = 1/(A11*A22-A12**2)
        Mag.append(mag)
    return np.array(Mag)

def comtimedelay(a, px, py, Solux, Soluy,betax,betay):
    T = []
    for thetax,thetay in zip(Solux,Soluy):
        r = np.sqrt((thetax-betax)**2 + (thetay-betay)**2)
        t = 0.5*r**2 - np.sum(a*np.log(r))
        T.append(t)
    return np.array(T)

def solvethetas(Rangex,Rangey,func,precision,prefix):
    solux, solustrx = [], []
    soluy, solustry = [], []
    for rangex in Rangex:
        for rangey in Rangey:
            s=fsolve(func,[rangex,rangey])
            # print(rangex*rad2arcsec, rangey*rad2arcsec, s[0]*rad2arcsec,s[1]*rad2arcsec,"solution")
            # print("funcres: ", func(s))
            # print(np.isnan(s[0]),s[0])
            ress = func(s)
            # print(ress)
            s = [i*rad2arcsec for i in s]
            if (not np.isnan(ress[0])) and (not np.isnan(ress[1])):
                if abs(ress[0])<precision and abs(ress[1]<precision):
                    if str(s[0])[:prefix] not in solustrx and str(s[1])[:prefix] not in solustry:
                        solustrx.append(str(s[0])[:prefix])
                        solustry.append(str(s[1])[:prefix])
                        solux.append(s[0])
                        soluy.append(s[1])
    return solux, soluy  

def checkMagTimedelay(a, px, py, solux,soluy, betax,betay):
    # for i in range(len(solux)):
    #     funcres = func([solux[i]*arcsec2rad,soluy[i]*arcsec2rad])
    #     print(funcres)
    Mag = commagnitude(a, px, py, np.array(solux)*arcsec2rad, np.array(soluy)*arcsec2rad)
    # print("Mag",Mag)
    T = comtimedelay(a, px, py, np.array(solux)*arcsec2rad, np.array(soluy)*arcsec2rad, betax,betay)
    # print("Time delay",T)
    # print("Mag maxindex {}, time delay maxindex {}".format(np.argmin(abs(Mag), axis=0),np.argmin(T, axis=0)))
    return np.argmax(abs(Mag), axis=0), np.argmax(T, axis=0), Mag, T

if __name__ == '__main__':

    # position of the N point mass
    # px = np.array([-1,1]).astype("float32") # in arcsec
    # py = np.array([0,0]).astype("float32")
    # thetaE = np.array([1,2]).astype("float32") # in arcsec

    px = np.array([-1,1]).astype("float32") # in arcsec
    py = np.array([0,0]).astype("float32")
    thetaE = np.array([1,2]).astype("float32") # in arcsec

    px *= arcsec2rad
    py *= arcsec2rad
    # Einstein radius of the N point mass
    
    thetaE *= arcsec2rad
    # a = thetaE**arcsec2rad # bug !!!
    a = thetaE**2


    # betaxy = np.array([1,0]).astype("float32") # in arcsec
    # betaxy *= arcsec2rad

    # betax, betay = betaxy[0], betaxy[1]
    
    # Question: How to find all solutions
    N = len(px)
    scale = 1
    Nscale = 10
    prefix = 7
    precision = 1e-11
    minx, maxx = - (np.max(thetaE)+np.max(np.abs(px))), (np.max(thetaE)+np.max(np.abs(px)))
    miny, maxy = - (np.max(thetaE)+np.max(np.abs(py))), (np.max(thetaE)+np.max(np.abs(py)))

    Rangex = np.linspace(scale*minx,scale*maxx,N*(Nscale+1))
    Rangey = np.linspace(scale*miny,scale*maxy,N*(Nscale+1))


    # betaxy = np.array([1,0]).astype("float32") # in arcsec
    # betaxy *= arcsec2rad
    # betax, betay = betaxy[0], betaxy[1]

    Numbetas = 1e1
    flag = False
    Betax,Betay = [],[]
    for betax in np.linspace(scale*minx,scale*maxx,Numbetas):
        for betay in np.linspace(scale*miny,scale*maxy,Numbetas):
            solux,soluy = solvethetas(Rangex,Rangey,func,precision,prefix)
            maxMagindex,maxTindex,_,_ = checkMagTimedelay(a, px, py, solux,soluy, betax,betay)
            # print("Mag maxindex {}, time delay maxindex {}".format(maxMagindex,maxTindex))
            if maxMagindex == maxTindex:
                flag = True
                Betay.append(betay)
                Betax.append(betax)

    if flag:
        print("Yes, there exist a image with max Magnification as well as max time delay")
    else:
        print("No, does not find any case")
    Betax = np.array(Betax)
    Betay = np.array(Betay)
    print("Betax: {}, Betay: {} in arcsec".format(Betax*rad2arcsec,Betay*rad2arcsec))
    for betax,betay in zip(Betax,Betay):
        solux,soluy = solvethetas(Rangex,Rangey,func,precision,prefix)
        maxMagindex,maxTindex, Mag , T = checkMagTimedelay(a, px, py, solux,soluy, betax,betay)
        print("Mag: {}, T: {}".format(Mag,T))

# Yes, there exist a image with max Magnification as well as max time delay
# Betax: [4.848136995860841e-06, 4.848136995860841e-06], Betay: [-7.541546438005753e-06, 7.541546438005753e-06]
# Mag: [ 1.08388316 -0.06154303 -0.00902427], T: [1.56257984e-09 1.42557201e-09 1.40557417e-09]
# Mag: [ 1.08388316 -0.06154303 -0.00902427], T: [1.41894733e-09 1.40230365e-09 1.41037294e-09]



    # for rangex in Rangex:
    #     for rangey in Rangey:
    #         s=fsolve(func,[rangex,rangey])
    #         # print(rangex*rad2arcsec, rangey*rad2arcsec, s[0]*rad2arcsec,s[1]*rad2arcsec,"solution")
    #         # print("funcres: ", func(s))
    #         # print(np.isnan(s[0]),s[0])
    #         ress = func(s)
    #         # print(ress)
    #         s = [i*rad2arcsec for i in s]
    #         if (not np.isnan(ress[0])) and (not np.isnan(ress[1])):
    #             if abs(ress[0])<precision and abs(ress[1]<precision):
    #                 if str(s[0])[:prefix] not in solustrx and str(s[1])[:prefix] not in solustry:
    #                     solustrx.append(str(s[0])[:prefix])
    #                     solustry.append(str(s[1])[:prefix])
    #                     solux.append(s[0])
    #                     soluy.append(s[1])

    # solux,soluy = solvethetas(Rangex,Rangey,func,precision,prefix)
    # maxMagindex,maxTindex = checkMagTimedelay(a, px, py, solux,soluy, betax,betay)
    # print("Mag maxindex {}, time delay maxindex {}".format(maxMagindex,maxTindex))

    # for i in range(len(solux)):
    #     funcres = func([solux[i]*arcsec2rad,soluy[i]*arcsec2rad])
    #     print(funcres)
    # # Mag = commagnitude(a, px, py, np.array(solux)*arcsec2rad, np.array(solux)*arcsec2rad) # such bug!!!
    # Mag = commagnitude(a, px, py, np.array(solux)*arcsec2rad, np.array(soluy)*arcsec2rad)
    # print("Mag",Mag)
    # T = comtimedelay(a, px, py, np.array(solux)*arcsec2rad, np.array(soluy)*arcsec2rad, betax,betay)
    # print("Time delay",T)
    # print("Mag maxindex {}, time delay maxindex {}".format(np.argmin(abs(Mag), axis=0),np.argmin(T, axis=0)))


