'''
Two point mass lens system
negtive mass ratio
with total mass scales to 1

lens1 position (X,0)
lens2 position (-X,0)

ref: Teraflop per second gravitational lensing
ray-shooting using graphics processing units
Alexander C. Thompson, et al.
'''

import numpy as np
import sys, time
from utils import *
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathos.pools import ProcessPool
from tqdm import tqdm
# from matplotlib import animation as ani
# from matplotlib.animation import FuncAnimation, writers
arcsec2rad = 1/3600/180*np.pi
rad2arcsec = 180/np.pi*3600


import dill                            #pip install dill --user
# filename = 'globalsave.pkl'
# dill.dump_session(filename)

# # and to load the session again:
# dill.load_session(filename)


'''

from tqdm import tqdm
import time
 
#total参数设置进度条的总长度
with tqdm(total=100) as pbar:
  for i in range(100):
    time.sleep(0.05)
    #每次更新进度条的长度
    pbar.update(1)
https://www.jb51.net/article/166648.htm
'''



class Onelens(object):
    """docstring for Onelens"""
    def __init__(self, mass, position, dis=0):
        super(Onelens, self).__init__()
        self.mass = mass
        self.pos = position
        self.dis = dis
    
class Nlenses(object):
    """docstring for Nlenses
    lenses: [(1, 0, 0), (-1, 0, 0)] mass, x, y
    """
    def __init__(self, masses, xs, ys):
        self.masses = masses
        self.xs = xs
        self.ys = ys
        self.betax, self.betay = None, None
        self.mag = None

    def inverse_ray_shooting(self, thetax=None, thetay=None):
        # Ri2s = (thetax-self.xs)**2 + (thetay-self.ys)**2  
        # Ri2s = []
        self.betax, self.betay = thetax.copy(), thetay.copy()
        for i in range(len(self.xs)):
            # Ri2s.append( (thetax-self.xs[i])**2 + (thetay-self.ys[i])**2 )
            Ri2s = (thetax-self.xs[i])**2 + (thetay-self.ys[i])**2
            self.betax -= self.masses[i]*(thetax - self.xs[i])/Ri2s
            self.betay -= self.masses[i]*(thetay - self.ys[i])/Ri2s

        # r1_2 = (thetax - self.lens1.pos[0])**2 + (thetay - self.lens1.pos[1])**2
        # r2_2 = (thetax - self.lens2.pos[0])**2 + (thetay - self.lens2.pos[1])**2
        # self.betax = thetax - self.mu1*(thetax - self.lens1.pos[0])/r1_2 - self.mu2*(thetax - self.lens2.pos[0])/r2_2
        # self.betay = thetay - self.mu1*(thetay - self.lens1.pos[1])/r1_2 - self.mu2*(thetay - self.lens2.pos[1])/r2_2
        
    def comp_mag_samez(self, thetax, thetay):
        # arcsec2rad = 1
        # x1,y1,x2,y2 = self.lens1.pos[0], self.lens1.pos[1], self.lens2.pos[0], self.lens2.pos[1]
        # x1,y1,x2,y2 = x1*arcsec2rad,y1*arcsec2rad,x2*arcsec2rad,y2*arcsec2rad
        # # print(x1,y1,x2,y2)
        # thetax *= arcsec2rad
        # thetay *= arcsec2rad

        A11, A12, A22 = 1, 0, 1
        for i in range(len(self.xs)):
            # Ri2s.append( (thetax-self.xs[i])**2 + (thetay-self.ys[i])**2 )
            Ri2s = (thetax-self.xs[i])**2 + (thetay-self.ys[i])**2
            # self.betax -= self.masses[i]*(thetax - self.xs[i])/Ri2s
            # self.betay -= self.masses[i]*(thetay - self.ys[i])/Ri2s
            A11 -= self.masses[i]*(1/Ri2s - 2*(thetax - self.xs[i])**2/Ri2s**2)
            A12 -= 2*self.masses[i]*( thetax - self.xs[i] ) * ( thetay - self.ys[i] )/Ri2s**2
            A22 -= self.masses[i]*(1/Ri2s - 2*(thetay - self.ys[i])**2/Ri2s**2)


        # R1 = np.sqrt((thetax-x1)**2+(thetay-y1)**2)
        # R2 = np.sqrt((thetax-x2)**2+(thetay-y2)**2)
        # A11 = 1 - self.mu1*(1/R1**2 - 2*(thetax-x1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetax-x2)**2/R2**4)
        # A12 = -2*self.mu1*(thetax-x1)*(thetay-y1)/R1**4 - 2*self.mu2*(thetax-x2)*(thetay-y2)/R2**4
        # A22 = 1 - self.mu1*(1/R1**2 - 2*(thetay-y1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetay-y2)**2/R2**4)


        detA = (A11*A22 - A12**2)
        self.mag = 1/np.abs(detA)

        print("Maximum Magnification: ",np.abs(self.mag).max())
        print("Minimum Magnification: ",np.abs(self.mag).min())


    def scatter_caustics(self):
        spescatter(self.betax,self.betay,s=0.01, xlabel = r"$\beta_x$/arcsec",ylabel = r"$\beta_y$/arcsec",title="Caustics",issqure=False)
    
    def img_mapping(self, inxdata, inydata, xlim, ylim, ImgSize,valarr = None, datatype=np.float64):
        IMG = np.zeros(ImgSize).astype(datatype)
        # IMG = np.ones(ImgSize)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        xdata = inxdata - xlim[0]
        ydata = inydata - ylim[0]
        # print(xdata.shape) #(250000,)
        # too slow:
        # print((len(xdata)))
        for i in range(len(xdata)):
            # x, y = xdata[i], ydata[i]
            try:
                # IMG[int(y*ratioy), int(x*ratiox)] += 1
                # print("1")
                IMG[int(ydata[i]*ratioy+0.5), int(0.5+xdata[i]*ratiox)] += valarr[i] #valarr
            except:
                pass
        return IMG

    def img_mapping_inone(self, inxdata, inydata, betainxdata, betainydata, xlim, ylim, ImgSize,valarr1 = None,valarr2 = None, datatype=np.float64):
        srcplaneIMG = np.zeros(ImgSize).astype(datatype)
        srcplaneIMG_withoutlens = np.zeros(ImgSize).astype(datatype)
        imgplaneIMG = np.zeros(ImgSize).astype(datatype)
        # IMG = np.ones(ImgSize)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        xdata = inxdata - xlim[0]# thetax
        ydata = inydata - ylim[0]# thetax
        betaxdata = betainxdata - xlim[0]
        betaydata = betainydata - ylim[0]


        self.inputs = [[srcplaneIMG,srcplaneIMG_withoutlens,imgplaneIMG],[((betaydata*ratioy+0.5).astype(np.int),(betaxdata*ratiox+0.5).astype(np.int))
        ,((ydata*ratioy+0.5).astype(np.int),(xdata*ratiox+0.5).astype(np.int)),((ydata*ratioy+0.5).astype(np.int),(xdata*ratiox+0.5).astype(np.int))],[valarr1,valarr2,valarr1]]
        
        pool = ProcessPool(nodes=3)
        results = pool.imap(lambda a: self.compIMGs(a), range(3))
        resultIMGs = list(results)
        srcplaneIMG,srcplaneIMG_withoutlens,imgplaneIMG = resultIMGs[0],resultIMGs[1],resultIMGs[2]


        srcplaneIMG /= srcplaneIMG_withoutlens
        return srcplaneIMG, imgplaneIMG, srcplaneIMG_withoutlens

    def compIMGs(self, a):  
        print("A new core for calculating images is started")
        resimg = self.inputs[0][a]
        y, x = self.inputs[1][a][0], self.inputs[1][a][1]
        val = self.inputs[2][a]
        for i in range(len(x)):
            try:
                resimg[y[i],x[i]] += val[i]
            except:
                pass
        return resimg





    def com_lightcurve(self,srcIMG, inX, inY, xlim, ylim, radius, npoints):
        # radius in pixel
        print("Generating lightcurve")
        ImgSize = srcIMG.shape
        lc = np.zeros(len(inX),)
        iny = inY - ylim[0]
        inx = inX - xlim[0]
        angs = np.linspace(0,2*np.pi,npoints)
        Radius = np.linspace(0,radius,npoints)

        cosangs, sinangs = np.cos(angs), np.sin(angs)
        rcos, rsin = np.array([],), np.array([],)
        for r in Radius:
            rcos, rsin = np.append(rcos, r*cosangs), np.append(rsin, r*sinangs)
        NUM = len(rcos)
        # print(NUM)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        inx, iny = inx*ratiox, iny*ratioy
        
        for i in range(len(inx)):
            Xs, Ys = (inx[i]+ rcos + 0.5).astype(np.int), (iny[i] + rsin + 0.5).astype(np.int)
            for x,y in zip(Xs, Ys):
                try:
                    lc[i] += srcIMG[y,x]
                except:
                    pass
            lc[i] /= NUM
        return lc

    def get_imgs_lessmem_v2(self, ImgSize, xlim, ylim, num,datatype=np.float64):
        srcplaneIMG = np.zeros(ImgSize).astype(datatype)
        srcplaneIMG_withoutlens = np.zeros(ImgSize).astype(datatype)
        imgplaneIMG = np.zeros(ImgSize).astype(datatype)
        # IMG = np.ones(ImgSize)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        incx = (xlim[1]-xlim[0])/num
        incy = (ylim[1]-ylim[0])/num


        with tqdm(total=num) as pbar:
            for thetax in np.linspace(xlim[0], xlim[1], num).astype(datatype):
                thetay = np.linspace(ylim[0],ylim[1],num).astype(datatype)
                thetax = np.ones(thetay.shape)*thetax
                betax, betay, MAG = self.ray_shoot_comp_mag_in_one(thetax, thetay)
                # print(mag, betax, betay)
                xdata = thetax - xlim[0]# thetax
                ydata = thetay - ylim[0]# thetax
                betaxdata = betax - xlim[0]
                betaydata = betay - ylim[0]
                # try:
                BJ = (betaydata*ratioy+0.5).astype(np.int)
                BI = (betaxdata*ratiox+0.5).astype(np.int)
                J = (ydata*ratioy+0.5).astype(np.int)
                I = (xdata*ratiox+0.5).astype(np.int)
                for idx in range(num):
                    bi, bj, i, j, mag = BI[idx], BJ[idx], I[idx], J[idx], MAG[idx]
                    if ((bi>=0 and bi<ImgSize[1]) and (bj>=0 and bj<ImgSize[0])):
                        srcplaneIMG[bj, bi] += mag
                    if ((i>=0 and i<ImgSize[1]) and (j>=0 and j<ImgSize[0])):
                        srcplaneIMG_withoutlens[j,i] += 1
                        imgplaneIMG[j,i] += mag
                pbar.update(1)
        srcplaneIMG /= srcplaneIMG_withoutlens
        return srcplaneIMG, imgplaneIMG, srcplaneIMG_withoutlens

    def get_imgs_lessmem(self, ImgSize, xlim, ylim, num,datatype=np.float64):
        srcplaneIMG = np.zeros(ImgSize).astype(datatype)
        srcplaneIMG_withoutlens = np.zeros(ImgSize).astype(datatype)
        imgplaneIMG = np.zeros(ImgSize).astype(datatype)
        # IMG = np.ones(ImgSize)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        incx = (xlim[1]-xlim[0])/num
        incy = (ylim[1]-ylim[0])/num


        with tqdm(total=num**2) as pbar:
            for thetax in np.linspace(xlim[0], xlim[1], num).astype(datatype):
                for thetay in np.linspace(ylim[0], ylim[1], num).astype(datatype):
                # for ix in range(num):
                #     thetax = xlim[0]+ix*incx
                #     for iy in range(num):
                #         thetay = ylim[0]+iy*incy
                        betax, betay, mag = self.ray_shoot_comp_mag_in_one(thetax, thetay)
                        # print(mag, betax, betay)
                        xdata = thetax - xlim[0]# thetax
                        ydata = thetay - ylim[0]# thetax
                        betaxdata = betax - xlim[0]
                        betaydata = betay - ylim[0]
                        # try:
                        bj = int(betaydata*ratioy+0.5)
                        bi = int(betaxdata*ratiox+0.5)
                        j = int(ydata*ratioy+0.5)
                        i = int(xdata*ratiox+0.5)
                        if ((bi>=0 and bi<ImgSize[1]) and (bj>=0 and bj<ImgSize[0])):
                            srcplaneIMG[bj, bi] += mag
                        if ((i>=0 and i<ImgSize[1]) and (j>=0 and j<ImgSize[0])):
                            srcplaneIMG_withoutlens[j,i] += 1
                            imgplaneIMG[j,i] += mag
                        # try:
                        #     srcplaneIMG[bj, bi] += mag
                        # except:
                        #     pass
                        # try:
                        #     srcplaneIMG_withoutlens[j,i] += 1
                        #     imgplaneIMG[j,i] += mag
                        # except:
                        #     pass
                        # except:
                        #     pass
                        pbar.update(1)

        srcplaneIMG /= srcplaneIMG_withoutlens
        return srcplaneIMG, imgplaneIMG, srcplaneIMG_withoutlens

    def ray_shoot_comp_mag_in_one(self, thetax, thetay):
        A11=1
        A12=0
        A22=1
        betax = thetax.copy()
        betay = thetay.copy()
        for i in range(len(self.masses)):
            Ri2 = (thetax - self.xs[i])**2 + (thetay - self.ys[i])**2
            betax -= self.masses[i]*(thetax - self.xs[i])/Ri2
            betay -= self.masses[i]*(thetay - self.ys[i])/Ri2
            # Ah Ah Ah, such a bug!!!
            # A11 -= self.masses[i]*( 1/Ri2 - 2*(thetax - self.xs[i])**2 )/Ri2**2
            # A22 -= self.masses[i]*( 1/Ri2 - 2*(thetay - self.ys[i])**2 )/Ri2**2
            A11 -= self.masses[i]*( 1/Ri2 - 2*(thetax - self.xs[i])**2/Ri2**2)
            A22 -= self.masses[i]*( 1/Ri2 - 2*(thetay - self.ys[i])**2 /Ri2**2)
            A12 -= 2*self.masses[i]*(thetax - self.xs[i])*(thetay - self.ys[i])/Ri2**2
        detA = A11*A22 - A12**2
        # print("detA: ", abs(detA))
        mag = 1/abs(detA)
        # print("mag: ", mag)

        # A11, A12, A22 = 1, 0, 1
        # for i in range(len(self.xs)):
        #     # Ri2s.append( (thetax-self.xs[i])**2 + (thetay-self.ys[i])**2 )
        #     Ri2s = (thetax-self.xs[i])**2 + (thetay-self.ys[i])**2
        #     # self.betax -= self.masses[i]*(thetax - self.xs[i])/Ri2s
        #     # self.betay -= self.masses[i]*(thetay - self.ys[i])/Ri2s
        #     A11 -= self.masses[i]*(1/Ri2s - 2*(thetax - self.xs[i])**2/Ri2s**2)
        #     A12 -= 2*self.masses[i]*( thetax - self.xs[i] ) * ( thetay - self.ys[i] )/Ri2s**2
        #     A22 -= self.masses[i]*(1/Ri2s - 2*(thetay - self.ys[i])**2/Ri2s**2)

        # R1 = np.sqrt((thetax-x1)**2+(thetay-y1)**2)
        # R2 = np.sqrt((thetax-x2)**2+(thetay-y2)**2)
        # A11 = 1 - self.mu1*(1/R1**2 - 2*(thetax-x1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetax-x2)**2/R2**4)
        # A12 = -2*self.mu1*(thetax-x1)*(thetay-y1)/R1**4 - 2*self.mu2*(thetax-x2)*(thetay-y2)/R2**4
        # A22 = 1 - self.mu1*(1/R1**2 - 2*(thetay-y1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetay-y2)**2/R2**4)


        # detA = (A11*A22 - A12**2)
        # self.mag = 1/np.abs(detA)

        return betax, betay, mag

        # float mag = this -> ray_shoot_comp_mag_in_one(thetax, thetay, betaxy);


def genimgs(xylim = 6, q = 0.5, posscale = 0.5, imgsz = 2048, num = 2000, method = "slow"):
    xlim, ylim = (-xylim,xylim), (-xylim,xylim)
    masses = [q/(1+q),1/(1+q)]
    xs = [1, -1]
    ys = [0,0]
    num = int(num)
    masses = np.array(masses)
    masses = masses/np.sum(masses)
    xs = np.array(xs)*posscale
    ys = np.array(ys)*posscale
    ImgSize = (imgsz,imgsz) # raw, colume
    datatype = np.float64
    # generate light rays in the image plane and creat the two lens system and run the inverse ray shooting code
    twolens = Nlenses(masses, xs, ys)
    if method == "fast":
        thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num, datatype = datatype)
        twolens.inverse_ray_shooting(thetax, thetay)
        twolens.comp_mag_samez(thetax, thetay)
        t0 = time.time()
        srcplaneIMG, imgplaneIMG, srcplaneIMG_withoutlens = twolens.img_mapping_inone(thetax, thetay,twolens.betax, twolens.betay, xlim, ylim, ImgSize,valarr1 = twolens.mag, valarr2 = np.ones([len(thetax),]), datatype=datatype)
        t1 = time.time()
    elif method == "slow":
        t0 = time.time()
        # srcplaneIMG, imgplaneIMG = twolens.get_imgs_lessmem(ImgSize, xlim, ylim, num, datatype = np.float64)
        srcplaneIMG, imgplaneIMG, srcplaneIMG_withoutlens = twolens.get_imgs_lessmem_v2(ImgSize, xlim, ylim, num, datatype = np.float64)
        t1 = time.time()
    # # mid fast method   
    else:
        thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num, datatype = datatype)
        twolens.inverse_ray_shooting(thetax, thetay)
        twolens.comp_mag_samez(thetax, thetay)
        t0 = time.time()
        srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize, valarr = twolens.mag) #valarr = np.ones([len(twolens.betax),])
        srcplaneIMG_withoutlens = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = np.ones([len(thetax),]))
        srcplaneIMG /= srcplaneIMG_withoutlens
        imgplaneIMG = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = twolens.mag)
        t1 = time.time()
    # print("time spent on inverse ray shooting: {}".format(t1-t0))
    return srcplaneIMG, imgplaneIMG

def com_lightcurve_v2(srcIMG, inX, inY, xlim, ylim, radius):
    # radius in pixel
    print("Generating lightcurve")
    ImgSize = srcIMG.shape
    lc = np.zeros(len(inX),)
    lc_withoutmag = np.zeros(len(inX),)
    iny = inY - ylim[0]
    inx = inX - xlim[0]
    ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
    inx, iny = inx*ratiox, iny*ratioy
    
    for i in range(len(inx)):
        # fot testx in range()
        # Xs, Ys = (inx[i]+ rcos + 0.5).astype(np.int), (iny[i] + rsin + 0.5).astype(np.int)
        basex, basey = (inx[i] + 0.5).astype(np.int), (iny[i] + 0.5).astype(np.int)
        # Xs, Ys = [],[]
        # for dx in range(radius):
        #     for dy in range(radius):
        # Xs = np.linspace(basex-radius, basex+radius, ).astype(np.int)
        Xs = np.arange(basex-radius, basex+radius+1, 1).astype(np.int)
        Ys = np.arange(basey-radius, basey+radius+1, 1).astype(np.int)
        num=0
        for x,y in zip(Xs, Ys):
            if (x-basex)**2+(y-basey)**2<=radius**2:
                try:
                    lc[i] += srcIMG[y,x]
                    num+=1
                    # lc_withoutmag[i] += srcIMG_nomag[y,x]
                except:
                    pass
        # lc[i] /= lc_withoutmag[i]
        lc[i] /= num
    return lc



if __name__ == '__main__':
    # plt.style.use('classic')
    # py3 .py 4.5 -0.2 1 2048 5000 0.35 fast/slow
    import sys
    xylim = float(sys.argv[1])
    # Nrows = int(sys.argv[2])
    q = float(sys.argv[2])
    posscale = float(sys.argv[3])
    imgsz = int(sys.argv[4])
    num = int(sys.argv[5])
    kscale = float(sys.argv[6])
    method = str(sys.argv[7])


    stupstr = "{}pms_{}q_{}X_{}rays".format(2,q,posscale,num)

    # python3 *.py 4.5 5 0.5 512 1000



    # xlim, ylim = (-4.5,4.5), (-4.5,4.5)
    xlim, ylim = (-xylim,xylim), (-xylim,xylim)
    ImgSize = (imgsz,imgsz) # raw, colume

    srcplaneIMG, imgplaneIMG =  genimgs(xylim = xylim, q = q, posscale = posscale, imgsz = imgsz, num = num, method = "slow")
    srcplaneIMG2, imgplaneIMG2 =  genimgs(xylim = xylim, q = -q, posscale = posscale, imgsz = imgsz, num = num, method = "slow")



    # generate light curve
    sourcesizeR_E = 0.03 # source size 0.03 R_E
    radius = int(sourcesizeR_E*ImgSize[1]/(xlim[1]-xlim[0]) + 0.5)
    k = 0.5
    kscale = 0.6 # source position for generate light curve
    B = np.linspace( 1, -1 , 3)*kscale
    KB = [(k,b) for b in B]
    cname = ["r","g","y","k","c","m","w"]
    COLLOR = cname[:len(KB)]
    Y, LC, LC2 = {}, {}, {}
    scale = 0.8
    # x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/2))#, endpoint=False
    x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/radius*3))#, endpoint=False
    for kb in KB:
        k,b = kb[0], kb[1]
        Y[kb] = k*x + b
    #     LC[kb] = twolens.com_lightcurve(srcplaneIMG, x, Y[kb], xlim, ylim, radius, npoints)
        LC[kb] = com_lightcurve_v2(srcplaneIMG, x, Y[kb], xlim, ylim, radius)
        LC2[kb] = com_lightcurve_v2(srcplaneIMG2, x, Y[kb], xlim, ylim, radius)





    '''
    # # plot critical lines and caustics
    timestr=time.ctime().replace(" ","")[3:-4]
    fontsize = 18
    cmap = plt.cm.get_cmap('inferno')#viridis, plasma
    font1 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : fontsize,
    }
    fig3 = plt.figure()
    fig3.set_size_inches(16,16)
    plt.subplots_adjust(left=0.04, top = 0.93, bottom = 0.05, right=0.96, hspace = 0, wspace = 0)
    plt.margins(0,0)

    grid = plt.GridSpec(12, 3, wspace=0.2, hspace=0)#, wspace=0.5, hspace=0.5
    plt.subplot(grid[0:3,0])
    # plt.subplot(131)

    imgplaneim = plt.imshow(np.log10(imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(imgplaneim)

    plt.subplot(grid[3:6,0])
    # plt.subplot(132)
    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    srcplaneim = plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(srcplaneim)

    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt],linewidth=0.8)#-ylim[0]
        cnt += 1


    #light curve
    axes = [1 for i in range(len(KB)*2)]
    cnt = 0

    for kb in KB:
    #     axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
        axes[cnt] = plt.subplot(grid[cnt*2:cnt*2+2,1:3])
    #     plt.tick_params(labelsize=fontsize)
        label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
        axes[cnt].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
        axes[cnt].plot([x[0],x[-1]],[0,0],linewidth=0.5,c="k")
        axes[cnt].plot([x[0],x[-1]],[-1,-1],linewidth=0.5,c="k")
        axes[cnt].set_ylim([-3, 5])
        axes[cnt].set_ylabel(r"$\log$ $\mu$")#
        
        axes[cnt].legend()#prop=font1

        cnt += 1

    ##########################################################################################################    

    plt.subplot(grid[6:9,0])
    # plt.subplot(131)
    # cmap = plt.cm.get_cmap('viridis')
    imgplaneim = plt.imshow(np.log10(imgplaneIMG2), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(imgplaneim)

    plt.subplot(grid[9:12,0])
    # plt.subplot(132)
    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    srcplaneim = plt.imshow(np.log10(srcplaneIMG2), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(srcplaneim)

    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt],linewidth=0.8)#-ylim[0]
        cnt += 1

    cnt = len(KB)

    for kb in KB:
    #     axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
        axes[cnt] = plt.subplot(grid[cnt*2:cnt*2+2,1:3])
    #     plt.tick_params(labelsize=fontsize)
        label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
        
        axes[cnt].plot(x,np.log10(LC2[kb]),color=COLLOR[cnt%len(KB)],label=label)
        axes[cnt].plot([x[0],x[-1]],[0,0],linewidth=0.5,c="k")
        axes[cnt].plot([x[0],x[-1]],[-1,-1],linewidth=0.5,c="k")
        axes[cnt].set_ylim([-3, 5])
        axes[cnt].set_ylabel(r"$\log$ $\mu$")#
          
        
        axes[cnt].legend()#prop=font1

        cnt += 1
    axes[cnt-1].set_xlabel(r"$t/t_E$",fontsize = fontsize)
    # fig3.suptitle('Main title') # or plt.suptitle('Main title')
    title = "Critical lines and Caustics (left) and Light curves (right) of 2 point mass lenses,\n"+r"$q = ${:.2f} (upper) and $q = ${:.2f} (bottom), $2X = ${:.2f}".format(q, -q,2*posscale)
    plt.suptitle(title, size=20)
    plt.subplots_adjust(hspace=.0)
    lightcurve_filename = "./resimgs/lightcurve_{}_{}.png".format(timestr, stupstr)

    '''

    # # plot critical lines and caustics
    timestr=time.ctime().replace(" ","")[3:-4]
    fontsize = 18
    cmap = plt.cm.get_cmap('inferno')#viridis, plasma
    font1 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : fontsize,
    }
    fig3 = plt.figure()
    fig3.set_size_inches(20,12)
    plt.subplots_adjust(left=0.04, top = 0.92, bottom = 0.05, right=0.96, hspace = 0, wspace = 0)
    # plt.subplots_adjust(hspace=.0)
    plt.margins(0,0)

    grid = plt.GridSpec(6, 3, wspace=0.15, hspace=0)#, wspace=0.5, hspace=0.5
    plt.subplot(grid[0:3,0])
    # plt.subplot(131)

    imgplaneim = plt.imshow(np.log10(imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])

    plt.subplot(grid[0:3,1])
    # plt.subplot(132)
    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    srcplaneim = plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(srcplaneim)
    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt],linewidth=0.8)#-ylim[0]
        cnt += 1
    #light curve
    axes = [1 for i in range(len(KB)*2)]
    cnt = 0

    for kb in KB:
    #     axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
        axes[cnt] = plt.subplot(grid[cnt,2])
    #     plt.tick_params(labelsize=fontsize)
        label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
        axes[cnt].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
        axes[cnt].plot([x[0],x[-1]],[0,0],linewidth=0.5,c="k")
        axes[cnt].plot([x[0],x[-1]],[-1,-1],linewidth=0.5,c="k")
        axes[cnt].set_ylim([-3, 5])
        axes[cnt].set_ylabel(r"$\log$ $\mu$")#
        
        
        axes[cnt].legend()#prop=font1
    #     axes[cnt].set_xlabel(r"$t/t_E$",fontsize = fontsize)
        
        if cnt == 0:
            pass
        cnt += 1

    ##########################################################################################################    

    plt.subplot(grid[3:6,0])
    # plt.subplot(131)
    # cmap = plt.cm.get_cmap('viridis')
    imgplaneim = plt.imshow(np.log10(imgplaneIMG2), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    plt.subplot(grid[3:6,1])
    # plt.subplot(132)
    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    srcplaneim = plt.imshow(np.log10(srcplaneIMG2), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(srcplaneim)

    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt],linewidth=0.8)#-ylim[0]
        cnt += 1


    #light curve
    # axes = [1 for i in range(len(KB))]
    cnt = len(KB)

    for kb in KB:
    #     axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
        axes[cnt] = plt.subplot(grid[cnt,2])
    #     plt.tick_params(labelsize=fontsize)
        label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
        
        axes[cnt].plot(x,np.log10(LC2[kb]),color=COLLOR[cnt%len(KB)],label=label)
        axes[cnt].plot([x[0],x[-1]],[0,0],c="k",linewidth=0.5)
        axes[cnt].plot([x[0],x[-1]],[-1,-1],c="k",linewidth=0.5)
        axes[cnt].set_ylim([-3, 5])
        axes[cnt].set_ylabel(r"$\log$ $\mu$")#
 
        
        axes[cnt].legend()#prop=font1
    #     axes[cnt].set_xlabel(r"$t/t_E$",fontsize = fontsize)
          
        cnt += 1
    axes[cnt-1].set_xlabel(r"$t/t_E$",fontsize = fontsize)
    # fig3.suptitle('Main title') # or plt.suptitle('Main title')
    title = "Critical lines (left) and Caustics (mid) and Light curves (right) of 2 point mass lenses,\n"+r"$q = ${:.2f} (upper) and $q = ${:.2f} (bottom), $2X = ${:.2f}".format(q, -q,2*posscale)
    plt.suptitle(title, size=18)
    lightcurve_filename = "./resimgs/lightcurve_{}_{}.png".format(timestr, stupstr)
    pkl_filename = "./resimgs/lightcurve_{}_{}.pkl".format(timestr, stupstr)


    fig3.savefig(lightcurve_filename, format='png', bbox_inches='tight', dpi=100, pad_inches = 0)

    HRlightcurve_filename = "./resimgs/HRlightcurve_{}_{}.png".format(timestr, stupstr)
    fig3.savefig(HRlightcurve_filename, format='png', bbox_inches='tight', dpi=300, pad_inches = 0)
    dill.dump_session(pkl_filename)