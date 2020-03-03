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


    def com_lightcurve_v2(self,srcIMG,srcIMG_nomag, inX, inY, xlim, ylim, radius):
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

    # masses = []
    # xs = []
    # ys = []
    # # Nrows = 5
    # minx = abs(Nrows//2)
    # for i in range(Nrows):
    #         masses += [(-1)**(i+k) for k in range(Nrows)]
    #         # masses += [(1)**(i+k) for k in range(Nrows)]
    #         xs += [k for k in range(-minx,Nrows-minx)]
    #         ys += [i-minx for k in range(Nrows)]
    # # posscale = 0.5 # 0.3 --> 1

    
    # masses = [-0.2, 1]
    
    masses = [q/(1+q),1/(1+q)]
    xs = [1, -1]
    ys = [0,0]
    # posscale = 0.35
    # X = 1, raynum = 4000(or 20000, the pattern is different)


    masses = np.array(masses)
    masses = masses/np.sum(masses)
    xs = np.array(xs)*posscale
    ys = np.array(ys)*posscale
    # print(masses, xs, ys)

    # srcxlim, srcylim = (-4,4), (-4, 4)
    # ImgSize = (512,512) # raw, colume
    ImgSize = (imgsz,imgsz) # raw, colume
    import time
    # num = 1000
    datatype = np.float64#np.single = np.float32
    
    twolens = Nlenses(masses, xs, ys)
    # print(twolens.masses, twolens.xs, twolens.ys)

    #https://blog.csdn.net/zj360202/article/details/78543141, struct.error: 'i' format requires -2147483648
    # File "/usr/lib/python2.7/pickle.py", line 494, in save_string
    # self.write(BINSTRING + pack("<i", n) + obj)
    # self.write(BINSTRING + pack("<i", n) + obj)  -->  self.write(BINSTRING + pack("<q", n) + obj)

    # # faster method
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
    print("time spent on inverse ray shooting: {}".format(t1-t0))

    # generate light curve
    # radius, npoints = 5, 40 # radius in pixel
    sourcesizeR_E = 0.03 # source size 0.03 R_E
    radius = int(sourcesizeR_E*ImgSize[1]/(xlim[1]-xlim[0]) + 0.5)
    # print(radius)# 12
    # radius = 15# radius in pixel
    k = 0.5
    # kscale = 0.6 # source position for generate light curve
    print("kscale",kscale)
    B = np.linspace( 1, -1 , 5)*kscale
    KB = [(k,b) for b in B]
    cname = ["r","g","y","k","c","m","w"]
    COLLOR = cname[:len(KB)]
    Y, LC = {}, {}
    scale = 0.7
    # x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/2))#, endpoint=False
    x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/radius*3))#, endpoint=False
    for kb in KB:
        k,b = kb[0], kb[1]
        Y[kb] = k*x + b
    #     LC[kb] = twolens.com_lightcurve(srcplaneIMG, x, Y[kb], xlim, ylim, radius, npoints)
        LC[kb] = twolens.com_lightcurve_v2(srcplaneIMG, srcplaneIMG_withoutlens, x, Y[kb], xlim, ylim, radius)

    # # plot critical lines and caustics
    timestr=time.ctime().replace(" ","")[3:-4]
    fontsize = 18
    font1 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : fontsize,
    }
    fig3 = plt.figure()
    fig3.set_size_inches(16,9)
    plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.2, wspace = 0.2)
    plt.margins(0,0)

    grid = plt.GridSpec(10, 3, wspace=0.2, hspace=0)#, wspace=0.5, hspace=0.5
    plt.subplot(grid[0:5,0])
    # plt.subplot(131)
    cmap = plt.cm.get_cmap('viridis')
    imgplaneim = plt.imshow(np.log10(imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(imgplaneim)
    title = "Critical lines of {} point masses, q={:.2f}, X = {:.2f}, (Log scale)".format(len(masses),q,posscale)
    # plt.title(title)
    plt.xlabel(r"$\xi_1(D_d\theta_E)$", font1)
    plt.ylabel(r"$\xi_2(D_d\theta_E)$", font1)
    # plt.tick_params(labelsize=fontsize)

    plt.subplot(grid[5:10,0])
    # plt.subplot(132)
    cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    srcplaneim = plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.colorbar()
    add_colorbar(srcplaneim)
    title = "Caustics of {} point masses(blue/red - smaller/larger than 0),\n, q={:.2f}, X = {:.2f}, (Log scale)".format(len(masses),q,posscale)
    for m, a, b in zip(masses , xs, ys):  
        if m < 0:
            c = "b"
        else:
            c = "r"
        plt.scatter(a,b,color=c,marker="x",s=0.1)
    # plt.title(title)

    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt],linewidth=0.8)#-ylim[0]
        cnt += 1
    plt.xlabel(r"$\eta_1(D_s\theta_E)$", font1)
    plt.ylabel(r"$\eta_2(D_s\theta_E)$", font1)
    # plt.tick_params(labelsize=fontsize)


    #light curve
    axes = [1 for i in range(len(KB))]
    cnt = 0

    for kb in KB:
    #     axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
        axes[cnt] = plt.subplot(grid[cnt*2:cnt*2+2,1:3])
    #     plt.tick_params(labelsize=fontsize)
        label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
        
    #     axes[cnt].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
    #     axes[cnt].set_ylim([-3, 4])
    #     axes[cnt].set_ylabel(r"$\log$ $\mu$")#
        
        axes[cnt].plot(x,(LC[kb]),color=COLLOR[cnt],label=label)
        axes[cnt].plot([x[0],x[-1]],[1,1],linewidth=0.5)
        axes[cnt].set_ylim([-1, 10])
        axes[cnt].set_ylabel(r"$\mu$")#,fontsize = fontsize
        
        
        axes[cnt].legend()#prop=font1
        axes[cnt].set_xlabel(r"$t/t_E$",fontsize = fontsize)
        
        if cnt == 0:
    #         pass
            title = r"Light Curve of a source star with radius {:.2f} $R_E$".format(sourcesizeR_E)+"\n"+r"2 point mass ratio: $q =$ {:.2f}, separation: $2X =$ {:.2f}".format(q,2*posscale)
            plt.title(title,fontsize=fontsize)
        
        cnt += 1
    plt.subplots_adjust(hspace=.0)
    lightcurve_filename = "./resimgs/lightcurve_{}_{}.png".format(timestr, stupstr)
    fig3.savefig(lightcurve_filename, format='png', bbox_inches='tight', dpi=100, pad_inches = 0)
    # twolens = Twolenses(lens1, lens2)
    # posscale = 0.2
    # lenses = [(1, -1*posscale, 1*posscale), (-1, 0*posscale, 1*posscale), (1, 1*posscale, 1*posscale),
    # (-1, -1*posscale, 0*posscale), (1, 0*posscale, 0*posscale), (-1, 1*posscale, 0*posscale), 
    # (1, -1*posscale, -1*posscale), (-1, 0*posscale, -1*posscale), (1, 1*posscale, -1*posscale)
    #  ]

    # masses = [1,-1,1,-1,1,-1,1,-1,1]# masses = [1,1,1,1,3,1,1,1,1]
    # xs = [-1,0,1,-1,0,1,-1,0,1]
    # ys = [1,1,1,0,0,0,-1,-1,-1]
    # posscale = 0.35 # 0.3 --> 1




    # # masses = [1,-0.3,1,-0.3,1,-0.3,1,-0.3,1]
    # masses = [1,0.3,1,0.3,1,0.3,1,0.3,1]
    # xs = [-1,0,1,-1,0,1,-1,0,1]
    # ys = [1,1,1,0,0,0,-1,-1,-1]
    # posscale = 0.85 # 0.3 --> 1

    # masses = [1,1]
    # xs = [-1, 1]
    # ys = [0,0]

    # N = 10
    # masses = np.random.uniform(0,1,N)
    # xs = np.random.uniform(-1,1,N)
    # ys = np.random.uniform(-1,1,N)
    # posscale = 1
    # np.random.shuffle(xs)
    # np.random.shuffle(ys)
    # np.random.shuffle(masses)

    

    # srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize, valarr = twolens.mag) #valarr = np.ones([len(twolens.betax),])
    # srcplaneIMG_withoutlens = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = np.ones([len(thetax),]))
    # srcplaneIMG /= srcplaneIMG_withoutlens
    # # srcplaneIMG += 1
    # # srcplaneIMG = np.log10(srcplaneIMG)
    # imgplaneIMG = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = twolens.mag)

    # where_smaller_than1 = np.zeros(srcplaneIMG.shape)
    # where_smaller_than1[:,:] = srcplaneIMG[:,:]
    # where_smaller_than1[where_smaller_than1>1] = 0

    # where_smaller_than1 = srcplaneIMG
    # where_smaller_than1[where_smaller_than1>1] = 0

    # srcplaneIMG += 1
    # print("time spent on img_mapping_inone: {}".format(t1-t0))
    # print("thetax datatype: ", thetax.dtype)
    # print("betax datatype: ", twolens.betax.dtype)
    # print("srcplaneIMG datatype: ", srcplaneIMG.dtype)
    # print("imgplaneIMG datatype: ", imgplaneIMG.dtype)

    # imgplaneIMG = np.log10(imgplaneIMG)
    # imgplaneIMG /= np.min(imgplaneIMG)
    
    #https://stackoverflow.com/questions/53101815/improve-color-contrast-in-matplotlib
    # try one of the continuous colormaps, like viridis or gist_ncar. 
    # If you don't need all colors to be different, 
    # try one of the repeating maps like flag or prism
    
    '''                
                    # fig = plt.figure()
                    # plt.subplot(122)
                    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
                    # plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
                    # #title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
                    # # title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
                    # plt.colorbar()
                    # # plt.scatter(xs, ys,  c='r',edgecolors='blue',marker='o',s=0)  
                    # # 1，plt.text(横坐标，纵坐标，‘显示文字’)
                    # # 2,  plt.annotate('文字',xy=(箭头坐标),xytext=(文字坐标),arrowprops=dict(facecolor='箭头颜色'))
                    # # for m, a, b in zip(masses , xs, ys):  
                    # #     # plt.text(a, b, "({:.2f},{:.2f},{:.2f})".format(a,b,m),ha='center', va='bottom', fontsize=5) 
                    # #     plt.text(a, b, "{:.2f}".format(m),ha='center', va='bottom', fontsize=4) 
                    # # title = "Caustics of {} point masses, (Log scale)".format(len(masses))
                    # title = "Caustics of {} point masses, q={:.2f}, X = {:.2f}, (Log scale)".format(len(masses),q,posscale)
                    # plt.title(title)
                    

                    # # fig = plt.figure()
                    # plt.subplot(121)
                    # cmap = plt.cm.get_cmap('viridis')
                    # plt.imshow(np.log10(imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
                    # #title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
                    # # title = "Two point mass lenses, img plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
                    # # title = "Critical lines of {} point masses (Log scale)".format(len(masses))
                    # plt.colorbar()
                    # title = "Critical lines of {} point masses, q={:.2f}, X = {:.2f}, (Log scale)".format(len(masses),q,posscale)
                    # plt.title(title)
                    

                    # # lensplaneIMG = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize)
                    # # if lensplaneIMG.min() == 0:
                    # #     lensplaneIMG[lensplaneIMG==0]=1
                    # # CausticsIMG = srcplaneIMG/lensplaneIMG

                    # # plt.figure()
                    # # plt.imshow(CausticsIMG, origin='lower',cmap=cmap)
                    # # plt.colorbar()

                    # # matplotlib.animation.FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)
                    # # ani.FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)

                    # # plt.tight_layout()

                    # fig.set_size_inches(16,9)


                    # # plt.subplots_adjust(top = 2, bottom = 2, right = 1, left = 1, hspace = 0, wspace = 0)
                    # # 有六个可选参数来控制子图布局。值均为0~1之间。其中left、bottom、right、top围成的区域就是子图的区域。wspace、hspace分别表示子图之间左右、上下的间距
                    # plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
                    # plt.margins(0,0)
                    # # caustics_filename = "./resimgs/neg_mass2point_samez/massratio_{:.3f}X1{:.2f}X2{:.2f}_caustics.png".format(massratio,lens1.pos[0],lens2.pos[0])
                    # # fig.savefig(caustics_filename, format='png', bbox_inches='tight', dpi=900, pad_inches = 0)#, transparent=True

                    # timestr=time.ctime().replace(" ","")[3:-4]
                    # caustics_filename = "./resimgs/foranndy/Caustics_{}_{}.png".format(timestr, stupstr)
                    # # fig.savefig(caustics_filename, format='png', bbox_inches='tight', dpi=600, pad_inches = 0)#, transparent=True

                    # # k,b = 0,0
                    # # scale = 0.5
                    
                    # # x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale))#ImgSize[1]
                    # # y = k*x + b 
                    # # lc = twolens.com_lightcurve(srcplaneIMG, x, y, xlim, ylim, )
                    # # plt.figure()
                    # # plt.plot(x, np.log10(lc))
                    # # print(np.log10(lc)[:10])
                    # # # plt.ylim((0,5))
                    # # plt.ylabel("Log Magnification")
                    # # plt.show()

                    # radius, npoints = 5, 40 # radius in pixel
                    # k = 0.5
                    # # B = np.array([0.5,0.4,0.3,0.2,0.1,0.01]) - 0.2
                    # # B = np.linspace( -minx, Nrows-minx , 5)*kscale
                    # B = np.linspace( -1, 1 , 5)*kscale
                    # KB = [(k,b) for b in B]
                    # #https://www.cnblogs.com/darkknightzh/p/6117528.html
                    # #b: blue g: green r: red c: cyan m: magenta y: yellow k: black w: white
                    # # COLLOR = ["r","g","y","k","c","m","w"]

                    # # cname = []
                    # # import matplotlib
                    # # for name, _ in matplotlib.colors.cnames.items():
                    #     # cname.append(name)

                    # # import matplotlib.colors as colors
                    # # import random
                    # # cname = list(colors._colors_full_map.keys())
                    # # random.shuffle(cname)

                    # cname = ["r","g","y","k","c","m","w"]
                    # COLLOR = cname[:len(KB)]

                    # # LABEL = ["label1","label2","label3","label4"]
                    # Y, LC = {}, {}
                    # scale = 0.9
                    # # print(xlim,"xlim")
                    # # print(ylim,"ylim")
                    # x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/2))#, endpoint=False
                    # # print(x)
                    # for kb in KB:
                    #     k,b = kb[0], kb[1]
                    #     Y[kb] = k*x + b
                    #     LC[kb] = twolens.com_lightcurve(srcplaneIMG, x, Y[kb], xlim, ylim, radius, npoints)
                    #     # print(x)

                    # fig = plt.figure()
                    # # fig, axs =  plt.subplots(len(KB), 2)
                    # # https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.subplots.html
                    # # fig, axes = plt.subplots(4, 2, sharex='col')
                    # fig.set_size_inches(16,9)
                    # plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
                    # plt.margins(0,0)

                    # plt.subplot(121)

                    # # plt.title("Source Place")
                    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
                    # plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
                    # # axes[0,0].imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
                    # # title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
                    # # title = "Caustics of {} point masses (Log scale)".format(len(masses))
                    # plt.colorbar()
                    # title = "Caustics of {} point masses(blue/red - smaller/larger than 0),\n, q={:.2f}, X = {:.2f}, (Log scale)".format(len(masses),q,posscale)
                    # for m, a, b in zip(masses , xs, ys):  
                    #     if m < 0:
                    #         c = "b"
                    #     else:
                    #         c = "r"
                    #     plt.scatter(a,b,color=c,marker="x")
                    # plt.title(title)
                    # # plt.colorbar()
                    # # plt.hold(True)
                    # cnt = 0
                    # for kb in KB:
                    #     plt.plot(x,Y[kb],color=COLLOR[cnt])#-ylim[0]
                    #     cnt += 1


                    # # plt.subplot(122)

                    # # plt.subplots(122, sharex=True)
                    # axes = [1 for i in range(len(KB))]
                    # cnt = 0
                    # # plt.subplot(len(KB),2)
                    # for kb in KB:
                    #     axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
                    #     label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
                    #     axes[cnt].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
                    #     # axes[cnt].ylim(0,5)
                    #     axes[cnt].set_ylim([-4, 5])
                    #     axes[cnt].legend()
                    #     if cnt == 0:
                    #         plt.title("Light Curve (Log scale)")
                    #     # axes[cnt,1].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
                    #     cnt += 1
                    # # axes[0].get_shared_x_axes().join(axes[0], axes[1], axes[2],axes[3])
                    # # axes[0].set_xticklabels([])
                    # plt.subplots_adjust(hspace=.0)
                    # # plt.subplot(122)
                    
                    
                    # # https://stackoverflow.com/questions/37737538/merge-matplotlib-subplots-with-shared-x-axis/37738851
                    
                    
                    # lightcurve_filename = "./resimgs/foranndy/lightcurve_{}_{}.png".format(timestr, stupstr)
    '''



    # lightcurve_filename = "./resimgs/neg_mass2point_samez/massratio_{:.3f}X1{:.2f}X2{:.2f}_lightcurve.png".format(massratio,lens1.pos[0],lens2.pos[0])
    # fig.savefig(lightcurve_filename, format='png', bbox_inches='tight', dpi=900, pad_inches = 0)#, transparent=True

    # fig, axs = plt.subplots(3, 1, sharex=True)
    # # Remove horizontal space between axes
    # fig.subplots_adjust(hspace=0)
    
    # # Plot each graph, and manually set the y tick values
    # axs[0].plot(t, s1)
    # axs[0].set_yticks(np.arange(-0.9, 1.0, 0.4))
    # axs[0].set_ylim(-1, 1)

    # axs[1].plot(t, s2)
    # axs[1].set_yticks(np.arange(0.1, 1.0, 0.2))
    # axs[1].set_ylim(0, 1)

    # axs[2].plot(t, s3)
    # axs[2].set_yticks(np.arange(-0.9, 1.0, 0.4))
    # axs[2].set_ylim(-1, 1)


    # fig=plt.figure()
    # ax1 = plt.subplot(211)
    # ax2 = plt.subplot(212)

    # ax1.plot(t,x)
    # ax2.plot(t,y)

    # ax1.get_shared_x_axes().join(ax1, ax2)
    # ax1.set_xticklabels([])
    # # ax2.autoscale() ## call autoscale if needed


    # fig = plt.figure()
    # plt.subplot(111)
    # cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    # plt.imshow(where_smaller_than1, origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # # title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    # # title = ""
    # title = "Region where magnification < 1, with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    # plt.title(title)
    # plt.colorbar()


    # plt.show()
    # if input("save imgs or not? (y/n)>>>: ")=="y":
    #     fig.savefig(caustics_filename, format='png', bbox_inches='tight', dpi=600, pad_inches = 0)#, transparent=True
    #     fig.savefig(lightcurve_filename, format='png', bbox_inches='tight', dpi=600, pad_inches = 0)#, transparent=True
    # if input("save imgs or not? (y/n)>>>: ")=="y":
    #     #     fig1.savefig(caustics_filename, format='png', bbox_inches='tight', dpi=300, pad_inches = 0)#, transparent=True
    #     #     fig2.savefig(lightcurve_filename, format='png', bbox_inches='tight', dpi=300, pad_inches = 0)#, transparent=True
    #     fig3.savefig(lightcurve_filename, format='png', bbox_inches='tight', dpi=100, pad_inches = 0)#, transparent=True