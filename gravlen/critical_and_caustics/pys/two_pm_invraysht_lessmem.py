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
from utils import genxy, spescatter
import matplotlib.pyplot as plt
from pathos.pools import ProcessPool
# from matplotlib import animation as ani
# from matplotlib.animation import FuncAnimation, writers
arcsec2rad = 1/3600/180*np.pi
rad2arcsec = 180/np.pi*3600

class Onelens(object):
    """docstring for Onelens"""
    def __init__(self, mass, position, dis=0):
        super(Onelens, self).__init__()
        self.mass = mass
        self.pos = position
        self.dis = dis
        

class Twolenses(object):
    """docstring for Twolenses"""
    def __init__(self, lens1, lens2):
        super(Twolenses, self).__init__()
        self.lens1, self.lens2 = lens1, lens2
        self.massratio = lens1.mass/lens2.mass
        self.mu1 = self.massratio/(1+self.massratio)
        self.mu2 = 1/(1+self.massratio)
        self.betax, self.betay = None, None
        self.mag = None
        # Erdl & Schneider 1993, equ (3), (4)
        if not self.lens2.dis == 0:
            M1 = self.mu1
            M2 = self.mu2
            D1s, D2, D2s, D1 = 1 - self.lens1.dis, self.lens2.dis, 1 - self.lens2.dis, self.lens1.dis
            self.m1 = M1/(M1+D2s*D1/(D1s*D2)*M2)
            self.m2 = 1 - self.m1
            self.beta = (D2 - D1)/(D1s*D2) # D12*Ds/D1s/D2 from Erdl & Schneider 1993

    
    def inverse_ray_shooting_diffz(self, thetax, thetay):
        # lens1.pos should be (0,0)
        r_2 = (thetax - self.lens1.pos[0])**2 + (thetay - self.lens1.pos[1])**2
        x2_x = thetax - self.m1*self.beta*thetax/r_2
        x2_y = thetay - self.m1*self.beta*thetay/r_2
        r2_2 = (x2_x - self.lens2.pos[0])**2 + (x2_y - self.lens2.pos[1])**2
        self.betax = thetax - self.m1*thetax/r_2 - self.m2*(x2_x - self.lens2.pos[0])/r2_2
        self.betay = thetay - self.m1*thetay/r_2 - self.m2*(x2_y - self.lens2.pos[1])/r2_2

    def inverse_ray_shooting(self, thetax, thetay):
        r1_2 = (thetax - self.lens1.pos[0])**2 + (thetay - self.lens1.pos[1])**2
        r2_2 = (thetax - self.lens2.pos[0])**2 + (thetay - self.lens2.pos[1])**2
        betax = thetax - self.mu1*(thetax - self.lens1.pos[0])/r1_2 - self.mu2*(thetax - self.lens2.pos[0])/r2_2
        betay = thetay - self.mu1*(thetay - self.lens1.pos[1])/r1_2 - self.mu2*(thetay - self.lens2.pos[1])/r2_2
        return betax, betay
    def comp_mag_samez(self, thetax, thetay):
        x1,y1,x2,y2 = self.lens1.pos[0], self.lens1.pos[1], self.lens2.pos[0], self.lens2.pos[1]
        R1 = np.sqrt((thetax-x1)**2+(thetay-y1)**2)
        R2 = np.sqrt((thetax-x2)**2+(thetay-y2)**2)
        A11 = 1 - self.mu1*(1/R1**2 - 2*(thetax-x1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetax-x2)**2/R2**4)
        A12 = -2*self.mu1*(thetax-x1)*(thetay-y1)/R1**4 - 2*self.mu2*(thetax-x2)*(thetay-y2)/R2**4
        A22 = 1 - self.mu1*(1/R1**2 - 2*(thetay-y1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetay-y2)**2/R2**4)
        detA = (A11*A22 - A12**2)
        # self.mag = 1/np.abs(detA)
        mag = 1/np.abs(detA)
        return mag

    # def main_plot(self, inxdata, inydata, xlim, ylim, ImgSize):
    #     ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
    #     xdata = inxdata - xlim[0]
    #     ydata = inydata - ylim[0]
    #     return round(ydata*ratioy), round(xdata*ratiox)

    def scatter_caustics(self):
        spescatter(self.betax,self.betay,s=0.01, xlabel = r"$\beta_x$/arcsec",ylabel = r"$\beta_y$/arcsec",title="Caustics",issqure=False)
    
    def img_mapping(self, inxdata, inydata, xlim, ylim, ImgSize,valarr = None):
        IMG = np.zeros(ImgSize)
        # IMG = np.ones(ImgSize)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        xdata = inxdata - xlim[0]
        ydata = inydata - ylim[0]
        # print(xdata.shape) #(250000,)
        # too slow:
        print((len(xdata)))
        for i in range(len(xdata)):
            # x, y = xdata[i], ydata[i]
            try:
                # IMG[int(y*ratioy), int(x*ratiox)] += 1
                IMG[round(ydata[i]*ratioy), round(xdata[i]*ratiox)] += valarr[i] #valarr
            except:
                pass
        # xdata *= ratiox
        # ydata *= ratioy
        # xdata, ydata = xdata.astype(np.int), ydata.astype(np.int)
        return IMG

    def img_mapping_inone(self, inxdata, inydata, betainxdata, betainydata, xlim, ylim, ImgSize,valarr1 = None,valarr2 = None):
        # srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize, valarr = twolens.mag) #valarr = np.ones([len(twolens.betax),])
        # srcplaneIMG_withoutlens = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = np.ones([len(thetax),]))
        # srcplaneIMG /= srcplaneIMG_withoutlens
        # imgplaneIMG = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = twolens.mag)
        #srcplaneIMG, imgplaneIMG = img_mapping_inone(thetax, thetay,twolens.betax, twolens.betay, xlim, ylim, ImgSize,valarr1 = twolens.mag, valarr2 = np.ones([len(thetax),]))
        srcplaneIMG = np.zeros(ImgSize)
        srcplaneIMG_withoutlens = np.zeros(ImgSize)
        imgplaneIMG = np.zeros(ImgSize)
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
        return srcplaneIMG, imgplaneIMG

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

    def get_imgs_lessmem(self,ImgSize,xlim,ylim,THETAX,THETAY):

        srcplaneIMG = np.zeros(ImgSize)
        srcplaneIMG_withoutlens = np.zeros(ImgSize)
        imgplaneIMG = np.zeros(ImgSize)

        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        cnt=1
        lenXY = len(THETAX)

        for thetax, thetay in zip(THETAX, THETAY):
            sys.stdout.write('done %d/%d\r' % (cnt, lenXY))
            cnt += 1

            betax, betay = self.inverse_ray_shooting(thetax, thetay)
            mag = self.comp_mag_samez(thetax, thetay)
        
            xdata = thetax - xlim[0]# thetax
            ydata = thetay - ylim[0]# thetax
            betaxdata = betax - xlim[0]
            betaydata = betay - ylim[0]

            try:
                srcplaneIMG[(betaydata*ratioy+0.5).astype(np.int), (betaxdata*ratiox+0.5).astype(np.int)]+=mag
            except:
                pass
            srcplaneIMG_withoutlens[(ydata*ratioy+0.5).astype(np.int), (xdata*ratiox+0.5).astype(np.int)]+=1
            imgplaneIMG[(ydata*ratioy+0.5).astype(np.int), (xdata*ratiox+0.5).astype(np.int)]+=mag


        # BETAX, BETAY = self.inverse_ray_shooting(THETAX, THETAY)
        # MAG = self.comp_mag_samez(THETAX, THETAY)
        # THETAX = ((THETAX - xlim[0])*ratiox+0.5).astype(np.int)
        # THETAY = ((THETAY - ylim[0])*ratioy+0.5).astype(np.int)
        # BETAX = ((BETAX - xlim[0])*ratiox+0.5).astype(np.int)
        # BETAY = ((BETAY - ylim[0])*ratioy+0.5).astype(np.int)
        

        # for thetax, thetay, betax, betay, mag in zip(THETAX, THETAY, BETAX, BETAY, MAG):
        #     sys.stdout.write('done %d/%d\r' % (cnt, lenXY))
        #     cnt += 1
        #     try:
        #         srcplaneIMG[betay, betax]+=mag
        #     except:
        #         pass
        #     srcplaneIMG_withoutlens[thetay, thetax]+=1
        #     imgplaneIMG[thetay, thetax]+=mag

        srcplaneIMG /= srcplaneIMG_withoutlens        
        return srcplaneIMG, imgplaneIMG

if __name__ == '__main__':

    import sys
    mass1 = float(sys.argv[1])
    lens1 = Onelens(mass1, (-0.3,0))
    lens2 = Onelens(1, (0.3,0))
    d1d2 = (0,0)
    lens1.dis = d1d2[0]
    lens2.dis = d1d2[1]
    massratio = lens1.mass/lens2.mass
    xlim, ylim = (-2.5,2.5), (-2.5,2.5)
    # srcxlim, srcylim = (-4,4), (-4, 4)
    ImgSize = (512,512) # raw, colume
    import time
    
    num = 1e4 # 1e4, 20:18:50-20:34:55=16minutes; VS: 80minute!!!
    twolens = Twolenses(lens1, lens2)

    THETAX, THETAY = genxy(xlim=xlim,ylim=ylim,num=num)

    srcplaneIMG, imgplaneIMG = twolens.get_imgs_lessmem(ImgSize,xlim,ylim,THETAX, THETAY)    

    fig = plt.figure()
    plt.subplot(122)
    cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    plt.title(title)
    plt.colorbar()


    plt.subplot(121)
    cmap = plt.cm.get_cmap('viridis')
    plt.imshow(np.log10(imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    #title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
    title = "Two point mass lenses, img plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    plt.title(title)
    plt.colorbar()

    fig.set_size_inches(16,9)


    plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
    plt.margins(0,0)

    radius, npoints = 5, 30 # radius in pixel
    k = 0.2
    B = np.array([0.5,0.4,0.3,0.2,0.1,0.01]) - 0.2
    KB = [(k,b) for b in B]

    cname = ["r","g","y","k","c","m","w"]
    COLLOR = cname[:len(KB)]

    # LABEL = ["label1","label2","label3","label4"]
    Y, LC = {}, {}
    scale = 0.9
    x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/2))#, endpoint=False
    # print(x)
    for kb in KB:
        k,b = kb[0], kb[1]
        Y[kb] = k*x + b
        LC[kb] = twolens.com_lightcurve(srcplaneIMG, x, Y[kb], xlim, ylim, radius, npoints)
        # print(x)

    fig = plt.figure()
    fig.set_size_inches(16,9)
    plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
    plt.margins(0,0)

    plt.subplot(121)

    cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # axes[0,0].imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    plt.title(title)
    plt.colorbar()


    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt])#-ylim[0]
        cnt += 1
    axes = [1 for i in range(len(KB))]
    cnt = 0

    for kb in KB:
        axes[cnt] = plt.subplot(len(KB), 2, cnt*2+2) #int("42"+str(cnt*2+2))
        label = "k={:.2f}, b={:.2f}".format(kb[0],kb[1])
        axes[cnt].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
        # axes[cnt].ylim(0,5)
        axes[cnt].set_ylim([0,5.5])
        axes[cnt].legend()
        if cnt == 0:
            plt.title("Light Curve (Log scale)")
        # axes[cnt,1].plot(x,np.log10(LC[kb]),color=COLLOR[cnt],label=label)
        cnt += 1
    plt.subplots_adjust(hspace=.0)

    plt.show()