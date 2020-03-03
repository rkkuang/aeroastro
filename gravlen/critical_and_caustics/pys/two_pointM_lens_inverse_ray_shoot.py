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
        self.betax = thetax - self.mu1*(thetax - self.lens1.pos[0])/r1_2 - self.mu2*(thetax - self.lens2.pos[0])/r2_2
        self.betay = thetay - self.mu1*(thetay - self.lens1.pos[1])/r1_2 - self.mu2*(thetay - self.lens2.pos[1])/r2_2
    def comp_mag_samez(self, thetax, thetay):
        arcsec2rad = 1
        x1,y1,x2,y2 = self.lens1.pos[0], self.lens1.pos[1], self.lens2.pos[0], self.lens2.pos[1]
        x1,y1,x2,y2 = x1*arcsec2rad,y1*arcsec2rad,x2*arcsec2rad,y2*arcsec2rad
        # print(x1,y1,x2,y2)
        thetax *= arcsec2rad
        thetay *= arcsec2rad
        R1 = np.sqrt((thetax-x1)**2+(thetay-y1)**2)
        R2 = np.sqrt((thetax-x2)**2+(thetay-y2)**2)
        A11 = 1 - self.mu1*(1/R1**2 - 2*(thetax-x1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetax-x2)**2/R2**4)
        A12 = -2*self.mu1*(thetax-x1)*(thetay-y1)/R1**4 - 2*self.mu2*(thetax-x2)*(thetay-y2)/R2**4
        A22 = 1 - self.mu1*(1/R1**2 - 2*(thetay-y1)**2/R1**4) - self.mu2*(1/R2**2 - 2*(thetay-y2)**2/R2**4)
        detA = (A11*A22 - A12**2)
        self.mag = 1/np.abs(detA)


        # detA = 1-(self.mu1/R1**2+self.mu2/R2**2)**2-16*self.mu1*self.mu2*x1*x2*(thetay-y1)*(thetay-y2)/(R1**4*R2**4)
        # self.mag = 1/np.abs(detA)

        # self.mag = np.abs(self.mag)
        # self.mag /= self.mag.min()

        print("Maximum Magnification: ",np.abs(self.mag).max())
        print("Minimum Magnification: ",np.abs(self.mag).min())


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

        # method 1, single core, time for num=5000: 80.57510709762573
        # for i in range(len(xdata)):
        #     try:
        #         srcplaneIMG[int(betaydata[i]*ratioy), int(betaxdata[i]*ratiox)] += valarr1[i] 
        #     except:
        #         pass
        #     try:
        #         srcplaneIMG_withoutlens[int(ydata[i]*ratioy), int(xdata[i]*ratiox)] += valarr2[i]#np.ones([len(twolens.betax),])
        #     except:
        #         pass
        #     try:
        #         imgplaneIMG[int(ydata[i]*ratioy), int(xdata[i]*ratiox)] += valarr1[i] #valarr
        #     except:
        #         pass

        # # method 2, three cores, time for num=5000: 42.20756411552429
        # self.inputs = [[srcplaneIMG,srcplaneIMG_withoutlens,imgplaneIMG],[((betaydata*ratioy).astype(np.int),(betaxdata*ratiox).astype(np.int))
        # ,((ydata*ratioy).astype(np.int),(xdata*ratiox).astype(np.int)),((ydata*ratioy).astype(np.int),(xdata*ratiox).astype(np.int))],[valarr1,valarr2,valarr1]]
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



if __name__ == '__main__':
    # plt.style.use('classic')


    # lens1 = Onelens(0.5, (0.5,0))
    # lens2 = Onelens(0.5, (-0.5,0))
    # xlim, ylim =(-1.5,1.5), (-1,1)
    # ImgSize = (684,1026) # raw, colume

    # lens1 = Onelens(0.5, (0.3,0))
    # lens2 = Onelens(0.5, (-0.3,0))
    # xlim, ylim =(-1.5,1.5), (-1.5,1.5)
    # ImgSize = (1026,1026) # raw, colume

    # lens1 = Onelens(0.5, (1,0))
    # lens2 = Onelens(0.5, (-1,0))
    # xlim, ylim =(-2,2), (-1,1)
    # ImgSize = (513,1026) # raw, colume

    # # same redshift:
    # lens1 = Onelens(1, (-1,0))
    # lens2 = Onelens(1, (1,0))
    # xlim, ylim =(-3.5,2), (-3,3)
    # ImgSize = (1026,1026) # raw, colume
    # twolens = Twolenses(lens1, lens2, 0)
    # thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=1000)
    # twolens.inverse_ray_shooting(thetax, thetay)

    # diff redshift
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
    num = 4000
    #https://blog.csdn.net/zj360202/article/details/78543141, struct.error: 'i' format requires -2147483648
    # File "/usr/lib/python2.7/pickle.py", line 494, in save_string
    # self.write(BINSTRING + pack("<i", n) + obj)
    # self.write(BINSTRING + pack("<i", n) + obj)  -->  self.write(BINSTRING + pack("<q", n) + obj)
    thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num)
    twolens = Twolenses(lens1, lens2)
    twolens.inverse_ray_shooting(thetax, thetay)
    twolens.comp_mag_samez(thetax, thetay)
    # srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize, valarr = twolens.mag) #valarr = np.ones([len(twolens.betax),])
    # srcplaneIMG_withoutlens = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = np.ones([len(thetax),]))
    # srcplaneIMG /= srcplaneIMG_withoutlens
    # # srcplaneIMG += 1
    # # srcplaneIMG = np.log10(srcplaneIMG)
    # imgplaneIMG = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize, valarr = twolens.mag)
    t0 = time.time()
    srcplaneIMG, imgplaneIMG = twolens.img_mapping_inone(thetax, thetay,twolens.betax, twolens.betay, xlim, ylim, ImgSize,valarr1 = twolens.mag, valarr2 = np.ones([len(thetax),]))
    t1 = time.time()

    where_smaller_than1 = np.zeros(srcplaneIMG.shape)
    where_smaller_than1[:,:] = srcplaneIMG[:,:]
    where_smaller_than1[where_smaller_than1>1] = 0

    # where_smaller_than1 = srcplaneIMG
    # where_smaller_than1[where_smaller_than1>1] = 0

    # srcplaneIMG += 1
    print("time spent on img_mapping_inone: {}".format(t1-t0))
    # imgplaneIMG = np.log10(imgplaneIMG)
    # imgplaneIMG /= np.min(imgplaneIMG)
    
    #https://stackoverflow.com/questions/53101815/improve-color-contrast-in-matplotlib
    # try one of the continuous colormaps, like viridis or gist_ncar. 
    # If you don't need all colors to be different, 
    # try one of the repeating maps like flag or prism
    

    fig = plt.figure()
    plt.subplot(122)
    cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    #title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
    title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    plt.title(title)
    plt.colorbar()

    # fig = plt.figure()
    plt.subplot(121)
    cmap = plt.cm.get_cmap('viridis')
    plt.imshow(np.log10(imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    #title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
    title = "Two point mass lenses, img plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    plt.title(title)
    plt.colorbar()

    # lensplaneIMG = twolens.img_mapping(thetax, thetay, xlim, ylim, ImgSize)
    # if lensplaneIMG.min() == 0:
    #     lensplaneIMG[lensplaneIMG==0]=1
    # CausticsIMG = srcplaneIMG/lensplaneIMG

    # plt.figure()
    # plt.imshow(CausticsIMG, origin='lower',cmap=cmap)
    # plt.colorbar()

    # matplotlib.animation.FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)
    # ani.FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)

    # plt.tight_layout()

    fig.set_size_inches(16,9)


    # plt.subplots_adjust(top = 2, bottom = 2, right = 1, left = 1, hspace = 0, wspace = 0)
    # 有六个可选参数来控制子图布局。值均为0~1之间。其中left、bottom、right、top围成的区域就是子图的区域。wspace、hspace分别表示子图之间左右、上下的间距
    plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
    plt.margins(0,0)
    # caustics_filename = "./resimgs/neg_mass2point_samez/massratio_{:.3f}X1{:.2f}X2{:.2f}_caustics.png".format(massratio,lens1.pos[0],lens2.pos[0])
    # fig.savefig(caustics_filename, format='png', bbox_inches='tight', dpi=900, pad_inches = 0)#, transparent=True

    # k,b = 0,0
    # scale = 0.5
    
    # x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale))#ImgSize[1]
    # y = k*x + b 
    # lc = twolens.com_lightcurve(srcplaneIMG, x, y, xlim, ylim, )
    # plt.figure()
    # plt.plot(x, np.log10(lc))
    # print(np.log10(lc)[:10])
    # # plt.ylim((0,5))
    # plt.ylabel("Log Magnification")
    # plt.show()

    radius, npoints = 5, 30 # radius in pixel
    k = 0.2
    B = np.array([0.5,0.4,0.3,0.2,0.1,0.01]) - 0.2
    KB = [(k,b) for b in B]
    #https://www.cnblogs.com/darkknightzh/p/6117528.html
    #b: blue g: green r: red c: cyan m: magenta y: yellow k: black w: white
    # COLLOR = ["r","g","y","k","c","m","w"]

    # cname = []
    # import matplotlib
    # for name, _ in matplotlib.colors.cnames.items():
        # cname.append(name)

    # import matplotlib.colors as colors
    # import random
    # cname = list(colors._colors_full_map.keys())
    # random.shuffle(cname)

    cname = ["r","g","y","k","c","m","w"]
    COLLOR = cname[:len(KB)]

    # LABEL = ["label1","label2","label3","label4"]
    Y, LC = {}, {}
    scale = 0.9
    # print(xlim,"xlim")
    # print(ylim,"ylim")
    x = np.linspace(xlim[0]*scale, xlim[1]*scale, int(ImgSize[1]*scale/2))#, endpoint=False
    # print(x)
    for kb in KB:
        k,b = kb[0], kb[1]
        Y[kb] = k*x + b
        LC[kb] = twolens.com_lightcurve(srcplaneIMG, x, Y[kb], xlim, ylim, radius, npoints)
        # print(x)

    fig = plt.figure()
    # fig, axs =  plt.subplots(len(KB), 2)
    # https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.subplots.html
    # fig, axes = plt.subplots(4, 2, sharex='col')
    fig.set_size_inches(16,9)
    plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
    plt.margins(0,0)

    plt.subplot(121)

    # plt.title("Source Place")
    cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
    plt.imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # axes[0,0].imshow(np.log10(srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
    plt.title(title)
    # plt.colorbar()
    # plt.hold(True)
    cnt = 0
    for kb in KB:
        plt.plot(x,Y[kb],color=COLLOR[cnt])#-ylim[0]
        cnt += 1


    # plt.subplot(122)

    # plt.subplots(122, sharex=True)
    axes = [1 for i in range(len(KB))]
    cnt = 0
    # plt.subplot(len(KB),2)
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
    # axes[0].get_shared_x_axes().join(axes[0], axes[1], axes[2],axes[3])
    # axes[0].set_xticklabels([])
    plt.subplots_adjust(hspace=.0)
    # plt.subplot(122)
    
    
    # https://stackoverflow.com/questions/37737538/merge-matplotlib-subplots-with-shared-x-axis/37738851
    
    
    

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
    plt.show()