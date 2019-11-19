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
# from matplotlib import animation as ani
# from matplotlib.animation import FuncAnimation, writers


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
        # Erdl & Schneider 1993, equ (3), (4)
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
    def scatter_caustics(self):
        spescatter(self.betax,self.betay,s=0.01 , xlabel = r"$\beta_x$/arcsec",ylabel = r"$\beta_y$/arcsec",title="Caustics",issqure=False)
    def img_mapping(self, xdata, ydata, xlim, ylim, ImgSize):
        IMG = np.zeros(ImgSize)
        ratiox, ratioy = (ImgSize[1]-1)/(xlim[1]-xlim[0]), (ImgSize[0]-1)/(ylim[1]-ylim[0])
        xdata -= xlim[0]
        ydata -= ylim[0]
        # print(xdata.shape) #(250000,)
        # too slow:
        for i in range(len(xdata)):
            # x, y = xdata[i], ydata[i]
            try:
                # IMG[int(y*ratioy), int(x*ratiox)] += 1
                IMG[int(ydata[i]*ratioy), int(xdata[i]*ratiox)] += 1
            except:
                pass
        # xdata *= ratiox
        # ydata *= ratioy
        # xdata, ydata = xdata.astype(np.int), ydata.astype(np.int)


        

        # IMG[(xdata*ratiox).astype(np.int),(ydata*ratioy).astype(np.int)]
        return IMG



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
    lens1 = Onelens(1, (0,0))
    lens2 = Onelens(1, (0,0))
    xlim, ylim =(-3.5,2), (-3,3)
    ImgSize = (1026,1026) # raw, colume
    twolens = Twolenses(lens1, lens2, 0.1)
    thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=1000)
    twolens.inverse_ray_shooting_diffz(thetax, thetay)


    
    #https://stackoverflow.com/questions/53101815/improve-color-contrast-in-matplotlib
    # try one of the continuous colormaps, like viridis or gist_ncar. 
    # If you don't need all colors to be different, 
    # try one of the repeating maps like flag or prism
    cmap = plt.cm.get_cmap('hot')# 'Paired', viridis, gist_ncar, 

    srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize)
    fig = plt.figure()
    plt.imshow(srcplaneIMG, origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    # plt.pcolor(srcplaneIMG)
    title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
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

    plt.show()

