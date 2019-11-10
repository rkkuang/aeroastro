'''
Draw the critical lines and caustics of two SIS lenses,
loacted at the same redshift



'''

import numpy as np 
import matplotlib.pyplot as plt
from two_pointM_lens import spescatter

class oneSIS():
    def __init__(self, sigma_v=200, z=1, D_ds=1, D_s=2, loc=(0.8,0)):
        self.sigma_v = sigma_v #in km/s
        self.z = z
        self.D_ds = D_ds
        self.D_s = D_s
        self.theta_E = 1.6*(self.sigma_v/200)**2*(self.D_ds/self.D_s) #in arcsec
        self.theta_E = self.theta_E/3600/180*np.pi # in rad
        # self.loc = loc #in srcsec, location at x axis of lens plane
        self.locx = loc[0]/3600/180*np.pi # !!!
        self.locy = loc[1]/3600/180*np.pi # !!!
class twoSISes():
    def __init__(self, lens1, lens2):
        self.lens1 = lens1
        self.lens2 = lens2
        self.detA = None
        self.mu = None
        self.criticalx = None
        self.criticaly = None
        self.causticsx = None
        self.causticsy = None
    def comp_magnification_diffz(self,x,y):
        a1, a2 = self.lens1.theta_E, self.lens2.theta_E
        rho1x, rho1y, rho2x, rho2y  = self.lens1.locx, self.lens1.locy, self.lens2.locx, self.lens2.locy
        r1 = np.sqrt((x-rho2x)**2+(y-rho2y)**2)
        b1 = a2*(1/r1 - (x-rho2x)**2/r1**3)
        b2 = a2*(x-rho2x)*(y-rho2y)/r1**3
        b3 = a2*(1/r1 - (y-rho2y)**2/r1**3)
        kx = x - a2/r1*(x-rho2x) - rho1x
        ky = y - a2/r1*(y-rho2y) - rho1y
        r2 = np.sqrt(kx**2 + ky**2)
        A11 = 1 - b1 - a1*( (1-b1)/r2 - kx/r2**3*( (1-b1)*kx + b2*ky ) )
        A12 = b2 - a1*( b2/r2 - kx/r2**3*( b2*kx + (1-b3)*ky ) )
        A21 = b2 - a1*( b2/r2 - ky/r2**3*( (1-b1)*kx + b2*ky ) )
        A22 = 1 - b3 - a1*( (1-b3)/r2 - ky/r2**3*( b2*kx + (1-b3)*ky ) )

        self.detA = np.multiply(A11,A22) - np.multiply(A12,A12)
        self.mu = 1/self.detA
        print(np.abs(self.mu).max())
        print(np.abs(self.mu).min())

    def comp_magnification_diffz_v2(self,x,y):
        a1, a2 = self.lens1.theta_E, self.lens2.theta_E
        rho1x, rho1y, rho2x, rho2y  = self.lens1.locx, self.lens1.locy, self.lens2.locx, self.lens2.locy
        r1 = np.sqrt((x-rho2x)**2+(y-rho2y)**2)
        # b1 = a2*(1/r1 - (x-rho2x)**2/r1**3)
        # b2 = a2*(x-rho2x)*(y-rho2y)/r1**3
        # b3 = a2*(1/r1 - (y-rho2y)**2/r1**3)
        # kx = x - a2/r1*(x-rho2x) - rho1x
        # ky = y - a2/r1*(y-rho2y) - rho1y
        # r2 = np.sqrt(kx**2 + ky**2)
        # A11 = 1 - b1 - a1*( (1-b1)/r2 - kx/r2**3*( (1-b1)*kx + b2*ky ) )
        # A12 = b2 - a1*( b2/r2 - kx/r2**3*( b2*kx + (1-b3)*ky ) )
        # A21 = b2 - a1*( b2/r2 - ky/r2**3*( (1-b1)*kx + b2*ky ) )
        # A22 = 1 - b3 - a1*( (1-b3)/r2 - ky/r2**3*( b2*kx + (1-b3)*ky ) )
        mx = x - a2*(x-rho2x)/r1
        my = y - a2*(y-rho2y)/r1
        r2 = np.sqrt((mx-rho1x)**2+(my-rho1y)**2)
        dmx_dx = 1 - a2*(1/r1 - (x-rho2x)**2/r1**3)

        dmx_dy = -a2*(-(x-rho2x)*y/r1**3)
        dmy_dx = -a2*(-(y-rho2y)*x/r1**3)
        dmy_dy = 1 - a2*(1/r1 - (y-rho2y)**2/r1**3)

        dr2_dx = 2 * ( (mx-rho1x)*dmx_dx + (my-rho1y)*dmy_dx)
        dr2_dy = 2 * ( (mx-rho1x)*dmx_dy + (my-rho1y)*dmy_dy)

        A11 = dmx_dx - a1*(dmx_dx/r2 - 0.5*(mx-rho1x)*dr2_dx/r2**3)
        A12 = dmx_dy - a1*(dmx_dy/r2 - 0.5*(mx-rho1x)*dr2_dy/r2**3)
        A21 = dmy_dx - a1*(dmy_dx/r2 - 0.5*(my-rho1y)*dr2_dx/r2**3)
        A22 = dmy_dy - a1*(dmy_dy/r2 - 0.5*(my-rho1y)*dr2_dy/r2**3)


        self.detA = np.multiply(A11,A22) - np.multiply(A12,A12)
        self.mu = 1/self.detA
        print(np.abs(self.mu).max())
        print(np.abs(self.mu).min())

    def comp_magnification(self,x,y):
        rho1x = self.lens1.locx
        rho2x = self.lens2.locx
        r1 = np.sqrt(np.multiply(x-rho1x,x-rho1x)+np.multiply(y,y))
        r2 = np.sqrt(np.multiply(x-rho2x,x-rho2x)+np.multiply(y,y))

        A11 = 1-(self.lens1.theta_E*(1/r1-(x-rho1x)**2/r1**3)
            + self.lens2.theta_E*(1/r2-(x-rho2x)**2/r2**3))

        A12 = -(-self.lens1.theta_E*(x-rho1x)*y/r1**3
             - self.lens2.theta_E*(x-rho2x)*y/r2**3)
        # A21 = A12
        A22 = 1-(self.lens1.theta_E*(1/r1-y**2/r1**3)
            + self.lens2.theta_E*(1/r2-y**2/r2**3))

        # self.detA = np.multiply(A11,A22) - np.multiply(A12,A12)
        self.detA = A11*A22-A12*A12
        # mu = 1/np.abs(detA)
        self.mu = 1/self.detA
        print(np.abs(self.mu).max())
        print(np.abs(self.mu).min())
        # print(self.detA.max())
        # print(self.detA.min())
        # input()
        # print(mu)
    def gen_critical_lines1(self,pltlens=False,issqure=True,anno = False, xlim=(-1,1),ylim=(-1,1),num=100,threshold = 0.985, magfunc=None, title="Critical lines",xylim = None):
        thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num)
        # self.comp_magnification(thetax, thetay)
        magfunc(thetax, thetay)
        # criticalx, criticaly = self.findcriticall0(thetax,thetay)
        self.criticalx, self.criticaly = self.findcritical(threshold,thetax,thetay)
        # print(criticalx)
        # print(criticaly)
        # title = "The Critical Curves of two SIS lens located at\n"+r"({:.1f},0),Ds/Dds={:.1f};({:.1f},0),Ds/Dds={:.1f} with $\theta_E=${:.1f};{:.1f} in arcsec".format(self.lens1.loc*180/np.pi*3600,self.lens1.D_s/self.lens1.D_ds,
        #     self.lens2.loc*180/np.pi*3600,self.lens2.D_s/self.lens2.D_ds,self.lens1.theta_E*180/np.pi*3600, self.lens2.theta_E*180/np.pi*3600)
        
        
        spescatter(self.criticalx*180/np.pi*3600,self.criticaly*180/np.pi*3600,r"$\theta_1$/arcsec",r"$\theta_2$/arcsec",title=title,issqure=issqure,xylim=xylim)
        if pltlens:
            markersize = 100  
            plt.scatter([self.lens1.locx*180/np.pi*3600],
                [self.lens1.locy*180/np.pi*3600],s=markersize,c="g",marker="x")
            plt.scatter([self.lens2.locx*180/np.pi*3600],
                [self.lens2.locy*180/np.pi*3600],s=markersize*(self.lens2.D_ds/self.lens1.D_ds),c="g",marker="x")



        # print((self.lens1.locx,self.lens1.locy))
        # print((self.lens2.locx,self.lens2.locy))

        # plt.annotate("Lens1 at Ds/Dds = {:.1f}\n (closer to source)".format(self.lens1.D_s/self.lens1.D_ds),
        #     xy=(self.lens1.locx*180/np.pi*3600,self.lens1.locy*180/np.pi*3600),
        #     xytext=(self.lens1.locx*180/np.pi*3600-abs(self.lens1.locx*180/np.pi*3600)/4,self.lens1.locy*180/np.pi*3600-abs(self.lens1.locx*180/np.pi*3600)/4)
        #     ,arrowprops=dict(facecolor='black', shrink=0.05),)#,textcoords='offset points'
        # plt.annotate("Lens2 at Ds/Dds = {:.1f}".format(self.lens2.D_s/self.lens2.D_ds),
        #     xy=(self.lens2.locx*180/np.pi*3600,self.lens2.locy*180/np.pi*3600),
        #     xytext=(self.lens2.locx*180/np.pi*3600-abs(self.lens2.locx*180/np.pi*3600)/4,self.lens2.locy*180/np.pi*3600+abs(self.lens2.locx*180/np.pi*3600)/6)
        #     ,arrowprops=dict(facecolor='black', shrink=0.05),)#,textcoords='offset points'
        
        #arrowstyle="<->", connectionstyle="arc3")
        if anno:
            arrowprops = dict(
                arrowstyle = "->",
                connectionstyle = "arc3")#dict(facecolor='black', shrink=0.05)

            offset = 50
            fontsize = 14
            plt.annotate("Lens1 at Ds/Dds = {:.1f}\n (closer to the source)".format(self.lens1.D_s/self.lens1.D_ds),
                xy=(self.lens1.locx*180/np.pi*3600,self.lens1.locy*180/np.pi*3600),
                xytext=(-1*offset, -1*offset)
                ,arrowprops=arrowprops,textcoords='offset points',color="green",fontsize=fontsize)#
            plt.annotate("Lens2 at Ds/Dds = {:.1f}".format(self.lens2.D_s/self.lens2.D_ds),
                xy=(self.lens2.locx*180/np.pi*3600,self.lens2.locy*180/np.pi*3600),
                xytext=(-1*offset, 1*offset)
                ,arrowprops=arrowprops,textcoords='offset points',color="green",fontsize=fontsize)

        

        # plt.text(self.lens2.locx,
        # self.lens2.locy-self.lens2.locx,
        # "Lens2 at D_ds = {}".format(self.lens2.D_ds),
        # fontsize=15,
        # verticalalignment="top",
        # horizontalalignment="right"
        # )
        '''
        annotate(s='str' ,xy=(x,y) ,xytext=(l1,l2) ,..)


        s 为注释文本内容 
        xy 为被注释的坐标点
        xytext 为注释文字的坐标位置
        '''

    def gen_caustic_lines_diffz_v2(self,issqure=False,xylim = None,title = "Caustics",xlabel=r"x (arcsec)",ylabel=r"y (arcsec)"):
        x = self.criticalx
        y = self.criticaly
        a1, a2 = self.lens1.theta_E, self.lens2.theta_E
        rho1x, rho1y, rho2x, rho2y  = self.lens1.locx, self.lens1.locy, self.lens2.locx, self.lens2.locy
        r1 = np.sqrt((x-rho2x)**2+(y-rho2y)**2)

        mx = x - a2*(x-rho2x)/r1
        my = y - a2*(y-rho2y)/r1
        r2 = np.sqrt((mx-rho1x)**2+(my-rho1y)**2)

        self.causticsx = mx - a1*(mx-rho1x)/r2
        self.causticsy = my - a1*(my-rho1y)/r2
        markersize = 100
        spescatter(self.causticsx*180/np.pi*3600,self.causticsy*180/np.pi*3600,xlabel=xlabel,ylabel=ylabel,
            title=title,issqure=issqure,xylim = xylim,c="r")
        plt.scatter([self.lens1.locx*180/np.pi*3600],
            [self.lens1.locy*180/np.pi*3600],s=markersize,c="g",marker="x")
        plt.scatter([self.lens2.locx*180/np.pi*3600],
            [self.lens2.locy*180/np.pi*3600],s=markersize*(self.lens2.D_ds/self.lens1.D_ds),c="g",marker="x")




    def gen_caustic_lines_diffz(self,issqure=False,xylim = None,title = "Caustics",xlabel=r"x (arcsec)",ylabel=r"y (arcsec)"):
        x = self.criticalx
        y = self.criticaly
        a1, a2 = self.lens1.theta_E, self.lens2.theta_E
        rho1x, rho1y, rho2x, rho2y  = self.lens1.locx, self.lens1.locy, self.lens2.locx, self.lens2.locy
        r1 = np.sqrt((x-rho2x)**2+(y-rho2y)**2)
        kx = x - a2/r1*(x-rho2x) - rho1x
        ky = y - a2/r1*(y-rho2y) - rho1y
        r2 = np.sqrt(kx**2 + ky**2)
        self.causticsx = x - a2/r1*(x-rho2x) - a1/r2*kx
        self.causticsy = y - a2/r1*(y-rho2y) - a1/r2*ky
        markersize = 100
        spescatter(self.causticsx*180/np.pi*3600,self.causticsy*180/np.pi*3600,xlabel=xlabel,ylabel=ylabel,
            title=title,issqure=issqure,xylim = xylim,c="r")
        plt.scatter([self.lens1.locx*180/np.pi*3600],
            [self.lens1.locy*180/np.pi*3600],s=markersize,c="g",marker="x")
        plt.scatter([self.lens2.locx*180/np.pi*3600],
            [self.lens2.locy*180/np.pi*3600],s=markersize*(self.lens2.D_ds/self.lens1.D_ds),c="g",marker="x")

        # arrowprops = dict(
        #     arrowstyle = "->",
        #     connectionstyle = "arc3")#dict(facecolor='black', shrink=0.05)

        # offset = 50
        # plt.annotate("Lens1 at Ds/Dds = {:.1f}\n (closer to source)".format(self.lens1.D_s/self.lens1.D_ds),
        #     xy=(self.lens1.locx*180/np.pi*3600,self.lens1.locy*180/np.pi*3600),
        #     xytext=(-1*offset, -1*offset)
        #     ,arrowprops=arrowprops,textcoords='offset points',color="red")#
        # plt.annotate("Lens2 at Ds/Dds = {:.1f}".format(self.lens2.D_s/self.lens2.D_ds),
        #     xy=(self.lens2.locx*180/np.pi*3600,self.lens2.locy*180/np.pi*3600),
        #     xytext=(-1*offset, 1*offset)
        #     ,arrowprops=arrowprops,textcoords='offset points',color="green")


    def gen_caustic_lines(self,issqure=True,xylim = None,title = "Caustics",xlabel=r"x (arcsec)",ylabel=r"y (arcsec)"):
        # beta = theta - alpha(theta)
        #beta = theta - thetaE*vectheta/|vectheta|
        # print(self.lens1.loc,self.lens2.loc)
        r1 = np.sqrt((self.criticalx-self.lens1.locx)**2+self.criticaly**2)
        r2 = np.sqrt((self.criticalx-self.lens2.locx)**2+self.criticaly**2)
        # print(r1)

        '''
        wrong codes !!!!!!!!: 
        程序里的2个bug，第一个方框多余，第二个方框是 alphay，当时写的时候不够仔细
        alphax = self.criticalx - (self.lens1.theta_E*(self.criticalx-self.lens1.loc)/r1+self.lens2.theta_E*(self.criticalx-self.lens2.loc)/r2)
        alphay = self.criticaly - (self.lens1.theta_E*self.criticaly/r1+self.lens2.theta_E*self.criticaly/r2)
        self.causticsx = self.criticalx - alphax
        self.causticsy = self.criticaly - alphax
        '''

        alphax =  (self.lens1.theta_E*(self.criticalx-self.lens1.locx)/r1+self.lens2.theta_E*(self.criticalx-self.lens2.locx)/r2)
        alphay =  (self.lens1.theta_E*self.criticaly/r1+self.lens2.theta_E*self.criticaly/r2)

        self.causticsx = self.criticalx - alphax
        # self.causticsy = self.criticaly - alphax  # wrong again!!!!!
        self.causticsy = self.criticaly - alphay
        spescatter(self.causticsx*180/np.pi*3600,self.causticsy*180/np.pi*3600,xlabel=xlabel,ylabel=ylabel,
            title=title,issqure=issqure,xylim = xylim,c="r")
        # plt.scatter([self.lens1.locx*180/np.pi*3600,self.lens2.locx*180/np.pi*3600],
            # [self.lens1.locy*180/np.pi*3600,self.lens2.locy*180/np.pi*3600],s=100,c="r",marker="x")
        markersize = 100
        plt.scatter([self.lens1.locx*180/np.pi*3600],
            [self.lens1.locy*180/np.pi*3600],s=markersize,c="g",marker="x")
        plt.scatter([self.lens2.locx*180/np.pi*3600],
            [self.lens2.locy*180/np.pi*3600],s=markersize*(self.lens2.D_ds/self.lens1.D_ds),c="g",marker="x")

    def findcritical(self,threshold,thetax,thetay):
        # threshold = self.mu.max()*threshold
        mu = np.abs(self.mu)
        criticalx = thetax[mu>threshold]
        criticaly = thetay[mu>threshold]
        return criticalx, criticaly
    def findcriticalg0(self,thetax,thetay):
        threshold = 0
        criticalx = thetax[self.mu>threshold]
        criticaly = thetay[self.mu>threshold]
        return criticalx, criticaly       
    def findcriticall0(self,thetax,thetay):
        threshold = 0
        mu = np.abs(self.mu)
        criticalx = thetax[mu<threshold]
        criticaly = thetay[mu<threshold]
        return criticalx, criticaly 

def genxy(xlim=(-1,1),ylim=(-1,1),num=100):
    x = np.linspace(xlim[0],xlim[1],num)
    y = np.linspace(ylim[0],ylim[1],num)
    X,Y = np.meshgrid(x,y)
    return X.reshape(1,-1), Y.reshape(1,-1)
    

if __name__ == '__main__':

    ###testing two SIS lenses at different redshift
    loc1 = (2.0,0) # in arcsec
    loc2 = (-2.0,0)
    #lens1 closer to the source
    lens1 = oneSIS(sigma_v=800, z=1, D_ds=1, D_s=4, loc=loc1)
    lens2 = oneSIS(sigma_v=500, z=1, D_ds=2, D_s=4, loc=loc2)
    lenses = twoSISes(lens1,lens2)
    maxloc = max([abs(loc1[0]), abs(loc1[1]), abs(loc2[0]), abs(loc2[1])])
    mx = max(lens1.theta_E+maxloc/3600/180*np.pi,lens2.theta_E+maxloc/3600/180*np.pi)*2.5
    xlim = (-mx,mx)
    ylim = xlim
    print("xlim in arcsec",xlim[0]*180/np.pi*3600,xlim[1]*180/np.pi*3600)
    # plt.subplot(121)
    crititle = "The Critical Curves(blue) of two SIS lens(green) located at different redshift\n"+r"loc: ({:.1f},{:.1f}),({:.1f},{:.1f}) and $\theta_E: ${:.1f},{:.1f} in arcsec".format(lenses.lens1.locx*180/np.pi*3600,
        lenses.lens1.locy*180/np.pi*3600,lenses.lens2.locx*180/np.pi*3600,lenses.lens2.locy*180/np.pi*3600,lenses.lens1.theta_E*180/np.pi*3600, lenses.lens2.theta_E*180/np.pi*3600)
    lenses.gen_critical_lines1(pltlens=True,anno = True,xlim=xlim,ylim=xlim,threshold = 3e3,num=1000,magfunc=lenses.comp_magnification_diffz_v2)
    # plt.subplot(122)
    cau_title = "\n and the corresponding Caustics(red) on source plane"
    lenses.gen_caustic_lines_diffz_v2(title=crititle+cau_title,xlabel=r"x (arcsec)",ylabel=r"y (arcsec)")
    plt.tight_layout()
    plt.show()

    # # # testing two SIS lenses at the same redshift
    # loc1 = 1. # in arcsec
    # loc2 = -1.
    # lens1 = oneSIS(sigma_v=220, z=1, D_ds=1, D_s=2, loc=(loc1,0))
    # lens2 = oneSIS(sigma_v=220, z=1, D_ds=1, D_s=2, loc=(loc2,0))
    # lenses = twoSISes(lens1,lens2)
    # mx = max(lens1.theta_E+max(abs(loc1),abs(loc2))/3600/180*np.pi,lens2.theta_E+max(abs(loc1),abs(loc2))/3600/180*np.pi)*2.5
    # xlim = (-mx,mx)
    # ylim = xlim
    # cri_title = "The Critical Curves(blue) of two SIS lens(green) located at\n"+r"({:.1f},0),Ds/Dds={:.1f};({:.1f},0),Ds/Dds={:.1f} with $\theta_E=${:.1f};{:.1f} in arcsec".format(lenses.lens1.locx*180/np.pi*3600,lenses.lens1.D_s/lenses.lens1.D_ds,
    #     lenses.lens2.locx*180/np.pi*3600,lenses.lens2.D_s/lenses.lens2.D_ds,lenses.lens1.theta_E*180/np.pi*3600, lenses.lens2.theta_E*180/np.pi*3600)
    # print("xlim in arcsec",xlim[0]*180/np.pi*3600,xlim[1]*180/np.pi*3600)
    # # plt.subplot(121)
    # cau_title = "\n"+r"and the corresponding Caustics(red) on source plane"
    # title = cri_title + cau_title
    # xylim = max(3,2*loc1)
    # lenses.gen_critical_lines1(xylim = xylim,issqure=True, xlim=xlim,ylim=xlim,threshold = 5e3,num=6000,magfunc=lenses.comp_magnification)
    # # plt.subplot(122)
    # lenses.gen_caustic_lines(xylim = xylim,issqure=True,title=title,xlabel=r"x (arcsec)",ylabel=r"y (arcsec)")
    # plt.grid()
    # plt.show()   