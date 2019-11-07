'''
Draw the critical lines and caustics of two SIS lenses,
loacted at the same redshift



'''

import numpy as np 
import matplotlib.pyplot as plt
from two_asym_lens import spescatter

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

    def comp_magnification(self,x,y):
        lens1.loc = lens1.locx
        lens2.loc = lens2.locx
        r1 = np.sqrt(np.multiply(x-lens1.loc,x-lens1.loc)+np.multiply(y,y))
        r2 = np.sqrt(np.multiply(x-lens2.loc,x-lens2.loc)+np.multiply(y,y))

        A11 = 1-(lens1.theta_E*(1/r1-(x-lens1.loc)**2/r1**3)
            + lens2.theta_E*(1/r2-(x-lens2.loc)**2/r2**3))

        A12 = -(-lens1.theta_E*(x-lens1.loc)*y/r1**3
             - lens2.theta_E*(x-lens2.loc)*y/r2**3)
        # A21 = A12
        A22 = 1-(lens1.theta_E*(1/r1-y**2/r1**3)
            + lens2.theta_E*(1/r2-y**2/r2**3))  

        self.detA = np.multiply(A11,A22) - np.multiply(A12,A12)
        # mu = 1/np.abs(detA)
        self.mu = 1/self.detA
        print(np.abs(self.mu).max())
        print(np.abs(self.mu).min())
        # print(self.detA.max())
        # print(self.detA.min())
        # input()
        # print(mu)
    def gen_critical_lines1(self,anno = False, xlim=(-1,1),ylim=(-1,1),num=100,threshold = 0.985, magfunc=None, title="Critical lines"):
        thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num)
        # self.comp_magnification(thetax, thetay)
        magfunc(thetax, thetay)
        # criticalx, criticaly = self.findcriticall0(thetax,thetay)
        self.criticalx, self.criticaly = self.findcritical(threshold,thetax,thetay)
        # print(criticalx)
        # print(criticaly)
        # title = "The Critical Curves of two SIS lens located at\n"+r"({:.1f},0),Ds/Dds={:.1f};({:.1f},0),Ds/Dds={:.1f} with $\theta_E=${:.1f};{:.1f} in arcsec".format(self.lens1.loc*180/np.pi*3600,self.lens1.D_s/self.lens1.D_ds,
        #     self.lens2.loc*180/np.pi*3600,self.lens2.D_s/self.lens2.D_ds,self.lens1.theta_E*180/np.pi*3600, self.lens2.theta_E*180/np.pi*3600)
        
        markersize = 100
        spescatter(self.criticalx*180/np.pi*3600,self.criticaly*180/np.pi*3600,r"$\theta_1$/arcsec",r"$\theta_2$/arcsec",title)
        plt.scatter([self.lens1.locx*180/np.pi*3600],
            [self.lens1.locy*180/np.pi*3600],s=markersize,c="r",marker="+")
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
                ,arrowprops=arrowprops,textcoords='offset points',color="red",fontsize=fontsize)#
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
    def gen_caustic_lines_multi(self,title=None):
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
        spescatter(self.causticsx*180/np.pi*3600,self.causticsy*180/np.pi*3600,r"$\beta_1$/arcsec",r"$\beta_2$/arcsec",title)
        plt.scatter([self.lens1.locx*180/np.pi*3600],
            [self.lens1.locy*180/np.pi*3600],s=markersize,c="r",marker="+")
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


    def gen_caustic_lines(self):
        # beta = theta - alpha(theta)
        #beta = theta - thetaE*vectheta/|vectheta|
        r1 = np.sqrt((self.criticalx-self.lens1.loc)**2+self.criticaly**2)
        r2 = np.sqrt((self.criticalx-self.lens2.loc)**2+self.criticaly**2)
        alphax = self.criticalx - (self.lens1.theta_E*(self.criticalx-self.lens1.loc)/r1+self.lens2.theta_E*(self.criticalx-self.lens2.loc)/r2)
        alphay = self.criticaly - (self.lens1.theta_E*self.criticaly/r1+self.lens2.theta_E*self.criticaly/r2)
        self.causticsx = self.criticalx - alphax
        self.causticsy = self.criticaly - alphax
        spescatter(self.causticsx*180/np.pi*3600,self.causticsy*180/np.pi*3600,r"$\beta_1$/arcsec",r"$\beta_2$/arcsec",
            "The corresponding Caustics on source plane")
        plt.scatter([self.lens1.locx*180/np.pi*3600,self.lens2.locx*180/np.pi*3600],
            [self.lens1.locy*180/np.pi*3600,self.lens2.locy*180/np.pi*3600],s=50,c="r",marker="+")

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

    # testing two SIS lenses at different redshift
    loc1 = (3.0,0) # in arcsec
    loc2 = (-3.0,0)
    #lens1 closer to the source
    lens1 = oneSIS(sigma_v=800, z=1, D_ds=1, D_s=4, loc=loc1)
    lens2 = oneSIS(sigma_v=500, z=1, D_ds=2, D_s=4, loc=loc2)
    lenses = twoSISes(lens1,lens2)
    maxloc = max([abs(loc1[0]), abs(loc1[1]), abs(loc2[0]), abs(loc2[1])])
    mx = max(lens1.theta_E+maxloc/3600/180*np.pi,lens2.theta_E+maxloc/3600/180*np.pi)*2.5
    xlim = (-mx,mx)
    ylim = xlim
    print("xlim in arcsec",xlim[0]*180/np.pi*3600,xlim[1]*180/np.pi*3600)
    plt.subplot(121)
    crititle = "The Critical Curves of two SIS lens located at different redshift\n"+r"loc: ({:.1f},{:.1f}),({:.1f},{:.1f}) and $\theta_E: ${:.1f},{:.1f} in arcsec".format(lenses.lens1.locx*180/np.pi*3600,
        lenses.lens1.locy*180/np.pi*3600,lenses.lens2.locx*180/np.pi*3600,lenses.lens2.locy*180/np.pi*3600,lenses.lens1.theta_E*180/np.pi*3600, lenses.lens2.theta_E*180/np.pi*3600)
    lenses.gen_critical_lines1(anno = True,xlim=xlim,ylim=xlim,threshold = 3e3,num=8000,magfunc=lenses.comp_magnification_diffz,title=crititle)
    plt.subplot(122)
    lenses.gen_caustic_lines_multi("The corresponding Caustics on source plane")
    plt.show()

    # # testing two SIS lenses at the same redshift
    # loc1 = .5 # in arcsec
    # loc2 = -3.5
    # lens1 = oneSIS(sigma_v=200, z=1, D_ds=1, D_s=2, loc=(loc1,0))
    # lens2 = oneSIS(sigma_v=200, z=1, D_ds=1, D_s=2, loc=(loc2,0))
    # lenses = twoSISes(lens1,lens2)
    # mx = max(lens1.theta_E+max(abs(loc1),abs(loc2))/3600/180*np.pi,lens2.theta_E+max(abs(loc1),abs(loc2))/3600/180*np.pi)*2.5
    # xlim = (-mx,mx)
    # ylim = xlim
    # print("xlim in arcsec",xlim[0]*180/np.pi*3600,xlim[1]*180/np.pi*3600)
    # plt.subplot(121)
    # lenses.gen_critical_lines1(xlim=xlim,ylim=xlim,threshold = 5e3,num=6000,magfunc=lenses.comp_magnification)
    # plt.subplot(122)
    # lenses.gen_caustic_lines()
    # plt.show()   