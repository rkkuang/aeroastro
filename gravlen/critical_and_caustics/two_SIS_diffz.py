'''
Draw the critical lines and caustics of two SIS lenses,
loacted at different redshift

ref: Two point mass lens model at different redshift:
Classification of the multiple deflection two point-mass gravitational lens models and
application of catastrophe theory in lensing.
https://ui.adsabs.harvard.edu/abs/1993A%26A...268..453E/abstract

'''

import numpy as np 
import matplotlib.pyplot as plt
from utils import spescatter

arcsec2rad = 1/3600/180*np.pi
rad2arcsec = 180/np.pi*3600

class oneSIS():
    def __init__(self, sigma_v=200, D=2, loc=(0.8,0)):
        self.sigma_v = sigma_v #in km/s
        self.D = D
        self.theta = 1.6*(self.sigma_v/200)**2 #in arcsec
        self.theta = self.theta*arcsec2rad # in rad
        # self.loc = loc #in srcsec, location at x axis of lens plane
        self.locx = loc[0]*arcsec2rad # !!!
        self.locy = loc[1]*arcsec2rad # !!!
class twoSISes():
    def __init__(self, lens1, lens2, D_s = 4):
        self.lens1 = lens1
        self.lens2 = lens2
        self.D_s = D_s
        self.detA = None
        self.mu = None
        self.criticalx = None
        self.criticaly = None
        self.causticsx = None
        self.causticsy = None
        self.D1s = self.D_s - self.lens1.D
        self.D2s = self.D_s - self.lens2.D
        self.D12 = abs(self.lens2.D - self.lens1.D)



    def comp_magnification_diffz_v2(self,x,y):
        # print("x",x)
        theta1, theta2 = self.lens1.theta, self.lens2.theta
        # print("theta1",theta1)
        rhox, rhoy = self.lens2.locx, self.lens2.locy
        # print("rhox",rhox)
        r1 = np.sqrt(x**2+y**2)
        # print(r1)
        m1 = self.D1s*theta1
        m2 = self.D2s*theta2
        beta = self.D12/(self.lens2.D*self.D1s)
        t1 = x - m1*beta*x/r1 - rhox
        t2 = y - m1*beta*y/r1 - rhoy
        r2 = np.sqrt(t1**2+t2**2)
        a = m1*(1/r1 - x**2/r1**3)
        # print("a",a)
        # print("m1",m1)
        b = m1*x*y/r1**3
        c = m1*(1/r1 - y**2/r1**3)
        dt1dx1 = 1 - a*beta
        dt1dx2 = b*beta
        dt2dx1 = dt1dx2
        dt2dx2 = 1 - c*beta

        A11 = 1 - a - m2*( dt1dx1/r2 - t1*(t1*dt1dx1 + t2*dt2dx1)/r2**3 )
        A12 = b + m2*( t1*(t1*dt1dx2 + t2*dt2dx2)/r2**3 )
        A21 = b + m2*( t2*(t1*dt1dx1 + t2*dt2dx1)/r2**3 )
        A22 = 1 - c - m2*( dt2dx2/r2 - t2*(t1*dt1dx2 + t2*dt2dx2)/r2**3 )


        self.detA = A11*A22 - A12*A21
        # print(self.detA)
        self.mu = 1/self.detA
        print(np.abs(self.mu).max())
        print(np.abs(self.mu).min())
        
    def gen_critical_lines1(self,pltlens=False,issqure=False,anno = False, xlim=(-1,1),ylim=(-1,1),num=100,threshold = 0.985, magfunc=None, title="Critical lines",xylim = None):
        thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num)
        # self.comp_magnification(thetax, thetay)
        magfunc(thetax, thetay)
        # criticalx, criticaly = self.findcriticall0(thetax,thetay)
        self.criticalx, self.criticaly = self.findcritical(threshold,thetax,thetay)
        # print(len(self.criticalx))
        spescatter(self.criticalx*rad2arcsec,self.criticaly*rad2arcsec,r"$\theta_1$/arcsec",r"$\theta_2$/arcsec",title=title,issqure=issqure,xylim=xylim)
        if pltlens:
            markersize = 100  
            plt.scatter([self.lens1.locx*rad2arcsec],
                [self.lens1.locy*rad2arcsec],s=markersize*(self.D2s/self.D1s),c="g",marker="x")
            plt.scatter([self.lens2.locx*rad2arcsec],
                [self.lens2.locy*rad2arcsec],s=markersize,c="g",marker="x")
        if anno:
            arrowprops = dict(
                arrowstyle = "->",
                connectionstyle = "arc3")#dict(facecolor='black', shrink=0.05)

            offset = 50
            fontsize = 14
            plt.annotate("Lens1 at Dd/Ds = {:.2f}\n (closer to the observer)".format(self.lens1.D/self.D_s),
                xy=(self.lens1.locx*rad2arcsec,self.lens1.locy*rad2arcsec),
                xytext=(-1*offset, -1*offset)
                ,arrowprops=arrowprops,textcoords='offset points',color="green",fontsize=fontsize)#
            plt.annotate("Lens2 at Dd/Ds = {:.2f}".format(self.lens2.D/self.D_s),
                xy=(self.lens2.locx*rad2arcsec,self.lens2.locy*rad2arcsec),
                xytext=(-1*offset, 1*offset)
                ,arrowprops=arrowprops,textcoords='offset points',color="green",fontsize=fontsize)


    def gen_caustic_lines_diffz_v2(self,pltlens=False,issqure=False,xylim = None,title = "Caustics",xlabel=r"x (arcsec)",ylabel=r"y (arcsec)"):
        theta1, theta2 = self.lens1.theta, self.lens2.theta
        rhox, rhoy = self.lens2.locx, self.lens2.locy
        x = self.criticalx
        y = self.criticaly
        r1 = np.sqrt(x**2+y**2)
        m1 = self.D1s*theta1
        m2 = self.D2s*theta2
        beta = self.D12/(self.lens2.D*self.D1s)
        t1 = x - m1*beta*x/r1 - rhox
        t2 = y - m1*beta*y/r1 - rhoy
        r2 = np.sqrt(t1**2+t2**2)

        self.causticsx = x - m1*x/r1 - m2*t1/r2
        # self.causticsy = y = m1*y/r1 - m2*t2/r2  # such a bug!!!
        self.causticsy = y - m1*y/r1 - m2*t2/r2
        markersize = 100
        spescatter(self.causticsx*rad2arcsec,self.causticsy*rad2arcsec,xlabel=xlabel,ylabel=ylabel,
            title=title,issqure=issqure,xylim = xylim,c="r")
        if pltlens:
            markersize = 100  
            plt.scatter([self.lens1.locx*rad2arcsec],
                [self.lens1.locy*rad2arcsec],s=markersize*(self.D2s/self.D1s),c="g",marker="x")
            plt.scatter([self.lens2.locx*rad2arcsec],
                [self.lens2.locy*rad2arcsec],s=markersize,c="g",marker="x")


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
    loc1 = (0,0) # in arcsec
    loc2 = (3,0)
    #lens1 closer to the source
    lens1 = oneSIS(sigma_v=200, D=0.5, loc=loc1)
    lens2 = oneSIS(sigma_v=200, D=0.5, loc=loc2)
    lenses = twoSISes(lens1,lens2, D_s = 2)
    # maxloc = max([abs(loc1[0]), abs(loc1[1]), abs(loc2[0]), abs(loc2[1])])
    # mx = max(lens1.theta+maxloc*arcsec2rad,lens2.theta+maxloc*arcsec2rad)*2.5
    # xlim = (-mx,mx)
    # ylim = xlim
    coe = 3
    xlim = (-lens1.theta*coe, abs(loc2[0])*arcsec2rad+lens2.theta*coe)
    ylim = (-lens1.theta*coe, lens2.theta*coe+abs(loc2[1])*arcsec2rad)

    print("xlim in arcsec",xlim[0]*rad2arcsec,xlim[1]*rad2arcsec)
    print("ylim in arcsec",ylim[0]*rad2arcsec,ylim[1]*rad2arcsec)
    plt.subplot(121)
    crititle = "The Critical Curves(blue) of two SIS lens(green) located at different redshift\n"+r"loc: ({:.1f},{:.1f}),({:.1f},{:.1f}) and $\sigma_v$: {},{} in $km/s$".format(lenses.lens1.locx*rad2arcsec,
        lenses.lens1.locy*rad2arcsec,lenses.lens2.locx*rad2arcsec,lenses.lens2.locy*rad2arcsec,lenses.lens1.sigma_v, lenses.lens2.sigma_v)
    cau_title = "\n and the corresponding Caustics(red) on source plane"
    lenses.gen_critical_lines1(issqure=True,pltlens=True,title=crititle+cau_title,anno = True,xlim=xlim,ylim=xlim,
        threshold = 500,num=8000,magfunc=lenses.comp_magnification_diffz_v2)
    plt.grid()
    lenses.gen_caustic_lines_diffz_v2(issqure=True,title=crititle+cau_title,xlabel=r"x (arcsec)",ylabel=r"y (arcsec)")
    plt.subplot(122)
    lenses.gen_caustic_lines_diffz_v2(issqure=True,title=cau_title,xlabel=r"x (arcsec)",ylabel=r"y (arcsec)")
    plt.show()

    # crititle = "The Critical Curves(blue) of two SIS lens(green) located at different redshift\n"+r"loc: ({:.1f},{:.1f}),({:.1f},{:.1f}) and $\sigma_v$: {},{} in $km/s$".format(lenses.lens1.locx*rad2arcsec,
    #     lenses.lens1.locy*rad2arcsec,lenses.lens2.locx*rad2arcsec,lenses.lens2.locy*rad2arcsec,lenses.lens1.sigma_v, lenses.lens2.sigma_v)
    # cau_title = "\n and the corresponding Caustics(red) on source plane"
    # lenses.gen_critical_lines1(issqure=True,pltlens=True,title=crititle+cau_title,anno = True,xlim=xlim,ylim=xlim,
    #     threshold = 500,num=8000,magfunc=lenses.comp_magnification_diffz_v2)
    # lenses.gen_caustic_lines_diffz_v2(issqure=True,title=crititle+cau_title,xlabel=r"x (arcsec)",ylabel=r"y (arcsec)")
    # plt.grid()
    # plt.show()