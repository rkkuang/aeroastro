'''
Draw the critical lines and caustics of two SIS lenses,
loacted at the same redshift



'''

import numpy as np 
import matplotlib.pyplot as plt
from two_asym_lens import spescatter

class oneSIS():
    def __init__(self, sigma_v=200, z=1, D_ds=1, D_s=2, loc=0.8):
        self.sigma_v = sigma_v #in km/s
        self.z = z
        self.D_ds = D_ds
        self.D_s = D_s
        self.theta_E = 1.6*(self.sigma_v/200)**2*(self.D_ds/self.D_s) #in arcsec
        self.theta_E = self.theta_E/3600/180*np.pi # in rad
        # self.loc = loc #in srcsec, location at x axis of lens plane
        self.loc = loc/3600/180*np.pi # !!!
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
    def comp_magnification(self,x,y):
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
    def gen_critical_lines1(self,xlim=(-1,1),ylim=(-1,1),num=100,threshold = 0.985):
        thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=num)
        self.comp_magnification(thetax, thetay)
        # criticalx, criticaly = self.findcriticall0(thetax,thetay)
        self.criticalx, self.criticaly = self.findcritical(threshold,thetax,thetay)
        # print(criticalx)
        # print(criticaly)
        # title = "The Critical Curves of two SIS lens located at\n"+r"({:.1f},0),Ds/Dds={:.1f};({:.1f},0),Ds/Dds={:.1f} with $\theta_E=${:.1f};{:.1f} in arcsec".format(self.lens1.loc*180/np.pi*3600,self.lens1.D_s/self.lens1.D_ds,
        #     self.lens2.loc*180/np.pi*3600,self.lens2.D_s/self.lens2.D_ds,self.lens1.theta_E*180/np.pi*3600, self.lens2.theta_E*180/np.pi*3600)
        title = "The Critical Curves of two SIS lens located at the same redshift\n"+r"({:.1f},0),({:.1f},0) with $\theta_E=${:.1f};{:.1f} in arcsec".format(self.lens1.loc*180/np.pi*3600,
            self.lens2.loc*180/np.pi*3600,self.lens1.theta_E*180/np.pi*3600, self.lens2.theta_E*180/np.pi*3600)

        spescatter(self.criticalx*180/np.pi*3600,self.criticaly*180/np.pi*3600,r"$\theta_1$/arcsec",r"$\theta_2$/arcsec",title)
        plt.scatter([self.lens1.loc*180/np.pi*3600,self.lens2.loc*180/np.pi*3600],[0,0],s=50,c="r",marker="+")
    def gen_caustic_lines(self):
        # beta = theta - alpha(theta)
        #beta = theta - thetaE*vectheta/|vectheta|
        r = np.sqrt(self.criticalx**2+self.criticaly**2)
        alphax = self.lens1.theta_E*self.criticalx/r+self.lens2.theta_E*self.criticalx/r
        alphay = self.lens1.theta_E*self.criticaly/r+self.lens2.theta_E*self.criticaly/r
        self.causticsx = self.criticalx - alphax
        self.causticsy = self.criticaly - alphax
        spescatter(self.causticsx*180/np.pi*3600,self.causticsy*180/np.pi*3600,r"$\beta_1$/arcsec",r"$\beta_2$/arcsec",
            "The corresponding Caustics on source plane")
        plt.scatter([self.lens1.loc*180/np.pi*3600,self.lens2.loc*180/np.pi*3600],[0,0],s=50,c="r",marker="+")

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
    loc1 = 1.5 # in arcsec
    loc2 = -1.5
    lens1 = oneSIS(sigma_v=300, z=1, D_ds=1, D_s=2, loc=loc1)
    lens2 = oneSIS(sigma_v=200, z=1, D_ds=1, D_s=2, loc=loc2)
    lenses = twoSISes(lens1,lens2)
    mx = max(lens1.theta_E+max(abs(loc1),abs(loc2))/3600/180*np.pi,lens2.theta_E+max(abs(loc1),abs(loc2))/3600/180*np.pi)*2.5
    xlim = (-mx,mx)
    ylim = xlim
    print("xlim in arcsec",xlim[0]*180/np.pi*3600,xlim[1]*180/np.pi*3600)
    plt.subplot(121)
    lenses.gen_critical_lines1(xlim=xlim,ylim=xlim,threshold = 5e3,num=6000)
    plt.subplot(122)
    lenses.gen_caustic_lines()
    plt.show()
