'''
ref: The two-point-mass lens - Detailed investigation of a special asymmetric gravitational lens
http://adsabs.harvard.edu/abs/1986A%26A...164..237S

by considering the lens equation as a 2-d "Lagrangian mapping"

two parameter (mass ratio and distance between the masses) family of models is well suited to study the behaviour of an asymmetrical gravitational lens.

notice in original paper there is a typo in eq (14)

Prof's homepage. https://astro.uni-bonn.de/de/m/peter
and his book: https://www.springer.com/gp/book/9783642540820#aboutAuthors


'''
import numpy as np 
import matplotlib.pyplot as plt
from utils import *




class GravLens():
    def __init__(self, massratio=0.5, X=0.1):
        '''
        massratio: mu1/mu2
        '''
        self.massratio = massratio
        self.X = X

    def gen_critical_lines(self, num = 100,title="title",xlabel="x",ylabel="y",xylim=2):
        if self.massratio == 0.5:
            if self.X**2 <= 0.125:
                # print("hello")
                ra = np.sqrt(0.5 - self.X**2 + np.sqrt(0.25 - 2*self.X**2)) # eq 13a
                rb = np.sqrt(0.5 - self.X**2 - np.sqrt(0.25 - 2*self.X**2)) # eq 13b
                rc = np.sqrt(-0.5 - self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13c
                re = np.sqrt(0.5 + self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13e
                # rc< r < rb, ra < r < re
                r_1 = np.linspace(rc, rb, num)
                r_2 = np.linspace(ra, re, num)
                r = np.concatenate((r_1, r_2),axis=0)
                specialr1 = np.array([re, 0, 0, 0])
                specialr2 = np.array([0, ra, rb, rc])
            elif self.X**2 <= 1:
                # print("hello")
                rc = np.sqrt(-0.5 - self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13c
                re = np.sqrt(0.5 + self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13e
                r = np.linspace(rc, re, num)

                specialr1 = np.array([re, 0])
                specialr2 = np.array([0, rc])
            else:
                rd = np.sqrt(0.5 + self.X**2 - np.sqrt(0.25 + 2*self.X**2)) # eq 13d
                re = np.sqrt(0.5 + self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13e
                r = np.linspace(rd, re, num)
                specialr1 = np.array([rd])
                specialr2 = np.array([0])

            cosphi2 = 1/(4*r**2*self.X**2)*(0.5+(r**2+self.X**2)**2-np.sqrt(0.25+2*(r**4+self.X**4)))
            # print(cosphi2)
            cosphi_pos = np.sqrt(cosphi2)
            phi1 = np.arccos(cosphi_pos)*180/np.pi # in degree
            phi2 = -phi1
            cosphi_neg = -cosphi_pos
            phi3 = np.arccos(cosphi_neg)*180/np.pi
            phi4 = -phi3
            # print(phi1[0:3])
            # print(phi2[0:3])
            # print(phi3[0:3])
            # print(phi4[0:3])
            # phi = np.concatenate((phi1,phi2,phi3,phi4),axis=1)
            # print(phi)
            rx1, ry1 = self.radi2car(r, phi1)
            rx2, ry2 = self.radi2car(r, phi2)
            rx3, ry3 = self.radi2car(r, phi3)
            rx4, ry4 = self.radi2car(r, phi4)



            self.r = np.concatenate((r,r,r,r),axis=0)
            # self.phi = np.concatenate((phi1,phi2,phi3,phi4),axis=0)
            self.r1 = np.concatenate((rx1,rx2,rx3,rx4),axis=0)
            self.r2 = np.concatenate((ry1,ry2,ry3,ry4),axis=0)

            #plot, r"$r_1=r\cos\phi$",r"$r_2=r\sin\phi$"
            spescatter(np.concatenate((self.r1,specialr1),axis=0),np.concatenate((self.r2,specialr2),axis=0),
                xlabel = xlabel,ylabel=ylabel,title=title,xylim=xylim)
            plt.scatter([-self.X,self.X],[0,0],s=100,c="g",marker="x")

    def gen_caustic_lines(self, num = 100,title="title",xlabel="x",ylabel="y",xylim=2):
        if self.massratio == 0.5:
            if not len(self.r1):
                print("Please generate critical lines first by calling gen_critical_lines()")
            else:
                if self.X**2 <= 0.125:
                    xe2 = self.X**2-0.5+0.125/self.X**2*(np.sqrt(1+8*self.X**2)-1)
                    xa2 = 0.125/self.X**2*(1-np.sqrt(1-8*self.X**2))-0.5-self.X**2
                    xb2 = 0.125/self.X**2*(1+np.sqrt(1+8*self.X**2))-0.5-self.X**2
                    xc2 = 0.125/self.X**2*((1+np.sqrt(1+8*self.X**2))**1.5+1-20*self.X**2-8*self.X**4)
                    specialx1 = np.array([np.sqrt(xe2),0,0,0])
                    specialx2 = np.array([0,-np.sqrt(xa2), -np.sqrt(xb2), -np.sqrt(xc2)])
                elif self.X**2 <= 1:
                    xe2 = self.X**2-0.5+0.125/self.X**2*(np.sqrt(1+8*self.X**2)-1)
                    xc2 = 0.125/self.X**2*((1+np.sqrt(1+8*self.X**2))**1.5+1-20*self.X**2-8*self.X**4)
                    specialx1 = np.array([np.sqrt(xe2),0])
                    specialx2 = np.array([0, -np.sqrt(xc2)])
                else:
                    xd2 = self.X**2-0.5-0.125/self.X**2*(np.sqrt(1+8*self.X**2)+1) # eq 15d
                    specialx1 = np.array([np.sqrt(xd2)])
                    specialx2 = np.array([0])

                # rcosphi = np.multiply(self.r,np.cos(self.phi*np.pi/180))
                # rsinphi = np.multiply(self.r,np.sin(self.phi*np.pi/180))
                rcosphi = self.r1
                rsinphi = self.r2
                r2 = np.multiply(self.r,self.r)
                # r2 = (np.multiply(rcosphi,rcosphi)+np.multiply(rsinphi,rsinphi))
                dem1 = r2+self.X**2-2*self.X*rcosphi
                dem2 = r2+self.X**2+2*self.X*rcosphi
                # notice in original paper there is a typo in eq (14):
                # self.x1 = rcosphi - 0.5*(rcosphi-self.X)/(dem1)-0.5*(rcosphi-self.X)/(dem2) # original paper eq (14)
                self.x1 = rcosphi - 0.5*(rcosphi-self.X)/(dem1)-0.5*(rcosphi+self.X)/(dem2) 
                self.x2 = rsinphi - 0.5*(rsinphi)/(dem1)-0.5*(rsinphi)/(dem2)

                #plot
                x1=np.concatenate((self.x1,specialx1),axis=0)
                x2=np.concatenate((self.x2,specialx2),axis=0)
                spescatter(x1,x2,
                    xlabel = xlabel,ylabel=ylabel, title=title,c="r",xylim = xylim)

                # r"The Caustic Lines of two point mass lens with $\mu_1=\mu_2$ and $X = {}$".format(self.X)
                # plt.scatter([-self.X,self.X],[0,0],s=50,c="r",marker="+")

    def radi2car(self, r, phi):
        x = np.multiply(r, np.cos(phi*np.pi/180))
        y = np.multiply(r, np.sin(phi*np.pi/180))
        return x,y



if __name__ == '__main__':
    mu = 0.5 # mass ratio
    X = .35 # position in arcsec
    glens = GravLens(mu, X)#8**(-0.5001)
    # plt.subplot(121)
    cri_title = r"The Critical Curves(blue) of two point mass(green) lens with $\mu_1=\mu_2$ and $X = {:.2f}$".format(glens.X)
    glens.gen_critical_lines(10000,title="title",xylim=max(1.5,2*X))
    # plt.subplot(122)
    cau_title = "\n"+r"and the corresponding Caustics(red) on source plane"
    final_title = cri_title+cau_title
    glens.gen_caustic_lines(title=final_title,xlabel=r"x (arcsec)",ylabel=r"y (arcsec)",xylim=max(1.5,2*X))
    plt.show()
    print("done")



