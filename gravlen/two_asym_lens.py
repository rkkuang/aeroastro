'''
ref: The two-point-mass lens - Detailed investigation of a special asymmetric gravitational lens
http://adsabs.harvard.edu/abs/1986A%26A...164..237S

by considering the lens equation as a 2-d "Lagrangian mapping"

two parameter (mass ratio and distance between the masses) family of models is well suited to study the behaviour of an asymmetrical gravitational lens.


'''
import numpy as np 
import matplotlib.pyplot as plt

# 设置横纵坐标的名称以及对应字体格式, https://blog.csdn.net/A_Z666666/article/details/81165123
font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 16,
         }

def spescatter(x,y,xlabel,ylabel,title):
    #plot
    plt.scatter(x,y,s=1,c="b")
    # plt.scatter([-self.X,self.X],[0,0],s=50,c="r",marker="+")
    plt.xlabel(xlabel,font2)
    plt.ylabel(ylabel,font2)
    plt.tick_params(labelsize=12)
    plt.grid()
    plt.title(title)
    # plt.axis('square')

class GravLens():
    def __init__(self, massratio=0.5, X=0.1,  ):
        '''
        massratio: mu1/mu2
        '''
        self.massratio = massratio
        self.X = X

    def gen_critical_lines(self, num = 100):
        if self.massratio == 0.5:
            if self.X**2 <= 0.125:
                print("hello")
                ra = np.sqrt(0.5 - self.X**2 + np.sqrt(0.25 - 2*self.X**2)) # eq 13a
                rb = np.sqrt(0.5 - self.X**2 - np.sqrt(0.25 - 2*self.X**2)) # eq 13b
                rc = np.sqrt(-0.5 - self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13c
                re = np.sqrt(0.5 + self.X**2 + np.sqrt(0.25 + 2*self.X**2)) # eq 13e
                # rc< r < rb, ra < r < re
                r_1 = np.linspace(rc, rb, num)
                r_2 = np.linspace(ra, re, num)
                r = np.append(r_1, r_2)
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
                x1, y1 = self.radi2car(r, phi1)
                x2, y2 = self.radi2car(r, phi2)
                x3, y3 = self.radi2car(r, phi3)
                x4, y4 = self.radi2car(r, phi4)

                specialr1 = np.array([re, 0, 0, 0])
                specialr2 = np.array([0, ra, rb, rc])

                self.r = np.concatenate((r,r,r,r),axis=0)
                self.phi = np.concatenate((phi1,phi2,phi3,phi4),axis=0)
                self.r1 = np.concatenate((x1,x2,x3,x4),axis=0)
                self.r2 = np.concatenate((y1,y2,y3,y4),axis=0)

                #plot
                spescatter(np.concatenate((self.r1,specialr1),axis=0),np.concatenate((self.r2,specialr2),axis=0),
                    r"$r_1=r\cos\phi$",r"$r_2=r\sin\phi$",r"The Critical Lines of two point mass lens with $\mu_1=\mu_2$ and $X = {:.2f}$".format(self.X))
                plt.scatter([-self.X,self.X],[0,0],s=50,c="r",marker="+")
    def gen_caustic_lines(self, num = 100):
        if self.massratio == 0.5:
            if self.X**2 <= 0.125:
                if not len(self.r1):
                    print("Please generate critical lines first by calling gen_critical_lines()")
                else:
                    # rcosphi = np.multiply(self.r,np.cos(self.phi*np.pi/180))
                    # rsinphi = np.multiply(self.r,np.sin(self.phi*np.pi/180))
                    rcosphi = self.r1
                    rsinphi = self.r2
                    r2 = np.multiply(self.r,self.r)
                    dem1 = r2+self.X**2-2*self.X*rcosphi
                    dem2 = r2+self.X**2+2*self.X*rcosphi
                    self.x1 = rcosphi - 0.5*(rcosphi-self.X)/(dem1)-0.5*(rcosphi-self.X)/(dem2)
                    self.x2 = rsinphi - 0.5*(rsinphi)/(dem1)-0.5*(rsinphi)/(dem2)
                    xe2 = self.X**2-0.5+0.125/self.X**2*(np.sqrt(1+8*self.X**2)-1)
                    xa2 = 0.125/self.X**2*(1-np.sqrt(1-8*self.X**2))-0.5-self.X**2
                    xb2 = 0.125/self.X**2*(1+np.sqrt(1+8*self.X**2))-0.5-self.X**2
                    xc2 = 0.125/self.X**2*((1+np.sqrt(1+8*self.X**2))**1.5+1-20*self.X**2-8*self.X**4)
                    specialx1 = np.array([np.sqrt(xe2),0,0,0])
                    specialx2 = np.array([0,-np.sqrt(xa2), -np.sqrt(xb2), -np.sqrt(xc2)])
                    #plot
                    spescatter(np.concatenate((self.x1,specialx1),axis=0),np.concatenate((self.x2,specialx2),axis=0),
                        r"$x_1$",r"$x_2$", "The corresponding Caustic Lines on source plane")
                    # r"The Caustic Lines of two point mass lens with $\mu_1=\mu_2$ and $X = {}$".format(self.X)
                    # plt.scatter([-self.X,self.X],[0,0],s=50,c="r",marker="+")

    def radi2car(self, r, phi):
        x = np.multiply(r, np.cos(phi*np.pi/180))
        y = np.multiply(r, np.sin(phi*np.pi/180))
        return x,y



if __name__ == '__main__':
    glens = GravLens(0.5, 8**(-0.5001))
    plt.subplot(121)
    glens.gen_critical_lines(10000)
    plt.subplot(122)
    glens.gen_caustic_lines()
    plt.show()



