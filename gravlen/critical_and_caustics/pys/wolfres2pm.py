import numpy as np
from utils import *

q = -0.01
X = 0.3


mu1 = q/(1+q)
mu2 = 1/(1+q)

def abc(mu1, mu2, X, r):
    a = 16*X**2*r**2*(r**4-mu2)
    b = 8*r*X*(mu1*mu2*r**2 - (r**2+4*X**2)*(r**4-mu2**2) )
    c = (r**2+4*X**2)*(r**4-mu2**2)-mu1**2*r**4-2*mu1*mu2*r**2*(r**2-4*X**2)
    insqrt = b**2 - 4*a*c
    return a,b,c,insqrt

def radi2car(r, phi):
    x = np.multiply(r, np.cos(phi*np.pi/180))
    y = np.multiply(r, np.sin(phi*np.pi/180))
    return x,y


num = 20000

# np.concatenate((r,r,r,r),axis=0)

rlim1 = 1.00291
rlim2 = 1.00669
rlim3 = 0.800514
rlim4 = 1.00672
case1r = np.concatenate((np.linspace(rlim1, rlim1, 1),np.linspace(rlim2, 1.1, num)),axis=0)
case2r = np.concatenate((np.linspace(rlim3, rlim1, 1),np.linspace(rlim2, rlim4, num)),axis=0)

a1,b1,c1,insqrt1 = abc(mu1, mu2, X, case1r)
cosphi_pos1 = 0.5*(-b1-np.sqrt(insqrt1))/a1

phi1_1 = np.arccos(cosphi_pos1)*180/np.pi # in degree
phi1_2 = -phi1_1

a2,b2,c2,insqrt2 = abc(mu1, mu2, X, case2r)
cosphi_pos2 = 0.5*(-b2-np.sqrt(insqrt2))/a2

phi2_1 = np.arccos(cosphi_pos2)*180/np.pi # in degree
phi2_2 = -phi2_1

r = np.concatenate((case1r,case1r,case2r,case2r),axis=0)
phi = np.concatenate((phi1_1, phi1_2, phi2_1, phi2_2),axis=0)

rx, ry = radi2car(r, phi)

rx -= X

spescatter(rx,ry,
    xlabel = "x",ylabel="y",title="title")
plt.scatter([-X,X],[0,0],s=50,c="g",marker="x")

plt.show()