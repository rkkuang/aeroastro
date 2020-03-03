'''
Draw the critical lines and caustics of the singular isotheral sphere model
lens plane: vevtheta = (theta1, theta2), r = norm(vectheta)

thetaE: Radius of the Einstein Ring
sigmav : dispersion velocity in km/s
thetaE = 1.6*(sigmav/200)**2*(D_ds/D_s) in arcsec
we can derive:
detA = theta_E*r**3-(theta_E**2+1)*r**2+2*theta_E*r-theta_E**2

solve detA = 0 ==> r, the critical lines

'''
import numpy as np 
import matplotlib.pyplot as plt
import sympy as sp
from two_asym_lens import spescatter
class oneSIS():
    def __init__(self, sigma_v, D_ds, D_s):
        self.sigma_v = sigma_v #in km/s
        self.D_ds = D_ds
        self.D_s = D_s
        self.theta_E = 1.6*(self.sigma_v/200)**2*(self.D_ds/self.D_s) #in arcsec
        self.theta_E = self.theta_E/3600/180*np.pi # in rad
        self.r = None
        self.theta1 = None
        self.theta2 = None
        self.beta1 = None
        self.beta2 = None
    def solve_r(self):
        r = sp.Symbol("r")
        f = self.theta_E*r**3-(self.theta_E**2+1)*r**2+2*self.theta_E*r-self.theta_E**2
        x = sp.solve(f)
        x = [complex(i).real for i in x]
        x.remove(max(x))
        self.r = x
    def gen_critical_lines(self, num = 3600):
        angs = np.linspace(0, 360, num)
        x1 = self.r[0]*np.cos(angs)
        y1 = self.r[0]*np.sin(angs)
        x2 = self.r[1]*np.cos(angs)
        y2 = self.r[1]*np.sin(angs)
        self.theta1 = np.concatenate((x1,x2),axis=0)*180/np.pi*3600
        self.theta2 = np.concatenate((y1,y2),axis=0)*180/np.pi*3600  
        spescatter(self.theta1,self.theta2,r"$\theta_1$/arcsec",r"$\theta_2$/arcsec",
            r"The Critical Curves of a SIS lens with $\sigma_v={:.1f}km/s$ and $Ds/Dds$={:.2f}".format(self.sigma_v, self.D_s/self.D_ds))
        # plt.xlim(min(x),max(x))
        # plt.ylim(min(y),max(y))
    def gen_caustic_lines(self): 
        R = np.sqrt(np.multiply(self.theta1,self.theta1)+np.multiply(self.theta2,self.theta2))
        k = self.theta_E/R
        self.beta1 = self.theta1 - np.multiply(k,self.theta1)
        self.beta2 = self.theta2 - np.multiply(k,self.theta2)
        spescatter(self.beta1,self.beta2,r"$\beta_1$/arcsec",r"$\beta_2$/arcsec",
            r"The corresponding Caustics on source plane")


if __name__ == '__main__':
    '''
    import sympy as sp    # 导入sympy包
    x = sp.Symbol('x')    # 定义符号变量
    f = x**3 - 3*x**2 + 3*x - 9/16   # 定义要求解的一元三次方程
    x = sp.solve(f)       # 调用solve函数求解方程
    '''
    lens = oneSIS(5000,1,2)#(sigma_v, D_ds, D_s), D_ds < D_s
    lens.solve_r()
    plt.subplot(121)
    lens.gen_critical_lines()
    plt.subplot(122)
    lens.gen_caustic_lines()
    plt.show()

