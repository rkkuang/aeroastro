'''
Two (N=2) galaxies are aligned perfectly with the Earth and a distant quasar. 
Each galaxy can be modelled as a singular isothermal sphere.
How many Einstein rings are formed as a result? 

In this python script I will check the situation of N = 2,

The distences from observer to lens2, from lens2 to lens1, from lens1 to source are d3, d2, d1, respectively,
and alpha1 is the deflection algle of Lens1, (not effective deflection angles)
alpha2 is the deflection algle of Lens2, 

As shown in geo.jpg
h1/h2 is the distance of the first/second turning point to the Lens1/Lens2

the idea is solving the equation(under small angle assumption):
for situation1:
    [1/d1, 1/d3; 1/d1, 1/d3-1/d2][h1, h2] = [alpha1+alpha2, alpha2]# in matlab manner
for situation2:
    [1/d1, 1/d3; 1/d1+1/d2, 1/d2][h3, h4] = [alpha1-alpha2, alpha1]

given a set of d1,d2,d3,alpha1,alpha2, we can solve out h1, h2 and h3, h4, respectively from above equations.

The physical meaningful solotion should satisfy: h1 > h2 > 0 ; h3, h4 > 0 ; and h1, h2, h3, h4 should be small
to validate the small angle assumption, e.g. 
h1 << d1
h2 << d3
(h1 - h2) << d2

h4 << d3
h3 << d1
(h3 + h4) << d2

if for a given d1,d2,d3,alpha1,alpha2, both (h1,h2) and (h3,h4) are exist, then there will be 2 Einstein Rings,
otherwise there will be only one.


How will your results generalize when you have N>2 galaxies?
if there are N galaxies, we can consider how many times j the light will cross the line of sight befor reach the observer:
if j = 0, there are 1 possibility,
if j = 1, there are C_{N-1}^{1} possibilities,
if j = 2, there are C_{N-1}^{2} possibilities,
...
if j = k, there are C_{N-1}^{k} possibilities,
, so depends on how those N galaxies are placed and the dispersion velocities, 
the maximun number of Einstein rings will be 2^{N-1}.


'''

# all d1, d2, d3 are normalized to d1, so I will fix d1 = 1
import numpy as np
alpha0 = 1.4 * np.pi / (3600*180)# 1.4 arcsec for a SIS model with dispersion 220km/s
# https://en.wikipedia.org/wiki/Velocity_dispersion


d1 = 1
num = 2*1e1
D2 = np.linspace(0.1, 4, num)
D3 = D2
# Soltype = np.zeros(D2.shape[0])
# 0 for no solution, 1 for has one solution for situation1, 2 for has one solution for situation2
# 3 for has both solution
Alpha1 = np.linspace(0.01*alpha0, 4*alpha0, num)
Alpha2 = Alpha1
only_sol_1 = []
only_sol_2 = []
both_sols = []
for d2 in D2:
    for d3 in D3:
        for alpha1 in Alpha1:
            for alpha2 in Alpha2:
                # [1/d1, 1/d3; 1/d1, 1/d3-1/d2][h1, h2] = [alpha1+alpha2, alpha2]
                h1h2 = np.linalg.solve(np.array([[1/d1, 1/d3],[1/d1, 1/d3-1/d2]]),np.array([alpha1+alpha2, alpha2]))
                # [1/d1, 1/d3; 1/d1+1/d2, 1/d2][h3, h4] = [alpha1-alpha2, alpha2]
                h3h4 = np.linalg.solve(np.array([[1/d1, 1/d3],[1/d1+1/d2, 1/d2]]),np.array([alpha1-alpha2, alpha1]))
                issol1 = h1h2[0] > h1h2[1] and h1h2[1] > 0
                issol2 = h3h4[0] > 0 and h3h4[1] > 0
                if issol1 and not issol2:
                    only_sol_1.append((d2,d3,alpha1,alpha2,h1h2[0],h1h2[1]))
                    # print("sol1")
                if not issol1 and issol2:
                    only_sol_2.append((d2,d3,alpha1,alpha2,h3h4[0],h3h4[1]))
                    # print("sol2")
                if issol1 and issol2:
                    both_sols.append((d2,d3,alpha1,alpha2,h1h2[0],h1h2[1],h3h4[0],h3h4[1]))
                    # print("both sols")
print("Find total {} solutions out of {} experiments: only situation1 {}, only situation2 {} and both situations {}"
    .format((len(only_sol_1)+len(only_sol_2)+len(both_sols)),len(D2)*len(D3)*len(Alpha1)*len(Alpha2),len(only_sol_1),len(only_sol_2),len(both_sols)))
np.save("resdata/only_sol_1.npy",only_sol_1)
np.save("resdata/only_sol_2.npy",only_sol_2)
np.save("resdata/both_sols.npy",both_sols)
for bsol in both_sols:
    # print(bsol)
    d2,d3,alpha1,alpha2,h1,h2,h3,h4 = bsol
    print("d2:{:.2f}, d3:{:.2f}, alpha1:{:.2f} arcsec, alpha2:{:.2f} arcsec, h1:{:.2e}, h2:{:.2e}, h3:{:.2e}, h4:{:.2e}"
        .format(d2,d3,alpha1*(3600*180)/np.pi,alpha2*(3600*180)/np.pi,h1h2[0],h1,h2,h3,h4))

# run: python3 n2_einstain_rings.py >> resdata/stdout.txt

# Find total 78215 solutions out of 160000 experiments: only situation1 59951, only situation2 15920 and both situations 2344

'''
>>> a = np.array([[1,1],[1,-1]])
>>> b = np.array([3,1])
>>> np.linalg.solve(a,b)
array([2., 1.])

2.numpy.save("filename.npy",a)

利用这种方法，保存文件的后缀名字一定会被置为.npy，这种格式最好只用

numpy.load("filename")来读取。


3.numpy.savetxt("filename.txt",a)

b =  numpy.loadtxt("filename.txt")

用于处理一维和二维数组
'''
