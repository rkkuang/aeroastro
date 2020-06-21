import numpy as np
import matplotlib.pyplot as plt
a = 0.8
b = 0.3

r = 0.6 # = sin(theta)
mu = (1-r**2)**0.5 # = cos(theta)

Gamma = 2*a/(3-a)
Lammb = 4*b/(5-b)

def f1(a, r):
    return 3/(3-a)*( 1-a*(1- (1-r**2)**0.5 ) )

def f2(b, r):
    return 5/(5-b)*( 1-b*(1- (1-r**2)**0.25 ) )
def model1(Gamma, Lammb, mu):
    return 1-Gamma*(1-1.5*mu)-Lammb*(1-1.25*mu**0.5)
def model2(a,b,r):
    return f1(a, r) + f2(b, r) - 1
def model3(Gamma, mu):
    return 1-Gamma*(1-1.5*mu)

print("model1: ", model1(Gamma, Lammb, mu))
print("model2: ", model2(a, b, r))


rs = np.linspace(0,1,1000)
mod1s = []
mod2s = []
mod3s = []

for r in rs:
    mu = (1-r**2)**0.5
    mod1s.append(model1(Gamma, Lammb, mu))
    mod2s.append(model2(a, b, r))
    mod2s.append(model3(Gamma, mu))

plt.plot(rs, mod1s, "-r", label="$S_\lambda(\mu), square root$")
plt.plot(rs, mod3s, "-g", label="$S_\lambda(\mu), linear$")
plt.plot(rs, mod2s, "-.b", label = "f(r)")
plt.xlabel("r")
plt.legend()
plt.text(0.6,1.35,"a = {}\nb = {}".format(a,b))


plt.show()