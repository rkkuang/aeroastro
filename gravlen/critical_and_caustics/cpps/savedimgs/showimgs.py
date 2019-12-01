import numpy as np
# from utils import genxy, spescatter
import matplotlib.pyplot as plt

srcplaneIMG = plt.imread("srcplaneIMG.png")
imgplaneIMG = plt.imread("imgplaneIMG.png")
xlim = (-2.5,2.5)
ylim = (-2.5,2.5)
fig = plt.figure()
plt.subplot(122)
cmap = plt.cm.get_cmap('viridis') # 'Paired', viridis, gist_ncar, 
plt.imshow((srcplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
# title = "Two point mass lenses, src plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
title="source plane"
plt.title(title)
plt.colorbar()


plt.subplot(121)
cmap = plt.cm.get_cmap('viridis')
plt.imshow((imgplaneIMG), origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
#title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta={:.1f}$".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0] ,twolens.beta)
# title = "Two point mass lenses, img plane (log scale), with mass ratio:\n"+r" $M_1/M_2=${:.3f} and $x_1=${:.1f}, $x_2 =${:.1f}, $D1={:.2f}$, $D2={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], d1d2[0], d1d2[1])
title = "Lens plane"
plt.title(title)
plt.colorbar()

fig.set_size_inches(16,9)


plt.subplots_adjust(left=0.04, top = 0.9, bottom = 0.05, right=0.96, hspace = 0.13, wspace = 0.13)
plt.margins(0,0)

plt.show()