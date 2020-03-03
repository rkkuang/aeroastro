from bi_inverse_raysht import *
import numpy as np
import sys
# s = 1 # separation between the two lenses in units of total ang. Einstein radii
# q = 0.5 # mass ratio: mass of the lens on the right divided by mass of the lens on the left
# rho = 0.01 # source radius in Einstein radii of the total mass.

s = float(sys.argv[1])
q = float(sys.argv[2])
rho = float(sys.argv[3])

xs = np.array([-q/(1+q), 1/(1+q)])*s
ys = np.array([0,0])

masses = np.array([1/(1+q),q/(1+q)]) # smaller mass is at origin
twolens = Nlenses(masses, xs, ys)

xlim, ylim = (-2,2), (-2,2)
pix = 2048
datatype = np.float32
ImgSize = (pix, pix)
raynum = int(sys.argv[4])#3000

srcplaneIMG_rayshootNomag, _, _ = twolens.get_imgs_lessmem_v3(ImgSize, xlim, ylim, raynum, datatype = datatype)

title="Ray{}s{}q{}rho{}".format(raynum, s, q, rho).replace(".","_")

# np.save("filename.npy",a) 
# b = np.load("filename.npy")

np.save(title+".npy",srcplaneIMG_rayshootNomag)
