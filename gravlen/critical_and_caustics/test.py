
import numpy as np
import matplotlib.pyplot as plt

ImgSize = (100,200)
IMG = np.ones(ImgSize)
IMG[10:20, 50:60]=0
print(IMG[99,99])
plt.imshow(IMG, origin='lower')
plt.colorbar()
plt.show()