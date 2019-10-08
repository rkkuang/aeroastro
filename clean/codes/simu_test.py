from telescope import Telescope
import numpy as np
import matplotlib.pyplot as plt

tele1 = Telescope([],[])
RC = (500,500)
site1 = (100,0,120,0.2)
site2 = (200,240,120,0.2)
site3 = (150,60,120,0.2)
tele1.fake_uvcover(RC,(site1,site2,site3))
tele1.plot(tele1.uvcover, "Fake uv coverage of 3 sites","u (pixel)", "v (pixel)")
tele1.dirty_beam()
tele1.plot(tele1.dirty_beam, "Dirty beam of 3 sites","x (pixel)", "y (pixel)")
plt.show()





