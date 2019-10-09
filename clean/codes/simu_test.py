from telescope import Telescope
import numpy as np
import matplotlib.pyplot as plt
from utils import plot, imgfft, readimg, imgifft, sample_visibility
import random

# tele1 = Telescope([],[])
# RC = (500,500)
# site1 = (100,0,120,0.2)
# site2 = (200,240,120,0.2)
# site3 = (150,60,120,0.2)
# tele1.fake_uvcover(RC,(site1,site2,site3))
# plot(tele1.uvcover, title = "Fake uv coverage of 3 sites",xlabel = "u (pixel)",ylabel =  "v (pixel)")
# tele1.gen_dirty_beam()
# plot(tele1.dirty_beam, title = "Dirty beam of 3 sites",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

# back_uvcover, mag = imgfft(tele1.dirty_beam)
# # plt.imshow(mag,cmap='gray')
# plot(mag, title = "FT of 3 sites dirty beam",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False)


# tele1.fake_uvcover2(RC, 30)
# plot(tele1.uvcover, title = "Fake uv coverage corresponding to a Gaussian dirty beam",xlabel = "u (pixel)", ylabel = "v (pixel)",)
# tele1.gen_dirty_beam()
# plot(tele1.dirty_beam, title = "Gaussian dirty beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

# back_uvcover, mag = imgfft(tele1.dirty_beam)
# # plt.imshow(mag,cmap='gray')
# plot(mag, title = "FT of Gaussian dirty beam",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False)

# plt.show()

##############################Telescope2######################################

# tele1 = Telescope([],[])
# RC = (512,512)
# tele1.fake_uvcover2(RC, 90)

# plt.subplot(231)
# plot(tele1.uvcover, title = "Fake uv coverage corresponding to a Gaussian dirty beam",xlabel = "u (pixel)", ylabel = "v (pixel)", colorbar = True)
# tele1.gen_dirty_beam()
# tele1.dirty_beam = np.fft.ifftshift(tele1.dirty_beam)
# plt.subplot(232)
# plot(tele1.dirty_beam, title = "Gaussian dirty beam(after np.fft.ifftshift)",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

# lenna = readimg("../imgs/Lenna.png")
# plt.subplot(233)
# plot(lenna, title = "Lenna img",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

# lenna_visibility, lenna_visibility_mag = imgfft(lenna)
# plt.subplot(234)
# plot(lenna_visibility_mag, title = "FT of Lenna img",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)

# measured_visibility = np.zeros((RC[0],RC[1],2))
# measured_visibility[:,:,0] = lenna_visibility[:,:,0] * tele1.uvcover
# measured_visibility[:,:,1] = lenna_visibility[:,:,1] * tele1.uvcover
# measured_lenna_visibility_mag = lenna_visibility_mag * tele1.uvcover
# plt.subplot(235)
# plot(measured_lenna_visibility_mag, title = "Measured visibility after sampling in uv plane",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)

# # dirty_map = imgifft(measured_visibility)
# # plot(dirty_map, title = "Dirty map",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

# lenna_map, lenna_map_mag = imgifft(measured_visibility)
# plt.subplot(236)
# plot(lenna_map_mag, title = "IFT of LennaFT",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

##############################Telescope1######################################
tele1 = Telescope([],[])
RC = (512,512)

# site0 = (10,0,300,0.1)
# site1 = (100,0,300,0.2)
# site2 = (200,240,300,0.2)
# site3 = (150,60,300,0.2)
# tele1.fake_uvcover(RC,(site1,site2,site3))

# sites = (
# (0.2,0,360,0.1),
# (0.4,0,360,0.1),
# (0.6,0,360,0.1),
# (0.8,0,360,0.1),
# (1,0,360,0.1),
# (1.2,0,360,0.1),
# (1.4,0,360,0.1),
# (1.6,0,360,0.1),
# (1.8,0,360,0.1),
# (2,0,360,0.1),
# (2.2,0,360,0.1),
# (2.4,0,360,0.1),
# (2.6,0,360,0.1),
# (2.8,0,360,0.1),
# (3,0,360,0.1),
# (3.2,0,360,0.1),
# (3.4,0,360,0.1),
# (3.6,0,360,0.1),
# (3.8,0,360,0.1),
# (5,0,360,0.1),
# (5.2,0,360,0.1),
# (5.4,0,360,0.1),
# (5.6,0,360,0.1),
# (5.8,0,360,0.1),
# (7,0,360,0.1),
# (8,0,360,0.1),
# (9,0,360,0.1),
# (10,0,300,0.1),
# (11,0,360,0.1),
# (12,0,360,0.1),
# (15,0,360,0.1),
# (20,0,360,0.1),
# (25,0,360,0.1),
# (30,0,360,0.1),
# (35,0,360,0.1),
# (40,0,360,0.1),
# (45,0,360,0.1),
# (50,0,360,0.1),
# (55,0,360,0.1),
# (60,0,360,0.1),
# (65,0,360,0.1),
# (70,0,360,0.1),
# (75,0,360,0.1),
# (80,0,360,0.1),
# (85,0,360,0.1),
# (90,0,360,0.1),
# (100,0,360,0.2),
# (200,240,360,0.2),
# (150,60,360,0.2)
#     )

numsites = 150
sites = []
dangle = 0.1
angle_arnge = 120
inner = 10
for i in range(int(numsites/3)):
    sites.append((random.random()*inner,360*random.random(),angle_arnge,dangle))
for i in range(numsites - int(numsites/3)):
    sites.append((random.random()*(RC[0]/2.1),360*random.random(),angle_arnge,dangle))

tele1.fake_uvcover(RC,sites)

plt.subplot(231)
plot(tele1.uvcover, title = "Fake uv coverage of {} sites".format(len(sites)),xlabel = "u (pixel)", ylabel = "v (pixel)", colorbar = True)
tele1.gen_dirty_beam()
tele1.dirty_beam = np.fft.ifftshift(tele1.dirty_beam)
plt.subplot(232)
plot(tele1.dirty_beam, title = "Dirty beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

lenna = readimg("../imgs/Lenna.png")
plt.subplot(233)
plot(lenna, title = "Lenna img (True brightness distribution)",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

lenna_visibility, lenna_visibility_mag = imgfft(lenna)
plt.subplot(234)
plot(lenna_visibility_mag, title = "FT of Lenna img (Theoretical visibility)",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)


measured_visibility, measured_lenna_visibility_mag = sample_visibility(lenna_visibility, lenna_visibility_mag, tele1.uvcover)

plt.subplot(235)
plot(measured_lenna_visibility_mag, title = "Measured visibility after sampling in uv plane",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)

# dirty_map = imgifft(measured_visibility)
# plot(dirty_map, title = "Dirty map",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

lenna_map, lenna_map_mag = imgifft(measured_visibility)
plt.subplot(236)
plot(lenna_map_mag, title = "Dirty map",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)



plt.show()



