from telescope import Telescope
import numpy as np
import matplotlib.pyplot as plt
from utils import plot, imgfft, readimg, imgifft, sample_visibility
import random
import scipy.optimize as opt
from cleaner import Cleaner

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
# tele1.dirty_beam = np.fft.ifftshift(tele1.dirty_beam) --> added into gen_dirty_beam() function
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

cl = Cleaner(lenna_map_mag, tele1.dirty_beam)
plt.figure()
plt.subplot(221)
# plt.matshow(cl.model_beam, cmap=plt.cm.gist_earth_r)
plot(cl.model_beam, title = "Model beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)
itertime = 500
loop_gain = 0.2
cl.clean(itertime = itertime, loop_gain = loop_gain)
# plt.imshow(cl.residual)
plt.subplot(222)
plot(cl.residual, title = "Residual after {} clean steps with loop gain {}".format(itertime, loop_gain),xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

plt.subplot(223)
plot(cl.model, title = "Final model",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

# plt.figure()
# plot(cl.dirty_beam, title = "Dirty beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

cl.add_residual()
plt.subplot(224)
plot(cl.cleandmap, title = "Cleaned map",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)


plt.show()



# worked: 
# from gaussfitter import *
# Create the gaussian data
# coresize = 512
# Xin, Yin = np.mgrid[0:coresize, 0:coresize]
# # data = gaussian(3, 100, 100, 20, 40)(Xin, Yin) + np.random.random(Xin.shape)
# # data = tele1.dirty_beam[206:206+coresize, 206: 206+coresize]
# data = tele1.dirty_beam
# # plt.figure()
# plt.matshow(data, cmap=plt.cm.gist_earth_r)

# params = fitgaussian(data)
# fit = gaussian(*params)
# fit_res = fit(*np.indices(data.shape))
# plt.contour(fit_res, cmap=plt.cm.copper)
# ax = plt.gca()
# (height, x, y, width_x, width_y) = params

# plt.text(0.95, 0.05, """
# x : %.1f
# y : %.1f
# width_x : %.1f
# width_y : %.1f""" %(x, y, width_x, width_y),
#         fontsize=16, horizontalalignment='right',
#         verticalalignment='bottom', transform=ax.transAxes)

# plt.matshow(fit_res, cmap=plt.cm.gist_earth_r)


# plt.show()

# from temp import twoD_Gaussian
# # 2d gaussian fitting 
# # Create x and y indices
# coresize = 100

# x = np.linspace(1, coresize, coresize)
# y = np.linspace(1, coresize, coresize)
# x, y = np.meshgrid(x, y)

# #create data
# # def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
# # data = twoD_Gaussian((x, y), 40000, 100, 100, 20, 40, 0, 0)
# data = tele1.dirty_beam[206:206+coresize, 206: 206+coresize]
# # print(data.shape)
# # input(">>>>>>>>>>>>>")
# # data = data.reshape(RC[0]*RC[1],1)
# data = data.reshape(-1,1)
# # print(data.shape)
# # input(">>>>>>>>>>>>>")
# # plot twoD_Gaussian data generated above
# plt.figure()
# plt.imshow(data.reshape(coresize,coresize))
# plt.colorbar()


# # plt.show()
# # input(">>>>>>>>>>>")

# # add some noise to the data and try to fit the data generated beforehand
# initial_guess = (40000,0,0,5,5,0,1)#(3,100,100,20,40,0,10)
# #amplitude, xo, yo, sigma_x, sigma_y, theta, offset

# # data_noisy = data + 0.2*np.random.normal(size=data.shape)
# data_noisy = data

# popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_noisy, p0=initial_guess)

# data_fitted = twoD_Gaussian((x, y), *popt)

# fig, ax = plt.subplots(1, 1)
# # ax.hold(True)
# ax.imshow(data_noisy.reshape(coresize,coresize), cmap=plt.cm.jet, origin='bottom',
#     extent=(x.min(), x.max(), y.min(), y.max()))
# ax.contour(x, y, data_fitted.reshape(coresize,coresize), 8, colors='w')
# plt.figure()
# plt.imshow(data_fitted.reshape(coresize,coresize))

# plt.show()




