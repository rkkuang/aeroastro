from telescope import Telescope
import numpy as np
import matplotlib.pyplot as plt
from utils import plot, imgfft, readimg, imgifft, sample_visibility, gen_point_source, img_conv
import random
import scipy.optimize as opt
from cleaner import Cleaner
from gaussfitter import gaussian


##############################Telescope1######################################
tele1 = Telescope([],[])
RC = (512,512)

numsites = 200
sites = []
dangle = 0.1
angle_arnge = 360
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

# lenna = readimg("../imgs/Lenna.png")
num_source = 30
point_source = gen_point_source(num_source, RC)

#####################################################
# point source conv with a known_kernel
kernelsize = (25,25)
Xin, Yin = np.mgrid[0:kernelsize[0], 0:kernelsize[1]]
# %gaussian(height, center_x, center_y, width_x, width_y)(Xin, Yin)
(height, center_x, center_y, width_x, width_y) = (1, kernelsize[0]//2, kernelsize[1]//2, 5, 5)
known_kernel = gaussian(height, center_x, center_y, width_x, width_y)(Xin, Yin)

point_source = img_conv(point_source, known_kernel)
#####################################################

noise = 0.01*np.random.normal(size=RC)

point_source += noise
plt.subplot(233)
plot(point_source, title = "{} point sources img cov with a gaussian kernel with some noise\n(True brightness distribution)".format(num_source),xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

lenna_visibility, lenna_visibility_mag = imgfft(point_source)
plt.subplot(234)
plot(lenna_visibility_mag, title = "FT of point source img (Theoretical visibility)",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)


measured_visibility, measured_lenna_visibility_mag = sample_visibility(lenna_visibility, lenna_visibility_mag, tele1.uvcover)

plt.subplot(235)
plot(measured_lenna_visibility_mag, title = "Measured visibility after sampling in uv plane",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)

# dirty_map = imgifft(measured_visibility)
# plot(dirty_map, title = "Dirty map",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

lenna_map, lenna_map_mag = imgifft(measured_visibility)
plt.subplot(236)
plot(lenna_map_mag, title = "Dirty map",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

########################################################################################################

#(coresize, height, center_x, center_y, width_x, width_y) = model_beam_para
model_beam_para = (512, 1, 256, 256, 5, 5)
cl = Cleaner(lenna_map_mag, tele1.dirty_beam, model_beam_para)

plt.figure()

plt.subplot(231)
plot(tele1.uvcover, title = "Fake uv coverage of {} sites".format(len(sites)),xlabel = "u (pixel)", ylabel = "v (pixel)", colorbar = True, cmap = plt.cm.jet)

plt.subplot(232)
plot(lenna_map_mag, title = "Dirty map (linear scale)",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False, cmap = plt.cm.jet)


plt.subplot(233)
# plt.matshow(cl.model_beam, cmap=plt.cm.gist_earth_r)
cl.gen_model_beam_byhand(coresize = 512, height=1, center_x=256, center_y=256, width_x=5, width_y=5)
# plot(cl.model_beam, title = "Model beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)
plot(point_source, title = "{} point sources img cov with a gaussian kernel with some noise\n(True brightness distribution)".format(num_source),xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False, cmap = plt.cm.jet)

itertime = 5000
loop_gain = 0.2
peak = 5000
# cl.clean(itertime = itertime, loop_gain = loop_gain, criteria = "max_itertime")

cl.clean(loop_gain = loop_gain, criteria = "peak", peak = peak)
# plt.imshow(cl.residual)
plt.subplot(234)
# plot(cl.residual, title = "Residual after {} clean steps with loop gain {} (linear scale)".format(itertime, loop_gain),xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True , islog = False)
plot(cl.residual, title = "Residual after peak < {} with loop gain {} (linear scale)".format(peak, loop_gain),xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False, cmap = plt.cm.jet)

plt.subplot(235)
plot(cl.model, title = "Final model",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, cmap = plt.cm.jet)

# plt.figure()
# plot(cl.dirty_beam, title = "Dirty beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

cl.add_residual()
plt.subplot(236)
plot(cl.cleandmap, title = "Cleaned map (linear scale)",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False, cmap = plt.cm.jet)


plt.show()