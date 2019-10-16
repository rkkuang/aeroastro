from telescope import Telescope
import numpy as np
import matplotlib.pyplot as plt
from utils import plot, imgfft, readimg, imgifft, sample_visibility, gen_point_source, img_conv, gen_site, plot_scatter_Baseline
import random
import scipy.optimize as opt
from cleaner import Cleaner
from gaussfitter import gaussian
# from matplotlib.ticker import FuncFormatter

##############################Telescope1######################################
# tele1 = Telescope([],[])
RC = (512,512)

LMT = gen_site("LMT", "-97:18:53", "18:59:06", -768713.9637, -5988541.7982, 2063275.9472, 15, 85, 560, 50)
PV = gen_site("PV", "-3:23:33.8", "37:03:58.2", 5088967.9, -301681.6, 3825015.8, 15, 85, 2900, 30)
ALMA50 = gen_site("ALMA50", "-67:45:11.4", "-23:01:09.4", 2225037.1851, -5441199.162, -2479303.4629, 15, 85, 110, 84.7)
ALMA50_1 = gen_site("ALMA50_1", "-67:45:11.4", "-23:01:09.4", -1827796.2, -5053406.8, 3428865.2, 15, 85, 110, 84.7)
ALMA50_12 = gen_site("ALMA50_1", "-67:45:11.4", "-23:01:09.4", -1826796.2, -5052406.8, 3429865.2, 15, 85, 110, 84.7)

ALMA50_2 = gen_site("ALMA50_2", "-67:45:11.4", "-23:01:09.4", 2227037.1851, -5442199.162, -2471303.4629, 15, 85, 110, 84.7)
ALMA50_3 = gen_site("ALMA50_3", "-67:45:11.4", "-23:01:09.4", 2228037.1851, -5443299.162, -2472303.4629, 15, 85, 110, 84.7)
ALMA50_4 = gen_site("ALMA50_4", "-67:45:11.4", "-23:01:09.4", 2229037.1851, -5444399.162, -2473303.4629, 15, 85, 110, 84.7)
ALMA50_5 = gen_site("ALMA50_5", "-67:45:11.4", "-23:01:09.4", 2224037.1851, -5445499.162, -2474303.4629, 15, 85, 110, 84.7)
ALMA50_6 = gen_site("ALMA50_6", "-67:45:11.4", "-23:01:09.4", 2223037.1851, -5446599.162, -2475303.4629, 15, 85, 110, 84.7)
ALMA50_7 = gen_site("ALMA50_7", "-67:45:11.4", "-23:01:09.4", 2222037.1851, -5447699.162, -2476303.4629, 15, 85, 110, 84.7)
ALMA50_8 = gen_site("ALMA50_8", "-67:45:11.4", "-23:01:09.4", 2221037.1851, -5448799.162, -2477303.4629, 15, 85, 110, 84.7)
ALMA50_9 = gen_site("ALMA50_9", "-67:45:11.4", "-23:01:09.4", 2220037.1851, -5449899.162, -2478303.4629, 15, 85, 110, 84.7)
ALMA50_10 = gen_site("ALMA50_10", "-67:45:11.4", "-23:01:09.4", 2230037.1851, -5440199.162, -2479303.4629, 15, 85, 110, 84.7)
ALMA50_11 = gen_site("ALMA50_11", "-67:45:11.4", "-23:01:09.4", 2231037.1851, -5441199.162, -2460303.4629, 15, 85, 110, 84.7)

#SMTO = gen_site("SMTO", "-109:52:19", "32:42:06", -1828796.2, -5054406.8, 3427865.2, 15, 85, 11900, 10)
SMT = gen_site("SMT", "-109:52:19", "32:42:06", -1828796.2, -5054406.8, 3427865.2, 15, 85, 11000, 10)
Hawaii8 = gen_site("Hawaii8", "-155:28:40.7", "19:49:27.4", -5464523.4, -2493147.08, 2150611.75, 15, 85, 4900, 20.8)
#PdBI = gen_site("PdBI", "05:54:28.5", "44:38:02.0", 4523998.4, 468045.24, 4460309.76, 15, 85, 1600, 36,7)
PDB = gen_site("PDB", "05:54:28.5", "44:38:02.0", 4523998.4, 468045.24, 4460309.76, 15, 85, 5200, 36.7)
SPT = gen_site("SPT", "-000:00:00.0", "-90:00:00", 0, 0, -6359587.3, 15, 85, 7300, 12)
GLT = gen_site("GLT", "72:35:46.4", "38:25:19.1", 1500692, -1191735, 6066409, 15, 85, 4744, 12)
CARMA8 = gen_site("CARMA8", "-118:08:30.3", "37:16:49.6", -2397431.3, -4482018.9, 3843524.5, 15, 85, 3500, 26.9)
SMA = gen_site("SMA", "-155:28:40.7", "19:49:27.4", -5464523.4, -2493147.08, 2150611.75, 15, 85, 4000, 20.8)
# Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT,CARMA8,SMA,PDB,ALMA50_1,ALMA50_2,ALMA50_3,ALMA50_4,ALMA50_5,ALMA50_6,ALMA50_7,ALMA50_8,ALMA50_9,ALMA50_10,ALMA50_11]
Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT,CARMA8,SMA,PDB,ALMA50_1,ALMA50_12]
# Array = [PV,ALMA50]


# Array = [LMT,PV,ALMA50,ALMA50_1,ALMA50_2,ALMA50_3,ALMA50_4,ALMA50_5,ALMA50_6,ALMA50_7,ALMA50_8,ALMA50_9,ALMA50_10,ALMA50_11]

# 用uvArray_Baseline 生成 UV 覆盖, 并用plot_scatter 画图(带legend)的测试程序
# Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT,CARMA8,SMA,PDB]
# Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT]
tele1 = Telescope([],[])
center = ALMA50#Array[0]
H = (0, 12, 1/2, 1/12, 1/3600)#表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每5秒采一次样)， 采完等待 10分钟，之后再采5分钟
Dec = 60 #Declination, 赤纬，度
freq = 227.297 * 10**9 # 1 GHz, 227297 MHz
UVW, Baselines_uv = tele1.uvArray_Baseline(Array=Array, center=[], H=H, Dec=Dec, freq=freq)
#uvArray_Baseline(self, Array=[], center=None, H=(0, 24, 1/12, 1/6,  6/3600), Del=60, freq=1*10**9):
# for site in resuvArray:
#     site["U"] = np.array(site["U"])/10**3
#     site["V"] = np.array(site["V"])/10**3
#     site["W"] = np.array(site["W"])/10**3

#Field Of View Size:         Right Ascension (arcseconds) 0.000204595   Declination (arcseconds) 0.000204595
FieldofViewSize_RightAscesion = 0.000204595
FieldofViewSize_Delination = 0.000204595
# 0.000204595 /3600 * pi /180 rad

Total_flux_density = 1 #Janskys

title1 = "UV coverage of {} sites at {:.1f} GHz\n".format(len(Array),freq/10**9)
title2 = r"Hour angle: ${}^h$ $\sim$ ${}^h$, Declination: ${}^\circ$".format(H[0],H[1],Dec)
atxHz = " at {:.1f} GHz".format(freq/10**9)
Map_center = r" Map center: RA: ${}^h$, Dec ${}^\circ$".format(H[0],Dec)

tele1.gen_uvimg_from_UVW(UVW, freq, RC, (FieldofViewSize_RightAscesion,FieldofViewSize_Delination))


# # temp
# numsites = 100
# sites = []
# dangle = 0.1
# angle_range = 30
# inner = 10
# for i in range(int(numsites/3)):
#     sites.append((random.random()*inner,360*random.random(),angle_range,dangle))
# for i in range(numsites - int(numsites/3)):
#     sites.append((random.random()*(RC[0]/2.1),360*random.random(),angle_range,dangle))
# tele1.fake_uvcover(RC,sites)

# func = lambda x,pos: "{:g}".format(x*1000)
# fmt = FuncFormatter(func)

# plt.colorbar(..., corlorbarformat=fmt)

fmt = r'Flux $(J_y)$'

fig = plt.figure()
plt.subplot(231)
plot_scatter_Baseline(Baselines_uv, title = title1+title2 ,dim=2,legend=False ,xlabel = r"u ($\lambda$)",zlabel = r"w ($\lambda$)",plotlabel = True, ylabel = r"v ($\lambda$)")
angle_per_pixel = FieldofViewSize_RightAscesion/RC[1]

xticks=(np.arange(0, 512, 100), np.arange(-256*angle_per_pixel, 256*angle_per_pixel, 100*angle_per_pixel))
# yticks=(np.arange(0, 512, 100), np.arange(0, 512, 100))
yticks = xticks
tele1.gen_dirty_beam()
tele1.dirty_beam = Total_flux_density*tele1.dirty_beam/np.sum(tele1.dirty_beam)
plt.subplot(232)
plot(tele1.dirty_beam, title = "Dirty beam\n"+atxHz+Map_center,xlabel = "Right Ascension (arcsecond)", origin='lower',corlorbarformat=fmt, cmap = plt.cm.jet, ylabel = "Relative Declination (arcsecond)", colorbar = True, islog = False,xticks=xticks,yticks=yticks)

num_source = 30
point_source = gen_point_source(num_source, RC)

# #####################################################
# # point source conv with a known_kernel
kernelsize = (25,25)
Xin, Yin = np.mgrid[0:kernelsize[0], 0:kernelsize[1]]
# %gaussian(height, center_x, center_y, width_x, width_y)(Xin, Yin)
(height, center_x, center_y, width_x, width_y) = (1, kernelsize[0]//2, kernelsize[1]//2, 5, 5)
known_kernel = gaussian(height, center_x, center_y, width_x, width_y)(Xin, Yin)


point_source = img_conv(point_source, known_kernel)
# #####################################################
noise = 0.001*np.random.normal(size=RC)

point_source += noise
point_source[point_source<0]=0
point_source = Total_flux_density*point_source/np.sum(point_source)

plt.subplot(233)
tbd_title1 = "{} point sources img cov with a gaussian kernel with some noise\n(True brightness distribution)".format(num_source)+atxHz+Map_center
plot(point_source,tbd_title1 ,xlabel = "Right Ascension (arcsecond)", ylabel = "Relative Declination (arcsecond)",corlorbarformat=fmt, origin='lower', colorbar = True, islog = False, cmap = plt.cm.jet,xticks=xticks,yticks=yticks)

lenna_visibility, lenna_visibility_mag = imgfft(point_source)
# plt.subplot(244)
# plot(lenna_visibility_mag, title = "FT of point source img (Theoretical visibility)",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)


measured_visibility, measured_lenna_visibility_mag = sample_visibility(lenna_visibility, lenna_visibility_mag, tele1.uvcover)

# plt.subplot(235)
# plot(measured_lenna_visibility_mag, title = "Measured visibility after sampling in uv plane",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)

lenna_map, lenna_map_mag = imgifft(measured_visibility)
lenna_map_mag = Total_flux_density*lenna_map_mag/np.sum(lenna_map_mag)
plt.subplot(234)
plot(lenna_map_mag, title = "Dirty map\n"+atxHz+Map_center,xlabel = "Right Ascension (arcsecond)", origin='lower',corlorbarformat=fmt, ylabel = "Relative Declination (arcsecond)", colorbar = True, islog = False, cmap = plt.cm.jet,xticks=xticks,yticks=yticks)


# ########################################################################################################

# #(coresize, height, center_x, center_y, width_x, width_y) = model_beam_para
model_beam_para = (512, 1, 256, 256, 5, 5)
cl = Cleaner(lenna_map_mag, tele1.dirty_beam, model_beam_para, Total_flux_density)


# plt.subplot(233)
# # plt.matshow(cl.model_beam, cmap=plt.cm.gist_earth_r)
cl.gen_model_beam_byhand(coresize = 512, height=1, center_x=256, center_y=256, width_x=5, width_y=5)
# # plot(cl.model_beam, title = "Model beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

itertime = 800
loop_gain = 0.15
peak = 9e-6
plt.subplot(235)
cl.clean(itertime = itertime, loop_gain = loop_gain, criteria = "max_itertime")


# cl.clean(loop_gain = loop_gain, criteria = "peak", peak = peak)
# plot(cl.residual, title = "Residual after peak < {} with loop gain {}\n".format(peak, loop_gain)+atxHz+Map_center,corlorbarformat=fmt, origin='lower',xlabel = "Right Ascension (arcsecond)", ylabel = "Relative Declination (arcsecond)", colorbar = True, islog = False, cmap = plt.cm.jet,xticks=xticks,yticks=yticks)

cl.add_residual()
plot(cl.residual, title = "Residual after {} clean steps with loop gain {}\n".format(itertime, loop_gain)+atxHz+Map_center,corlorbarformat=fmt, origin='lower',xlabel = "Right Ascension (arcsecond)", ylabel = "Relative Declination (arcsecond)", colorbar = True, islog = False, cmap = plt.cm.jet,xticks=xticks,yticks=yticks)



plt.subplot(236)
plot(cl.cleandmap, islog = False, title = "Cleaned map\n"+atxHz+Map_center,xlabel = "Right Ascension (arcsecond)",  origin='lower',corlorbarformat=fmt,ylabel = "Relative Declination (arcsecond)", colorbar = True, cmap = plt.cm.jet,xticks=xticks,yticks=yticks)

# plt.tight_layout()
# fig.tight_layout()
# plt.savefig('clean_{}_sources_{}_steps.png'.format(num_source,itertime),dpi=900)
plt.show()