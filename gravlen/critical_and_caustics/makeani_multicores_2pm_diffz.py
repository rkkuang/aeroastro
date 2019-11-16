from two_pointM_lens_inverse_ray_shoot import *
import numpy as np
from utils import genxy
import matplotlib.pyplot as plt
# from matplotlib import animation as ani
from matplotlib.animation import FuncAnimation, writers
import sys
from pathos.pools import ProcessPool
import time

cmap = plt.cm.get_cmap('hot')

testNum = 60
# lens1massRange = np.linspace(1, -0.7 , testNum)
betaRange1 = np.linspace(0, 0.25 , 15, endpoint=False)
betaRange2 = np.linspace(0.25, 0.5 , 15, endpoint=False)
betaRange3 = np.linspace(0.5, 0.75 , 15, endpoint=False)
betaRange4 = np.linspace(0.75, 1 , 15, endpoint=False)
frameinputs = [betaRange1,betaRange2,betaRange3,betaRange4]
ALLbeta = np.concatenate((betaRange1,betaRange2,betaRange3,betaRange4),axis=0)


# massratio = mass1/mass2 change from 2 to -2, 
# mass2 fixed to 1, then only need to change mass1 from 2 to -2
lens1 = Onelens(2, (0,0))
lens2 = Onelens(1, (0.6,0))



massratio = lens1.mass/lens2.mass
testinfo = "x1_{:.1f}x2_{:.1f}massratio_{:.1f}".format(lens1.pos[0],lens2.pos[0],massratio)
mp4name = "./resimgs/2pointmass_diffz/"+testinfo+'.mp4'
htmlname = "./resimgs/2pointmass_diffz/"+testinfo+'.html'
gifname = "./resimgs/2pointmass_diffz/"+testinfo+'.gif'


xlim, ylim =(-0.5,2.5), (-1.5,1.5)
ImgSize = (1026,1026) # raw, colume
thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=3000)
# fig, axs = plt.figure()
fig, axs = plt.subplots(1, 1,figsize=(8,6))

cntIM = 0
def compIMGs(a):
    global cntIM
    
    print("A new core for calculating images is started")
    betaR = frameinputs[a]
    subIMGs = {}
    for beta in betaR:
        sys.stdout.write('done %d/%d\r' % (cntIM, testNum))
        cntIM += 1

        twolens = Twolenses(lens1, lens2, beta=beta)
        twolens.inverse_ray_shooting_diffz(thetax, thetay)
        srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize)
        subIMGs[beta] = srcplaneIMG
    return subIMGs

cnt = 0
def update(beta):
    global cnt
    sys.stdout.write('done %d/%d\r' % (cnt, testNum))
    cnt += 1
    global IMGs
    # title = r"Two point mass lenses, with mass ratio $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}".format(lens1mass, lens1.pos[0], lens2.pos[0])
    title = "Two point mass lenses, with mass ratio:\n"+r" $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}, $\beta=D12*Ds/(D1s*D2)={:.2f}$".format(massratio, lens1.pos[0], lens2.pos[0], beta)
    fig.suptitle(title)
    axs.clear()
    axs.imshow(IMGs[beta], origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])

pool = ProcessPool(nodes=4)
results = pool.imap(lambda a: compIMGs(a), range(4))
resultIMGs = list(results)


IMGs = {}
for img in resultIMGs:
    IMGs.update(img)
print(len(IMGs),"len IMGs")


ani = FuncAnimation(fig, update, frames = ALLbeta,blit = False, interval=1000/2, save_count=300)

# save it into html5 format (need ffmpeg)
print('Begin saving html5\n')
with open(htmlname, 'w') as f:
    # f.write('<!DOCTYPE html> <html> <head> <meta charset="UTF-8"> <title>Test</title> </head> <body> ')
    f.write(ani.to_html5_video())
print('Finished.\n')

# save it into gif (need imagemagick)
# cnt = 0
# print('Begin saving gif')
# ani.save(gifname, writer='imagemagick', fps=5, dpi=240)
# print('Finished.\n')

# save it into mp4 (need ffmpeg)
# fps : number, optional
# Frames per second in the movie. Defaults to None, which will use the animation's specified interval to set the frames per second.
# https://matplotlib.org/api/_as_gen/matplotlib.animation.Animation.html#matplotlib.animation.Animation.save
cnt = 0
print('Begin saving mp4')
FFMpegWriter = writers['ffmpeg']
writer = FFMpegWriter(fps=5, metadata=dict(title='None', artist='None', comment="None"), bitrate=9600)
ani.save(mp4name, writer=writer,dpi=240)
print('Finished.')

plt.show()