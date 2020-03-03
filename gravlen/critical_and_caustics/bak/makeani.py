from two_pointM_lens_inverse_ray_shoot import *
import numpy as np
from utils import genxy
import matplotlib.pyplot as plt
# from matplotlib import animation as ani
from matplotlib.animation import FuncAnimation, writers
import sys

#ref: https://littleround.cn/2019/01/04/Python%E5%88%B6%E4%BD%9C%E5%8A%A8%E6%80%81%E5%9B%BE-matplotlib.animation/

cmap = plt.cm.get_cmap('hot')

testNum = 50
lens1massRange = np.linspace(1, -0.7 , testNum)

lens2 = Onelens(1, (0.3,0))
lens1 = Onelens(1, (-0.3,0))

testinfo = "x1_{:.1f}x2_{:.1f}".format(lens1.pos[0],lens2.pos[0])
mp4name = testinfo+'.mp4'
htmlname = testinfo+'.html'


xlim, ylim =(-2,2), (-1.5,1.5)
ImgSize = (770,1026) # raw, colume
thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=3000)
# fig, axs = plt.figure()
fig, axs = plt.subplots(1, 1,figsize=(8,6))

cnt = 0
def update(lens1mass):
    global cnt
    sys.stdout.write('done %d/%d\r' % (cnt, testNum))
    cnt += 1

    lens1.mass = lens1mass
    twolens = Twolenses(lens1, lens2)
    twolens.inverse_ray_shooting(thetax, thetay)
    srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize)
    
    title = r"Two point mass lenses, with mass ratio $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}".format(twolens.massratio, twolens.lens1.pos[0], twolens.lens2.pos[0])
    fig.suptitle(title)
    
    axs.clear()
    axs.imshow(srcplaneIMG, origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])


ani = FuncAnimation(fig, update, frames = lens1massRange,blit = False, interval=1000, save_count=300)
#, interval=1000.0/12, save_count=100

# https://matplotlib.org/api/animation_api.html
# fig 进行动画绘制的figure
# func 自定义动画函数，即传入刚定义的函数animate
# frames 动画长度，一次循环包含的帧数
# init_func 自定义开始帧，即传入刚定义的函数init
# interval 更新频率，以ms计
# blit 选择更新所有点，还是仅更新产生变化的点。应选择True，但mac用户请选择False，否则无法显示动画

# save it into html5 format (need ffmpeg)
print('Begin saving html5')
with open(htmlname, 'w') as f:
    # f.write('<!DOCTYPE html> <html> <head> <meta charset="UTF-8"> <title>Test</title> </head> <body> ')
    f.write(ani.to_html5_video())
    # f.write('</body> </html>')
print('Finished.')

# # save it into gif (need imagemagick)
# print('Begin saving gif')
# ani.save('test.gif', writer='imagemagick', fps=5, dpi=240)
# print('Finished.')

# save it into mp4 (need ffmpeg)
# fps : number, optional
# Frames per second in the movie. Defaults to None, which will use the animation's specified interval to set the frames per second.
# https://matplotlib.org/api/_as_gen/matplotlib.animation.Animation.html#matplotlib.animation.Animation.save
print('Begin saving mp4')
FFMpegWriter = writers['ffmpeg']
writer = FFMpegWriter(fps=5, metadata=dict(title='None', artist='None', comment="None"), bitrate=9600)
ani.save(mp4name, writer=writer,dpi=240)
print('Finished.')



plt.show()