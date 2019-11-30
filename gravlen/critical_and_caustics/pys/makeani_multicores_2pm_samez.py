# not working multiprocess.pool.MaybeEncodingError: 
# Error sending result: 
# '[<matplotlib.animation.FuncAnimation object at 0x7f58be20c2d0>]'. 
# Reason: 'TypeError("can't pickle _tkinter.tkapp objects")'
# so compute all IMGs first then make animation


from two_pointM_lens_inverse_ray_shoot import *
import numpy as np
from utils import genxy
import matplotlib.pyplot as plt
# from matplotlib import animation as ani
from matplotlib.animation import FuncAnimation, writers
import sys
from pathos.pools import ProcessPool
import time

# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
# from multiprocessing.pool import Pool
# import pathos.multiprocessing as mp

# #ref: https://littleround.cn/2019/01/04/Python%E5%88%B6%E4%BD%9C%E5%8A%A8%E6%80%81%E5%9B%BE-matplotlib.animation/
# import dill
# try:
#     import copy_reg
# except:
#     import copyreg as copy_reg
# import types

# def _pickle_method(method):
#     attached_object = method.im_self or method.im_class
#     func_name = method.im_func.func_name

#     if func_name.startswith('__'):
#         func_name = filter(lambda method_name: method_name.startswith('_') and method_name.endswith(func_name), dir(attached_object))[0]

#     return (getattr, (attached_object, func_name))

# copy_reg.pickle(types.MethodType, _pickle_method)

cmap = plt.cm.get_cmap('hot')

testNum = 60
# lens1massRange = np.linspace(1, -0.7 , testNum)
lens1massRange1 = np.linspace(1, 0.2 , 10, endpoint=False)
lens1massRange2 = np.linspace(0.2, 0 , 20, endpoint=False)
lens1massRange3 = np.linspace(0, -0.2 , 20, endpoint=False)
lens1massRange4 = np.linspace(-0.2, -0.75 , 10, endpoint=False)
frameinputs = [lens1massRange1,lens1massRange2,lens1massRange3,lens1massRange4]
ALLmass = np.concatenate((lens1massRange1,lens1massRange2,lens1massRange3,lens1massRange4),axis=0)


lens2 = Onelens(1, (1,0))
lens1 = Onelens(1, (-1,0))

testinfo = "x1_{:.1f}x2_{:.1f}".format(lens1.pos[0],lens2.pos[0])
mp4name = "./resimgs/2pointmass_samez/"+testinfo+'.mp4'
htmlname = "./resimgs/2pointmass_samez/"+testinfo+'.html'
gifname = "./resimgs/2pointmass_samez/"+testinfo+'.gif'


xlim, ylim =(-3.5,2), (-3,3)
ImgSize = (1026,1026) # raw, colume
thetax, thetay = genxy(xlim=xlim,ylim=ylim,num=30)
# fig, axs = plt.figure()
fig, axs = plt.subplots(1, 1,figsize=(8,6))

cntIM = 0
def compIMGs(a):
    global cntIM
    sys.stdout.write('done %d/%d\r' % (cntIM, testNum))
    print("A new core for calculating images is started")
    lens1massR = frameinputs[a]
    subIMGs = {}
    for lens1mass in lens1massR:
        lens1.mass = lens1mass
        twolens = Twolenses(lens1, lens2)
        twolens.inverse_ray_shooting(thetax, thetay)
        srcplaneIMG = twolens.img_mapping(twolens.betax, twolens.betay, xlim, ylim, ImgSize)
        # IMGs.append(srcplaneIMG)
        subIMGs[lens1mass] = srcplaneIMG
    return subIMGs

cnt = 0
def update(lens1mass):
    global cnt
    sys.stdout.write('done %d/%d\r' % (cnt, testNum))
    cnt += 1
    global IMGs
    title = r"Two point mass lenses, with mass ratio $\mu_1/\mu_2=${:.2f} and $x_1=${:.1f}, $x_2 =${:.1f}".format(lens1mass, lens1.pos[0], lens2.pos[0])
    fig.suptitle(title)
    axs.clear()
    axs.imshow(IMGs[lens1mass], origin='lower',cmap=cmap, extent=[xlim[0],xlim[1],ylim[0],ylim[1]])

# def partani(a):
#     frameinput = frameinputs[a]
#     return FuncAnimation(fig, update, frames = frameinput,blit = False, interval=1000, save_count=300)

# ani = FuncAnimation(fig, update, frames = lens1massRange,blit = False, interval=1000, save_count=300)

# anis = []
# for frameinput in frameinputs:
#     future = executor.submit(partani,frameinput)
#     res = future.done()
#     anis.append(res)

# #多进程
# #AttributeError: Can't pickle local object 'FuncAnimation.__init__.<locals>.<lambda>'
# #https://stackoverflow.com/questions/8804830/python-multiprocessing-picklingerror-cant-pickle-type-function
# pool = ProcessPoolExecutor(max_workers=4)
# resultanis = list(pool.map(partani,frameinputs))
# print('主进程 done')

#https://xbuba.com/questions/53049044
#ThreadPoolExecutor可以工作，但ProcessPoolExecutor没有

# pool = Pool(processes=4)#,  maxtasksperchild=100
# resultanis = []
# for x in pool.imap_unordered(partani,frameinputs):
#     resultanis.append(x)


# with ProcessPoolExecutor(max_workers = 4) as executor:
#     print("Exec")
#     resultanis = executor.map(lambda a: partani(a), frameinputs)

#https://matthewrocklin.com/blog/work/2013/12/05/Parallelism-and-Serialization
#https://pathos.readthedocs.io/en/latest/pathos.html
# p = mp.Pool(4)
# resultIMGs = p.map(lambda a: compIMGs(a), range(4))
# resultIMGs = resultIMGs.get()
# # p.close()
# # p.join()

pool = ProcessPool(nodes=4)
# results = pool.amap(lambda a: compIMGs(a), range(4))
# while not results.ready():
#     time.sleep(5); print(".", end=' ')
# resultIMGs = results.get()
results = pool.imap(lambda a: compIMGs(a), range(4))
resultIMGs = list(results)


# print(type(resultIMGs))
# print(len(resultIMGs[0]))
# input()
IMGs = {}
for img in resultIMGs:
    IMGs.update(img)
print(len(IMGs),"len IMGs")


ani = FuncAnimation(fig, update, frames = ALLmass,blit = False, interval=1000/2, save_count=300)

# input("continue")


# # 多线程
# pool = ThreadPoolExecutor(max_workers=4)
# resultanis = list(pool.map(partani,frameinputs))
# print('主线程结束')



#, interval=1000.0/12, save_count=100

# https://matplotlib.org/api/animation_api.html
# fig 进行动画绘制的figure
# func 自定义动画函数，即传入刚定义的函数animate
# frames 动画长度，一次循环包含的帧数
# init_func 自定义开始帧，即传入刚定义的函数init
# interval 更新频率，以ms计
# blit 选择更新所有点，还是仅更新产生变化的点。应选择True，但mac用户请选择False，否则无法显示动画

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