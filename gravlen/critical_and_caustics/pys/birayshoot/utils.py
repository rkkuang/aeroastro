import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
# 设置横纵坐标的名称以及对应字体格式, https://blog.csdn.net/A_Z666666/article/details/81165123
font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 16,
         }

def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)


def genxy(xlim=(-1,1),ylim=(-1,1),num=100,datatype = np.float64):
    #.astype(np.float32) -- 3.6G -> 8G, Maximum Magnification:  266846.88
    #.astype(np.float64) -- 3.6G -> 9G, Maximum Magnification:  328266.3319540504
    x = np.linspace(xlim[0],xlim[1],num).astype(datatype)#np.single, np.double
    y = np.linspace(ylim[0],ylim[1],num).astype(datatype)
    X,Y = np.meshgrid(x,y)
    # return X.reshape(1,-1), Y.reshape(1,-1) # shape: (1,10)
    return X.reshape(1,-1)[0], Y.reshape(1,-1)[0] # shape: (3,)

def spescatter(x,y,xlabel="x",ylabel="y",title="title",c="b",s=1,issqure=True,xylim=None):
    #plot
    plt.scatter(x,y,s=s,c=c)
    # plt.scatter([-self.X,self.X],[0,0],s=50,c="r",marker="+")
    plt.xlabel(xlabel,font2)
    plt.ylabel(ylabel,font2)
    plt.tick_params(labelsize=12)
    plt.grid()
    plt.title(title)
    if issqure:
        plt.axis('square')
    if xylim:
        plt.ylim(-xylim,xylim)
        plt.xlim(-xylim,xylim)

def genD1D2(D1, num):
    # generate a list of (D1, D2) while D2 > D1
    D2nplist = np.linspace(D1, 1, num)
    res = []
    for d2 in D2nplist:
        res.append((D1, d2))
    return res