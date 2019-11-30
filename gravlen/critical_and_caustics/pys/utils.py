import numpy as np
import matplotlib.pyplot as plt
# 设置横纵坐标的名称以及对应字体格式, https://blog.csdn.net/A_Z666666/article/details/81165123
font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 16,
         }

def genxy(xlim=(-1,1),ylim=(-1,1),num=100):
    x = np.linspace(xlim[0],xlim[1],num)
    y = np.linspace(ylim[0],ylim[1],num)
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