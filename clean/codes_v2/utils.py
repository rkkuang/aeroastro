import matplotlib.pyplot as plt
from skimage import exposure
import cv2
import numpy as np
import scipy.optimize as opt
import random
from scipy import signal
from mpl_toolkits.mplot3d import Axes3D  # 空间三维画图
from matplotlib.ticker import ScalarFormatter,FuncFormatter
from matplotlib.ticker import StrMethodFormatter
# plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places
# plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places

# https://www.cnblogs.com/qqhfeng/p/5567539.html
# cbar = plt.colorbar(gci)  
# cbar.set_label('$T_B(K)$',fontdict=font)  
# cbar.set_ticks(np.linspace(160,300,8))  
# cbar.set_ticklabels( ('160', '180', '200', '220', '240',  '260',  '280',  '300'))  

yHeight = 1.05

# 设置一下cmap, https://blog.csdn.net/qq_28485501/article/details/82656614
# plt.imshow(train_images[0], cmap='binary')
# plt.imshow(train_images[0], cmap='Greys_r')
# plt.imshow(train_images[0], cmap='Greys')
def mjrFormatter(x, pos):
    return "{0:.1f}".format(x)
    # return 
# from matplotlib.ticker import ScalarFormatter

class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here

def plot(whichimg, title = "title", xlabel="x", ylabel="y",origin = None, colorbar = False, islog = False,corlorbarformat=None, cmap = 'gray',xticks=(np.arange(0, 512, 100), np.arange(0, 512, 100)), yticks=(np.arange(0, 512, 100), np.arange(0, 512, 100))):
    # plt.figure()
    # plt.yticks(np.arange(0, 1024, 100), np.arange(10000, 11024, 100))
    # #第一个参数表示原来的坐标范围，100是每隔100个点标出一次
    # #第二个参数表示将展示的坐标范围替换为新的范围，同样每隔100个点标出一次
    # plt.xticks(np.arange(0, 2000, 500), np.arange(0, 50000, 500)) 
    # #同理将x轴的表示范围由（0,2000）扩展到（0,50000）每隔500个点标出一次
    

    # https://stackoverflow.com/questions/12750355/python-matplotlib-figure-title-overlaps-axes-label-when-using-twiny
    # Python Matplotlib figure title overlaps axes label when using twiny

    # https://cloud.tencent.com/developer/ask/190276
    # 在Matplotlib中更改Colorbar的缩放/单位


    # figure_1 = plt.figure()
    # ax = figure_1.add_subplot(111)
    if islog:
        # whichimg = exposure.rescale_intensity(whichimg,in_range=(0,255))
        # print(min(whichimg.all()))
        whichimg = exposure.adjust_log(whichimg)
    yticks1,yticks2 = yticks[0],yticks[1]#(np.arange(0, 512, 100), np.arange(0, 512, 100))
    xticks1,xticks2 = xticks[0],xticks[1]
    #extent=extent, origin='lower', (-middle,middle,-middle,middle)
    #https://stackoverflow.com/questions/9382664/python-matplotlib-imshow-custom-tickmarks
    # extent = (-1,1,-1,1)#(xticks[0],xticks[:-1],yticks[0],yticks[:-1])
    extent = (xticks2[0],xticks2[-1],yticks2[0],yticks2[-1])
    # print(extent)
    # input(">>>>>>>>>>>>")
    plt.imshow(whichimg, cmap = cmap,  origin='lower',extent=extent )#cmap = plt.cm.gray_r
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, y=yHeight)
    
    # ax.set_yticks(yticks1,yticks2)
    # ax.set_xticks(xticks1,xticks2)
    


    # plt.ticklabel_format(axis='x', style='sci')
    # plt.ticklabel_format(axis='y', style='sci')
    # #http://www.cocoachina.com/articles/83405
    # #https://blog.csdn.net/junxinwoxin/article/details/86610914
    # ax.yaxis.major.formatter.set_powerlimits((0,1))## 将坐标轴的base number设置为一位。
    # ax.xaxis.major.formatter.set_powerlimits((0,1))
    ax = plt.gca()
    # yfmt = ScalarFormatterForceFormat()
    # yfmt.set_powerlimits((0,0))
    # ax.yaxis.set_major_formatter(yfmt)
    # ax.xaxis.set_major_formatter(yfmt)

    # plt.yticks(yticks1,yticks2)
    # plt.xticks(xticks1,xticks2)

    # plt.yticks(yticks1,yticks2)
    # plt.xticks(xticks1,xticks2)

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    # plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    # ****** https://www.jianshu.com/p/929a3e711f0e
    # 使用python中的matplotlib进行绘图分析数据
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
    ax.yaxis.set_major_formatter(xfmt)
    ax.xaxis.set_major_formatter(xfmt)

    # plt.ticklabel_format(axis='y', style='sci')
    # ax.yaxis.major.formatter.set_powerlimits((0,0))

    # ax.yaxis.set_major_formatter(FuncFormatter(mjrFormatter))
    # ax.xaxis.set_major_formatter(FuncFormatter(mjrFormatter))

    # ax.xaxis.set_major_formatter(ScalarFormatter())
    # ax.yaxis.set_major_formatter(ScalarFormatter())

    #https://blog.csdn.net/wuzlun/article/details/80053277

    if colorbar:
        # plt.colorbar()
        cbar = plt.colorbar()
        cbar.set_label(corlorbarformat)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
def plot_scatter(data, title = "title", xlabel="x",origin = None, ylabel="y",zlabel = 'z',dim = 3, plotlabel = False, colorbar = False, islog = False, cmap = 'gray'):
    # plt.figure()
    if islog:
        # whichimg = exposure.rescale_intensity(whichimg,in_range=(0,255))
        # print(min(whichimg.all()))
        whichimg = exposure.adjust_log(whichimg)
    if not plotlabel:
        ax = plt.gca()
        plt.scatter(data[0],data[1],s=1,c="b")#cmap = plt.cm.gray_r
        plt.axis('square')
        plt.ticklabel_format(axis='x', style='sci')
        plt.ticklabel_format(axis='y', style='sci')
        #http://www.cocoachina.com/articles/83405
        #https://blog.csdn.net/junxinwoxin/article/details/86610914
        ax.yaxis.major.formatter.set_powerlimits((0,1))## 将坐标轴的base number设置为一位。
        ax.xaxis.major.formatter.set_powerlimits((0,1))
    else:
        if dim == 2:
            for site in data:
                ax = plt.gca()
                plt.scatter(site["U"],site["V"],s=1,marker='o',label=site["name"])
                plt.ticklabel_format(axis='x', style='sci')
                plt.ticklabel_format(axis='y', style='sci')
                #http://www.cocoachina.com/articles/83405
                #https://blog.csdn.net/junxinwoxin/article/details/86610914
                ax.yaxis.major.formatter.set_powerlimits((0,1))## 将坐标轴的base number设置为一位。
                ax.xaxis.major.formatter.set_powerlimits((0,1))

                plt.axis('square')
        else:
            fig = plt.figure()
            ax = Axes3D(fig)
            for site in data:
                #https://blog.csdn.net/qq_41149269/article/details/81774026
                ax.scatter(site["U"],site["V"],site["W"],s=2,marker='o',label=site["name"])
                # plt.zlabel(zlabel)
            ax.set_zlabel(zlabel)


        plt.legend()#loc = 'upper right'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if colorbar:
        plt.colorbar()

def plot_scatter_Baseline(data, title = "title",legend=True, xlabel="x", ylabel="y",zlabel = 'z',dim = 3, plotlabel = False, colorbar = False, islog = False, cmap = 'gray'):
    if dim == 2:
        for key,value in data.items():
            ax = plt.gca()
            plt.scatter(value["U"],value["V"],s=1,marker='o',label=key)
            # plt.ticklabel_format(axis='x', style='sci')
            # plt.ticklabel_format(axis='y', style='sci')
            # #http://www.cocoachina.com/articles/83405
            # #https://blog.csdn.net/junxinwoxin/article/details/86610914
            # ax.yaxis.major.formatter.set_powerlimits((0,1))## 将坐标轴的base number设置为一位。
            # ax.xaxis.major.formatter.set_powerlimits((0,1))

            xfmt = ScalarFormatter(useMathText=True)
            xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
            ax.yaxis.set_major_formatter(xfmt)
            ax.xaxis.set_major_formatter(xfmt)

            plt.axis('square')
    else:
        fig = plt.figure()
        ax = Axes3D(fig)
        for site in data:
            #https://blog.csdn.net/qq_41149269/article/details/81774026
            ax.scatter(site["U"],site["V"],site["W"],s=2,marker='o',label=site["name"])
            # plt.zlabel(zlabel)
        ax.set_zlabel(zlabel)

    if legend:
        plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, y=yHeight)
    if colorbar:
        plt.colorbar()


# cv2.IMREAD_COLOR：加载彩色图像。图像的任何透明度都将被忽略。这是默认标志。 flags=1
# cv2.IMREAD_GRAYSCALE：以灰度模式加载图像 flags=0
# cv2.IMREAD_UNCHANGED：加载包含Alpha通道的图像 flags=-1

def readimg(imgpath, flag = 0):
    return cv2.imread(imgpath, flag)

def saveimg(img, savename = "imgname.png"):
    cv2.imwrite(savename, img)

# def imgfft(img):
#     # fft = np.fft.fft2(img)
#     # 将空间域转化为频率域
#     # shift = np.fft.fftshift(fft)# --> 将低频部分移动到图像中心
#     # 由于这里直接生成了 UV 覆盖，低频本来就在图像中心，所以这一步就不做了
#     # return shift

#     dft = cv2.dft(np.float32(img),flags=cv2.DFT_COMPLEX_OUTPUT)
#     #plt.imshow(20*np.log(cv2.magnitude(dft[:,:,0],dft[:,:,1])),cmap='gray')
#     dft_shift = np.fft.fftshift(dft)# --> 将低频部分移动到图像中心
#     #plt.imshow(20*np.log(cv2.magnitude(dft_shift[:,:,0],dft_shift[:,:,1])),cmap='gray')

#     mag = 20*np.log(cv2.magnitude(dft_shift[:,:,0],dft_shift[:,:,1]))
#     # mag = cv2.magnitude(dft_shift[:,:,0],dft_shift[:,:,1])# UV覆盖会扩展为对称形状
#     # real = dft_shift[:,:,1]
#     return dft_shift, mag

# def imgifft(FTimg):
#     ishift = np.fft.ifftshift(FTimg)# 将低频部分从中心移动回到左上角
#     ifft = np.fft.ifft2(ishift) #将频率域转化回空间域
#     dirty_beam = np.real(ifft)
#     dirty_beam = np.fft.ifftshift(dirty_beam)
#     # print(dirty_beam.shape)
#     # print(len(dirty_beam.shape))
#     if len(dirty_beam.shape) == 3:
#         # mag = cv2.magnitude(dirty_beam[:,:,0],dirty_beam[:,:,1])# UV覆盖会扩展为对称形状
#         mag = 20*np.log(cv2.magnitude(dirty_beam[:,:,0],dirty_beam[:,:,1]))
#     else:
#         mag = dirty_beam
#     return dirty_beam, mag

def sample_visibility(lenna_visibility, lenna_visibility_mag, uvcover):
    measured_visibility = np.zeros((uvcover.shape[0],uvcover.shape[1],2))
    measured_visibility[:,:,0] = lenna_visibility[:,:,0] * uvcover
    measured_visibility[:,:,1] = lenna_visibility[:,:,1] * uvcover
    measured_lenna_visibility_mag = lenna_visibility_mag * uvcover
    return measured_visibility, measured_lenna_visibility_mag

def imgfft(img):
    dft = cv2.dft(np.float32(img),flags=cv2.DFT_COMPLEX_OUTPUT)
    dft_shift = np.fft.fftshift(dft)# --> 将低频部分移动到图像中心
    mag = 20*np.log(cv2.magnitude(dft_shift[:,:,0],dft_shift[:,:,1]))
    return dft_shift, mag

def imgifft(FTimg):
    idft_shift = np.fft.ifftshift(FTimg)# 将低频部分从中心移动回到左上角
    idft = cv2.idft(idft_shift)
    # idft = np.fft.ifftshift(idft)# to be deleted
    if len(idft.shape) == 3:
        mag = cv2.magnitude(idft[:,:,0],idft[:,:,1])
    else:
        mag = cv2.magnitude(idft,idft)
    return idft, mag

def gen_point_source(numofpoints, imgsize):
    resimg = np.zeros(imgsize)
    for i in range(numofpoints):
        r = int((0.05+0.95*random.random())*imgsize[0]/1.2)
        c = int((0.05+0.95*random.random())*imgsize[0]/1.2)
        resimg[r,c] = 1
    return resimg
def img_conv(img, kernel,mode = "same",boundary='symm',fillvalue=0):
    return signal.convolve2d(img,kernel,mode=mode,boundary=boundary,fillvalue=fillvalue)
def gen_site(name, East_Longitude, Latitude , X_position, Y_position, Z_position , Lower_Elevation, Upper_Elevation, SEFD, Diameter):
    return {"name":name,"East_Longitude":East_Longitude,"Latitude":Latitude,"X_position":X_position,"Y_position":Y_position,"Z_position":Z_position,"Lower_Elevation":Lower_Elevation,"Upper_Elevation":Upper_Elevation,"SEFD":SEFD,"Diameter":Diameter}

def XYZ2uvw(site,Hdec):
    # Hdec = (H,Dec)
    XYZ = np.array(site)
    H = Hdec[0]*np.pi/180
    Dec = Hdec[1]*np.pi/180
    HdecMatrix = np.array([[np.sin(H),np.cos(H),0],
        [-np.sin(Dec)*np.cos(H),np.sin(Dec)*np.sin(H),np.cos(Dec)],
        [np.cos(Dec)*np.cos(H),-np.cos(Dec)*np.sin(H),np.sin(Dec)]]
        )
    return np.dot(HdecMatrix, XYZ) # or np.dot(HdecMatrix, XYZ.transpose()), the same result



#https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
#python 2 dimensional gaussian fitting






if __name__ == '__main__':
    # lenna = readimg("../imgs/Lenna.png")
    # print(lenna.shape)
    # # cv2.imshow("Lenna", lenna)
    # # cv2.waitKey(0)
    # # cv2.destroyAllWindows()
    # # plt.figure()
    # # plt.imshow(lenna,cmap='gray')
    # # saveimg(lenna, "../imgs/Lennagray.png")
    # plot(lenna, title = "Lenna img",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

    # lenna_visibility, lenna_visibility_mag = imgfft(lenna)
    # plot(lenna_visibility_mag, title = "FT of Lenna img",xlabel = "u (pixel)", ylabel = "v (pixel)", islog = False, colorbar = True)
    
    # ifftlenna, ifftlenna_mag = imgifft(lenna_visibility)
    # # plt.figure()
    # # plt.imshow(ifftlenna_mag,cmap='gray')
    # plot(ifftlenna_mag, title = "IFT of LennaFT",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True)

    # plt.show()
    pass