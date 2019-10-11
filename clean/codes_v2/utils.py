import matplotlib.pyplot as plt
from skimage import exposure
import cv2
import numpy as np
import scipy.optimize as opt
import random
from scipy import signal

# 设置一下cmap, https://blog.csdn.net/qq_28485501/article/details/82656614
# plt.imshow(train_images[0], cmap='binary')
# plt.imshow(train_images[0], cmap='Greys_r')
# plt.imshow(train_images[0], cmap='Greys')
def plot(whichimg, title = "title", xlabel="x", ylabel="y", colorbar = False, islog = False, cmap = 'gray'):
    # plt.figure()
    if islog:
        # whichimg = exposure.rescale_intensity(whichimg,in_range=(0,255))
        # print(min(whichimg.all()))
        whichimg = exposure.adjust_log(whichimg)
    plt.imshow(whichimg, cmap = cmap)#cmap = plt.cm.gray_r
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
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