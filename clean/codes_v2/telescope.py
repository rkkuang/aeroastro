import numpy as np
import matplotlib.pyplot as plt
import cv2
# from skimage import exposure
from utils import plot, imgifft

class Telescope():
    uvcover = None
    dirty_beam = None
    def __init__(self, poslist, directionlist):
        self.poslist = poslist
        self.directionlist = directionlist
    def uvcoverage(self, dt, t0, t1, freq):
        pass
    def fake_uvcover(self, RC, teles):
        # RC is the size of generate uv coverage image, teles are sites in different radiu and angle_position
        Row = RC[0]
        Col = RC[1]
        self.uvcover = np.zeros((Row,Col))#in pixel 
        for tele in teles:
            radius = tele[0] 
            start_angle = tele[1]/180*np.pi
            end_angle = (tele[1]+tele[2])/180*np.pi
            dtheta = tele[3]/180*np.pi
            for i in np.arange(start_angle,end_angle,dtheta):
                self.uvcover[int(Row/2-radius*np.sin(i)),int(Col/2+radius*np.cos(i))]=1
        # return self.uvcover
    def fake_uvcover2(self, RC, radius):
        # RC is the size of generate uv coverage image, teles are sites in different radiu and angle_position
        Row = RC[0]
        Col = RC[1]
        self.uvcover = np.zeros((Row,Col))#in pixel 
        for i in range(Row):
            for j in range(Col):
                if (Row/2-i)**2+(Col/2-j)**2<radius**2:
                    self.uvcover[i,j]=1
        # return self.uvcover

    # def plot(self, whichimg, title = "title", xlabel="x", ylabel="y", colorbar = False, islog = False):
    #     plt.figure()
    #     if islog:
    #         whichimg = exposure.rescale_intensity(whichimg,in_range=(0,100))
    #         # print(min(whichimg.all()))
    #         whichimg = exposure.adjust_log(whichimg)
    #     plt.imshow(whichimg, cmap = plt.cm.gray_r)#cmap = plt.cm.gray_r
    #     plt.xlabel(xlabel)
    #     plt.ylabel(ylabel)
    #     plt.title(title)
    #     if colorbar:
    #         plt.colorbar()

    def gen_dirty_beam(self):
        # generate dirty beam from uv coverage
        #https://blog.csdn.net/giffordy/article/details/92838671
        # f_ishift = np.fft.ifftshift(self.uvcover)
        # self.dirty_beam = cv2.idft(f_ishift)

        # self.dirty_beam = cv2.idft(self.uvcover)
        # self.dirty_beam = np.fft.ifftshift(self.dirty_beam)

        #print(img_back.shape)
        #input(">>>>>>>>")
        #img_back = cv2.magnitude(img_back[:,:,0],img_back[:,:,1])

        # https://blog.csdn.net/a13602955218/article/details/84448075

        # fft = np.fft.fft2(img)
        # 将空间域转化为频率域

        #shift = np.fft.fftshift(fft) --> 将低频部分移动到图像中心
        # 由于这里直接生成了 UV 覆盖，低频本来就在图像中心，所以这一步就不做了

        # worked code: now moved to utils
        # ishift = np.fft.ifftshift(self.uvcover)# 将低频部分从中心移动回到左上角
        # ifft = np.fft.ifft2(ishift) #将频率域转化回空间域
        # self.dirty_beam = np.real(ifft)
        # self.dirty_beam = np.fft.ifftshift(self.dirty_beam)
        # # self.dirty_beam = cv2.normalize(self.dirty_beam)
        # # self.dirty_beam = 20*np.log(np.abs(self.dirty_beam))

        # # idft_shift = np.fft.ifftshift(self.uvcover)
        # # idft = cv2.idft(idft_shift)
        # # self.dirty_beam = cv2.magnitude(idft[:,:,0],idft[:,:,1])

        _, self.dirty_beam = imgifft(self.uvcover)
        self.dirty_beam = np.fft.ifftshift(self.dirty_beam)
        # return self.dirty_beam


if __name__ == "__main__":
    # poslist = [(45,120),(30,10)]
    # directionlist = []
    # dt = 100 # seconds
    # t0 = 0
    # t1 = 10**4
    # freq = 100*10**6# 100 MHz
    # tele1 = Telescope(poslist,directionlist)
    # uvcover = tele1.uvcoverage(dt, t0, t1, freq)

    tele1 = Telescope([],[])
    RC = (500,500)
    site1 = (100,0,120,0.2)
    site2 = (200,240,120,0.2)
    site3 = (150,60,120,0.2)
    tele1.fake_uvcover(RC,(site1,site2,site3))
    plot(tele1.uvcover, title = "Fake uv coverage of 3 sites",xlabel = "u (pixel)",ylabel =  "v (pixel)")
    tele1.gen_dirty_beam()
    plot(tele1.dirty_beam, title = "Dirty beam of 3 sites",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

    tele1.fake_uvcover2(RC, 20)
    plot(tele1.uvcover, title = "Fake uv coverage corresponding to a Gaussian dirty beam",xlabel = "u (pixel)", ylabel = "v (pixel)",)
    tele1.gen_dirty_beam()
    plot(tele1.dirty_beam, title = "Gaussian dirty beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

    plt.show()