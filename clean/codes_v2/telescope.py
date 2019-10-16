import numpy as np
import matplotlib.pyplot as plt
import cv2
# from skimage import exposure
from utils import plot, imgifft, XYZ2uvw, plot_scatter, gen_site, plot_scatter_Baseline

SPEEDOFLIGHT = 3*10**8 # m/s

class Telescope():
    uvcover = None
    dirty_beam = None
    earth_omega = 360/24/3600 #degree per second
    #0.0002 - 512 pixel
    def __init__(self, poslist, directionlist, angle_per_pixel = 1.0/3600, freq = 10**9):
        self.poslist = poslist
        self.directionlist = directionlist
        self.angle_per_pixel = angle_per_pixel # 1.0/3600 --> one arcsec per pixel
        self.freq = freq # 10**9 --> 1 GHz
    def uvcoverage(self, dt, t0, t1):
        pass
    def uvXYZ(self, sites=[(-768713.9637, -5988541.7982, 2063275.9472),(5088967.9, -301681.6, 3825015.8)], center=(-768713.9637, -5988541.7982, 2063275.9472), H=(0, 24, 1/12, 1/6,  6/3600), Del=60, freq=1*10**9):
        # sites: [(X1,Y1,Z1), (Z2,Y2,Z2) in meter
        # center, center site's (X,Y,Z) in meter
        # H: (H0,H1,lH,dH) in hour (15 degrees per hour), e.g. (6.1, 7.1, 1/12, 1/6) 表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每6秒采一次样)， 采完等待 10分钟，之后再采5分钟
        # Dec: Declination, in degree
        # freq: frequency in Hz
        u = []
        v = []
        w = []
        h = H[0]
        wavelength = SPEEDOFLIGHT/freq
        for site in sites:
            # site = (site[i]-center[i] for i in range(3)) # wrong! can not code like this
            #<generator object <genexpr> at 0x7ff25457be60>
            # site = (i/wavelength for i in site)
            site[0],site[1],site[2] = site[0]-center[0],site[1]-center[1],site[2]-center[2]
            site[0],site[1],site[2] = site[0]/wavelength,site[1]/wavelength,site[2]/wavelength
        while h<H[1]:
            for site in sites:
                temph = h
                for i in np.arange(0,H[2],H[4]):
                    uvw = XYZ2uvw(site,(temph*15, Dec))
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                    temph = h + i# temph += i, is wrong
            h += (H[2]+H[3])
        return (u,v,w)
    def uvArray(self, Array=[], center=None, H=(0, 24, 1/12, 1/6,  6/3600), Del=60, freq=1*10**9):
        # Array: [ALMA, PDB], ALMA is a dict contains information of ALMA
        # center, center site dict
        # H: (H0,H1,lH,dH) in hour (15 degrees per hour), e.g. (6.1, 7.1, 1/12, 1/6) 表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每6秒采一次样)， 采完等待 10分钟，之后再采5分钟
        # Dec: Declination, in degree
        # freq: frequency in Hz
        # this is wrong, should use baseline
        u = []
        v = []
        w = []
        h = H[0]
        wavelength = SPEEDOFLIGHT/freq
        for site in Array:
            site["X_position"],site["Y_position"],site["Z_position"] = site["X_position"]-center["X_position"],site["Y_position"]-center["Y_position"],site["Z_position"]-center["Z_position"]
            site["X_position"],site["Y_position"],site["Z_position"] = site["X_position"]/wavelength, site["Y_position"]/wavelength, site["Z_position"]/wavelength
            site["U"] = []
            site["V"] = []
            site["W"] = []
        while h<H[1]:
            for site in Array:
                temph = h
                for i in np.arange(0,H[2],H[4]):
                    uvw = XYZ2uvw((site["X_position"],site["Y_position"],site["Z_position"]),(temph*15, Dec))
                    u.append(uvw[0]), site["U"].append(uvw[0])
                    v.append(uvw[1]), site["V"].append(uvw[1])
                    w.append(uvw[2]), site["W"].append(uvw[2])
                    temph = h + i# temph += i, is wrong
            h += (H[2]+H[3])
        return (u,v,w), Array
    def uvArray_Baseline(self, Array=[], center=None, H=(0, 24, 1/12, 1/6,  6/3600), Dec=60, freq=1*10**9):
        # Array: [ALMA, PDB], ALMA is a dict contains information of ALMA
        # center, center site dict
        # H: (H0,H1,lH,dH) in hour (15 degrees per hour), e.g. (6.1, 7.1, 1/12, 1/6) 表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每6秒采一次样)， 采完等待 10分钟，之后再采5分钟
        # Dec: Declination, in degree
        # freq: frequency in Hz
        u = []
        v = []
        w = []
        
        Baselines_uv = {}
        wavelength = SPEEDOFLIGHT/freq
        ploted = []
        for center in Array:
            h = H[0]
            ploted.append(center)
            for site in Array:
                if site not in ploted:
                    Baselines_uv[center["name"]+"-"+site["name"]] = {}
                    Baselines_uv[center["name"]+"-"+site["name"]]["U"] = []
                    Baselines_uv[center["name"]+"-"+site["name"]]["V"] = []
                    Baselines_uv[center["name"]+"-"+site["name"]]["W"] = []
            while h<H[1]:
                for site in Array:
                    if site not in ploted:
                        temph = h
                        x,y,z = site["X_position"]-center["X_position"],site["Y_position"]-center["Y_position"],site["Z_position"]-center["Z_position"]
                        x,y,z = x/wavelength, y/wavelength, z/wavelength
                        x2,y2,z2 = -x,-y,-z
                        for i in np.arange(0,H[2],H[4]):
                            # site1 - site2
                            uvw = XYZ2uvw((x,y,z),(temph*15, Dec))
                            u.append(uvw[0]), Baselines_uv[center["name"]+"-"+site["name"]]["U"].append(uvw[0])
                            v.append(uvw[1]), Baselines_uv[center["name"]+"-"+site["name"]]["V"].append(uvw[1])
                            w.append(uvw[2]), Baselines_uv[center["name"]+"-"+site["name"]]["W"].append(uvw[2])
                            # site2 - site1
                            uvw = XYZ2uvw((x2,y2,z2),(temph*15, Dec))
                            u.append(uvw[0]), Baselines_uv[center["name"]+"-"+site["name"]]["U"].append(uvw[0])
                            v.append(uvw[1]), Baselines_uv[center["name"]+"-"+site["name"]]["V"].append(uvw[1])
                            w.append(uvw[2]), Baselines_uv[center["name"]+"-"+site["name"]]["W"].append(uvw[2])
                            temph = h + i# temph += i, is wrong
                h += (H[2]+H[3])
        return (u,v,w), Baselines_uv
    def gen_uvimg_from_UVW(self, UVW, freq, RC, FOV):
        uvimg = np.zeros(RC)
        wavelength = SPEEDOFLIGHT/ freq
        U = np.array(UVW[0])/wavelength
        V = np.array(UVW[1])/wavelength
        # scaleu = FOV[0]/3600*np.pi/180
        # scalev = FOV[1]/3600*np.pi/180
        scaleu = np.max((np.max(np.abs(U)),np.max(np.abs(V))))*1.05
        scalev = scaleu
        U = U/scaleu
        V = V/scalev
        r = (RC[0]/2*(1-V))
        c = (RC[1]/2*(1+U))
        flag = 1
        # uvimg[int(r),int(v)]=1
        uvimg[r.astype(np.int),c.astype(np.int)]=1
        # for i in range(len(r)):
        #     try:
        #         uvimg[int(r[i]),int(c[i])] = 1
        #     except:
        #         flag = 0

        # for i in range(len(U)):
        #     u,v = U[i],V[i]
        #     try:
        #         # uvimg[int(u/scaleu*RC[0]+RC[0]/2),int(v/scalev*RC[1]+RC[1]/2)] = 1
        #         # uvimg[int(RC[0]/2 - v/scalev*RC[0]),int(u/scaleu*RC[1]+RC[1]/2)] = 1

        #     except:
        #         flag = 0
        if not flag:
            print("uvcover out of FOV")
        self.uvcover = uvimg
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

    # tele1 = Telescope([],[])
    # RC = (500,500)
    # site1 = (100,0,120,0.2)
    # site2 = (200,240,120,0.2)
    # site3 = (150,60,120,0.2)
    # tele1.fake_uvcover(RC,(site1,site2,site3))
    # plot(tele1.uvcover, title = "Fake uv coverage of 3 sites",xlabel = "u (pixel)",ylabel =  "v (pixel)")
    # tele1.gen_dirty_beam()
    # plot(tele1.dirty_beam, title = "Dirty beam of 3 sites",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

    # tele1.fake_uvcover2(RC, 20)
    # plot(tele1.uvcover, title = "Fake uv coverage corresponding to a Gaussian dirty beam",xlabel = "u (pixel)", ylabel = "v (pixel)",)
    # tele1.gen_dirty_beam()
    # plot(tele1.dirty_beam, title = "Gaussian dirty beam",xlabel = "x (pixel)", ylabel = "y (pixel)", colorbar = True, islog = False)

    # plt.show()

    LMT = gen_site("LMT", "-97:18:53", "18:59:06", -768713.9637, -5988541.7982, 2063275.9472, 15, 85, 560, 50)
    PV = gen_site("PV", "-3:23:33.8", "37:03:58.2", 5088967.9, -301681.6, 3825015.8, 15, 85, 2900, 30)
    ALMA50 = gen_site("ALMA50", "-67:45:11.4", "-23:01:09.4", 2225037.1851, -5441199.162, -2479303.4629, 15, 85, 110, 84.7)
    #SMTO = gen_site("SMTO", "-109:52:19", "32:42:06", -1828796.2, -5054406.8, 3427865.2, 15, 85, 11900, 10)
    SMT = gen_site("SMT", "-109:52:19", "32:42:06", -1828796.2, -5054406.8, 3427865.2, 15, 85, 11000, 10)
    Hawaii8 = gen_site("Hawaii8", "-155:28:40.7", "19:49:27.4", -5464523.4, -2493147.08, 2150611.75, 15, 85, 4900, 20.8)
    #PdBI = gen_site("PdBI", "05:54:28.5", "44:38:02.0", 4523998.4, 468045.24, 4460309.76, 15, 85, 1600, 36,7)
    PDB = gen_site("PDB", "05:54:28.5", "44:38:02.0", 4523998.4, 468045.24, 4460309.76, 15, 85, 5200, 36.7)
    SPT = gen_site("SPT", "-000:00:00.0", "-90:00:00", 0, 0, -6359587.3, 15, 85, 7300, 12)
    GLT = gen_site("GLT", "72:35:46.4", "38:25:19.1", 1500692, -1191735, 6066409, 15, 85, 4744, 12)
    CARMA8 = gen_site("CARMA8", "-118:08:30.3", "37:16:49.6", -2397431.3, -4482018.9, 3843524.5, 15, 85, 3500, 26.9)
    SMA = gen_site("SMA", "-155:28:40.7", "19:49:27.4", -5464523.4, -2493147.08, 2150611.75, 15, 85, 4000, 20.8)
    Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT,CARMA8,SMA,PDB]


    # 用uvXYZ 生成 UV 覆盖, 并用plot_scatter 画图的测试程序
    # uvXYZ(self, sites, center, H, delta, lambda)
    # tele1 = Telescope([],[])
    # sites = []
    # for site in Array:
    #     XYZ = [site["X_position"],site["Y_position"],site["Z_position"]]
    #     sites.append(XYZ)
    # center = sites[0]
    # H = (0, 2, 1/12, 1/6, 5/3600)#表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每5秒采一次样)， 采完等待 10分钟，之后再采5分钟
    # Dec = 60 #Declination, 赤纬，度
    # freq = 227.297 * 10**9 # 1 GHz, 227297 MHz
    # UVW = tele1.uvXYZ(sites, center, H, Dec, freq)
    # x = np.array(UVW[0])/10**3
    # y = np.array(UVW[1])/10**3
    # plot_scatter((x,y), title = "UV coverage",xlabel = r"u ($k\lambda$)", ylabel = r"v ($k\lambda$)")
    # plt.show()

    # 用uvArray 生成 UV 覆盖, 并用plot_scatter 画图(带legend)的测试程序
    # uvArray(self, Array, center, H, delta, lambda)

    # The frequency range available to ALMA is divided into different receiver bands. 
    # Data can only be taken in one band at a time. These bands range from band 3, 
    # starting at 84 GHz, to band 10, ending at ~950 GHz. For comparison, 
    # a frequency of 300 GHz translates to a wavelength of approximately 1mm.

    # ============================================================
    # tele1 = Telescope([],[])
    # center = ALMA50#Array[0]
    # H = (0, 6, 1/3, 1/2, 1/3600)#表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每5秒采一次样)， 采完等待 10分钟，之后再采5分钟
    # Dec = 60 #Declination, 赤纬，度
    # freq = 227.297 * 10**9 # 1 GHz, 227297 MHz
    # UVW,resuvArray = tele1.uvArray(Array, center, H, Dec, freq)
    # for site in resuvArray:
    #     site["U"] = np.array(site["U"])/10**3
    #     site["V"] = np.array(site["V"])/10**3
    #     site["W"] = np.array(site["W"])/10**3

    # title1 = "UV coverage with center site: {} at {:.1f} GHz\n".format(center["name"],freq/10**9)
    # title2 = r"Hour angle: ${}^h$ $\sim$ ${}^h$, Declination: ${}^\circ$".format(H[0],H[1],Dec)
    # plot_scatter(resuvArray, title = title1+title2 ,dim=2 ,xlabel = r"u ($k\lambda$)",zlabel = r"w ($k\lambda$)",plotlabel = True, ylabel = r"v ($k\lambda$)")
    # # x = np.array(UVW[0])/10**3
    # # y = np.array(UVW[1])/10**3
    # # plot_scatter((x,y), title = "UV coverage",xlabel = r"u ($k\lambda$)", ylabel = r"v ($k\lambda$)")
    # plt.show()

    # 用uvArray_Baseline 生成 UV 覆盖, 并用plot_scatter 画图(带legend)的测试程序
    # Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT,CARMA8,SMA,PDB]
    Array = [LMT,PV,ALMA50,SMT,Hawaii8,GLT]
    tele1 = Telescope([],[])
    center = ALMA50#Array[0]
    H = (0, 6, 1/3, 1/2, 1/3600)#表示时角范围从 6.1 到 7.1, 每次采样 5分钟(这5分钟里每5秒采一次样)， 采完等待 10分钟，之后再采5分钟
    Dec = 60 #Declination, 赤纬，度
    freq = 227.297 * 10**9 # 1 GHz, 227297 MHz
    UVW, Baselines_uv = tele1.uvArray_Baseline(Array, [], H, Dec, freq)
    # for site in resuvArray:
    #     site["U"] = np.array(site["U"])/10**3
    #     site["V"] = np.array(site["V"])/10**3
    #     site["W"] = np.array(site["W"])/10**3
    
    tele1.gen_uvimg_from_UVW(UVW, freq, (512,512), (0.000204595,0.000204595))

    title1 = "UV coverage of {} sites at {:.1f} GHz\n".format(len(Array),freq/10**9)
    title2 = r"Hour angle: ${}^h$ $\sim$ ${}^h$, Declination: ${}^\circ$".format(H[0],H[1],Dec)
    plot_scatter_Baseline(Baselines_uv, title = title1+title2 ,dim=2 ,xlabel = r"u ($\lambda$)",zlabel = r"w ($\lambda$)",plotlabel = True, ylabel = r"v ($\lambda$)")
    # x = np.array(UVW[0])/10**3
    # y = np.array(UVW[1])/10**3
    # plot_scatter((x,y), title = "UV coverage",xlabel = r"u ($k\lambda$)", ylabel = r"v ($k\lambda$)")
    # print(np.max(tele1.uvcover))
    plt.figure()

    plot(tele1.uvcover, title = "UV coverage image of {} sites".format(len(Array)),xlabel = "u (pixel)", ylabel = "v (pixel)", colorbar = True)
    plt.show()