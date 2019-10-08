import numpy as np
import matplotlib.pyplot as plt

class Telescope():
    uvcover = None
    def __init__(self, poslist, directionlist):
        self.poslist = poslist
        self.directionlist = directionlist
    def uvcoverage(self, dt, t0, t1, freq):
        pass
    def fake_uvcover(self, RC, teles):
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
        return self.uvcover
    def uvplot(self,title):
        plt.imshow(self.uvcover, cmap = plt.cm.gray_r)
        plt.xlabel("u (pixel)")
        plt.ylabel("v (pixel)")
        plt.title(title)
        plt.show()


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
    tele1.uvplot("Fake uv coverage of 3 sites")