import numpy as np
import matplotlib.pyplot as plt

class Cleaner():
    residual = None
    model_beam = None
    cleand_map = None
    brightness = None
    uvcover = None
    def __init__(self, dirtymap):
        self.dirtymap = dirtymap
    # def plot(self, whichimg, title, scale = "log",  save = False, colorbar = False):
    #     plt.imshow(whichimg)
    #     plt.xlabel(xlabel)
    #     plt.ylabel(ylabel)
    #     plt.title(title)
    #     if colorbar:
    #         plt.colorbar()



