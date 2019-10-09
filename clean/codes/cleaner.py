import numpy as np
import matplotlib.pyplot as plt

class Cleaner():
    residual = None
    model_beam = None
    cleand_map = None
    # brightness = None
    # uvcover = None
    def __init__(self, dirtymap):
        self.dirtymap = dirtymap
    def clean(self, itertime = 10, gamma = 0.1, criteria = "max_itertime", minrms = 10):
        #self.dirty_map
        self.residual = self.dirtymap
        if criteria == "max_itertime":
            for i in range(itertime):
                #clean_once(self.residual)
                pass
        elif criteria == "rms":
            while "XXX" < minrms:
                pass
        else:
            pass
        # return self.residual, self.model_beam, self.cleand_map
        # do not need return
    def clean_once(self):
        # uodate self.residual
        # do not need return


if __name__ == '__main__':

    cl = Cleaner(lenna_map_mag)
    cl.clean()