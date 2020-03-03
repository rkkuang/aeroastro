import numpy as np
import matplotlib.pyplot as plt
from gaussfitter import fitbeam, gaussian

class Cleaner():
    residual = None
    model_beam = None
    cleandmap = None
    # brightness = None
    # uvcover = None
    def __init__(self, dirty_map, dirty_beam, model_beam_para, Total_flux_density):
        self.dirty_map = dirty_map
        self.dirty_beam = dirty_beam
        self.residual = dirty_map
        self.Total_flux_density=Total_flux_density
        self.total_subed_flux = 0
        self.total_model_flux = 0
        # self.model_beam = np.zeros(self.dirty_beam.shape)
        # model_beam is usually an elliptical Gaussian fitted to the central lobe of the dirty beam
        (coresize, height, center_x, center_y, width_x, width_y) = model_beam_para
        # height = self.find_peak(self.dirty_beam)[0]
        self.gen_model_beam_byhand(coresize, height, center_x, center_y, width_x, width_y)
        # 

        self.dirty_beam = self.normalize_img(self.dirty_beam)
        # self.gen_model_beam_from_dirty_beam()
        self.model_beam = self.normalize_img(self.model_beam)
        self.model = np.zeros((self.dirty_map.shape[0],self.dirty_map.shape[1]))
    def clean(self, itertime = 10, loop_gain = 0.1, criteria = "max_itertime", minrms = 10, peak = 1):
        #self.dirty_map
        # self.residual = self.dirty_map
        if criteria == "max_itertime":
            for i in range(itertime):
                self.clean_once(loop_gain = loop_gain)
        elif criteria == "rms":
            while "XXX" < minrms:
                pass
        elif criteria == "peak":
            while self.find_peak(self.residual)[0] > peak:
                self.clean_once(loop_gain = loop_gain)
        else:
            pass
        # return self.residual, self.model_beam, self.cleand_map
        # do not need return
        # self.residual = self.residual * (self.Total_flux_density)/np.sum(self.residual)
    def clean_once(self, loop_gain = 0.1):
        # uodate self.residual
        # do not need return

        #1. Find the strength and position of the peak (i.e., of the greatest absolute intensity) in the dirty image
        peak_value, position = self.find_peak(self.residual)
        #2. Subtract from the dirty image, at the position of the peak, the dirty beam B multiplied by the peak strength and a damping factor
        tobesub_residual = loop_gain*peak_value*self.gen_tobe_sub(position, self.dirty_beam)

        tempresidual = self.residual - tobesub_residual
        tempresidual[tempresidual<0] = 0
        self.total_subed_flux += np.sum(self.residual-tempresidual)
        self.residual = tempresidual


        tobeadd_model = loop_gain*peak_value*self.gen_tobe_sub(position, self.model_beam)
        self.model += tobeadd_model
        
        self.total_model_flux += np.sum(tobeadd_model)

        # self.residual = self.residual * (self.Total_flux_density - self.total_subed_flux)/np.sum(self.residual)
        # self.model = self.model * (self.total_subed_flux)/np.sum(self.model)

    # clean normalize with total flux density
    # def clean_once(self, loop_gain = 0.1):
    #     # uodate self.residual
    #     # do not need return

    #     #1. Find the strength and position of the peak (i.e., of the greatest absolute intensity) in the dirty image
    #     peak_value, position = self.find_peak(self.residual)
    #     #2. Subtract from the dirty image, at the position of the peak, the dirty beam B multiplied by the peak strength and a damping factor
    #     tobesub_residual = loop_gain*peak_value*self.gen_tobe_sub(position, self.dirty_beam)
    #     # self.residual -= tobesub_residual
    #     # self.residual[self.residual<0] = 0

    #     tempresidual = self.residual - tobesub_residual
    #     # 
    #     tempresidual[tempresidual<0] = 0
    #     self.total_subed_flux += np.sum(self.residual-tempresidual)
    #     self.residual = tempresidual


    #     tobeadd_model = loop_gain*peak_value*self.gen_tobe_sub(position, self.model_beam)
    #     self.model += tobeadd_model
        
    #     self.total_model_flux += np.sum(tobeadd_model)

    #     self.residual = self.residual * (self.Total_flux_density - self.total_subed_flux)/np.sum(self.residual)
    #     self.model = self.model * (self.total_subed_flux)/np.sum(self.model)

        # 3. Go to (1) unless any remaining peak is below some user-specified level. The search for peaks may be constrained to specified areas of the image, called `CLEAN' windows.
        # 4. Convolve the accumulated point source model with an idealized `CLEAN' beam (usually an elliptical Gaussian fitted to the central lobe of the dirty beam).
        # 5. Add the residuals of the dirty image to the `CLEAN' image.
    def add_residual(self):
        # cl.residual = cl.residual * (Total_flux_density - tele1.total_subed_flux)/np.sum(cl.residual)
        self.cleandmap = self.residual + self.model
        self.cleandmap = self.cleandmap * self.Total_flux_density / np.sum(self.cleandmap)
        print("restored total flux: {} Jy".format(np.sum(self.cleandmap)))
    def find_peak(self, beam_or_map, clean_window = [(0,0),(512,512)]):
        #https://thispointer.com/find-max-value-its-index-in-numpy-array-numpy-amax/
        # beam_or_map[0:512,0:512]
        R,C = beam_or_map.shape
        clean_window = [(0,0),(R,C)]
        peak_value = np.amax(beam_or_map[clean_window[0][0]:clean_window[1][0], clean_window[0][1]:clean_window[1][1]])
        res=np.where(beam_or_map==peak_value)
        # return peak_value, list(zip(res[0],res[1]))[0] #(273, 116)
        # print(peak_value)
        # print(list(zip(res[0],res[1])))
        # input(">>>>>>>>>>>>>>>>>>")
        return peak_value, list(zip(res[0],res[1]))[0] #(273, 116)
    def gen_tobe_sub(self, position, dirty_or_model_beam):
        R,C = dirty_or_model_beam.shape
        rpos = position[0]
        cpos = position[1]
        enlarged = np.zeros((3*R, 3*C))
        rcorner = int(R+rpos-R/2)
        ccorner = int(C+cpos-C/2)
        enlarged[rcorner:rcorner+R, ccorner:ccorner+C] = dirty_or_model_beam
        return enlarged[R:2*R,C:2*C]
    def gen_model_beam_from_dirty_beam(self):
        #model_beam is usually an elliptical Gaussian fitted to the central lobe of the dirty beam
        #fitbeam(coresize, data):
        self.model_beam = fitbeam(self.dirty_beam.shape[0], self.dirty_beam)
    def gen_model_beam_byhand(self, coresize = 512, height=10, center_x=256, center_y=256, width_x=5, width_y=10):
        #gaussian(height, center_x, center_y, width_x, width_y):
        Xin, Yin = np.mgrid[0:coresize, 0:coresize]
        # self.model_beam = gaussian(self.find_peak(self.dirty_beam)[0], center_x, center_y, width_x, width_y)(Xin, Yin)
        self.model_beam = gaussian(height, center_x, center_y, width_x, width_y)(Xin, Yin)
    def normalize_img(self, dirty_or_model_beam):
        peak_value, _ = self.find_peak(dirty_or_model_beam)
        return dirty_or_model_beam/peak_value

if __name__ == '__main__':

    # cl = Cleaner(lenna_map_mag, tele1.dirty_beam)
    # cl.clean()

    # cl = Cleaner(lenna_map_mag, tele1.dirty_beam)
    # plt.matshow(cl.model_beam, cmap=plt.cm.gist_earth_r)

    


    # plt.show()

    pass