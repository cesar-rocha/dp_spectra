import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

class drake_passage_map(object):
    """ a class that represents the Drake Passage Basemap object """
    
    def __init__(
        self,
        lonmin=-75,
        lonmax=-50,
        latmin=-67,
        latmax=-50,
        lat_1=-55,
        lon_0=-75,
        proj="lcc",
        linewidth=.5,
        fontsize=18.,
        fig_label="c",
        fig_title="AlTiKa"):

            self.lonmin = lonmin
            self.lonmax = lonmax
            self.latmin = latmin
            self.latmax = latmax
            self.lat_1 = lat_1
            self.lon_0 = lon_0
            self.proj = proj
            self.lw = linewidth
            self.fs = fontsize
            self.label= fig_label
            self.title = fig_title

            self.m = Basemap(llcrnrlon=self.lonmin,llcrnrlat=self.latmin,\
                    urcrnrlon=self.lonmax, urcrnrlat=self.latmax,
                    rsphere=(6378137.00,6356752.3142),
                    resolution='i',area_thresh=1000.,projection=self.proj,
                    lat_1=self.lat_1,lon_0=self.lon_0)


    def draw_par_mer(self):
        """ draw parallel and meridians """
        self.m.drawparallels(np.arange(self.latmin, self.latmax, 4),
            labels=[1, 0, 0, 0], linewidth=self.lw, fontsize=self.fs)
        self.m.drawmeridians(np.arange(self.lonmin, self.lonmax+4, 7),
            labels=[0, 0, 0, 1], linewidth=self.lw, fontsize=self.fs)

    def set_label(self,pos=(1650212,1485371)):
        plt.text(pos[0],pos[1], self.label,size=32)

    def set_title(self,pos=(1400212,1495371)):
        plt.text(pos[0],pos[1], self.title, size=25, rotation=0.,
            ha="center", va="center",bbox = dict(boxstyle="round",ec='k',fc='w'))



