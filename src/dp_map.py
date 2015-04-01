import numpy as np
from mpl_toolkits.basemap import Basemap

class drake_passage_map(object):

    
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
        fontsize=25.
        ):

            self.lonmin = lonmin
            self.lonmax = lonmax
            self.latmin = latmin
            self.latmax = latmax
            self.lat_1 = lat_1
            self.lon_0 = lon_0
            self.proj = proj
            self.lw = linewidth
            self.fs = fontsize
            
            self.m = Basemap(llcrnrlon=self.lonmin,llcrnrlat=self.latmin,\
                    urcrnrlon=self.lonmax, urcrnrlat=self.latmax,
                    rsphere=(6378137.00,6356752.3142),
                    resolution='i',area_thresh=1000.,projection=self.proj,
                    lat_1=self.lat_1,lon_0=self.lon_0)


    def draw_par_mer():



