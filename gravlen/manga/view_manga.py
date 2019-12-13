from marvin.tools.maps import Maps
import marvin
import sys
import matplotlib.pyplot as plt
plateifu_ID = sys.argv[1]
# maps = Maps(plateifu='8086-6101')
# #maps = Maps(plateifu=plateifu_ID)
# print(maps)
# # get an emission line map
# haflux = maps.emline_gflux_ha_6564
# values = haflux.value
# ivar = haflux.ivar
# mask = haflux.mask
# haflux.plot()

# my_cube = marvin.tools.Cube('8086-6101')
# flux = my_cube.flux
# spectrum = flux[:, 1, 2]
# spectrum.plot(show_std=True)


my_cube = marvin.tools.Maps('7443-12703')
ha = my_cube.emline_gflux_ha_6564
fig, ax = ha.plot()


plt.show()
