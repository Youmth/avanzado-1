from _3DHR_Utilities import *

'''Bulk reconstruction.

This code takes all the images withing the path that have a sequential name (hard coded) 
and reconstructs all of them one by one.'''


dx = 3.5 #micras
dy = 3.5

lambda_ = 0.633 #micras

region = 3

focus_points = 5

path = "C:/Users/dacv0/Desktop/holograms/"

grab = []

for i in range(1, 1165):
    try:
        grab.append(read(path + f'pacienteSCD_definitivo_hipoxia{i}.bmp'))
    except:
        continue

for i in range(len(grab)):
    holo = compensate(grab[i], dx, dy, lambda_, region, step=0.5, depth=3)
    holo, _, _ = prop_focus(holo, lambda_, dx, dy, -200, 200, focus_points, 40)

    save(holo, f'pacienteSCD_definitivo_hipoxia_RE_{i}', path=path+'resultados2/', ext='bmp', out_amp = False)