from _3DHR_Utilities import *

dx = 3.5 #micras
dy = 3.5

lambda_ = 0.633 #micras

region = 3

path = "C:/Users/dacv0/Desktop/holograms/"

grab = []

for i in range(1, 1165):
    try:
        grab.append(read(path + f'pacienteSCD_definitivo_hipoxia{i}.bmp'))
    except:
        continue

for i in range(len(grab)):
    holo = compensate(grab[i], dx, dy, lambda_, region, step=0.5, depth=3)

    save(holo, f'pacienteSCD_definitivo_hipoxia_RE_{i}', path=path+'resultados/', ext='bmp', out_amp = False)