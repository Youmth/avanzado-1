from _3DHR_Utilities import *

path = "C:/Users/dacv0/Desktop/holograms/resultados/"

grab = []

for i in range(1, 1165):
    try:
        grab.append(read(path + f'pacienteSCD_definitivo_hipoxia_RE_{i}.bmp'))
    except:
        continue

for i in range(len(grab)):
    
    plt.imsave(path+ 'Unwrapped/' + f'pacienteSCD_definitivo_hipoxia_RE_{i}.bmp', grab[i], cmap='gray')