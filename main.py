import Utilities as ut
import simulationOffaxis as sof
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

from _3DHR_Utilities import *

def propagate(U: np.array,
              z: float,
              lmbda: float,
              dx:float,
              dy:float,
              scale_factor = 1) -> np.ndarray:
  ''' Function to propagate a scalar wave field with its angular spectrum
  Inputs:
    U: matrix containing the input image array
    z: propagated distance
    lmbda: wavelength
    dx,dy: pixel pitch of the digital sensor in each direction
    scale_factor: magnification of the optical system

  Outputs:
    Uz_complex: Propagated complex field
  '''

  # Getting the shape of the input field
  N, M = np.shape(U)[:2]

  #Creating a meshgrid with the correct spatial coordinates
  mm, nn = np.linspace(-M/2, M/2-1, M), np.linspace(-N/2, N/2-1, N)

  m, n = np.meshgrid(mm, nn)

  #angular spectrum of the image at the origin
  A = np.fft.fft2(U)

  kernel = np.exp(-1j * 2 * np.pi * z * scale_factor *
                  np.sqrt((1/(lmbda**2)) - (m/(M*dx))**2 - (n/(N*dy))**2 + 0j))


  #optical field at z
  Uz_complex = np.fft.ifft2(A*np.fft.fftshift(kernel))

  return Uz_complex

def normalize(x: np.ndarray, scale: float) -> np.ndarray:
    '''Normalize every value of an array to the 0-scale interval.'''

    return scale*(x-np.min(x))/(np.max(x)-np.min(x))

def sphere_phase_shift(input_field, radius, pos, n, lambda_, dxy, scale_factor, n0=1):
    '''Produces a spherical phase shift in the input field
   
    The input field is affected by a pure-phase object in the shape of a sphere,
    locally inducing a phase shift without affecting its amplitude

    radius: real radius of the sphere in length units
    pos: real position of the sphere in the image in length units relative to the
    center of the image
    n: Refractive index of the sphere
    z: real distance from the sphere to the image in length units
    lambda_: wavelength of light in length units
    dx, dy: pixel pitches in length units
    scale_factor: scale factor induced by microscope magnification, affects z
    n0: refractive index of the liquid medium of the sample, set to air by default
    '''
    N, M = input_field.shape

    angle = np.zeros((M, N))

    radius = radius*scale_factor/dxy

    for x in range(N):
      for y in range(M):
      
        cx, cy = pos[0]*dxy+N//2, pos[1]*dxy+M//2

        if (x-cx)**2 + (y-cy)**2 >= radius**2:
          angle[y, x] = 0
        else:
          # The phase shift corresponding to a change in index of refraction is equal to 2pi(OPL)/lambda where
          # OPL=optical path length, and is equal to (n2-n1)*d where d is the distance, in this case
          # d corresponds to the width of a sphere in the direction of z. Naturally, the scale factor introduces
          # a difference between the real size of the object and its apparent size in the camera, changing the 
          # aparent optical path difference, which is corrected by dividing by the scale factor

          angle[y, x] = 2*np.pi*(n-n0)/lambda_ * (radius - np.sqrt((x-cx)**2 + (y-cy)**2))*dxy/scale_factor

    complex = np.abs(input_field)*np.exp(1j*(np.angle(input_field)+angle))

    return complex

def sphere_sample(input_field, radii, xys, zs, ns, lambda_, dxy, scale_factor, n0 = 1, final_z=None):
  '''Simulates a number of pure-phase spherical samples.

  The spheres can all have different 3D positions, radii and refraction indices

  Position format is (x, y, z), input field is assumed to be in z=0

  The final image of the sample will be propagated to a z position equal to final_prop, this position
  is equal to the position of the last sphere on the array by default

  There is not a restriction of order for the z positions of the spheres or the final propagation
  '''

  N = len(radii)

  for n in range(N):
    if n==0:
      input_field = propagate(input_field, zs[n], lambda_, dxy, dxy, scale_factor)
      input_field = sphere_phase_shift(input_field, radii[n], xys[n], ns[n], lambda_, dxy, scale_factor, n0)
    else:
      input_field = propagate(input_field, zs[n]-zs[n-1], lambda_, dxy, dxy, scale_factor)
      input_field = sphere_phase_shift(input_field, radii[n], xys[n], ns[n], lambda_, dxy, scale_factor, n0)     

  if final_z==None:
    return input_field
 
  input_field = propagate(input_field, final_z-zs[-1], lambda_, dxy, dxy, scale_factor)

  return input_field


# Lines to load the image
# image = ut.imageRead('data/1.png')
# image = image. resize((512, 512))
# img_array = np.array(image)
# image = np.where(img_array < 128, 0, 1)
# image = Image.fromarray((image * 255).astype(np.uint8))

# ut.imageShow(image, 'Image')
# width, height = ut.imgInfo(image)

input_field = np.ones((1000, 1000))
radii = [3, 5]

xys = [(-60, 0), (60, 0)]

zs = [100, 250]

ns = [1.1, 1.1]

lambda_ = 0.633

dxy = 3.5

scale_factor = 40

final_z = 400


field = sphere_sample(input_field, radii, xys, zs, ns, lambda_, dxy, scale_factor, final_z = final_z)

# Define parameters
wavelength = 0.000633
dxy = 0.0065
distance = 0

fx_temp = 200
fy_temp = 200

#distance = np.arange(-50,50,1)
print(distance)

ref_Wave = sof.referenceWave(fx_temp, fy_temp, wavelength, dxy, 1000, 1000)
#obj_wave = sof.objectWave(image, wavelength, distance, dxy, method='AngularSpectrum', show=True)
hologram = sof.simuHolo(ref_Wave, field, show=True)

fftHolo = ut.ft(np.abs(hologram))
int_fft = ut.intensity(fftHolo, True)
ut.imageShow(int_fft, 'FFT')
