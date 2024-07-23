# -*- coding: utf-8 -*-
"""
Title-->            Utility script
Authors-->          Raul Castaneda, Cartlos Trujillo, Ana Doblas
Date-->             10/08/2023
Universities-->     EAFIT University (Applied Optics Group)
                    UMASS (Optical Imaging Research Laboratory)
Links-->
"""

import math
import numpy as np
import Utilities as ut


#  Function to simulate the reference beam
def referenceWave(fx_temp, fy_temp, wavelength, dxy, width, height):
    # inputs:
    # fx_temp - pixel location +1 term in x-axis
    # fy_temp - pixel location -1 term in y-axis
    # wavelength - wavelength used to simulate the hologram
    # dxy - pixel pitch
    # width -
    # height -

    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY

    fx_0 = height / 2
    fy_0 = width / 2
    k = (2 * math.pi) / wavelength
    theta_x = math.asin((fx_0 - fx_temp) * wavelength / (height * dxy))
    theta_y = math.asin((fy_0 - fy_temp) * wavelength / (width * dxy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dxy) + (math.sin(theta_y) * Y * dxy)))

    return ref_wave


#  Function to simulate the object wave
def objectWave(image, wavelength, distance, dxy, method, show):
    # inputs:
    # image - image to create the complex object field
    # wavelength - wavelength used to simulate the hologram
    # distance - distance for propagation the optical field
    # method - propagator (Angular Spectrum of Fresnel)
    # image = np.array(image)/255
    objectWave = np.exp(1j*np.array(image))
    # objectWave = image
    if method == 'AngularSpectrum':
        objectWave = angularSpectrum(objectWave, wavelength, distance, dxy)
    elif method == 'Fresnel':
        objectWave = fresnel(objectWave, wavelength, distance, dxy)
    else:
        objectWave

    if show == True:
        objePropagate = ut.amplitude(objectWave, False)
        ut.imageShow(objePropagate, 'propagate object')

    return objectWave


#  Function to simulate the hologram
def simuHolo(ref_wave, object, show):
    hologram = np.power((ref_wave + object), 2)

    if show == True:
        hologram_show = ut.intensity(hologram, False)
        ut.imageShow(hologram_show, 'Hologram')

    return hologram


# Function to propagate an optical field using the Fresnel approach
def fresnel(field, wavelength, z, dxy):
    # Inputs:
    # field - complex field
    # wavelength - wavelength
    # z - propagation distance
    # dxy - sampling pitches
    field = np.array(field)
    M, N = field.shape
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')

    dxout = (wavelength * z) / (M * dxy)
    dyout = (wavelength * z) / (N * dxy)

    z_phase = np.exp(1j * 2 * math.pi * z / wavelength) / (1j * wavelength * z)

    out_phase = np.exp((1j * math.pi / (wavelength * z)) * (np.power(X * dxout, 2) + np.power(Y * dyout, 2)))
    in_phase = np.exp((1j * math.pi / (wavelength * z)) * (np.power(X * dxy, 2) + np.power(Y * dxy, 2)))

    tmp = (field * in_phase)
    tmp = np.fft.fftshift(tmp)
    tmp = np.fft.fft2(tmp)
    tmp = np.fft.fftshift(tmp)

    out = z_phase * out_phase * dxy * dxy * tmp

    return out


# Function to propagate an optical field using the Angular Spectrum approach
def angularSpectrum(field, wavelength, z, dxy):
    # Inputs:
    # field - complex field
    # wavelength - wavelength
    # z - propagation distance
    # dxy - sampling pitches
    field = np.array(field)
    M, N = field.shape
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')

    dfx = 1 / (dxy * M)
    dfy = 1 / (dxy * N)

    field_spec = np.fft.fftshift(field)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)

    kernel = np.power(1 / wavelength, 2) - (np.power(X * dfx, 2) + np.power(Y * dfy, 2))
    phase = np.exp(1j * z * 2 * math.pi * np.sqrt(kernel))

    tmp = field_spec * phase
    out = np.fft.ifftshift(tmp)
    out = np.fft.ifft2(out)
    out = np.fft.ifftshift(out)

    return out


# Function to detemrine the z critical value
def z_critical(field, wavelength, dxy):
    field = np.array(field)
    M, N = field.shape
    z = 0

    return z