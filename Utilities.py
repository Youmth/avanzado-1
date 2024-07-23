# -*- coding: utf-8 -*-
"""
Title-->            Utility script
Authors-->          Raul Castaneda, Cartlos Trujillo, Ana Doblas
Date-->             10/08/2023
Universities-->     EAFIT University (Applied Optics Group)
                    UMASS (Optical Imaging Research Laboratory)
Links-->
"""

# import lybraries
import numpy as np
from PIL import Image, ImageOps
from matplotlib import pyplot as plt


# Function to read an image file from the disk
def imageRead(namefile):
    # inputs:
    # namefile - direction image to read
    imagen = Image.open(namefile)
    loadImage = ImageOps.grayscale(imagen)

    return loadImage


# Function to display an Image
def imageShow(image, title):
    # inputs:
    # image - The image to show
    # title - Title of the displayed image
    plt.imshow(image, cmap='gray'), plt.title(title)
    plt.show()

    return


# Function to compute the amplitude of a given complex field
def amplitude(complexField, log):
    # inputs:
    # complexField - The input complex field to compute the amplitude
    # log - boolean variable to determine if a log representation is applied
    out = np.abs(complexField)

    if log == True:
        out = 20 * np.log(out)

    return out


# Function to compute the Intensity of a given complex field
def intensity(complexField, log):
    # inputs:
    # complexField - The input complex field to compute the intensity
    # log - boolean variable to determine if a log representation is applied
    out = np.abs(complexField)
    out = out * out

    if log == True:
        out = 20 * np.log(out)
        out[out == np.inf] = 0
        out[out == -np.inf] = 0

    return out


# Function to compute the phase of a given complex field
def phase(complexField):
    # inputs:
    # complexField - The input complex field to compute the phase
    out = np.angle(complexField)

    return out


# Function to compute the Fourier Transform
def ft(field):
    # inputs:
    # field - The input to compute the Fourier Transform
    ft = np.fft.fft2(field)
    ft = np.fft.ifftshift(ft)

    return ft


# Function to get image information
def imgInfo(img):
    # inputs:
    # img - The input img to get the information
    width, height = img.size
    print(f"Image size: {width} x {height} pixels")

    return width, height


