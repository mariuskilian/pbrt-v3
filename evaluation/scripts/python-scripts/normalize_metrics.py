import os
import sys
import pyexr
import numpy as np
from PIL import Image
from matplotlib import cm

# imag = pyexr.open("crown.exr").get()
# imag /= np.max(imag)
# imag = np.where(imag<=0.0031308,12.92 * imag, 1.055*(imag**(1/2.4)) - 0.055)
# viridis = cm.get_cmap('viridis', 255)
# imag = viridis(imag[:,:,0])
# imag = (imag * 255).astype(np.uint8)
# Image.fromarray(imag).save("crown_3.png")

def normalize_all(name):
    images = []
    filepaths = []
    maxvalue = 0
    for file in os.listdir("output"):
        if file.endswith(".exr") and file.startswith(name):
            filepath = os.path.join("output", file)
            image = pyexr.open(filepath).get()
            maxvalue = max(maxvalue, np.max(image))
            images.append(image)
            filepaths.append((fp := filepath.split("/"))[0] + "/nrm_" + fp[1][:-3] + "png")
    for i in range(len(images)):
        images[i] /= maxvalue
        images[i] = np.where(images[i]<=0.0031308,12.92 * images[i], 1.055*(images[i]**(1/2.4)) - 0.055)
        viridis = cm.get_cmap('viridis', 255)
        images[i] = viridis(images[i][:,:,0])
        images[i] = (images[i] * 255).astype(np.uint8)
        Image.fromarray(images[i]).save(filepaths[i])


def exec():
    name = str(sys.argv[1])
    normalize_all(name)

exec()