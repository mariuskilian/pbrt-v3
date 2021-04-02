import os
import sys
import imageio
import pyexr
import numpy as np
from PIL import Image
        
def exec():
    filepaths = []
    savepaths = []
    images = []
    maxvalues = []
    # Prep variable
    filelist = os.listdir("output")
    for file in filelist:
        if file.endswith(".exr"):
            filepath = os.path.join("output", file)
            savepath = sys.argv[0][:-len("exr2png.py")] + "../../plots/renders/"
            image = pyexr.open(filepath).get()
            images.append(image)
            maxvalues.append(np.max(image))
            filepaths.append(filepath)
            savepaths.append(savepath + file[:-3] + "png")
    for i in range(len(images)):
        #images[i] *= 16 / maxvalues[i]
        images[i] = np.where(images[i]<=0.0031308,12.92 * images[i], 1.055*(images[i]**(1/2.4)) - 0.055)
        images[i] = np.clip(images[i], 0, 1)
        images[i] = (images[i] * 255).astype(np.uint8)
        Image.fromarray(images[i]).save(savepaths[i])
    
exec()
