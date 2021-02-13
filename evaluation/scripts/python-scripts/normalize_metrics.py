import os
import sys
import pyexr
import numpy as np
from PIL import Image
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# imag = pyexr.open("crown.exr").get()
# imag /= np.max(imag)
# imag = np.where(imag<=0.0031308,12.92 * imag, 1.055*(imag**(1/2.4)) - 0.055)
# viridis = cm.get_cmap('viridis', 255)
# imag = viridis(imag[:,:,0])
# imag = (imag * 255).astype(np.uint8)
# Image.fromarray(imag).save("crown_3.png")

images = []
filepaths = []
maxvalue = 0

order={"embree":0, "bvh":10, "bvh-bfs":20, "kdtree":30, "octree":40, "octree-bfs":50}

def read_and_save_input(prefix):
    global images, filepaths, maxvalue
    # Prep variable
    filelist = os.listdir("output")
    for file in filelist:
        if file.endswith(".exr") and file.startswith(prefix):
            filepath = os.path.join("output", file)
            image = pyexr.open(filepath).get()
            maxvalue = max(maxvalue, np.max(image))
            images.append(image)
            fp = filepath.split("/")
            filepaths.append(fp[0] + "/nrm_" + fp[1][:-3] + "png")
    # Sort by order described in order dictionary above
    filepaths, images = zip(*[[fp,img] for fp,img in sorted(zip(filepaths,images),
        key=lambda pair: order[pair[0].split(prefix)[-1].strip('-').split(".png")[0]])])
    filepaths = list(filepaths)
    images = list(images)

def normalize_all():
    global images, filepaths, maxvalue
    for i in range(len(images)):
        images[i] /= maxvalue
        images[i] = np.where(images[i]<=0.0031308,12.92 * images[i], 1.055*(images[i]**(1/2.4)) - 0.055)
        viridis = cm.get_cmap('viridis', 255)
        images[i] = viridis(images[i][:,:,0])
        images[i] = (images[i] * 255).astype(np.uint8)
        Image.fromarray(images[i]).save(filepaths[i])

def make_plot(prefix):
    fig=plt.figure(figsize=(16, 4))
    # Set Figure title
    t = prefix.split('=')[-1].split('-')[0]
    title = "Number of "
    if t == "primitives": title += "Primitives"
    elif t == "nodes": title += "Nodes"
    elif t == "leafnodes": title += "Leaf Nodes"
    title += " Intersected per Pixel"
    fig.suptitle(title, fontsize=16)
    # Add all images as subplots
    num_images = len(images)
    for i in range(1, num_images + 1):
        fig.add_subplot(1, num_images + 1, i)
        plt.imshow(images[i-1])
        plt.axis('off')
        t = filepaths[i-1].split("/")[-1].split(prefix + '-')[-1].split(".png")[0]
        if t == "bvh": title = "Original BVH"
        elif t == "octree": title = "Original Octree"
        elif t == "bvh-bfs": title = "Compressed BVH"
        elif t == "octree-bfs": title = "Compressed Octree"
        elif t == "embree": title = "Embree BVH"
        plt.title(title)
    # Add placeholder image to show colorbar
    fig.add_subplot(1, num_images, num_images)
    plt.imshow(np.ones((1, 1, 3)), vmin = 0, vmax = maxvalue)
    plt.axis('off')
    plt.colorbar()
    plt.savefig(filepaths[0].split("/nrm_")[0] + "/plot_" + prefix + ".png")
        

def exec():
    prefix = str(sys.argv[1])
    read_and_save_input(prefix)
    normalize_all()
    make_plot(prefix)

exec()
