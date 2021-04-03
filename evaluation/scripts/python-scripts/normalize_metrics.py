import os
import sys
import pyexr
import numpy as np
from PIL import Image
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import argparse

plots_path = ""
images = []
filepaths = []
maxvalue = 0

order={"embree":0, "bvh":10, "bvh-bfs":20, "kdtree":30, "octree":40, "octree-bfs":50}

def read_and_save_input(prefix):
    global plots_path, images, filepaths, maxvalue
    # Prep variable
    filelist = os.listdir("output")
    for file in filelist:
        if file.endswith(".exr") and file.startswith(prefix):
            filepath = os.path.join("output", file)
            image = pyexr.open(filepath).get()
            maxvalue = max(maxvalue, np.max(image))
            images.append(image)
            filepaths.append("output/" + file[:-3] + "png")
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
        #Image.fromarray(images[i]).save(filepaths[i])

def determine_paths(scene, notitle=False, shorttitle=False):
    global plots_path
    fullsavepath = [sys.argv[0].rstrip("normalize_metrics.py") + "../../plots"]
    fullsavepath.append("visualizations")
    fullsavepath.append(scene)
    test_name = str(os.getcwd()).split('/')[-1][4:]
    savepath = ""
    for path in fullsavepath:
        savepath += path + "/"
        if not os.path.exists(savepath): os.mkdir(savepath)
    plots_path = savepath
    savepath += test_name
    if shorttitle: savepath += "_shorttitle"
    if notitle: savepath += "_notitle"
    savepath += ".pdf"
    return savepath

def make_plot(prefix, scene, savepath, notitle=False, shorttitle=False):
    fig, axes = plt.subplots(nrows=1, ncols=len(images), figsize=(3*(len(images) + 1), 5))
    # Set Figure title
    t = prefix.split('=')[-1].split('-')[0]
    title = "Number of "
    if t == "primitives": title += "Primitives"
    elif t == "nodes": title += "Nodes"
    elif t == "leafnodes": title += "Leaf Nodes"
    title += " Intersected per Pixel for "
    if shorttitle:
        title += "different Scenes"
    else:
        title += "Scene \"" + scene.capitalize() + "\""
    if not notitle: fig.suptitle(title, fontsize=16)
    # add images to subplot
    for i in range(len(axes.flat)):
        ax = axes.flat[i]
        ax.set_axis_off()
        im = ax.imshow(images[i], cmap='viridis',
                    vmin=0, vmax=maxvalue)
        t = filepaths[i].split("/")[-1].split(prefix + '-')[-1].split(".png")[0]
        if t == "bvh": title = "Basic BVH"
        elif t == "octree": title = "Basic Octree"
        elif t == "bvh-bfs": title = "Quantized BVH"
        elif t == "octree-bfs": title = "1-Bit Octree"
        elif t == "embree": title = "Embree BVH"
        if not notitle: ax.title.set_text(title)
    # Adjust subplot dimensions etc.
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                        wspace=0.02, hspace=0.02)
    # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
    cb_ax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    fig.colorbar(im, cax=cb_ax)
    # cbar.set_ticks(np.arange(0,int(maxvalue),int(maxvalue/4)))
    # save plot
    plt.savefig(savepath, bbox_inches='tight', dpi=600)
        

def exec():
    parser = argparse.ArgumentParser(description="Process render data from pbrt")
    parser.add_argument("--notitle", action='store_true', default=False)
    parser.add_argument("--shorttitle", action='store_true', default=False)
    args, _ = parser.parse_known_args()

    scene = str(sys.argv[1])
    integrator = str(sys.argv[2])

    print("\nnormalize_metrics.py:\n\tScene: ", scene)
    print("\tType:", integrator.split('=')[1])

    prefix = scene + '-' + integrator
    savepath = determine_paths(scene)
    print("\tPath:", savepath)
    print("Converting Images to PNG...")
    read_and_save_input(prefix)
    print("Normalizing and Colormapping Images...")
    normalize_all()
    print("Creating Plot...")
    make_plot(prefix, scene, savepath, args.notitle, args.shorttitle)
    print("Done!\n")

exec()
