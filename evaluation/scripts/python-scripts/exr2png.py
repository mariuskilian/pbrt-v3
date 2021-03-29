import os
import sys
import pyexr
from PIL import Image

order={"embree":0, "bvh":10, "bvh-bfs":20, "kdtree":30, "octree":40, "octree-bfs":50}
        
def exec():
    savepaths = []
    images = []
    # Prep variable
    filelist = os.listdir("output")
    for file in filelist:
        if file.endswith(".exr"):
            filepath = os.path.join("output", file)
            savepath = sys.argv[0][:-len("exr2png.py")] + "../../plots/renders/"
            image = pyexr.open(filepath).get()
            images.append(image)
            savepaths.append(savepath + file[:-3] + "png")
    for i in range(len(savepaths)):
        im = Image.fromarray(images[i])
        im.save(savepaths[i])
    
exec()
