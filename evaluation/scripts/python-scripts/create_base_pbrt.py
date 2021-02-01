import os
import sys

def create_file(scene, accelerator):
    path = "/Users/marius/Documents/BachelorArbeit/pbrt-v3/scenes/" + scene 
    if os.path.exists(path + "/" + scene + ".pbrt"):
        f = open(path + "/eval_base.pbrt", "w")
        f.write("Accelerator \"" + accelerator + "\"\n\n")
        f.write("Include \"" + scene + ".pbrt\"\n")
        f.close()
    
def exec():
    scene = str(sys.argv[1])
    accelerator = str(sys.argv[2])
    create_file(scene, accelerator)

exec()