import os
import sys

def remove_file(scene, accelerator):
    path = "/Users/marius/Documents/BachelorArbeit/pbrt-v3/scenes/" + scene 
    if os.path.exists(path + "/eval_base.pbrt"):
        os.remove(path + "/eval_base.pbrt")

def exec():
    scene = str(sys.argv[1])
    accelerator = str(sys.argv[2])
    remove_file(scene, accelerator)

exec()