import os
import sys

def remove_file(scene):
    path = str(sys.argv[0]).split("remove_base_pbrt.py")[0]
    path += "../../../scenes/" + scene 
    if os.path.exists(path + "/eval_base.pbrt"):
        os.remove(path + "/eval_base.pbrt")

def exec():
    scene = str(sys.argv[1])
    remove_file(scene)

exec()