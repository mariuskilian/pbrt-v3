import os
import sys

def create_file(scene, accelerator, extra_args=[]):
    path = str(sys.argv[0]).split("create_base_pbrt.py")[0]
    path += "../../../scenes/" + scene 
    print(path)
    if os.path.exists(path + "/" + scene + ".pbrt"):
        f = open(path + "/eval_base.pbrt", "w")
        f.write("Accelerator \"" + accelerator + "\"\n\n")
        for line in extra_args:
            f.write(line + "\n")
        f.write("Include \"" + scene + ".pbrt\"\n")
        f.close()
    
def exec():
    scene = str(sys.argv[1])
    accelerator = str(sys.argv[2])
    extra_args = str(sys.argv[3]).split(';')
    create_file(scene, accelerator, extra_args)

exec()