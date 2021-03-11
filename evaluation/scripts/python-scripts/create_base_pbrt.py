import os
import sys

def create_file(integrator, scene, accelerator, extra_args=[]):
    path = str(sys.argv[0]).split("create_base_pbrt.py")[0]
    path += "../../../scenes/" + scene 
    if os.path.exists(path + "/" + scene + ".pbrt"):
        f = open(path + "/eval_base.pbrt", "w")
        if integrator.startswith("metric"):
            f.write("Sampler \"sobol\" \"integer pixelsamples\" 1\n\n")
            intgr = integrator.split('=')
            f.write("Integrator \"" + intgr[0] + "\"\n")
            f.write("\"string metric\" \"" + intgr[1] + "\"\n\n")
        else:
            f.write("Sampler \"sobol\" \"integer pixelsamples\" 4\n\n")
            if integrator == "path":
                f.write("Integrator \"path\" \"integer maxdepth\" 20\n\n")
            else:
                f.write("Integrator \"" + integrator + "\"\n\n")
        f.write("Accelerator \"" + accelerator + "\"\n\n")
        for line in extra_args:
            if line != "":
                var_type = line.split(':')[0]
                var_info = line.split(':')[1].split('=')
                f.write("\"" + var_type + " " + var_info[0] + "\" " + var_info[1] + "\n")
        f.write("Include \"" + scene + ".pbrt\"\n\n")
        f.close()
    
def exec():
    integrator = str(sys.argv[1])
    scene = str(sys.argv[2])
    accelerator = str(sys.argv[3])
    extra_args = str(sys.argv[4]).split(';')[1:]
    create_file(integrator, scene, accelerator, extra_args)

exec()