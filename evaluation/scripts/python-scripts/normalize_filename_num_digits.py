import sys
import os
import re

def normalize_names(scene):
    filelist = os.listdir(scene)
    max_magnitude = 1
    for f in filelist:
        magnitude = len(re.search(r"=(\d+).log", f)[1])
        max_magnitude = max(magnitude, max_magnitude)
    filenames = []
    for f in filelist:
        m = re.search(r"(.*=)(\d+).log", f)
        prefix = ""
        for _ in range(max_magnitude-len(m[2])): prefix += "0"
        filenames.append(m[1] + prefix + m[2] + ".log")
    for i in range(len(filenames)):
        wd = os.getcwd() + "/" + scene + "/"
        old_name = wd + filelist[i]
        new_name = wd + filenames[i]
        os.rename(old_name, new_name)


def exec():
    scene = sys.argv[1]
    normalize_names(scene)

exec()