import matplotlib.pyplot as plt
import re
import sys
import os

# order={"embree":0, "bvh":10, "bvh-bfs":20, "kdtree":30, "octree":40, "octree-bfs":50}

def get_prof(path, key):
    # extract value of flattened profile key in seconds from file path
    with open(path) as f:
        # Skip to Profile (flattened) section in file
        for line in f:
            if "Profile (flattened)" in line:
                break
        # Find line with specified key
        for line in f:
            if key in line:
                # find & convert (   0:00:00.00) substring to seconds
                m = re.search(r"\(\s*(\d*):(\d*):(\d*).(\d*)\)", line)
                return 60 * 60 * int(m[1]) + 60 * int(m[2]) + int(m[3]) + int(m[4]) / 100

def get_stat(path, category):
    with open(path) as f:
        for line in f:
            if category in line:
                break
        for line in f:
            if "Total #" in line:
                m = re.search(r"\d+", line)
                return int(m[0])

def get_dist(path, category):
    with open(path) as f:
        for line in f:
            if category in line: break
        for line in f:
            if "Distribution" in line:
                m = re.search(r"(\d+.\d\d\d) avg \[range (\d+.?\d*) - (\d+.?\d*)\]", line)
                return float(m[1]), int(m[2]), int(m[3])

def get_yLabel(tp, stat):
    # Time
    if tp == "prof" or tp == "time": return "Execution Time (s)"
    # Stats
    ylabel = ""
    if tp == "dist": ylabel += "Average "
    if tp == "stat" or tp  == "dist":
        ylabel += "Number of "
        if "leaf" in stat: ylabel += "Leaf "
        if "node" in stat: ylabel += "Node "
        else: ylabel += stat.capitalize() + " "
        ylabel += "Intersections"
    if tp == "dist": ylabel += " per Ray"
    return ylabel

def get_xLabel(accellist, filelist):
    if all(accel == accellist[0] for accel in accellist):
        if "maxprims=" in filelist[0]:
            xlabel = "Primitive Threshold"
            xitems = [re.search(r"maxprims=(\d+)", file)[1] for file in filelist]
        elif "multthresh=" in filelist[0]:
            xlabel = "Multiplication Threshold"
            xitems = [re.search(r"multthresh=(\d.\d)", file)[1] for file in filelist]
        else:
            xlabel = sys.argv[3]
            xitems = [re.search(r"=?(\d.\d)", file)[1] for file in filelist]
        return xlabel, xitems
    else: return "Acceleration Structure", accellist

def exec():
    # Format input params
    scene = sys.argv[1]
    tp, stat = sys.argv[2], ""
    if ':' in tp:
        _ = tp.split(':')
        tp = _[0]
        stat = _[1]
    # Determine key (only needed for stat and dist types)
    if tp == "stat" or tp == "dist":
        key = " - Intersects - "
        if "node" in stat:
            key += "Node - "
            if "leaf" in stat: key += "Leaf"
            else: key += "Total"
        else: key += stat.capitalize()
    # Sort filelist alphabetically
    filelist = os.listdir(scene)
    accellist = []
    for file in filelist:
        if file.endswith(".log"):
            if "embree" in file: accellist.append("embree")
            elif "bvh" in file: accellist.append("bvh")
            elif "octree" in file: accellist.append("octree")
            if "-bfs" in file: accellist[-1] += "-bfs"
    filelist, accellist = zip(*[[f,accel] for f,accel in sorted(zip(filelist, accellist))])
    filelist = list(filelist)
    accellist = list(accellist)
    # Retrieve stat for every different file
    statlist = []
    for file in filelist:
        fp = scene + "/" + file
        if (tp == "prof"):
            statlist.append(get_prof(fp, "Accelerator::Intersect()"))
        elif (tp == "stat"):
            statlist.append(get_stat(fp, key))
        elif (tp == "dist"):
            statlist.append(get_dist(fp, key)[0])
    plt.ylabel(get_yLabel(tp, stat))
    xLabel, xItems = get_xLabel(accellist, filelist)
    plt.xlabel(xLabel)
    plt.bar(xItems, statlist)
    plt.show()

exec()

## System arguments format (in that order):
##   scene: crown, killeroo, etc.
##   type: prof, time, stat:<name>, dist:<name>
##      name: primitive, chunk, leafnode, node

# def plot_prof():
#     key = "Accelerator::Intersect"
#     bvh_time = get_prof("path-bvh.log", key)
#     kdtree_time = get_prof("kdtree.log", key)
#     plt.ylabel("Execution Time (s)")
#     plt.bar(["BVH", "KdTree"], [bvh_time, kdtree_time])
#     plt.show()