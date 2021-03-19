import matplotlib.pyplot as plt
import re
import sys
import os
from textwrap import wrap
import argparse

order={"Embree BVH":0, "Basic BVH":10, "Quantized BVH":20, "k-d Tree":30, "Basic Octree":40, "1-Bit Octree":50}

def get_nIsects(path):
    with open(path) as f:
        for line in f:
            if "Regular ray intersection tests" in line:
                return int(re.search(r"(\d+)", line)[1])

def get_prof(path, key):
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

def get_mem(path, key):
    with open(path) as f:
        for line in f:
            if "  Memory" in line: break
        for line in f:
            if key in line:
                m = re.search(r"\s*(\d+.\d\d)\s*((?:Gi|Mi|k)B)", line)
                if m == None: return None
                if str(m[2]).startswith("Gi"):
                    return 1000 * float(m[1])
                elif str(m[2]).startswith("Mi"):
                    return float(m[1])
                elif str(m[2]).startswith('k'):
                    return float(m[1]) / 1000

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


def get_info(scene, accellist, filelist, tp, stat):
    savepath = sys.argv[0].rstrip("plot_data.py") + "/../../plots/"
    test_name = str(os.getcwd()).split('/')[-1]
    test_type = ""
    if test_name.endswith("_stats") or test_name.endswith("_time"):
        test_type = test_name.split('_')[-1] + '_'
        test_name = test_name[:-len(test_name.split('_')[-1])].rstrip('_')
    savepath += test_name + '_' + scene + '_' + test_type + tp

    title = ""

    # x label / x items
    if all(accel == accellist[0] for accel in accellist):
        title += accellist[0] + ": "
        if "maxprims=" in filelist[0]:
            xlabel = "Primitive Threshold"
            xitems = [str(int(re.search(r"maxprims=(\d+)", file)[1])) for file in filelist]
        elif "multthresh=" in filelist[0]:
            xlabel = "Multiplication Threshold"
            xitems = [str(float(re.search(r"multthresh=(\d.\d)", file)[1])) for file in filelist]
        elif "relkeys=" in filelist[0]:
            xlabel = "Quantization Key Root Node"
            xitems = ["Parent Node" if "=true" in f else "Chunk Root Node" for f in filelist]
        else:
            xlabel = sys.argv[3]
            xitems = [re.search(r"=?(\d.\d)", file)[1] for file in filelist]
    else:
        xlabel = "Acceleration Structure"
        xitems = accellist

    # y label
    ylabel = ""
    # Time
    if tp == "prof" or tp == "time": ylabel += "Intersection Time per Ray (μs)"
    # Stats
    if tp == "dist": ylabel += "Average "
    if tp == "stat" or tp  == "dist":
        savepath += '=' + stat
        ylabel += "Number of "
        if "leaf" in stat: ylabel += "Leaf "
        if "node" in stat: ylabel += "Node "
        else: ylabel += stat.capitalize() + " "
        ylabel += "Intersections"
    if tp == "dist": ylabel += " per Ray"
    if tp == "mem":
        ylabel += "Size of "
        if stat == "":
            ylabel += "Total Structure"
        else:
            savepath += '=' + stat
            ylabel += stat.capitalize()
        ylabel += " (MiB)"
    if tp == "memprof":
        ylabel += "Execution Time (s) * Size of "
        if stat == "":
            ylabel += "Total Structure"
        else:
            savepath += '=' + stat
            ylabel += stat.capitalize()
        ylabel += " (MiB)"


    filler = " for each " if " per " in ylabel else " per "
    title += ylabel + filler + xlabel + " for Scene \"" + scene.capitalize() + "\""
    title = re.sub(r"\s\([^()]*\)", "", title)

    savepath += '.pdf'
    return title, xlabel, ylabel, xitems, savepath

def plot(title, xlabel, ylabel, xitems, savepath, statlist):
    plt.figure(figsize=(5, 5))
    plt.suptitle("\n".join(wrap(title, 55)), y=1)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.bar(xitems, statlist)
    plt.xticks(rotation=45)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.savefig(savepath, bbox_inches='tight', dpi=600)

def exec():
    # Format input params
    scene = sys.argv[1]
    tp, stat = sys.argv[2], ""
    if ':' in tp:
        _ = tp.split(':')
        tp = _[0]
        stat = _[1]

    # Determine key (only needed for stat and dist and mem types)
    if tp == "stat" or tp == "dist":
        key = " - Intersects - "
        if "node" in stat:
            key += "Node - "
            if "leaf" in stat: key += "Leaf"
            else: key += "Total"
        else: key += stat.capitalize()
    elif tp == "mem" or tp == "memprof":
        if stat == "topology": key = " topology"
        else: key = " tree"

    # Sort filelist alphabetically
    _filelist = os.listdir(scene)
    filelist = []
    accellist = []
    for file in _filelist:
        if file.endswith(".log"):
            filelist.append(file)
            if "embree" in file: accellist.append("Embree BVH")
            elif "bvh" in file:
                if "-bfs" in file: accellist.append("Quantized BVH")
                else: accellist.append("Basic BVH")
            elif "octree" in file:
                if "-bfs" in file: accellist.append("1-Bit Octree")
                else: accellist.append("Basic Octree")
    if all(accel == accellist[0] for accel in accellist):
        filelist, accellist = zip(*[[f,accel] for f,accel in sorted(zip(filelist, accellist))])
    else:
        filelist, accellist = zip(*[[f,accel] for f,accel in sorted(zip(filelist, accellist),
                key = lambda pair: order[pair[1]])])
    filelist = list(filelist)
    accellist = list(accellist)

    # Retrieve stat for every different file
    statlist = []
    for file in filelist:
        fp = scene + "/" + file
        if (tp == "prof"):
            time = get_prof(fp, "Accelerator::Intersect()")
            nIsects = get_nIsects(fp)
            if time == None or nIsects == None:
                statlist.append(None)
            else:
                statlist.append(1000 * 1000 * time / nIsects)
        elif (tp == "stat"): statlist.append(get_stat(fp, key))
        elif (tp == "dist"):
            value = get_dist(fp, key)
            if value == None: statlist.append(None)
            else: statlist.append(value[0])
        elif (tp == "mem"): statlist.append(get_mem(fp, key))
        elif (tp == "memprof"):
            mem = get_mem(fp, key)
            if mem == None: statlist.append(None)
            else: statlist.append(mem * get_prof(fp, "Accelerator::Intersect()"))

    title, xlabel, ylabel, xitems, savepath = get_info(scene, accellist, filelist, tp, stat)
    # embree doesnt have some stats
    nostat_ids = []
    for i in range(len(statlist)):
        if statlist[i] == None: nostat_ids.append(i-len(nostat_ids))
    for idx in nostat_ids:
        del xitems[idx]
        del statlist[idx]
    # :(
    
    parser = argparse.ArgumentParser(description="Process render data from pbrt")
    parser.add_argument("--plot", action='store_true', default=False)
    parser.add_argument("--tex", action='store_true', default=False)
    args, _ = parser.parse_known_args()

    if args.plot:
        plot(title, xlabel, ylabel, xitems, savepath, statlist)

exec()

## System arguments format (in that order):
##   scene: crown, killeroo, etc.
##   type: prof, time, stat:<name>, dist:<name>, mem[:<memtype>]
##      name: primitive, chunk, leafnode, node
##      memtype: topology
