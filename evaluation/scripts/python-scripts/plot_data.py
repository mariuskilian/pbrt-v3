import matplotlib.pyplot as plt
import re
import sys
import os
from textwrap import wrap
import argparse

OCTREE = "Basic Octree"
OCTREEBFS = "1-Bit Octree"
BVH = "Basic BVH"
BVHBFS = "Quantized BVH"
EMBREE = "Embree BVH"
KDTREE = "k-d Tree"

order={EMBREE:10, BVH:20, BVHBFS:30, KDTREE:40, OCTREE:50, OCTREEBFS:60}

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

def get_accelInfo(path, key, stat):
    with open(path) as f:
        for line in f:
            if key in line:
                if stat == "duplprims":
                    # Duplicate Primitives                          11891548 /      3540215 (3.36x)
                    m = re.search(r"\s*\d+\s*/\s*\d+\s*\((\d+\.\d*)x\)", line)
                    return float(m[1])
                elif stat == "chunks":
                    # Chunks - # Total                                                43523
                    m = re.search(r"Total\s*(\d+)", line)
                    return int(m[1])
                elif "chunkfill" in stat:
                    # Chunks - Fill %                                                 31.570 avg [range 12.500000 - 100.000000]
                    m = re.search(r"\s*(\d+\.\d*)\s*avg\s*\[range\s*(\d+\.\d*)\s*-\s*(\d+\.\d*)", line)
                    if stat == "chunkfill": return float(m[1])
                    elif stat == "chunkfillmin": return float(m[2])
                    elif stat == "chunkfillmax": return float(m[3])

def get_info(scenes, accellist, filelist, tp, stat):
    fullsavepath = [sys.argv[0].rstrip("plot_data.py") + "../../plots/"]
    test_name = str(os.getcwd()).split('/')[-1]
    test_type = ""
    if test_name.endswith("_stats") or test_name.endswith("_time"):
        test_type = test_name.split('_')[-1] + '_'
        test_name = test_name[:-len(test_name.split('_')[-1])].rstrip('_')
    fullsavepath.append(test_name)
    scenenames = ""
    for scene in scenes: scenenames += scene + ':'
    fullsavepath.append(scenenames[:-1])
    savepath = test_type + tp
    if stat != "": savepath += '=' + stat

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
        elif "CHUNK_SIZE" in filelist[0]:
            xlabel = "Chunk Size (Bytes)"
            xitems = [str(int(re.search(r"CHUNK_SIZE=(\d+)", file)[1])) for file in filelist]
        elif "BF_SIZE=" in filelist[0]:
            xlabel = "Bitfield Data Type"
            xitems = ["uint" + str(int(re.search(r"BF_SIZE=(\d+)", file)[1])) for file in filelist]
        else:
            xlabel = sys.argv[3]
            xitems = [re.search(r"=?(\d.\d)", file)[1] for file in filelist]
    else:
        xlabel = "Acceleration Structure"
        xitems = accellist
    if len(scenes) > 1:
        xitems = [xitems[i] + "\n(" + filelist[i].split('/')[0].capitalize() + ')' for i in range(len(xitems))]

    # y label
    ylabel = ""
    # Time
    if tp == "prof" or tp == "time": ylabel += "Intersection Time per Ray (μs)"
    # Stats
    if tp == "dist": ylabel += "Average "
    if tp == "stat" or tp  == "dist":
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
            ylabel += stat.capitalize()
        ylabel += " (MiB)"
    if tp == "memprof":
        ylabel += "Execution Time (s) * Size of "
        if stat == "":
            ylabel += "Total Structure"
        else:
            ylabel += stat.capitalize()
        ylabel += " (MiB)"
    if tp == "accel":
        if stat == "chunks":
            ylabel += "Total Number of Chunks"
        elif stat == "chunkfill":
            ylabel += "Average Chunk Fill Percentage (%)"
        elif stat == "chunkfillmin":
            ylabel += "Minimum Chunk Fill Percentage (%)"
        elif stat == "duplprims":
            ylabel += "Primitive Multiplication Factor"

    filler = " for each " if " per " in ylabel else " per "
    title += ylabel + filler + xlabel + " for Scene"
    if len(scenes) > 1: title += 's'
    title += " "
    for scene in scenes: title += '\"' + scene.capitalize() + "\", "
    title = title[:-2]
    title = re.sub(r"\s\([^()]*\)", "", title)

    savepath += '.pdf'
    fullsavepath.append(savepath)
    return title, xlabel, ylabel, xitems, fullsavepath

def plot_conf(title, xlabel, ylabel, xitems, savepath, statlist):
    plt.figure(figsize=(5, 5))
    plt.suptitle("\n".join(wrap(title, 55)), y=1)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.bar(xitems, statlist)
    if len(statlist) > 5: plt.xticks(rotation=45)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

def plot(title, xlabel, ylabel, xitems, savepath, statlist):
    plot_conf(title, xlabel, ylabel, xitems, savepath, statlist)
    plt.savefig(savepath, bbox_inches='tight', dpi=600)

def show(title, xlabel, ylabel, xitems, savepath, statlist):
    plot_conf(title, xlabel, ylabel, xitems, savepath, statlist)
    plt.show()

def exec():
    # Format input params
    scenes = sys.argv[1].split(',')
    tp, stat = sys.argv[2], ""
    if ':' in tp:
        _ = tp.split(':')
        tp = _[0]
        stat = _[1]

    print("\nplot_data.py:\n\tProcessing scenes: ", scenes)
    print("\tType:", tp, stat)

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
    elif tp == "accel":
        if "chunkfill" in stat: key = "Chunks - Fill %"
        elif stat == "chunks": key = "Chunks - # Total"
        elif stat == "duplprims": key = "Duplicate Primitives"

    # Sort filelist alphabetically
    filelist = []
    accellist = []
    for scene in scenes:
        _filelist = os.listdir(scene)
        for file in _filelist:
            if file.endswith(".log"):
                filelist.append(scene + '/' + file)
                if "embree" in file: accellist.append(EMBREE)
                elif "bvh" in file:
                    if "-bfs" in file: accellist.append(BVHBFS)
                    else: accellist.append(BVH)
                elif "octree" in file:
                    if "-bfs" in file: accellist.append(OCTREEBFS)
                    else: accellist.append(OCTREE)
    if all(accel == accellist[0] for accel in accellist):
        filelist, accellist = zip(*[[f,accel] for f,accel in sorted(zip(filelist, accellist))])
    else:
        def sort(pair):
            res = ""
            if len(scenes) > 1: res += pair[0].split('/')[0]
            res += str(order[pair[1]])
            return res
        filelist, accellist = zip(*[[f,accel] for f,accel in sorted(zip(filelist, accellist), key = lambda pair: sort(pair))])
    filelist = list(filelist)
    accellist = list(accellist)

    # Retrieve stat for every different file
    statlist = []
    for i in range(len(filelist)):
        fp = filelist[i] #filepath
        if tp == "prof":
            time = get_prof(fp, "Accelerator::Intersect()")
            nIsects = get_nIsects(fp)
            if time == None or nIsects == None:
                statlist.append(None)
            else:
                statlist.append(1000 * 1000 * time / nIsects)
        elif tp == "stat": statlist.append(get_stat(fp, key))
        elif tp == "dist":
            value = get_dist(fp, key)
            if value == None: statlist.append(None)
            else: statlist.append(value[0])
        elif "mem" in tp:
            fullkey = ""
            if accellist[i] == OCTREEBFS: fullkey += "Octree-BFS"
            elif accellist[i] == OCTREE: fullkey += "Octree"
            elif accellist[i] == BVHBFS: fullkey += "BVH-BFS"
            elif accellist[i] == BVH: fullkey += "BVH"
            fullkey += key
            mem = get_mem(fp, fullkey)
            if tp == "memprof":
                if mem == None: statlist.append(None)
                else: statlist.append(mem * get_prof(fp, "Accelerator::Intersect()"))
            else: statlist.append(mem)
        elif tp == "accel":
            accelinfo = get_accelInfo(fp, key, stat)
            statlist.append(accelinfo)

    title, xlabel, ylabel, xitems, fullsavepath = get_info(scenes, accellist, filelist, tp, stat)
    savepath = ""
    for path in fullsavepath[:-1]:
        savepath += path + "/"
        if not os.path.exists(savepath): os.mkdir(savepath)
    savepath += fullsavepath[-1]
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
    parser.add_argument("--show", action='store_true', default=False)
    parser.add_argument("--tex", action='store_true', default=False)
    args, _ = parser.parse_known_args()

    if args.plot or args.show:
        print("Plot information:")
        print("\tTitle:\t", title)
        print("\txLabel:\t", xlabel)
        print("\tyLabel:\t", ylabel)
        print("\txItems:\t", xitems)
        print("\tStats:\t", statlist)
        print("\tFiles:\t", filelist)
        print("\tAccels:\t", accellist)
        if args.plot:
            print("\tPath:\t" + savepath)
            plot(title, xlabel, ylabel, xitems, savepath, statlist)
        if args.show:
            show(title, xlabel, ylabel, xitems, savepath, statlist)
        print('\n')

exec()

## System arguments format (in that order):
##   scene: crown, killeroo, etc.
##   type: prof, time, stat:<name>, dist:<name>, mem[:<memtype>], accel:<acceldata>
##      name: primitive, chunk, leafnode, node
##      memtype: topology
##      acceldata: chunkfill, chunkfillmin, chunks, duplprims
