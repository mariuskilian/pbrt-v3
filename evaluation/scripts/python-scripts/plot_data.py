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

ECOSYS="Ecosys"
CROWN="Crown"
MEASUREONE="Measure-one"
HAIR="Hair"
LANDSCAPE="Landscape"

order={EMBREE:10, BVH:20, BVHBFS:30, KDTREE:40, OCTREE:50, OCTREEBFS:60}
sceneorder={ECOSYS:10, CROWN:20, MEASUREONE:30, HAIR:40, LANDSCAPE:50}

# === GET VALUES FROM LOG FILES ===

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

def get_totalProf(path, key):
    with open(path) as f:
        # Skip to Profile (flattened) section in file
        for line in f:
            if "Profile" in line:
                break
        # Find line with specified key
        time = 0
        num_stats = 0
        for line in f:
            if key in line:
                # find & convert (   0:00:00.00) substring to seconds
                m = re.search(r"\(\s*(\d*):(\d*):(\d*).(\d*)\)", line)
                time += 60 * 60 * int(m[1]) + 60 * int(m[2]) + int(m[3]) + int(m[4]) / 100
                num_stats += 1
            if num_stats == 2: return time

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
                if   "duplprim"     in stat:
                    # Duplicate Primitives                          11891548 /      3540215 (3.36x)
                    m = re.search(r"\s*\d+\s*/\s*\d+\s*\((\d+\.\d*)x\)", line)
                    return float(m[1])
                elif "primitive" in stat or "leafnode" in stat or "node" in stat or "chunklayer" in stat:
                    # Primitives - # Total                                         11891548
                    # Nodes - # Leaf                                                 769455
                    # Nodes - # Total (incl. implicit root node)                     879377
                    # Chunks - # Layers (excl. root chunk)                                6
                    m = re.search(r"(\d+)", line)
                    return int(m[1])
                elif stat == "chunk" or stat == "chunks":
                    # Chunks - # Total                                                43523
                    m = re.search(r"Total\s*(\d+)", line)
                    return int(m[1])
                elif "chunkfill" in stat:
                    # Chunks - Fill %                                                 31.570 avg [range 12.500000 - 100.000000]
                    m = re.search(r"\s*(\d+\.\d*)\s*avg\s*\[range\s*(\d+\.\d*)\s*-\s*(\d+\.\d*)", line)
                    if stat == "chunkfill": return float(m[1])
                    elif stat == "chunkfillmin": return float(m[2])
                    elif stat == "chunkfillmax": return float(m[3])
                elif key == "CHUNK_SIZE=":
                    m = re.search(r"CHUNK_SIZE=(\d+)", line)
                    return int(m[1])
                elif key == "Chosen Accelerator: ":
                    m = re.search(r"(BVH|Octree)", line)
                    return m[1]
                

# =========

# === GET META INFORMATION FOR PLOTS === 

def get_info(scenes, accellist, filelist, tp, stat):
    fullsavepath = [sys.argv[0].rstrip("plot_data.py") + "../../plots/"]
    test_name = str(os.getcwd()).split('/')[-1]
    test_type = ""
    if test_name.endswith("_stats") or test_name.endswith("_time"):
        test_type = test_name.split('_')[-1] + '_'
        test_name = test_name[:-len(test_name.split('_')[-1])].rstrip('_')
    fullsavepath.append(test_name)
    scenenames = ""
    for scene in sorted(scenes): scenenames += scene + ':'
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
            xlabel = "Acceleration Structure"
            xitems = accellist
    else:
        if all("Octree" in accel for accel in accellist):
            title += "Octree: "
        if all("BVH" in accel and "Embree" not in accel for accel in accellist):
            title += "BVH: "

        xlabel = "Acceleration Structure"
        xitems = accellist
    if len(scenes) > 1:
        xitems = [xitems[i] + "\n(" + filelist[i].split('/')[0].capitalize() + ')' for i in range(len(xitems))]

    # y label
    ylabel = ""
    # Time
    if tp == "prof" or tp == "time":
        if stat == "total": ylabel += "Intersection "
        else: ylabel += "Traversal "
        ylabel += "Time per Ray (μs)"
    # Stats
    if tp == "dist": ylabel += "Average "
    if tp == "stat" or tp  == "dist":
        ylabel += "Number of "
        if "leaf" in stat: ylabel += "Leaf "
        if "node" in stat: ylabel += "Node "
        else: ylabel += stat.capitalize() + " "
        ylabel += "Intersections"
    if tp == "dist": ylabel += " per Ray"
    if tp == "mem" or tp == "memcomp":
        ylabel += "Size of "
        if stat == "":
            ylabel += "Total Structure"
        else:
            ylabel += stat.capitalize()
        ylabel += " (MB)"
    if tp == "memprof":
        ylabel += "Execution Time (s) * Size of "
        if stat == "":
            ylabel += "Total Structure"
        else:
            ylabel += stat.capitalize()
        ylabel += " (MiB)"
    if tp == "accel":
        if   "duplprim"     in stat: ylabel = "Primitive Multiplication Factor"
        elif "primitive"    in stat: ylabel = "Number of Total Primitives"
        elif "leafnode"     in stat: ylabel = "Number of Leaf Nodes"
        elif "node"         in stat: ylabel = "Number of Total Nodes"
        elif "chunkfillmin" in stat: ylabel = "Minimum Chunk Fill Percentage (%)"
        elif "chunkfill"    in stat: ylabel = "Average Chunk Fill Percentage (%)"
        elif "chunklayer"   in stat: ylabel = "Depth of Chunk Tree"
        elif "chunk"        in stat: ylabel = "Total Number of Chunks"

    filler = " for each " if " per " in ylabel else " per "
    title += ylabel + filler + xlabel
    if tp == "memcomp":
        title = title.replace("Size of ", "Size Compression of ")
    if not tp.endswith("comp"):
        if len(scenes) <= 3:
            title += " for Scene"
            if len(scenes) > 1: title += 's'
            title += " "
            for scene in scenes: title += '\"' + scene.capitalize() + "\", "
            title = title[:-2]
        else:
            title += " for different Scenes"
    title = re.sub(r"\s\([^()]*\)", "", title)
    
    if tp.endswith("comp") or len(scenes) > 1: xlabel = "Scene"

    savepath += '.pdf'
    fullsavepath.append(savepath)
    return title, xlabel, ylabel, xitems, fullsavepath

# =========

# === PLOT CONFIGURATIONS ===

def bar_plot(ax, data, colors=None, total_width=0.8, single_width=1, legend=True):
    # Check if colors where provided, otherwhise use the default color cycle
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # Number of bars per group
    n_bars = len(data)
    # The width of a single bar
    bar_width = total_width / n_bars
    # List containing handles for the drawn bars, used for the legend
    bars = []
    # Iterate over all data
    for i, (_, values) in enumerate(data.items()):
        # The offset in x direction of that bar
        x_offset = (i - n_bars / 2) * bar_width + bar_width / 2
        # Draw a bar for every value of that type
        for x, y in enumerate(values):
            bar = ax.bar(x + x_offset, y, width=bar_width * single_width, color=colors[i % len(colors)])
        # Add a handle to the last drawn bar, which we'll need for the legend
        bars.append(bar[0])
    # Draw legend if we need
    if legend:
        ax.legend(bars, data.keys())

def plot_conf(ax, xitems, statlist):
    plt.bar(xitems, statlist)
    if len(statlist) > 5: plt.xticks(rotation=45)

def comp_plot_conf(ax, xitems, statlist):
    statlists = [statlist[i::2] for i in range(2)]
    scenes = {x.split('\n')[1] for x in xitems}
    legend = [xitems[i].split('\n')[0] for i in range(len(scenes))]
    xitems = [x.split('\n')[1][1:-1] for x in xitems[::2]]
    colrs = ['tab:blue', 'midnightblue']
    alpha = [0.6, 1]
    for i in range(2):
        plt.bar(xitems, statlists[i], color=colrs[i], alpha = alpha[i], label = legend[i])
    plt.legend()
    if len(statlists) > 5: plt.xticks(rotation=45)

def scenes_plot_conf(ax, xitems, statlist):
    dupl_scenes = [x.split('\n')[1][1:-1] for x in xitems]
    scenes = []
    for scene in dupl_scenes:
        if scene not in scenes: scenes.append(scene)
    numbars = int((len(xitems) + 1)/len(scenes))
    statlists = [statlist[i::numbars] for i in range(numbars)]
    legend = [xitems[i].split('\n')[0] for i in range(numbars)]
    data = {legend[i]:statlists[i] for i in range(len(legend))}
    bar_plot(ax, data)
    xitems = [scene for scene in scenes]
    plt.bar(xitems, len(xitems) * [0])


def visualize(title, xlabel, ylabel, xitems, savepath, statlist, plottype):
    # Configure based on plot type
    if plottype == "":
        plot = plot_conf
    elif plottype == "compare":
        plot = comp_plot_conf
    elif plottype == "scenes":
        plot = scenes_plot_conf
    fig = plt.figure(figsize=(5, 5))
    ax = fig.subplots()
    plot(ax, xitems, statlist)
    # Meta
    # plt.suptitle("\n".join(wrap(title, 55)), y=1)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

# === MAIN FUNCTION ===

def exec():
    # Format input params
    scenes = sys.argv[1].split(',')
    scenes.sort(key=lambda scene: sceneorder[scene.capitalize()])
    tp, stat = sys.argv[2], ""
    if ':' in tp:
        _ = tp.split(':')
        tp = _[0]
        stat = _[1]

    print("\nplot_data.py:\n\tProcessing scenes: ", scenes)
    print("\tType:", tp, stat)

    # Determine key (only needed for stat and dist and mem types)
    plottype = ""
    if tp == "stat" or tp == "dist":
        key = " - Intersects - "
        if "node" in stat:
            key += "Node - "
            if "leaf" in stat: key += "Leaf"
            else: key += "Total"
        else: key += stat.capitalize()
    elif "mem" in tp:
        if stat == "topology": key = " topology"
        else: key = " tree"
        if tp == "memcomp": plottype = "compare"
    elif tp == "accel":
        if   "duplprim"     in stat: key = "Duplicate Primitives"
        elif "primitive"    in stat: key = "Primitives - # Total"
        elif "leafnode"     in stat: key = "Nodes - # Leaf"
        elif "node"         in stat: key = "Nodes - # Total"
        elif "chunkfill"    in stat: key = "Chunks - Fill %"
        elif "chunklayer"   in stat: key = "Chunks - # Layers"
        elif "chunk"        in stat: key = "Chunks - # Total"
    if len(scenes) > 1: plottype = "scenes"

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
        def sort(file):
            return (sceneorder[file.split('/')[0].capitalize()], file.split('/')[1])
        filelist.sort(key = lambda f: sort(f))
    else:
        def sort(pair):
            return (sceneorder[pair[0].split('/')[0].capitalize()], order[pair[1]])
        filelist, accellist = zip(*[[f,accel] for f,accel in sorted(zip(filelist, accellist), key = lambda pair: sort(pair))])
        filelist = list(filelist)
        accellist = list(accellist)

    # Retrieve stat for every different file
    statlist = []
    for i in range(len(filelist)):
        fp = filelist[i] #filepath
        if tp == "prof":
            if stat == "total":
                time = get_totalProf(fp, "Accelerator::Intersect()")
            else:
                time = get_prof(fp, "Accelerator::Intersect()")
            nIsects = get_nIsects(fp)
            if time == None or nIsects == None:
                statlist.append(None)
            else:
                if time == "total":
                    statlist.append(1000 * 1000 * time / nIsects)
                else:
                    statlist.append(1000 * 1000 * time / nIsects)
        elif tp == "stat": statlist.append(get_stat(fp, key))
        elif tp == "dist":
            # value = get_dist(fp, key) messed up through instancing!!
            value = get_stat(fp, key)
            if value == None: statlist.append(None)
            else:
                numRays = get_nIsects(fp)
                statlist.append(value / numRays)
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
            if stat == "chunkfill":
                numnodes = get_accelInfo(fp, "Nodes - # Total", "node")
                numchunks = get_accelInfo(fp, "Chunks - # Total", "chunk")
                chunksize = get_accelInfo(fp, "CHUNK_SIZE=", "")
                accel = get_accelInfo(fp, "Chosen Accelerator: ", "")
                if accel == "BVH":
                    chunksize = max(chunksize, 64)
                    numnodesperchunk = 8 * (chunksize - 40)
                if accel == "Octree":
                    numnodesperchunk = 8 * (chunksize - 8)
                avgchunkfill = 100 * numnodes / (numchunks * numnodesperchunk)
                statlist.append(avgchunkfill)
            else:
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
    parser.add_argument("--debug", action='store_true', default=False)
    parser.add_argument("--tex", action='store_true', default=False)
    args, _ = parser.parse_known_args()

    if args.plot or args.show or args.debug:
        print("Plot information:")
        print("\tTitle:\t", title)
        print("\txLabel:\t", xlabel)
        print("\tyLabel:\t", ylabel)
        print("\txItems:\t", xitems)
        print("\tStats:\t", statlist)
        print("\tFiles:\t", filelist)
        print("\tAccels:\t", accellist)
        visualize(title, xlabel, ylabel, xitems, savepath, statlist, plottype)
        if args.plot:
            print("\tPath:\t" + savepath)
            plt.savefig(savepath, bbox_inches='tight', dpi=600)
        if args.show: plt.show()
        print('\n')

exec()

# ==========

## System arguments format (in that order):
##   scene: crown, killeroo, etc.
##   type: prof, time, stat:<name>, dist:<name>, <mem, memprof, memcomp>[:<memtype>], accel:<acceldata>,
##      name: primitive, chunk, leafnode, node
##      memtype: topology
##      acceldata: chunkfill, chunkfillmin, chunks, duplprims
