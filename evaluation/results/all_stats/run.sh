set -e

SCENES=$1

INTGR="path"
SOURCE=../../..
BUILD=$SOURCE/build/Evaluation
PYSCRIPTS=../../scripts/python-scripts
RUNPRE="./../../scripts/start_eval.sh $INTGR"

for arg in "$@"; do
    if [[ $arg == "-n="* ]]; then
        NPIXELSAMPLES=$arg
    fi
done

if [ ! -d "output" ]
then
    mkdir output
fi

if ! [[ $* == *--skip-render* ]]; then
    cmake -S $SOURCE -B $BUILD -DCOUNT_STATS=True
    make -C $BUILD -j

    SCENE_LIST=($(echo $SCENES | tr "," "\n"))
    for i in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
        SCENE=${SCENE_LIST[$i]}

        if [ ! -d $SCENE ]; then
            mkdir $SCENE
        fi
        RUN="$RUNPRE $SCENE $BUILD"
        $RUN "bvh" -n=1
        $RUN "bvh-bfs" -n=1
        $RUN "octree" -n=1
    done
fi

PLOT="python3 $PYSCRIPTS/plot_data.py $SCENES"
$PLOT stat:primitive --plot
$PLOT stat:leafnode --plot
$PLOT stat:node --plot
$PLOT dist:primitive --plot
$PLOT dist:leafnode --plot
$PLOT dist:node --plot
$PLOT accel:primitive --plot
$PLOT accel:leafnode --plot
$PLOT accel:node --plot