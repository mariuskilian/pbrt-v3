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
    cmake -S $SOURCE -B $BUILD
    make -C $BUILD -j

    SCENE_LIST=($(echo $SCENES | tr "," "\n"))
    for i in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
        SCENE=${SCENE_LIST[$i]}

        if [ ! -d $SCENE ]; then
            mkdir $SCENE
        fi
        RUN="$RUNPRE $SCENE $BUILD"
        $RUN "bvh" $NPIXELSAMPLES
        $RUN "bvh-bfs" $NPIXELSAMPLES
        $RUN "octree" $NPIXELSAMPLES
        $RUN "octree-bfs" $NPIXELSAMPLES
        $RUN "embree" $NPIXELSAMPLES
    done
fi

PLOT="python3 $PYSCRIPTS/plot_data.py $SCENES"
$PLOT prof --plot
$PLOT prof:total --plot
