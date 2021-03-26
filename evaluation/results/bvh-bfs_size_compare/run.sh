set -e

SCENES=$1

INTGR="path"

RUN="./../../scripts/start_eval.sh $INTGR"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation

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
    COUNT_STATS="-DCOUNT_STATS=False" # Because we want to compare time
    cmake -S $SOURCE -B $BUILD
    make -C $BUILD -j

    SCENE_LIST=($(echo $SCENES | tr "," "\n"))
    for i in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
        SCENE=${SCENE_LIST[$i]}

        if [ ! -d $SCENE ]; then
            mkdir $SCENE
        fi

        $RUN $SCENE $BUILD bvh $NPIXELSAMPLES
        $RUN $SCENE $BUILD bvh-bfs $NPIXELSAMPLES
    done
fi

PLOT="python3 $PYSCRIPTS/plot_data.py $SCENES"
$PLOT prof --plot
$PLOT memcomp:topology --plot