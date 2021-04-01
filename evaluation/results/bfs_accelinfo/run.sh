set -e

SCENES=$1

INTGR="path"

RUN="./../../scripts/start_eval.sh $INTGR"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation

NPIXELSAMPLES="-n=1"

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

        $RUN $SCENE $BUILD bvh-bfs $NPIXELSAMPLES
        $RUN $SCENE $BUILD octree-bfs $NPIXELSAMPLES
    done
fi

PLOT="python3 $PYSCRIPTS/plot_data.py $SCENES"
$PLOT accel:duplprim --plot
$PLOT accel:primitive --plot
$PLOT accel:leafnode --plot
$PLOT accel:node --plot
$PLOT accel:chunkfillmin --plot
$PLOT accel:chunkfill --plot
$PLOT accel:chunklayer --plot
$PLOT accel:chunk --plot