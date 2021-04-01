set -e

SCENE=$1

INTGR="metric=nodes"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation
RUN="./../../scripts/start_eval.sh $INTGR $SCENE $BUILD"

for arg in "$@"; do
    if [[ $arg == "-n="* ]]; then
        NPIXELSAMPLES=$arg
    fi
done

if [ ! -d $SCENE ]
then
    mkdir $SCENE
fi

if [ ! -d "output" ]
then
    mkdir output
fi

if ! [[ $* == *--skip-render* ]]; then
    cmake -S $SOURCE -B $BUILD
    make -C $BUILD -j

    $RUN "bvh" $NPIXELSAMPLES
    $RUN "bvh-bfs" $NPIXELSAMPLES
    $RUN "octree" $NPIXELSAMPLES
fi

python3 $PYSCRIPTS/normalize_metrics.py $SCENE-$INTGR