set -e

SCENE=$1

INTGR="path"
SOURCE=../../..
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
    cmake -S $SOURCE -B $BUILD -DCOUNT_STATS=True -DREL_KEYS=True -DBF_SIZE=64 -DCHUNK_SIZE=64
    make -C $BUILD -j

    $RUN "bvh" $NPIXELSAMPLES
    $RUN "bvh-bfs" $NPIXELSAMPLES
    $RUN "octree" $NPIXELSAMPLES
    $RUN "octree-bfs" $NPIXELSAMPLES
    $RUN "embree" $NPIXELSAMPLES
fi