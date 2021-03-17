set -e

SCENE=$1

INTGR="path"
SOURCE=../../..
BUILD=$SOURCE/build/Evaluation
RUN="./../../scripts/start_eval.sh $INTGR $SCENE $BUILD"

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

    $RUN "bvh"
    $RUN "bvh-bfs"
    $RUN "octree"
    $RUN "octree-bfs"
    $RUN "embree"
fi