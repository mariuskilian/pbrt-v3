set -e

SCENE=$1

INTGR="metric=primitives"
RUN="./../../scripts/start_eval.sh $INTGR $SCENE"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation

if [ ! -d $SCENE ]
then
    mkdir $SCENE
fi

if [ ! -d "output" ]
then
    mkdir output
fi

if ! getopts p flag; then
    cmake -S $SOURCE -B $BUILD -DCOUNT_STATS=True -DREL_KEYS=True -DBF_SIZE=64 -DCHUNK_SIZE=64
    make -C $BUILD -j

    $RUN "bvh"
    $RUN "bvh-bfs"
    $RUN "octree-bfs"
    $RUN "embree"
fi

python3 $PYSCRIPTS/normalize_metrics.py $SCENE-$INTGR