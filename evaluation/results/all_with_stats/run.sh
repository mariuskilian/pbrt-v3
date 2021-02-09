set -e

SCENE=$1
RUN=./../../scripts/start_eval.sh

if [ ! -d $SCENE ]
then
    mkdir $SCENE
fi

if [ ! -d "output" ]
then
    mkdir output
fi

cmake -S ../../.. -B ../../../build -DCOUNT_STATS=True -DREL_KEYS=True -DBF_SIZE=64 -DCHUNK_SIZE=64
make -C ../../../build -j

$RUN $SCENE "bvh"
$RUN $SCENE "bvh-bfs"
$RUN $SCENE "octree"
$RUN $SCENE "octree-bfs"
$RUN $SCENE "kdtree"