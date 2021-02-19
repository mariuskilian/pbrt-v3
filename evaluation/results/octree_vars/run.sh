set -e

SCENE=$1

INTGR="path"
RUN="./../../scripts/start_eval.sh $INTGR $SCENE"
SOURCE=../../..
BUILD=$SOURCE/build/Evaluation

if [ ! -d $SCENE ]
then
    mkdir $SCENE
fi

if [ ! -d "output" ]
then
    mkdir output
fi

COUNT_STATS="-DCOUNT_STATS=False" # Because we want to compare time

cmake -S $SOURCE -B $BUILD $COUNT_STATS
make -C $BUILD -j

$RUN "octree" float:multthresh=1
$RUN "octree" float:multthresh=1.5
$RUN "octree" float:multthresh=2
$RUN "octree" float:multthresh=2.5
$RUN "octree" float:multthresh=3
$RUN "octree" float:multthresh=3.5
$RUN "octree" float:multthresh=4
$RUN "octree" float:multthresh=4.5
$RUN "octree" float:multthresh=5
$RUN "octree" float:multthresh=5.5
$RUN "octree" float:multthresh=6
$RUN "octree" float:multthresh=6.5
$RUN "octree" float:multthresh=7
$RUN "octree" float:multthresh=7.5
$RUN "octree" float:multthresh=8