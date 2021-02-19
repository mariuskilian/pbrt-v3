set -e

SCENE=$1

INTGR="path"
RUN="./../../scripts/start_eval.sh $INTGR $SCENE octree integer:maxprims=1"
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

$RUN float:multthresh=1.5
$RUN float:multthresh=2
$RUN float:multthresh=2.5
$RUN float:multthresh=3
$RUN float:multthresh=3.5
$RUN float:multthresh=4
$RUN float:multthresh=4.5
$RUN float:multthresh=5
$RUN float:multthresh=5.5