set -e

SCENE=$1

INTGR="path"
ACCEL="octree"
RUN="./../../scripts/start_eval.sh $INTGR $SCENE $ACCEL"
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

COUNT_STATS="-DCOUNT_STATS=False" # Because we want to compare time

cmake -S $SOURCE -B $BUILD $COUNT_STATS
make -C $BUILD -j

$RUN float:multthresh=1.5
$RUN float:multthresh=2.0
$RUN float:multthresh=2.5
$RUN float:multthresh=3.0
$RUN float:multthresh=3.5
$RUN float:multthresh=4.0
$RUN float:multthresh=4.5
$RUN float:multthresh=5.0
$RUN float:multthresh=5.5

python3 $PYSCRIPTS/plot_data.py $SCENE $ACCEL prof