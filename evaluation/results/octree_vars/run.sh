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
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.3 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.3 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.3 float:prm_thresh=0.9
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.6 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.6 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.6 float:prm_thresh=0.9
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.9 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.9 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=8 float:vol_thresh=0.9 float:prm_thresh=0.9

$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.3 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.3 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.3 float:prm_thresh=0.9
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.6 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.6 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.6 float:prm_thresh=0.9
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.9 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.9 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=16 float:vol_thresh=0.9 float:prm_thresh=0.9

$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.3 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.3 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.3 float:prm_thresh=0.9
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.6 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.6 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.6 float:prm_thresh=0.9
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.9 float:prm_thresh=0.3
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.9 float:prm_thresh=0.6
$RUN "octree" integer:maxdepth=32 float:vol_thresh=0.9 float:prm_thresh=0.9