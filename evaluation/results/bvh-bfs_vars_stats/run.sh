set -e

SCENE=$1

INTGR="path"
ACCEL="bvh-bfs"
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

COUNT_STATS="-DCOUNT_STATS=True"

cmake -S $SOURCE -B $BUILD $COUNT_STATS "-DREL_KEYS=True"
make -C $BUILD -j
$RUN "string:relkeys=true"

cmake -S $SOURCE -B $BUILD $COUNT_STATS "-DREL_KEYS=False"
make -C $BUILD -j
$RUN "string:relkeys=false"

python3 $PYSCRIPTS/plot_data.py $SCENE stat:primitive
python3 $PYSCRIPTS/plot_data.py $SCENE stat:node
python3 $PYSCRIPTS/plot_data.py $SCENE stat:leafnode

python3 $PYSCRIPTS/plot_data.py $SCENE dist:primitive
python3 $PYSCRIPTS/plot_data.py $SCENE dist:node
python3 $PYSCRIPTS/plot_data.py $SCENE dist:leafnode