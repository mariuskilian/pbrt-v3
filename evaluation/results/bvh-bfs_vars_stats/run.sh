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

if ! [[ $* == *--skip-render* ]]; then
    COUNT_STATS="-DCOUNT_STATS=True"

    BUILD1="$BUILD-1"
    cmake -S $SOURCE -B $BUILD1 $COUNT_STATS "-DREL_KEYS=True"
    make -C $BUILD1 -j
    $RUN $BUILD1 "string:relkeys=true"

    BUILD2="$BUILD-2"
    cmake -S $SOURCE -B $BUILD2 $COUNT_STATS "-DREL_KEYS=False"
    make -C $BUILD2 -j
    $RUN $BUILD2 "string:relkeys=false"
fi

python3 $PYSCRIPTS/plot_data.py $SCENE stat:primitive
python3 $PYSCRIPTS/plot_data.py $SCENE stat:node
python3 $PYSCRIPTS/plot_data.py $SCENE stat:leafnode

python3 $PYSCRIPTS/plot_data.py $SCENE dist:primitive
python3 $PYSCRIPTS/plot_data.py $SCENE dist:node
python3 $PYSCRIPTS/plot_data.py $SCENE dist:leafnode