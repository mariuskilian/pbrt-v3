set -e

SCENE=$1

INTGR="path"
ACCEL="bvh-bfs"
RUN="./../../scripts/start_eval.sh $INTGR $SCENE"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation

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
    COUNT_STATS="-DCOUNT_STATS=False" # Because we want to compare time

    BUILD1="$BUILD-1"
    cmake -S $SOURCE -B $BUILD1 $COUNT_STATS "-DREL_KEYS=True"
    make -C $BUILD1 -j
    $RUN $BUILD1 $ACCEL "string:relkeys=true" $NPIXELSAMPLES

    BUILD2="$BUILD-2"
    cmake -S $SOURCE -B $BUILD2 $COUNT_STATS "-DREL_KEYS=False"
    make -C $BUILD2 -j
    $RUN $BUILD2 $ACCEL "string:relkeys=false" $NPIXELSAMPLES
fi

python3 $PYSCRIPTS/plot_data.py $SCENE prof --plot