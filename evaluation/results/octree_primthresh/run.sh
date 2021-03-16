set -e

SCENE=$1
MIN=$2
MAX=$3
STEP=$4

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

if ! [[ $* == *--skip-render* ]]; then
    COUNT_STATS="-DCOUNT_STATS=False" # Because we want to compare time

    cmake -S $SOURCE -B $BUILD $COUNT_STATS
    make -C $BUILD -j

    for i in $(seq $MIN $STEP $MAX); do
        $RUN integer:maxprims=$i
    done
fi

python3 $PYSCRIPTS/normalize_filename_num_digits.py $SCENE
python3 $PYSCRIPTS/plot_data.py $SCENE prof
python3 $PYSCRIPTS/plot_data.py $SCENE mem:topology