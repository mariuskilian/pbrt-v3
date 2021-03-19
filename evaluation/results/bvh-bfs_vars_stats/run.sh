set -e

SCENES=$1

INTGR="path"
ACCEL="bvh-bfs"
RUN="./../../scripts/start_eval.sh $INTGR"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation

for arg in "$@"; do
    if [[ $arg == "-n="* ]]; then
        NPIXELSAMPLES=$arg
    fi
done

if [ ! -d "output" ]
then
    mkdir output
fi

if ! [[ $* == *--skip-render* ]]; then
    COUNT_STATS="-DCOUNT_STATS=True"
    BUILD1="$BUILD-1"
    BUILD2="$BUILD-2"

    SCENE_LIST=($(echo $SCENES | tr "," "\n"))
    for i in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
        SCENE=${SCENE_LIST[$i]}
        
        if [ ! -d $SCENE ]; then
            mkdir $SCENE
        fi

        cmake -S $SOURCE -B $BUILD1 $COUNT_STATS "-DREL_KEYS=True"
        make -C $BUILD1 -j
        $RUN $SCENE $BUILD1 $ACCEL "string:relkeys=true" $NPIXELSAMPLES

        cmake -S $SOURCE -B $BUILD2 $COUNT_STATS "-DREL_KEYS=False"
        make -C $BUILD2 -j
        $RUN $SCENE $BUILD2 $ACCEL "string:relkeys=false" $NPIXELSAMPLES
    done
fi

python3 $PYSCRIPTS/plot_data.py $SCENE dist:primitive --plot