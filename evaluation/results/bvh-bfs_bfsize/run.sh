set -e

SCENES=$1

INTGR="path"
RUN="./../../scripts/start_eval.sh $INTGR"
SOURCE=../../..
BUILD=$SOURCE/build/Evaluation
PYSCRIPTS=../../scripts/python-scripts
BUILD1="$BUILD-1"
BUILD2="$BUILD-2"
BUILD3="$BUILD-3"
BUILD4="$BUILD-4"
BUILD5="$BUILD-5"
BUILD6="$BUILD-6"

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

    SCENE_LIST=($(echo $SCENES | tr "," "\n"))
    for i in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
        SCENE=${SCENE_LIST[$i]}

        if [ ! -d $SCENE ]
        then
            mkdir $SCENE
        fi
        COUNT_STATS="COUNT_STATS=False"

        BFSIZE="BF_SIZE=8"
        cmake -S $SOURCE -B $BUILD1 -D$COUNT_STATS -D$BFSIZE
        make -C $BUILD1 -j
        FILENAME_POSTFIX="integer:$BFSIZE"
        $RUN $SCENE $BUILD1 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        BFSIZE="BF_SIZE=16"
        cmake -S $SOURCE -B $BUILD2 -D$COUNT_STATS -D$BFSIZE
        make -C $BUILD2 -j
        FILENAME_POSTFIX="integer:$BFSIZE"
        $RUN $SCENE $BUILD2 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        BFSIZE="BF_SIZE=32"
        cmake -S $SOURCE -B $BUILD3 -D$COUNT_STATS -D$BFSIZE
        make -C $BUILD3 -j
        FILENAME_POSTFIX="integer:$BFSIZE"
        $RUN $SCENE $BUILD3 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        BFSIZE="BF_SIZE=64"
        cmake -S $SOURCE -B $BUILD4 -D$COUNT_STATS -D$BFSIZE
        make -C $BUILD4 -j
        FILENAME_POSTFIX="integer:$BFSIZE"
        $RUN $SCENE $BUILD4 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        # Scripts to have per scene plots should go here
    done
fi

python3 $PYSCRIPTS/normalize_filename_num_digits.py $SCENES
#Scripts to have plot show all scenes at once should go here
PLOT_DATA="python3 $PYSCRIPTS/plot_data.py $SCENES"
$PLOT_DATA prof --plot
$PLOT_DATA mem:topology --plot