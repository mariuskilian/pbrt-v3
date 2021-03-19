set -e

SCENES=$1

INTGR="path"
RUN="./../../scripts/start_eval.sh $INTGR"
SOURCE=../../..
BUILD=$SOURCE/build/Evaluation
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

        CHUNKSIZE="CHUNK_SIZE=16"
        cmake -S $SOURCE -B $BUILD1 -D$COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
        make -C $BUILD1 -j
        FILENAME_POSTFIX="integer:$CHUNKSIZE"
        $RUN $SCENE $BUILD1 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        CHUNKSIZE="CHUNK_SIZE=32"
        cmake -S $SOURCE -B $BUILD1 -D$COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
        make -C $BUILD1 -j
        FILENAME_POSTFIX="integer:$CHUNKSIZE"
        $RUN $SCENE $BUILD2 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        CHUNKSIZE="CHUNK_SIZE=64"
        cmake -S $SOURCE -B $BUILD1 -D$COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
        make -C $BUILD2 -j
        FILENAME_POSTFIX="integer:$CHUNKSIZE"
        $RUN $SCENE $BUILD3 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        CHUNKSIZE="CHUNK_SIZE=128"
        cmake -S $SOURCE -B $BUILD2 $COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
        make -C $BUILD3 -j
        FILENAME_POSTFIX="integer:$CHUNKSIZE"
        $RUN $SCENE $BUILD4 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

        # Scripts to have per scene plots should go here
    done
fi

#Scripts to have plot show all scenes at once should go here