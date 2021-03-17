set -e

SCENE=$1

INTGR="path"
RUN="./../../scripts/start_eval.sh $INTGR $SCENE"
SOURCE=../../..
BUILD=$SOURCE/build/Evaluation
BUILD1="$BUILD-1"
BUILD2="$BUILD-2"
BUILD2="$BUILD-3"
BUILD2="$BUILD-4"

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
    COUNT_STATS="COUNT_STATS=True"

    BFSIZE="BF_SIZE=32"
    CHUNKSIZE="CHUNK_SIZE=64"
    cmake -S $SOURCE -B $BUILD1 -D$COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
    make -C $BUILD1 -j
    FILENAME_POSTFIX="integer:$BFSIZE integer:$CHUNKSIZE"
    $RUN $BUILD1 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES
    $RUN $BUILD1 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

    BFSIZE="BF_SIZE=32"
    CHUNKSIZE="CHUNK_SIZE=128"
    cmake -S $SOURCE -B $BUILD2 $COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
    make -C $BUILD2 -j
    FILENAME_POSTFIX="integer:$BFSIZE integer:$CHUNKSIZE"
    $RUN $BUILD2 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES
    $RUN $BUILD2 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

    BFSIZE="BF_SIZE=64"
    CHUNKSIZE="CHUNK_SIZE=64"
    cmake -S $SOURCE -B $BUILD3 $COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
    make -C $BUILD3 -j
    FILENAME_POSTFIX="integer:$BFSIZE integer:$CHUNKSIZE"
    $RUN $BUILD3 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES
    $RUN $BUILD3 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES

    BFSIZE="BF_SIZE=64"
    CHUNKSIZE="CHUNK_SIZE=128"
    cmake -S $SOURCE -B $BUILD4 $COUNT_STATS -D$BFSIZE -D$CHUNKSIZE
    make -C $BUILD4 -j
    FILENAME_POSTFIX="integer:$BFSIZE integer:$CHUNKSIZE"
    $RUN $BUILD4 "octree-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES
    $RUN $BUILD4 "bvh-bfs" $FILENAME_POSTFIX $NPIXELSAMPLES
fi