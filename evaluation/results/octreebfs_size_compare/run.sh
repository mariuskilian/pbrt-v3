set -e

INTGR="path"

RUN="./../../scripts/start_eval.sh $INTGR"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation
NPIXELSAMPLES=-n=1

SCENE1=crown
SCENE2=measure-one

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

    cmake -S $SOURCE -B $BUILD
    make -C $BUILD1 -j

    $RUN $SCENE1 $BUILD octree $NPIXELSAMPLES
    $RUN $SCENE2 $BUILD octree $NPIXELSAMPLES
    $RUN $SCENE1 $BUILD octree-bfs $NPIXELSAMPLES
    $RUN $SCENE2 $BUILD octree-bfs $NPIXELSAMPLES
fi

python3 $PYSCRIPTS/plot_data.py $SCENE mem --plot
python3 $PYSCRIPTS/plot_data.py $SCENE mem:topology --plot