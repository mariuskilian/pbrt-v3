set -e

SCENES=$1

INTGR="metric=primitives"
SOURCE=../../..
PYSCRIPTS=../../scripts/python-scripts
BUILD=$SOURCE/build/Evaluation
RUN="./../../scripts/start_eval.sh $INTGR"

NPIXELSAMPLES="-n=1"


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

        cmake -S $SOURCE -B $BUILD
        make -C $BUILD -j

        $RUN $SCENE $BUILD embree $NPIXELSAMPLES
        $RUN $SCENE $BUILD bvh $NPIXELSAMPLES
        $RUN $SCENE $BUILD bvh-bfs $NPIXELSAMPLES
        $RUN $SCENE $BUILD octree $NPIXELSAMPLES
    done
fi

SCENE_LIST=($(echo $SCENES | tr "," "\n"))
for i in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
    SCENE=${SCENE_LIST[$i]}
    python3 $PYSCRIPTS/normalize_metrics.py $SCENE $INTGR
    python3 $PYSCRIPTS/normalize_metrics.py $SCENE $INTGR --notitle
    python3 $PYSCRIPTS/normalize_metrics.py $SCENE $INTGR --shorttitle
done
