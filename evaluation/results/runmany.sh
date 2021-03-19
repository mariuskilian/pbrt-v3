set -e

TESTS=$1
SCENES=$2

TEST_LIST=($(echo $TESTS | tr ":" "\n"))
SCENE_LIST=($(echo $SCENES | tr ":" "\n"))
for i in $(seq 0 1 $((${#TEST_LIST[@]} - 1))); do
    cd ${TEST_LIST[$i]}
    for j in $(seq 0 1 $((${#SCENE_LIST[@]} - 1))); do
        SCENE=${SCENE_LIST[$j]}
        echo $PWD
        ./run.sh $SCENE $3 $4 $5 $6 $7 $8 $9
    done
    cd ..
done