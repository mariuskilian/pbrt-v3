# Relative paths. This file is called from ../results/[case]/run.sh
#   so all relative paths must start from there
SCRIPTS_PATH=../../scripts
SRC_PATH=../../..

# Parameter names
SCENE=$1
ACCEL=$2

# Find any extra parameters
EXTRA_ARGS="# Custom extra arguments;"
i=1;
for arg in "$@" 
do
    if [ $i -gt 2 ]
    then
        EXTRA_ARGS+="$arg;"
    fi
    i=$((i + 1));
done

# Create base .pbrt file which sets Scene and Accelerator (and extra arguments)
python3 $SCRIPTS_PATH/python-scripts/create_base_pbrt.py $SCENE $ACCEL "$EXTRA_ARGS"
# Run pbrt and log the results
./$SRC_PATH/build/pbrt --outfile output/$SCENE-$ACCEL.exr $SRC_PATH/scenes/$SCENE/eval_base.pbrt | tee $SCENE/$ACCEL.log
# Remove previously created base .pbrt file
python3 $SCRIPTS_PATH/python-scripts/remove_base_pbrt.py $SCENE