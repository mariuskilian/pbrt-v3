# Relative paths. This file is called from ../results/[case]/run.sh
#   so all relative paths must start from there
SCRIPTS_PATH=../../scripts
SRC_PATH=../../..

# Parameter names
INTGR=$1
SCENE=$2
ACCEL=$3

# Find any extra parameters
EXTRA_ARGS="# Custom extra arguments;"
FILENAME_POSTFIX=""
i=1;
for arg in "$@" 
do
    if [ $i -gt 3 ]
    then
        EXTRA_ARGS+="$arg;"
        EXTRA_ARGS_LIST=($(echo $arg | tr ":" "\n"))
        FILENAME_POSTFIX+="_${EXTRA_ARGS_LIST[1]}"
    fi
    i=$((i + 1));
done

# Create base .pbrt file which sets Scene and Accelerator (and extra arguments)
python3 $SCRIPTS_PATH/python-scripts/create_base_pbrt.py $INTGR $SCENE $ACCEL "$EXTRA_ARGS"
# Run pbrt and log the results
./$SRC_PATH/build/Evaluation/pbrt --outfile output/$SCENE-$INTGR-$ACCEL$FILENAME_POSTFIX.exr $SRC_PATH/scenes/$SCENE/eval_base.pbrt | tee $SCENE/$INTGR-$ACCEL$FILENAME_POSTFIX.log
# Remove previously created base .pbrt file
# python3 $SCRIPTS_PATH/python-scripts/remove_base_pbrt.py $SCENE