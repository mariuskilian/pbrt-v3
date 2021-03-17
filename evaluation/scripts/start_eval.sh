# Relative paths. This file is called from ../results/[case]/run.sh
#   so all relative paths must start from there
SCRIPTS_PATH=../../scripts
SRC_PATH=../../..

# Parameter names
INTGR=$1
SCENE=$2
BUILD=$3
ACCEL=$4

# Find any extra parameters
EXTRA_ARGS="#CustomExtraArguments;"
FILENAME_POSTFIX=""
i=1;
for arg in "$@" 
do
    if [ $i -gt 3 ]
    then
        if [[ $arg == "-n="* ]]; then
            NPIXELSAMPLES="--npixelsamples ${arg:3}"
        else
            EXTRA_ARGS+="$arg;"
            EXTRA_ARGS_LIST=($(echo $arg | tr ":" "\n"))
            FILENAME_POSTFIX+="_${EXTRA_ARGS_LIST[1]}"
        fi
    fi
    i=$((i + 1));
done

# Create base .pbrt file which sets Scene and Accelerator (and extra arguments)
python3 $SCRIPTS_PATH/python-scripts/create_base_pbrt.py $INTGR $SCENE $ACCEL "$EXTRA_ARGS" $NPIXELSAMPLES
# Run pbrt and log the results
./$BUILD/pbrt --outfile output/$SCENE-$INTGR-$ACCEL$FILENAME_POSTFIX.exr --nthreads 8 $SRC_PATH/scenes/$SCENE/eval_base.pbrt | tee $SCENE/$INTGR-$ACCEL$FILENAME_POSTFIX.log
# Remove previously created base .pbrt file
# python3 $SCRIPTS_PATH/python-scripts/remove_base_pbrt.py $SCENE