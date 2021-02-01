# Parameter names
SCENE="crown"
ACCEL="bvh-chunk-bfs"
# Constant 
# Create base .pbrt file which sets Scene and Accelerator
python3 python-scripts/create_base_pbrt.py $SCENE $ACCEL
./../build/RelWithDebInfo/pbrt ../scenes/
# Remove previously created base .pbrt file
python3 python-scripts/remove_base_pbrt.py $SCENE $ACCEL