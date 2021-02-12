To run start_eval use the following arguments:

start_eval <integrator> <scene> <accelerator> (<extra_args>)

scene = {crown, killeroo, ...}
accelerator = {embree, bvh, bvh-bfs, octree, octree-bfs}    
integrator = {path, metric=<metric_name>} with
    metric_name as {primitives, nodes, leafnodes, time}
extra_args need to have form: <type>:<name>=<value> with
    type as {integer, float, string, ...},
    name as the name given in code,
    value as the value the variable should have 

Examples:
./start_eval crown octree-bfs metric=primitives string:maxprims=15 float:vol_thresh=0.8
./start_eval killeroo bvh-bfs path