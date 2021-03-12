set -e

SCENE=$1

cd ../bvh-bfs_vars_stats
./run.sh $SCENE

cd ../bvh-bfs_vars_time
./run.sh $SCENE