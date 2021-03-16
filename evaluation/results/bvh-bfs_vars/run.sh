set -e

SCENE=$1

cd ../bvh-bfs_vars_stats
./run.sh $SCENE $2 $3 $4 $5 $6

cd ../bvh-bfs_vars_time
./run.sh $SCENE $2 $3 $4 $5 $6