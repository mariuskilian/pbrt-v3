set -e

cd ../vis_leafnodes
./run.sh $1 $2 $3 $4 $5 $6

cd ../vis_nodes
./run.sh $1 $2 $3 $4 $5 $6

cd ../vis_primitives
./run.sh $1 $2 $3 $4 $5 $6