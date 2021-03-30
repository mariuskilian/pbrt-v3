set -e

cd ../all_stats
./run.sh $1 $2 $3 $4 $5 $6

cd ../all_time
./run.sh $1 $2 $3 $4 $5 $6