# How to Run?: bash merge_trajectories.sh $1 $2 $3 $4 $5 $6
# $1: prmtop, $2: trajectory files
cpptraj -p $1 -y $2 -y $3 -y $4 -y $5 -x $6
