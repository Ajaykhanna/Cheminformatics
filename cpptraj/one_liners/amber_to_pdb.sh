# How to Run?: bash amber_to_pdb.sh amber.prmtop amber.rst7 amber.pdb
cpptraj -p $1 -y $2 -x $3
quit
