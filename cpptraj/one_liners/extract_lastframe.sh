# How to Run?: bash extract_lastframe.sh
cpptraj -p $1 -y $2 -ya 'lastframe' -x $3
quit
