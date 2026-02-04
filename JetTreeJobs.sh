#!/bin/sh

echo "Running Jobs"
echo "Host: $(hostname)"
echo "Start time: $(date)"

# Arguments passed from Condor
INPUTFILE=$1
OUTPUTFILE=$2
RVALS=$3
REMOVEELECTRONS=$4
NHITCUT=$5

# Enter the eic-shell environment
/gpfs/mnt/gpfs02/eic/ddesouza/epic-dd4hep/eic-shell << EOF

# Load EPIC environment
source /opt/detector/epic-main/bin/thisepic.sh

echo "Environment loaded"
echo "Running ROOT macro..."

# Run your ROOT macro (modify macro name as needed)
root -l -b -q 'src/JetTreesRecluster.C("$INPUTFILE","$OUTPUTFILE",std::vector<float>{$RVALS},$REMOVEELECTRONS,$NHITCUT)'

echo "ROOT job finished"

exit
EOF

echo "End time: $(date)"

