#!/bin/bash -e
# Badri Adhikari, 5-21-2017
# The main script for making DNCON2 contact prediction after feature files are generated


ROOT=$(dirname $0)

dir_scripts=$ROOT
dir_config=$ROOT/../model-config-n-weights

fileX=$1
fileRR=$2
file_stg2_tmp=$3

if [ ! -e "$fileX" ]; then
	echo "ERROR! Input feature file does not exist!"
	echo "Usage: $0 <input-feature-file> <output-RR-file-name> <intermediate-feature-file-name>"
	exit
fi

if [ -z $fileRR ]; then
	echo "ERROR! Output RR file name not supplied!"
	echo "Usage: $0 <input-feature-file> <output-RR-file-name> <intermediate-feature-file-name>"
	exit
fi

if [ -z $file_stg2_tmp ]; then
	echo "ERROR! intermediate-feature-file-name not supplied!"
	echo "Usage: $0 <input-feature-file> <output-RR-file-name> <intermediate-feature-file-name>"
	exit
fi

# Predict using input feature file and stage1 models
if [ -f $file_stg2_tmp ]; then
	echo "Stage2 file already present.. using it.."
else
	echo "# Predict stage1 coevolution-based features (and prepare stage2 feature file).."
	for dist_thres in '60' '75' '80' '85' '10' ; do
		echo ""
		echo "Running prediction using coevo-${dist_thres}A .. "
		python $dir_scripts/cnn-predict-and-append-to-X.py $dir_config "stage1-${dist_thres}A.hdf5" "${dist_thres}A" $fileX $file_stg2_tmp
	done 
fi

echo ""
echo "All features in final X:"
grep "#" $file_stg2_tmp

# Predict stage2 with new/updated feature file
# Need to write output to a .tmp file and later rename the file because, if the file
# is detected and downloaded while it is being written, incomplete results will be transferred.
echo ""
echo "Predict stage2.."
python $dir_scripts/cnn-predict-stage2.py $dir_config $file_stg2_tmp $fileRR.tmp
mv $fileRR.tmp $fileRR
