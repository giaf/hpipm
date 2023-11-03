#! /usr/bin/bash

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]
then
	echo "Script is being sourced"
else
	echo "ERROR: Script is a subshell"
	echo "To affect your current shell enviroment source this script with:"
	echo "source env.sh"
	exit
fi

# check that this file is run
export ENV_RUN=true

# if hpipm folder not specified assume parent of this folder
HPIPM_MAIN_FOLDER=${HPIPM_MAIN_FOLDER:-"$(pwd)/../.."}
export HPIPM_MAIN_FOLDER
echo
echo "HPIPM_MAIN_FOLDER=$HPIPM_MAIN_FOLDER"

# if blasfeo folder not specified assume alongside the parent of this folder
BLASFEO_MAIN_FOLDER=${BLASFEO_MAIN_FOLDER:-"$(pwd)/../../../blasfeo"}
export BLASFEO_MAIN_FOLDER
echo
echo "BLASFEO_MAIN_FOLDER=$BLASFEO_MAIN_FOLDER"

# automatically export LD_LIBRARY_PATH or DYLD_LIBRARY_PATH in linux or macos
UNAME_S=$(uname -s)
echo $UNAME_S
if [ "$UNAME_S" = "Linux" ]
then
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HPIPM_MAIN_FOLDER/lib:$BLASFEO_MAIN_FOLDER/lib
	echo
	echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
else
	if [ "$UNAME_S" = "Darwin" ]
	then
		export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$HPIPM_MAIN_FOLDER/lib:$BLASFEO_MAIN_FOLDER/lib
		echo
		echo "DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH"
	else
		echo "LD_LIBRARY_PATH not set automatically for this OS. Try to it by hand in the correct way for your OS."
	fi
fi
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HPIPM_MAIN_FOLDER/lib:$BLASFEO_MAIN_FOLDER/lib
#echo
#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

