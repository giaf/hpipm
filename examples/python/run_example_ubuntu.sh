# BLASFEO root folder path or installation path
BLASFEO_PATH=/opt/blasfeo
# HPIPM root folder path or installation path
HPIPM_PATH=/opt/hpipm

#export path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BLASFEO_PATH/lib:$HPIPM_PATH/lib

# run example
python3 getting_started.py
