#!/usr/bin/env bash
#SBATCH --time=1-12:00:00
#SBATCH --partition=genx
#SBATCH --mem 64g
#SBATCH -c 24


# Load bash config.
source "${HOME}"'/.bashrc'
if [ $? != 0 ]; then
    echo 'Failed sourcing '"${HOME}"'/.bashrc'
    exit 1
fi

# Start up anaconda.
conda activate 'charm_altsp'
if [ $? != 0 ]; then
    echo 'Failed to activate conda environment.'
    exit 1
fi

# Change directories.
SRC_DIR='.'
cd "${SRC_DIR}"
if [ $? != 0 ]; then
    echo 'Failed changing directories to '"${SRC_DIR}"
    exit 1
fi

# Run train script.
/usr/bin/time -v python run_jem_tests.py --data data-V7 --event-type $1
echo $?

# Deactivate conda.
conda deactivate

