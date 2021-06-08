#PBS -lselect=1:ncpus=1:mem=64gb
#PBS -lwalltime=72:00:00
#PBS -J 0-49

# subjob's index within the array

## All subjobs run independently of one another


# Load modules for any applications

module load anaconda3/personal

# Change to the submission directory

cd $PBS_O_WORKDIR

# Run program, passing the index of this subjob within the array

python simulate_cmd.py 2021-03-18 $PBS_ARRAY_INDEX 50


