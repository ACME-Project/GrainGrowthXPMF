# Python script to automagically make directory trees and appropriate 
# sbatch files for a set of experiments

# Use structure is as follows: python make_sbatch /target/directory/root/ trunk* 

import glob, os, sys, string, shutil

# place should start and end with the  "/" character, and be the absolute directory. 
# name is the experiment name and name.log is the log file name
def makeBatch(place, job_name, nodes, num_steps, inc, exe_name="q_sparse.out", start=False, start_file="file.000000.dat", dependency=None):
	os.chdir(place)
	sb = "#SBATCH "
	name = os.getlogin()[:4]
	
	batch_name = 'sbatch.sh'
	if(start):
		batch_name = 'start.sh'
	bfile = open(batch_name,'w')
	
	bfile.write("#!/bin/sh" + '\n' + '\n')
	if(start):
		bfile.write(sb + "--job-name=" + job_name + 'S \n')
	else:
		bfile.write(sb + "--job-name=" + job_name + 'R \n')
	bfile.write(sb + "--account=" + name + '\n')
	
	if(start):
		bfile.write(sb + "--partition small" + '\n')
		bfile.write(sb + "--time 1440" + '\n')
		bfile.write(sb + "--nodes 1" + '\n')
		bfile.write(sb + "--ntasks 1" + '\n')
		bfile.write(sb + "--overcommit" + '\n' + '\n')
	else:
		if(nodes <= 64):
			bfile.write(sb + "--partition small" + '\n')
			bfile.write(sb + "--time 1440" + '\n')
		elif(nodes <= 512):
			bfile.write(sb + "--partition medium" + '\n')
			bfile.write(sb + "--time 720" + '\n')
		elif(nodes <= 2048):
			bfile.write(sb + "--partition large" + '\n')
			bfile.write(sb + "--time 360" + '\n')
		elif(nodes <= 4096):
			bfile.write(sb + "--partition verylarge" + '\n')
			bfile.write(sb + "--time 360" + '\n')

		bfile.write(sb + "--nodes " + str(nodes) + '\n')
		bfile.write(sb + "--ntasks " + str(nodes*16) + '\n')
		bfile.write(sb + "--overcommit" + '\n' + '\n')
		
	if dependency is not None:
		bfile.write(sb + "--dependency=afterok:" + str(dependency) + '\n' + '\n')

	bfile.write(sb + "-o " + place + job_name + ".log" + '\n')
	bfile.write(sb + "-D " + place + '\n' + '\n')
	
	if(start):
		bfile.write("srun --unbuffered " + place + exe_name + " --example 2 " + place + start_file)
	else:
		bfile.write("srun --unbuffered " + place + exe_name +" "+ place + start_file +" "+ str(num_steps) + " " + str(inc))

	bfile.close()

name = os.getlogin()
group = name[:4]
root_place = "/gpfs/u/home/" + group + "/" + name + "/scratch/" #all experiments are performed in the scratch directory

exp_name = str(sys.argv[1]) #The name of the experiment set, as determined by the directory targeted
os.chdir(root_place + exp_name)
directory_name_trunk = str(sys.argv[2]) #pattern name for trunk directory of simulation microstructure data. 
if(directory_name_trunk[-1]!='*'):
	directory_name_trunk += '*'
exp_set = glob.glob(directory_name_trunk)
num_steps = 100000
inc = 250
nodes = 64

exe_name="q_sparse.out"

exp_dir = root_place + exp_name	
		
for ex in exp_set:
	makeBatch((exp_dir + ex + "/" + ), ex, nodes, num_steps, inc, exe_name, start=True)
	makeBatch((exp_dir + ex + "/" + ), ex, nodes, num_steps, inc, exe_name, start=False) 
		




#!/bin/sh

#SBATCH --job-name=TwoXSmall
#SBATCH --account=GGST
#SBATCH --mail-user=cristd2@rpi.edu
#SBATCH --mail-type=END

#SBATCH --partition medium
#SBATCH --time 720
#SBATCH --nodes 512
#SBATCH --ntasks 8192
#SBATCH --overcommit

#SBATCH -o /gpfs/u/home/GGST/GGSTcrdv/scratch/taylor/xFactor/smallt/twoX/xTwo.log
#SBATCH -D /gpfs/u/home/GGST/GGSTcrdv/scratch/taylor/xFactor/smallt/twoX

#srun --unbuffered /gpfs/u/home/GGST/GGSTcrdv/scratch/taylor/xFactor/smallt/twoX/q_GGt.out --example 2 /gpfs/u/home/GGST/GGSTcrdv/scratch/taylor/xFactor/smallt/twoX/file.000000.dat
#srun --unbuffered /gpfs/u/home/GGST/GGSTcrdv/scratch/taylor/xFactor/smallt/twoX/q_GGt.out /gpfs/u/home/GGST/GGSTcrdv/scratch/taylor/xFactor/smallt/twoX/file.200000.dat 400000 5000

