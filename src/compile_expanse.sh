# module load cpu/0.17.3b gcc/10.2.0/npcyll4 openmpi/4.1.1/ygduf2r
# make yes-AGENT
# make g++_openmpi -j 4
# cp ./lmp_g++_openmpi ~/bin/lmp_body_latest  

module load cpu/0.17.3b
module load gcc/10.2.0/npcyll4
module load mvapich2/2.3.7/iyjtn3x
make yes-AGENT
make -j4 g++_mpich
cp lmp_g++_mpich ~/bin/lmp_body_latest

# module load cpu/0.15.4
# module load intel/19.1.1.217
# module load mvapich2/2.3.6
# make -j 4 intel_cpu_mpich
# cp lmp_intel_cpu_mpich ~/bin/lmp_body_latest


