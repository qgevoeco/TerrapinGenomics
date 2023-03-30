# INSTALLATION
cd /scratch/sbw0033/TerrapinGenomics/Data/Demography/GADMA #Get into the home directory

module load python/anaconda/3.9.13

python3 -m virtualenv gadma_env
source gadma_env/bin/activate

pip install gadma #Install GADMA
pip install git+https://bitbucket.org/simongravel/moments.git #Installs moments as well

# Test whether it works
gadma --test

#Leave the environment
deactivate
