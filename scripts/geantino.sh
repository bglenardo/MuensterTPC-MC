# Shell script to run the MS-TPC Monte Carlo using a Decay0 macro

# Your bashrc command aliases to load root and geant4
load_geant4
load_root6

echo 'Loaded ROOT6 and Geant4.10'

# Specify macro/file-name and directory to save files in
filename='geantino'
filedir='/Users/christianwittweg/Desktop/0vECEC/'

echo 'Filename is '$filename

cd $G4WORKDIR/MuensterTPC-MC

#Run the MC, change verbosity for more info with '-v 1'
runstring='./MuensterTPC-MC -f ./macros/src_'$filename'.mac -o '$filename'.root -n 1 -v 1'

$runstring

# Move files and go to filedir
mv *.root $filedir

echo 'Moving ROOT-file to the working directory'

#cd $filedir
