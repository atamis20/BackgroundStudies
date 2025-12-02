# BackgroundStudies
A guide for running background simulations and analyzing them using irt-2.1b branch of EICrecon as a base

Follow all instructions for getting 2.1b branch set up, following instructions here, using 2.1b instead of 2.1a for branch:
https://github.com/eic/EICrecon/tree/irt-2.1a/irt-sandbox
copied here:

# Change to your local working directory;
cd your-working-directory

# Install eic-shell script and a docker container image known to work;
curl -L https://github.com/eic/eic-shell/raw/main/install.sh     | bash -s --  --version 25.07.0-stable --no-cvmfs

# Run 'eic-shell';
./eic-shell

# Use git branch irt-2.1b for all repositories;
export branch="irt-2.1b"

# Download EICrecon;
git clone -b ${branch} https://github.com/eic/EICrecon.git

# Run a sandbox environment script;
. EICrecon/irt-sandbox/environ.sh

# Install EDM4eic;
git clone -b ${branch} https://github.com/eic/EDM4eic.git

cmake -S EDM4eic -B EDM4eic/build -DBUILD_DATA_MODEL=ON -DCMAKE_INSTALL_PREFIX=$EIC_SHELL_PREFIX -Wno-dev

cmake --build EDM4eic/build -j8

cmake --install EDM4eic/build

# Install IRT;
git clone -b ${branch} https://github.com/eic/irt.git

cmake -S irt -B irt/build -DCMAKE_BUILD_TYPE=Debug -DDELPHES=OFF -DEVALUATION=OFF -DCMAKE_INSTALL_PREFIX=$EIC_SHELL_PREFIX -Wno-dev

cmake --build irt/build -j8

cmake --install irt/build

# Install epic;
git clone -b  ${branch} https://github.com/eic/epic.git
cmake -S epic -B epic/build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$EIC_SHELL_PREFIX -DWITH_IRT=YES -Wno-dev

cmake --build epic/build -j8

cmake --install epic/build

# Install EICrecon;
cmake -S EICrecon -B EICrecon/build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_FIND_DEBUG_MODE=OFF -DCMAKE_INSTALL_PREFIX=$EIC_SHELL_PREFIX -Wno-dev

cmake --build EICrecon/build -j8

cmake --install EICrecon/build



# Analysis Code
Download contents of this directory into the irt-sandbox directory

#Enter eic-shell and source environ

cd your-working-directory

./eic-shell

. EICrecon/irt-sandbox/environ.sh

## Source files for beam gas are stored in jlab, similar files also exist for 18GeV and 100GeV, just change all instances of 10 or 275.  Can get information about rate on the wiki

https://wiki.bnl.gov/EPIC/index.php?title=Electron_Beam_Gas

https://wiki.bnl.gov/EPIC/index.php?title=Hadron_Beam_Gas


# Run Simuilation for electron beam+gas
npsim --runType run --compactFile $EIC_SHELL_PREFIX/share/epic/epic_tracking_and_pfrich.xml --outputFile ./beam_gas_e_10.edm4hep.tracking_and_pfrich.root --part.userParticleHandler= --inputFiles root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/electron/beam_gas_ep_10GeV_foam_emin10keV_10Mevt_vtx_cs_info.hepmc3.tree.root -N 1000000

# Run Simulation for proton beam+gas
npsim --runType run --compactFile $EIC_SHELL_PREFIX/share/epic/epic_tracking_and_pfrich.xml --outputFile ./beam_gas_proton_275.edm4hep.tracking_and_pfrich.root --part.userParticleHandler= --inputFiles root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/275GeV/pythia8.306-1.0_ProtonBeamGas_275GeV_run001.hepmc3.tree.root -N 100000

I Generated DIS files locally, using code which was quickly based off the pythia example of writing to hepmc, and can be used in the examples folder of a pythia instillation with 'make mymain_WriteToHepMC', or you can copy all of the settings into your hepmc writing file of choice. 

# Run Simulation for DIS
npsim --runType run --compactFile $EIC_SHELL_PREFIX/share/epic/epic_tracking_and_pfrich.xml --outputFile ./DIS_275_10.edm4hep.tracking_and_pfrich.root --part.userParticleHandler= --inputFiles ./PYTHIA8_eic_50000Events_100GeV_decays_off.hepmc -N 20000


# Make Directory to store output plots
mkdir Plots

mkdir Plots/bg_elec Plots/bg_proton Plots/DIS Plots/Total

# Run Analysis Code
root -l 'BackgroundAnalysis.C()'

#No input will default to the file names generated in the previous steps, alternatively
root -l 'BackgroundAnalysis.C("DISFileName","bg_elecFileName","bg_protonFileName");

# Changing Input Files
Luminosity of beam, background rates, and DIS cross section are all set at the top of the file, and should be changed when using differently defined hepmcs as inputs.









