# MuensterTPC-MC
The MuensterTPC Monte Carlo (MC) simulation was developed in order to built and test a small dual phase xenon TPC (height: 17cm, diameter: 8cm) at Muenster. This includes simulations with optical photons, external and internal radioactive sources. Some simulation results and plots can be found here: [bachelor thesis Althüser (2015)](https://www.uni-muenster.de/imperia/md/content/physik_kp/agweinheimer/theses/bachelor_lutz_althueser.pdf). This software stands under the BSD license (see [LICENSE](./LICENSE.md)).

Some analysis scripts are available in another repository: [MC-Analyzer](https://github.com/l-althueser/MC-Analyzer)

## Installation Prerequisites
This software can only be used with `GEANT4 10.02.p01 64bit` in Linux or Windows Subsystem for Linux. 

## Installing the MuensterTPC-Simulation
There are two ways to install this simulation. Option one is to use the GitHub interface to download and extract the master branch to a specific directory ([download link](https://github.com/l-althueser/MuensterTPC-Simulation/archive/master.zip)). Option two will follow the real GitHub way of doing this ..  

First `cd` to the folder you want the Analyzer to be installed (this should be your `geant4_workdir`). Therefore we assume that Geant4 is correctly installed. Then run:
```
cd $G4WORKDIR
git clone https://github.com/l-althueser/MuensterTPC-Simulation.git
```
Second you have to compile the Geant4 code and link correctly to the generated binary file. (Note: You can always rename `MuensterTPC-Simulation.cc` to `MuensterTPC-*.cc` and the makefile will still work.)
```
cd MuensterTPC-MC
make clean
make
make link
```
Now you should be able to run the simulation in `interactive` mode by typing:
```
./MuensterTPC-MC
```

## MuensterTPC-Simulation Tutorial
This section assumes that the MuensterTPC-Simulation is installed and all prerequisits are given. 

### Overview of the Geant4 Simulation
![Overview_01](/drawings/muensterTPCsim_overview_01.png)
![Overview_02](/drawings/muensterTPCsim_overview_02.png)

### Usage
The simulation offers the possibility to use some arguments in order to adjust every run time parameter.
```
./MuensterTPC-MC -p <custom_preinit.mac> -f <source_definition.mac> -o <outputfilename> -n <number_of_events> -v <verbositie_level> -i
```
* `-p <custom_preinit.mac>`: A default `preinit.mac` will be used if no custom file is given.
* `-f <source_definition.mac>`: This parameter has to be specified if `-i` is not set.
* `-o <outputfilename>`: The output file name will be `events.root` or `<source_definition>.root` if not specified.
* `-n <number_of_events>`: Has to be specified if `-i` is not set.
* `-v <verbositie_level>`: The verbosity level is `0` per default.
* `-i`: This activates the `interactive` mode in a Qt window.

### Simple `opticalphoton` simulation
```
./MuensterTPC-MC -f ./macros/src_optPhot_DP_S1.mac -o optPhot_S1_1e5.root -n 100000
```

### Advanced custom simulation
There are two options to confine the generation of the primary particle vertexes: 
* confine into a specific region, for example inside a cylinder or a cube  
* confine into a detector volume;
You can confine primary particles wherever you like using the following command lines: 
```
/Xe/gun/type Volume
/Xe/gun/shape Cylinder
/Xe/gun/halfz 135 mm
/Xe/gun/radius 50 mm
/Xe/gun/center 0 0 -84.5 mm
```
or as a point source:
```
/Xe/gun/type Point
/Xe/gun/center 0 -160 -84.5 mm
/Xe/gun/angtype iso
```
For example, in the first snippet of code you are confining particles in a volume with a cylindrical shape centred in (0.,0.,-84.5) and with height of 270mm and radius of 50mm. Primary particles will be then generated uniformly inside that volume.  
You can also be more specific by confining the generation volume to a specific detector volume. Let's confine for example a generation into the LXe (see below for volume names). In that case, we have to add to the previous code the following line: 
```
/Xe/gun/confine LXe
```
Pay attention to define the dimension of the confinement volume slightly larger than the detector volume. Once you have confined the generation, you can choose the type of primary particle. You can select a geantino, a neutron, etc…. by typing 
```
/Xe/gun/energy 6.98 eV
/Xe/gun/particle opticalphoton
```
or:
```
/Xe/gun/energy 0 keV
/Xe/gun/particle ion
/Xe/gun/ion 27 60 0 0
```
Alternatively, you can give an energy different from 0 or you can extract the energy spectrum from a file. Let's now suppose that we want to simulate neutrons, in this case we substitute the previous snippet of code with the following one:
```
/Xe/gun/particle neutron
/Xe/gun/energytype Spectrum
/Xe/gun/energyspectrum macros/xenon1t/spectra/neutron/238U.dat
```
In this case we are generating neutrons with energy that follows the energy spectrum defined in the `238U.dat` file.

### The output file/file format
You can simply view the generated simulation data with any version of [ROOT](https://root.cern.ch/). Just type ..
```
root
new TBrowser
```
.. and navigate to the generated `events.root` file. Now you can do the normal click and drag routine of ROOT.  

The output file has an specific file format which is described in the following.
#### Top directory
| Name | type | description |  
| --- | --- | --- |
| G4VERSION | TName<string> | version of geant4 |  
| MC_TAG | TName<string> | name of the simulation toolkit ("muensterTPC") |  

#### TDirectory::events
| Name | type | description |  
| --- | --- | --- |
| nbevents | TParameter<int> | number of simulated events |  

#### TDirectory::events/events
| Name | type | description |  
| --- | --- | --- |
| eventid | int | event number |
| ntpmthits | int | |
| nbpmthits | int | |
| pmthits | int | |
| etot | float | total G4 energy deposit in this event |
| nsteps | int | number of G4 steps |
| trackid  | int | track ID |
| type  | string | particle type |
| parentid  | int | track ID of parent |
| parenttype  | string | particle type of parent |
| creaproc  | string | process that created this particle |
| edproc  | string | process for this particular energy deposit |
| xp  | vector<float> | x coordinate of energy deposit (mm) |
| yp  | vector<float> | y coordinate of energy deposit (mm) |
| zp  | vector<float> | z coordinate of energy deposit (mm) |
| ed  | vector<float> | energy deposit (keV) |
| time  | vector<float> | timestamp of the current particle/trackid |
| type_pri  | string | particle type of primary  |
| e_pri  | vector<float> | energy of primary (keV) |
| xp_pri  | vector<float> | x coordinate of primary particle (mm) |
| yp_pri  | vector<float> | y coordinate of primary particle (mm) |
| zp_pri  | vector<float> | z coordinate of primary particle (mm) |

### Detector geometry
You can use the `interactive` mode to determine every volume name. This are the most recent ones:
* LXe
* GXe (with LXe as mother volume)
* CopperRings*
* Pmt*
* PTFE*
* GXeGridMesh*
* LXeGridMesh*
* ...  

### Sensitive detectors
Two sensitive detectors are defined in the code: 
* muensterTPCLXeSensitiveDetector  
* muensterTPCPmtSensitiveDetector
 
Both detectors are created in each simulation (as well as the corresponding hits collections) but the filled hits depends on the particle type. For example, only simulating optical photons you can fill the PmtHitsCollection.

### Using decay0
Currently there are macros for 0vKb+ and 2vKb+ decays which use events generated with DECAY0. These files in the `macros/events` folder contain 10000 events each. Example scripts how to run the corresponding MC can be found in the `scripts` folder. Example ROOT files can be found in `example_spectra_decay0`.
