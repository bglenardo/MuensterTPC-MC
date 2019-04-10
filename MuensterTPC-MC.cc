/******************************************************************
 * muensterTPCsim
 * 
 * Simulations of the Muenster dual phase xenon TPC
 * 
 * @author Lutz Althueser
 *
 * @comment This is the main class of the Muenster dual phase xenon TPC.
 *					It includes ..
 *						- argument handling for interactive or batch mode
 *						- naming of the output data file
 *						- enable/disable Multithreading
 *
 *					The simulation is ready for ..
 * 						- dual and single phase simulations (change LXe -> GXe)
 * 						  see DetectorMessenger.cc
 *						- different source types (shape and particle)
 *						  see ParticleSourceMessenger.cc
 *						- different PhysicsLists are available (incl. Cerenkov) 
 *							see PhysicsList.cc
 ******************************************************************/
// include C++ classes
#include <iostream>
#include <string>
#include <sstream>
#include <unistd.h>
#include <sys/time.h>

// include GEANT4 classes
//#ifdef G4MULTITHREADED
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif
#include <G4UImanager.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>

//include Muenster TPC simulation classes
#include "muensterTPCDetectorConstruction.hh"
#include "muensterTPCPhysicsList.hh"
#include "muensterTPCPrimaryGeneratorAction.hh"
#include "muensterTPCAnalysisManager.hh"
#include "muensterTPCStackingAction.hh"
#include "muensterTPCRunAction.hh"
#include "muensterTPCEventAction.hh"
#include "muensterTPCActionInitialization.hh"

void usage();
inline bool fileexists (const std::string& name);
inline bool fileexists (const char* name);

//******************************************************************/
// call of the main function with optional arguments 
//******************************************************************/
int main(int argc, char **argv)
{
	G4cout 	<< "===========================================" 	<< G4endl;
	G4cout	<< "Starting MuensterTPC-Simulation ..."			<< G4endl;
	G4cout 	<< "===========================================" 	<< G4endl;
	
	// get local time for the output file time stamp
	time_t time_un;
	tm *time_now;
	time_un = time(0);
	time_now = localtime(&time_un);
	std::stringstream hTimeStamp;
	hTimeStamp << time_now->tm_year+1900 << "-" << time_now->tm_mon+1 
		     << "-" << time_now->tm_mday << "_" << time_now->tm_hour
		     << "-" << time_now->tm_min << "-" << time_now->tm_sec;

	// switches
	int c = 0;
	
	bool bInteractive = false;
	bool bPreInitFromFile = false;
	bool bMacroFile = false;
	bool bDataFilename = false;
	bool bVerbosities = false;
	int iVerbosities = 0;
	int iNbEventsToSimulate = 0;
	std::string hPreInitFilename, hMacroFilename, hDataFilename;
	std::stringstream hStream;
	
	// parse switches
	// p: preinit with custom file
	// f: macro file (init source type)
	// o: output file name
	// n: number of events to simulate
	// i: interactive session
	// v: turn on debug verbosities
	if ( argc == 1 ) { bInteractive = true; }
	while((c = getopt(argc,argv,"p:f:o:n:v:i")) != -1) {
		switch(c)	{
			case 'p':
				bPreInitFromFile = true;
				if (fileexists(optarg))
					hPreInitFilename = optarg;
				else
					{ G4cout << "File '" << optarg << "' not found!" << G4endl; usage(); }
				break;

			case 'f':
				bMacroFile = true;
				if (fileexists(optarg))
					hMacroFilename = optarg;
				else
					{ G4cout << "File '" << optarg << "' not found!" << G4endl; usage(); }
				break;

			case 'o':
				bDataFilename = true;
				hDataFilename = optarg;
				break;

			case 'n':
				hStream.str(optarg);
				hStream.clear();
				hStream >> iNbEventsToSimulate;
				break;
				
			case 'v':
				bVerbosities = true;
				hStream.str(optarg);
				hStream.clear();
				hStream >> iVerbosities;
				break;
				
			case 'i':
				bInteractive = true;
				break;

			default:
				usage();
		}
	}
	
	G4UIExecutive* ui = 0;
	if (bInteractive) {	ui = new G4UIExecutive(argc, argv, "Qt"); }
	else if (!bMacroFile) { usage(); }
	
	std::stringstream DatafileName;
	DatafileName.clear();
 	size_t found;
	if (bDataFilename) {
		found=hDataFilename.find_last_of("/");
 		if (found!=std::string::npos)
			DatafileName << hDataFilename.substr(0,found+1) << hTimeStamp.str() << "_" << hDataFilename.substr(found+1);
		else
			DatafileName << hTimeStamp.str() << "_" << hDataFilename;
	}
	else {
		if (bMacroFile) {
			std::string hMacroFilenameshort = hMacroFilename;
			found=hMacroFilenameshort.find_last_of(".mac");
			if (found!=std::string::npos) { hMacroFilenameshort=hMacroFilenameshort.substr(0,found-3); }
			found=hMacroFilenameshort.find_last_of("/");
			if (found!=std::string::npos)
				DatafileName << hTimeStamp.str() << "_" << hMacroFilenameshort.substr(found+1) << ".root";
			else
				DatafileName << hTimeStamp.str() << "_" << hMacroFilenameshort << ".root";
		}
		else {
			DatafileName << hTimeStamp.str() << "_" << "events" << ".root";
		}
	}

	// create the run manager
	//#ifdef G4MULTITHREADED
	//	G4MTRunManager *pRunManager = new G4MTRunManager;
	//	pRunManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
	//#else
		G4RunManager *pRunManager = new G4RunManager;
	//#endif
	
	// set user-defined initialization classes
	pRunManager->SetUserInitialization(new muensterTPCDetectorConstruction);
	pRunManager->SetUserInitialization(new muensterTPCPhysicsList);
	
	// create the primary generator action
        G4cout << "Creating the primary generator action..." << G4endl;
	muensterTPCPrimaryGeneratorAction *pPrimaryGeneratorAction = new muensterTPCPrimaryGeneratorAction();
	pRunManager->SetUserInitialization(new muensterTPCActionInitialization(DatafileName.str(), pPrimaryGeneratorAction));

	// start visualization and ui manager
	G4VisManager* pVisManager = new G4VisExecutive;
	pVisManager->Initialize();
	G4UImanager* pUImanager = G4UImanager::GetUIpointer();
	
	G4String hCommand;
	if(bPreInitFromFile) {
		hCommand = "/control/execute " + hPreInitFilename;
		pUImanager->ApplyCommand(hCommand);
	}
	else if(bVerbosities) {
		hStream.str("");
		hStream.clear();
		hStream << "/control/verbose " << iVerbosities;
		pUImanager->ApplyCommand(hStream.str());
		pUImanager->ApplyCommand("/control/saveHistory");
		hStream.str("");
		hStream.clear();
		hStream << "/run/verbose " << iVerbosities;
		pUImanager->ApplyCommand(hStream.str());
		hStream.str("");
		hStream.clear();
		hStream << "/event/verbose " << iVerbosities;
		pUImanager->ApplyCommand(hStream.str());
		hStream.str("");
		hStream.clear();
		hStream << "/tracking/verbose " << iVerbosities;
		pUImanager->ApplyCommand(hStream.str());
		pUImanager->ApplyCommand("/control/execute macros/preinit_debug.mac");
	}
	else {
		pUImanager->ApplyCommand("/control/execute macros/preinit.mac");
	}
	
	pRunManager->Initialize();
	
 	if ( bInteractive ) { pUImanager->ApplyCommand("/control/execute macros/vis.mac"); }
	
	// run time parameter settings
	if(bMacroFile) {
		hCommand = "/control/execute macros/setup_optical_Muenster.mac";
		pUImanager->ApplyCommand(hCommand);
		hCommand = "/control/execute " + hMacroFilename;
		pUImanager->ApplyCommand(hCommand);
	}
		
	if(iNbEventsToSimulate)	{
		hStream.str("");
		hStream.clear();
		hStream << "/run/beamOn " << iNbEventsToSimulate;
		pUImanager->ApplyCommand(hStream.str());
	}
	
	if ( bInteractive ) { 
			// start ui and delete it when the user closes the ui
		ui->SessionStart();  	
		delete ui;
	}
  
	// delete all created objects
	if(bInteractive) delete pVisManager;
	delete pRunManager;
	
	G4cout 	<< "Datafile: " <<  DatafileName.str()				<< G4endl;
	G4cout 	<< "===========================================" 	<< G4endl;
	G4cout	<< "Simulation successfully completed :)"			<< G4endl;
	G4cout 	<< "===========================================" 	<< G4endl;
	return 0;	
}

void usage() {
  exit(0);
}

inline bool fileexists (const std::string& name) {
	return fileexists(name.c_str());
}

inline bool fileexists (const char* name) {
    std::ifstream f(name);
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

