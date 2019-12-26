/******************************************************************
 * muensterTPCsim
 * 
 * Simulations of the Muenster TPC
 * 
 * @author Lutz Althüser
 * @date   2015-04-14
 * @update 2016-01-31
 *
 * @comment - added timestamps
 *					- added 'progress line'
 *					- added some ascii art :)
 ******************************************************************/
#include <G4Event.hh>
#include <string>
#include <stdio.h>
#include <sstream>
#include <unistd.h>
#include <sys/time.h>
#include "muensterTPCEventAction.hh"

muensterTPCEventAction::muensterTPCEventAction(muensterTPCAnalysisManager *pAnalysisManager) {
	m_pAnalysisManager = pAnalysisManager;
	time_un = time(0);
	time_now = localtime(&time_un);

}

muensterTPCEventAction::~muensterTPCEventAction() {
}

void muensterTPCEventAction::BeginOfEventAction(const G4Event *pEvent) {
	if(pEvent->GetEventID() == 0)
	{
		m_iNbEventsToSimulate = m_pAnalysisManager->GetNbEventsToSimulate();
		starttime.str(std::string());
		starttime << time_now->tm_year+1900 << "-" << time_now->tm_mon+1 
		     << "-" << time_now->tm_mday << " " << time_now->tm_hour
		     << "-" << time_now->tm_min << "-" << time_now->tm_sec; 
		starttimeunix = mktime ( time_now ); // for UTC time use "timegm"
		G4cout << "================================================================" << G4endl;
		G4cout << "================================================================" << G4endl;
		G4cout << "  ___                            ___	" << G4endl;
		G4cout << " ||  |   ___         ___        ||  | " << G4endl;
		G4cout << " || _|__/  _\\_______/  _\\_______|| _|  SIMULATION IN PROGRESS!" << G4endl;
		G4cout << " ||(___(  (________(  (_________||((_)	 " << G4endl;
		G4cout << " ||  |  \\___/       \\___/       ||  |   PLEASE DO NOT SWITCH  " << G4endl;
		G4cout << " ||  |         ___              ||  |        OFF THE PC!!       " << G4endl;
		G4cout << " || _|________/  _\\_____________|| _|	" << G4endl;
		G4cout << " ||(_________(  (_______________||((_)	" << G4endl;
		G4cout << " ||  |        \\___/             ||  |    .-=-. ..		 " << G4endl;
		G4cout << " ||  |                          ||  |   ,|`~'|   take		" << G4endl;
		G4cout << " ||  |                          ||  |   `|   |       some 		   " << G4endl;
		G4cout << " ||  |                          ||  |     `~'            coffee!" << G4endl;
		G4cout << "================================================================" << G4endl;		
		G4cout << "================================================================" << G4endl;
		G4cout << "Data file stamp: " << starttime.str() << G4endl;
		G4cout << "================================================================" << G4endl;
	}
	if ( (pEvent->GetEventID() % 1000 == 0) || (pEvent->GetEventID() % int(m_iNbEventsToSimulate/5.+1.) == 0) )
	{
		time_un = time(0);
		time_now = localtime(&time_un);
		currtime.str(std::string());
		currtime << time_now->tm_year+1900 << "-" << time_now->tm_mon+1 
		     << "-" << time_now->tm_mday << " " << time_now->tm_hour
		     << "-" << time_now->tm_min << "-" << time_now->tm_sec;
		currtimeunix = time_un - starttimeunix;
		if ( (pEvent->GetEventID() > 0) && ((pEvent->GetEventID()/currtimeunix) > 0) ) {
			G4cout << "if ( (pEvent->GetEventID() > 0) && ((pEvent->GetEventID()/currtimeunix) > 0) ) {" << G4endl;
			endtimeunix = starttimeunix + m_iNbEventsToSimulate / (pEvent->GetEventID()/currtimeunix);
			time_end = localtime(&endtimeunix);
			strftime(endtime, sizeof(endtime), "%Y-%m-%d %H-%M-%S", time_end);
			G4cout 	<< "\r" << currtime.str() << " || Event " << pEvent->GetEventID()
							<< " / " << m_iNbEventsToSimulate << " || E/s " << pEvent->GetEventID()/currtimeunix 
							<< " || ETA: " << endtime << " " << std::flush;
		}
		else if (pEvent->GetEventID() == 0) {
			G4cout << currtime.str() << " || Start of first event || " << currtimeunix << "s initialization" << G4endl;
		}
	}
	if(m_pAnalysisManager)
		m_pAnalysisManager->BeginOfEvent(pEvent);
}

void muensterTPCEventAction::EndOfEventAction(const G4Event *pEvent) {
	if(pEvent->GetEventID() == ( m_iNbEventsToSimulate - 1 ) )
	{
		time_un = time(0);
		time_now = localtime(&time_un);
		currtime.str(std::string());
		currtime << time_now->tm_year+1900 << "-" << time_now->tm_mon+1 
		     << "-" << time_now->tm_mday << " " << time_now->tm_hour
		     << "-" << time_now->tm_min << "-" << time_now->tm_sec;
		currtimeunix = time_un - starttimeunix;
		G4cout << G4endl << currtime.str() << " || End of last event || " << currtimeunix << "s running time" << G4endl;
		G4cout << "================================================================" << G4endl;
	}	

	if(m_pAnalysisManager)
		m_pAnalysisManager->EndOfEvent(pEvent);
}


