{

	TFile * infile = new TFile("2019-12-26_17-39-11_Be7_neutrino_ccapture_to_133keV_level.root");
	TDirectory * dir = (TDirectory *)infile->Get("events");
	TTree * evts = (TTree *)dir->Get("events");
	evts->Print();


}
