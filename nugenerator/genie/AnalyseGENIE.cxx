// Do make to compile 
// Do AnalyseGENIE -f ghep_file.root -n events nbr to analyse

//____________________________________________________________________________
/*!

\program gtestEventLoop

\brief   Example event loop. Use that as a template for your analysis code.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
		 STFC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
		 For the full text of the license visit http://copyright.genie-mc.org
		 or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>
#include <TH1.h>
#include <TGraph.h>

#include "/data/pmartins/GENIE/v284/R-2_8_4/src/EVGCore/EventRecord.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/GHEP/GHepParticle.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/Ntuple/NtpMCFormat.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/Ntuple/NtpMCTreeHeader.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/Ntuple/NtpMCEventRecord.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/Messenger/Messenger.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/PDG/PDGCodes.h"
#include "/data/pmartins/GENIE/v284/R-2_8_4/src/Utils/CmdLnArgParser.h"

using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
int    gOptNEvt;
string gOptInpFilename;
string gOptOutFilename;
string gOptFluxFile;
string gOptSplineFile;
//___________________________________________________________________
int main(int argc, char ** argv)
{
	GetCommandLineArgs (argc, argv);
	// Open the ROOT file (XXX.ghep.root) (-f argument)
	// and get the TTree & its header
	TTree *           tree = 0;
	NtpMCTreeHeader * thdr = 0;
 	TFile file(gOptInpFilename.c_str(),"READ");
 	tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
 	thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
 	if(!tree) return 1;
 	NtpMCEventRecord * mcrec = 0;
 	tree->SetBranchAddress("gmcrec", &mcrec);
 	// Get the nbr of evts to analyse (-n argument)
 	int nev = (gOptNEvt > 0) ?
 	TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
 	(int) tree->GetEntries();
 	// Create a ROOT file (-o argument) with all the relevent histos
 	TFile outfile(gOptOutFilename.c_str(),"RECREATE");

  	TH1D * h_pion_nrj    = new TH1D("pion_nrj",  "pion_nrj",    500,0,5);   // Pion energy
  	TH1D * h_pion_mom    = new TH1D("pion_mom",  "pion_mom",    500,0,5);   // Pion mnomentum
  	TH1D * h_pion_ang    = new TH1D("pion_ang",  "pion_ang",    180,0,180); // Pion angle
  	TH1D * h_Q2          = new TH1D("Q2",        "Q2",          100,0,1);   // Q2
	
	// Get the splines used (converted in root using gspl2root)
	// This is needed to normalize your distribution
	TFile splinefile(gOptSplineFile.c_str(),"READ");
	TGraph* XS_graph;
	gDirectory->GetObject("nu_mu_C12/coh_cc",XS_graph);

	double nbrcoh=0;

	// Loop over all events
	for(int i = 0; i < nev; i++) {
		// get next tree entry
		tree->GetEntry(i);
		// get the GENIE event
		EventRecord &  event = *(mcrec->event);
		// Print out all events infos
		// LOG("myAnalysis", pNOTICE) << event;

		// Initialize the particles info for each event
		double pion_nrj=0;
		double pion_momentum=0;
		double pion_angle=0;
		double Q2=0;

		// Select only CC coherent interactions
		Interaction * in = event.Summary();
		const ProcessInfo & proc = in->ProcInfo();

		if(proc.IsCoherent() && proc.IsWeakCC()){
			nbrcoh++;

			GHepParticle * neu = event.Probe();
			const TLorentzVector & neu4vec = *(neu->P4());
			double total_xs = XS_graph->Eval(neu->E());

		 	TLorentzVector Pion4vec ;
		 	TLorentzVector Muon4vec ;

			GHepParticle * p = 0;
			TIter event_iter(&event);
			//LOG("myAnalysis", pNOTICE) <<"---------------------------------" ;
			//LOG("myAnalysis", pNOTICE) << "Event number  : " << i  <<"Neutrino E = "<< neu->E() <<" GeV , coh XS = "<<XS_graph->Eval(neu->E()); 
			//LOG("myAnalysis", pNOTICE) <<"---------------------------------" ;
			// Loop over all particles in this event
			while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))){
		 		if(p->Pdg() == 211){ // we have the pion
		 			Pion4vec = *(p->P4());
		 			double momentum = sqrt(pow(Pion4vec.Px(),2) + pow(Pion4vec.Py(),2) + pow(Pion4vec.Pz(),2));
		 			//LOG("myAnalysis", pNOTICE)  << "Got a : " << p->Name() << "("<<p->Pdg()<<")"<<" with mom = " << momentum << " GeV"; 
		 			pion_momentum = momentum;
		 			pion_nrj = p->E();
		 			pion_angle = 180*acos((neu4vec.Px() * Pion4vec.Px() + neu4vec.Py() * Pion4vec.Py() + neu4vec.Pz() * Pion4vec.Pz()) / (momentum * sqrt(pow(neu4vec.Px(),2) + pow(neu4vec.Py(),2) + pow(neu4vec.Pz(),2))))/3.14;
		 		}
		 		if(p->Pdg() == 13){ // we have the muon
		 			Muon4vec = *(p->P4());
		 		}
		 		TLorentzVector q = (neu4vec - Muon4vec);

		 		Q2 =  -1 * q.M2();
			}// end loop over particles	

			// Fill the variables weighted by the total xs from the spline
			h_pion_nrj->Fill( pion_nrj,total_xs);
			h_pion_mom->Fill( pion_momentum,total_xs);
			h_pion_ang->Fill( pion_angle, total_xs);
			h_Q2->Fill( Q2,total_xs);
		}
		mcrec->Clear(); // clear current mc event record
	}//end loop over events

	// Get the flux used - not always needed
	//TFile fluxfile(gOptFluxFile.c_str(),"READ");
	//TH1D * fluxhisto = (*TH1D)fluxfile->Get("spectrum");

	// Scale to take binning and nbr generated events into account
	double rate = 10 *  1 /nbrcoh;
	h_pion_nrj->Scale(rate/0.01);
	h_pion_mom->Scale(rate/0.01);
	h_pion_ang->Scale(rate);
	h_Q2->Scale(rate/0.01);
// close input GHEP event file
outfile.Write();
outfile.Close();

file.Close();
LOG("myAnalysis", pNOTICE)  << "Done!";
return 0;
}

//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
	LOG("myAnalysis", pINFO) << "Parsing commad line arguments";
	CmdLnArgParser parser(argc,argv);

  	// get GENIE event sample
	if( parser.OptionExists('i') ) {
		LOG("myAnalysis", pINFO) << "Reading event sample filename";
		gOptInpFilename = parser.ArgAsString('i');
	} else {
		LOG("myAnalysis", pFATAL) << "Unspecified input filename - Exiting";
		exit(1);
	}

  	// number of events to analyse
	if( parser.OptionExists('n') ) {
		LOG("myAnalysis", pINFO) << "Reading number of events to analyze";
		gOptNEvt = parser.ArgAsInt('n');
	} else {
		LOG("myAnalysis", pINFO)<< "Unspecified number of events to analyze - Use all";
		gOptNEvt = -1;
	}

	// Output root file
	if( parser.OptionExists('o') ) {
		LOG("myAnalysis", pINFO) << "Creating output ROOT file";
		gOptOutFilename = parser.ArgAsString('o');
	} else {
		LOG("myAnalysis", pFATAL) << "Unspecified output filename";
	}

	// Flux file used in root file
	if( parser.OptionExists('f') ) {
		LOG("myAnalysis", pINFO) << "Reading flux file";
		gOptFluxFile = parser.ArgAsString('f');
	} else {
		LOG("myAnalysis", pFATAL) << "Unspecified flux filename";
	}

	// Splines used in root file
	if( parser.OptionExists('s') ) {
		LOG("myAnalysis", pINFO) << "Reading splines file";
		gOptSplineFile = parser.ArgAsString('s');
	} else {
		LOG("myAnalysis", pFATAL) << "Unspecified splines filename";
	}
}
//_________________________________________________________________________________
