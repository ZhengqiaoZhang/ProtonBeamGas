// Create pythia eP collisions and save to HepMC
// based on example/main36.cc

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "Pythia8Plugins/HepMC3.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "eicBeamShape.h"


#include<iostream>
using std::cout;
using std::cerr;
using std::endl;

#include<string>
using std::string;

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using namespace Pythia8;


int main(int argc, char* argv[])
{
	if(argc != 8)
	{
		cerr << "Wrong number of arguments" << endl;
		cerr << "program.exe steer configuration hadronE leptonE xangle out.hist.root out.hepmc" << endl;
		exit(EXIT_FAILURE);
	}

	const int config = atoi(argv[2]);
	const double hadE = atof(argv[3]);
	const double lepE = atof(argv[4]);
	const double xing = atof(argv[5]);
	const char* rootOut = argv[6];
	const char* hepmcOut = argv[7];

	cout << "Steering File = " << argv[1] << endl;
	if(config == 1) cout << "Configuration = High Divergence" << endl;
	if(config == 2) cout << "Configuration = High Acceptance" << endl; 
	if(config == 3) cout << "Ion Running" << endl;
	cout << "Hadron Energy = " << hadE << endl;
	cout << "Lepton Energy = " << lepE << endl;
	cout << "Beam Crossing Angle (Radians) = " << xing << endl;
	cout << "Root Output = " << rootOut << endl;
	cout << "HepMC Output = " << hepmcOut << endl;

	// Open Root File
	TFile *ofile = TFile::Open(rootOut,"recreate");

	// Histos
	//
	//
	TH1D *eCM = new TH1D("eCM","Modified - Nominal CM Energy",1000,-1.0,1.0);

	TH2D *pXY1 = new TH2D("pXY1","Hadron Beam Py Vs Px",1000,6.5,7.2,1000,-0.1,0.1);
	TH2D *pXZProd1 = new TH2D("pXZProd1","Hadron Beam Px Vs Vertex z",500,-500.,500.,1000,-10.0,10.0);
	TH2D *pYZProd1 = new TH2D("pYZProd1","Hadron Beam Py Vs Vertex z",500,-500.,500.,1000,-10.0,10.0);
	TH1D *pZ1 = new TH1D("pZ1","Hadron Beam Pz",500,265,280.0);

	TH1D *vtxX = new TH1D("vtxX","Vertex x;[mm]",500,-5.0,5.0);
	TH1D *vtxY = new TH1D("vtxY","Vertex y;[mm]",500,-0.1,0.1);
	TH1D *vtxZ = new TH1D("vtxZ","Vertex z;[mm]",500,-5500.0,5500.0);

	TH2D *vtxYvsX = new TH2D("vtxYvsX","Vertex Y vs X;X [mm];Y [mm]",500,-5.0,5.0,500,-0.5,0.5);
	TH2D *vtxXvsZ = new TH2D("vtxXvsZ","Vertex X vs Z;Z [mm];X [mm]",500,-5500.0,550.0,5000,-5.0,5.0);
	TH2D *vtxYvsZ = new TH2D("vtxYvsZ","Vertex Y vs Z;Z [mm];Y [mm]",500,-5500.0,500.0,5000,-5.0,5.0);

	// Particle Quantities
	TH1D *partPtHist_e = new TH1D("partPt_e","Final State Particle Pt per event",500,0.,50.);
	TH1D *partPtHist = new TH1D("partPt","Final State Particle Pt",500,0.,50.);
	TH1D *partEtaHist_e = new TH1D("partEta_e","Final State Particle Eta per event",400,-10.,10.);
	TH1D *partEtaHist = new TH1D("partEta","Final State Particle Eta",400,-10.,10.);
	TH1D *partPhiHist_e = new TH1D("partPhi_e","Final State Particle Phi per event",100,-1.0*TMath::Pi(),TMath::Pi());
	TH1D *partPhiHist = new TH1D("partPhi","Final State Particle Phi",100,-1.0*TMath::Pi(),TMath::Pi());

	TH2D *partPtVsEtaHist = new TH2D("partPtVsEta","Final State Particle Pt Vs Eta",400,-10.,10.,500,0.,50.);
	TH2D *partPhiVsEtaHist = new TH2D("partPhiVsEta","Final State Particle Phi Vs Eta",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

	TH2D *partStatusVsEtaHist = new TH2D("partStatusVsEta","Particle Status Code Vs Eta",400,-10.,10.,200,0.,200.);

	TH2D *vtxZVsEtaHist = new TH2D("vtxZVsEta","vtx Z Vs Eta",400,-10.,10.,500,-5500,5500.);

	// Hi Pt
	TH1D *partPtHiHist_e = new TH1D("partPtHi_e","Final State Particle Pt per event(Pt > 1 GeV)",500,0.,50.);
	TH1D *partEtaHiHist_e = new TH1D("partEtaHi_e","Final State Particle Eta per event(Pt > 1 GeV)",400,-10.,10.);
	TH1D *partPhiHiHist_e = new TH1D("partPhiHi_e","Final State Particle Phi per event (Pt > 1 GeV)",100,-1.0*TMath::Pi(),TMath::Pi());
	TH1D *partPtHiHist = new TH1D("partPtHi","Final State Particle Pt (Pt > 1 GeV)",500,0.,50.);
	TH1D *partEtaHiHist = new TH1D("partEtaHi","Final State Particle Eta (Pt > 1 GeV)",400,-10.,10.);
	TH1D *partPhiHiHist = new TH1D("partPhiHi","Final State Particle Phi (Pt > 1 GeV)",100,-1.0*TMath::Pi(),TMath::Pi());

	TH2D *partPhiVsEtaHiHist = new TH2D("partPhiVsEtaHi","Final State Particle Phi Vs Eta (Pt > 1 GeV)",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

	TH2D *partStatusVsEtaHiHist = new TH2D("partStatusVsEtaHi","Particle Status Code Vs Eta (Pt > 1 GeV)",400,-10.,10.,200,0.,200.);




	//TH1D *ppPx = new TH1D("ppPx","ppPx",


	// Interface for conversion from Pythia8::Event to HepMC event.
	HepMC3::Pythia8ToHepMC3 topHepMC;

	// Specify file where HepMC events will be stored.
	//HepMC3::WriterAscii ascii_io(hepmcOut); // Write in HepMC3 Format
	HepMC3::WriterAsciiHepMC2 ascii_io(hepmcOut); // Write in HepMC2 Format for Delphes

	// Set Up Pythia Event
	Pythia8::Pythia p8;
	Pythia8::Event &event = p8.event;

	// A class to generate beam parameters according to own parametrization.
	BeamShapePtr myBeamShape = make_shared<eicBeamShape>(config,hadE,lepE,xing);

	// Hand pointer to Pythia.
	// If you comment this out you get internal Gaussian-style implementation.
	p8.setBeamShapePtr(myBeamShape);

	// Read in Steering File & Make Other Settings
	p8.readFile(argv[1]);
	p8.readString("Main:timesAllowErrors = 10000"); // allow more errors, eP is brittle

	// Initialize Pythia
	p8.init();

	// Record Nominal COM Energy
	double eCMnom = p8.info.eCM();


	// Run
	int nevents = p8.mode("Main:numberOfEvents");

	//double rate = 31.45/nevents; ///rate due to the beam gas pressure_1 275GeV
	//double rate = 30.74/nevents; ///rate due to the beam gas pressure_1 100GeV
	double rate = 30.96/nevents; ///rate due to the beam gas pressure_1 41GeV
	double scale_e = 1.0/nevents;

	if(nevents > 1000000)
	{
		cout << "Limit nEvents to < 10000 due to size of HepMC file" << endl;
		exit(EXIT_FAILURE);
	}
	for(int ev=0; ev<nevents; ev++)
	{
		if(!p8.next()) continue;
		// p8.event.list();

		// Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
		//
		//
		//p8.event.list();
		double eCMnow = p8.info.eCM();

		eCM->Fill(eCMnow - eCMnom);

		pXY1->Fill(p8.process[1].px(),p8.process[1].py());
		pXZProd1->Fill(p8.process[1].zProd(),p8.process[1].px());
		pYZProd1->Fill(p8.process[1].zProd(),p8.process[1].py());
		pZ1->Fill(p8.process[1].pz());


		vtxX->Fill(p8.process[0].xProd());
		vtxY->Fill(p8.process[0].yProd());
		vtxZ->Fill(-1.0*p8.process[0].zProd());

		vtxYvsX->Fill(p8.process[0].xProd(),p8.process[0].yProd());
		vtxXvsZ->Fill(-1.0*p8.process[0].zProd(),p8.process[0].xProd());
		vtxYvsZ->Fill(-1.0*p8.process[0].zProd(),p8.process[0].yProd());



		for(int i=0; i < p8.event.size(); i++)
		{
			bool partFin = p8.event[i].isFinal();
			double partPt = p8.event[i].pT();
			double partEta = p8.event[i].eta();
			double partPhi = p8.event[i].phi();

			double px = p8.event[i].px();
			double py = p8.event[i].py();
			double pz = p8.event[i].pz();
			double E = p8.event[i].e();

			if(partFin && partEta>-10.0 && partEta<10.0 )
			{
				partPtHist->Fill(partPt);
				partPhiHist->Fill(partPhi);
				partEtaHist->Fill(partEta);
				partPtHist_e->Fill(partPt);
                                partPhiHist_e->Fill(partPhi);
                                partEtaHist_e->Fill(partEta);

				vtxZVsEtaHist->Fill(partEta,-1.0*p8.process[0].zProd());

				partPtVsEtaHist->Fill(partEta,partPt);
				partPhiVsEtaHist->Fill(partEta,partPhi);

				partStatusVsEtaHist->Fill(partEta,p8.event[i].status());

				if(partPt > 1.0)
				{
					partPtHiHist->Fill(partPt);
					partPhiHiHist->Fill(partPhi);
					partEtaHiHist->Fill(partEta);
		
					 partPtHiHist_e->Fill(partPt);
                                        partPhiHiHist_e->Fill(partPhi);
                                        partEtaHiHist_e->Fill(partEta);

					partPhiVsEtaHiHist->Fill(partEta,partPhi);

					partStatusVsEtaHiHist->Fill(partEta,p8.event[i].status());
				}
			}

		}


		// Construct new empty HepMC event and fill it.
		HepMC3::GenEvent hepmcevt;
		topHepMC.fill_next_event( p8, &hepmcevt );

		// Write the HepMC event to file.
		ascii_io.write_event(hepmcevt);
	}

	// List Statistics  
	p8.stat();

	eCM->Scale(rate);

	pXY1->Scale(rate); 
	pXZProd1->Scale(rate); 
	pYZProd1->Scale(rate); 
	pZ1->Scale(rate); 

	vtxX->Scale(rate); 
	vtxY->Scale(rate); 
	vtxZ->Scale(rate); 

	vtxYvsX->Scale(rate); 
	vtxXvsZ->Scale(rate); 
	vtxYvsZ->Scale(rate);

	partPtHist->Scale(rate);
	partEtaHist->Scale(rate);
	partPhiHist->Scale(rate);

	partPtHist_e->Scale(scale_e);
        partEtaHist_e->Scale(scale_e);
        partPhiHist_e->Scale(scale_e);


	partPtVsEtaHist->Scale(rate);
	partPhiVsEtaHist->Scale(rate); 

	partStatusVsEtaHist->Scale(rate);

	vtxZVsEtaHist->Scale(rate); 

	partPtHiHist->Scale(rate);
	partEtaHiHist->Scale(rate);
	partPhiHiHist->Scale(rate);

	  partPtHiHist_e->Scale(scale_e);
        partEtaHiHist_e->Scale(scale_e);
        partPhiHiHist_e->Scale(scale_e);

	partPhiVsEtaHiHist->Scale(rate); 

	partStatusVsEtaHiHist->Scale(rate); 



	// Write and Close Root File
	ofile->Write();
	ofile->Close();
	// Done.
	return 0;


}
