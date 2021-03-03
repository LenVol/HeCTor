#include "vector"
#include "G4ios.hh"
#include "TTree.h"
#include "TFile.h"
#include "TH3D.h"
#include "globals.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4EmCalculator.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "OrganicMaterial.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Analysis.hh"
#include "Randomize.hh"

Analysis* Analysis::theAnalysis=NULL;
Analysis::~Analysis(){theAnalysis=NULL;}


Analysis::Analysis(G4String theName){
  theAnalysis  = this;

  theGenerator = PrimaryGeneratorAction::GetInstance();
  theDetector  = DetectorConstruction::GetInstance();

  // Check if the file name already exists in order not to overwrite large simulation data with nothing
  string file = Form("%s_%.0f.root",theName.data(),theDetector->theAngle);
  ifstream f(file.c_str());
  if(f.good()){
	char cont;
	cout << "File " << file << " exists. Are you sure you want to continue? Type y or n." << endl;
	cin >> cont;
	if (cont!='y'){
		cout << "Simulation will be aborted!" << endl;
		exit(1);
	}
  }

  //I use the root environment to handle data

  f1 = new TFile(file.c_str(),"recreate");

  t = new TTree("phase","PS");

  //Data for debugging and particle imaging
  t->Branch("x0",&x0,"x0/D"); // Real starting point
  t->Branch("y0",&y0,"y0/D"); // This contains a gaussian sampling around the set pencil beam position
  t->Branch("z0",&z0,"z0/D");
  t->Branch("y0TPS",&y0TPS,"y0TPS/D"); // Pencil beam position
  t->Branch("z0TPS",&z0TPS,"z0TPS/D"); 
  
  t->Branch("px0",&px0,"px0/D"); //initial direction (no gaussian sampling, divergence of initial beam after BAMS ignored)
  t->Branch("py0",&py0,"py0/D");
  t->Branch("pz0",&pz0,"pz0/D");

  t->Branch("Einit",&Einit,"Einit/D"); // inital energy after BAMS (to debug rifi spectrum)

  t->Branch("x1",&x1,"x1/D"); //Exit position measured on ideal scorer 
  t->Branch("y1",&y1,"y1/D");
  t->Branch("z1",&z1,"z1/D");

  t->Branch("px1",&px1,"px1/D"); //Exit direction measured with ideal scorer
  t->Branch("py1",&py1,"py1/D");
  t->Branch("pz1",&pz1,"pz1/D");

  t->Branch("Estop",&Estop,"Estop/D"); // energy after crossing the phantom measured on ideal scorer
  t->Branch("Id",&Id,"Id/I"); // Particle ID

  // Data relevant to paper plots:
  calHist = new TH1D("calHist", "Energy deposit in calorimeter",120, 0., 120.); // 1mm (sampled in 1mm for later rebinning before plotting)
  calHistBirks = new TH1D("calHistBirks", "Energy deposit in calorimeter",120, 0., 120.); // 1mm (same as above)

  // 3D histogram of dose overlapping with CT scan
  doseHist = new TH3D("doseHist","3D Dose distribution",theDetector->NbinsX,-theDetector->NbinsX*theDetector->halfX+theDetector->shift.x(),theDetector->NbinsX*theDetector->halfX+theDetector->shift.x(),
		theDetector->NbinsY,-theDetector->NbinsY*theDetector->halfY+theDetector->shift.y(),theDetector->NbinsY*theDetector->halfY+theDetector->shift.y(),
		theDetector->NbinsZ,-theDetector->NbinsZ*theDetector->halfZ+theDetector->shift.z(), theDetector->NbinsZ*theDetector->halfZ+theDetector->shift.z());
}



void Analysis::analyseHit(G4Step* aStep, G4String theName) //Analyse hit on the ideal scorer
{
  
  theSteppingAction = SteppingAction::GetInstance(); 
  f1->cd();
  if(theName=="sd2"){
    x0  = theGenerator->x0;
    y0  = theGenerator->y0;
    z0  = theGenerator->z0;
    y0TPS  = theGenerator->y0TPS;
    z0TPS  = theGenerator->z0TPS;
    px0  = theGenerator->px0;
    py0  = theGenerator->py0;
    pz0  = theGenerator->pz0;
    Einit = theGenerator->Einit;    
    
    x1  = aStep->GetPostStepPoint()->GetPosition()[0];
    y1  = aStep->GetPostStepPoint()->GetPosition()[1];
    z1  = aStep->GetPostStepPoint()->GetPosition()[2];
    
    px1 = aStep->GetPostStepPoint()->GetMomentumDirection()[0];
    py1 = aStep->GetPostStepPoint()->GetMomentumDirection()[1];
    pz1 = aStep->GetPostStepPoint()->GetMomentumDirection()[2];
    
    Id    = aStep->GetTrack()->GetTrackID();
    Estop = aStep->GetPostStepPoint()->GetKineticEnergy();
    t->Fill();  
  }
}

void Analysis::calorimeterHit(G4Step *aStep){ //Analise a hit in the calorimeter
	Edep = aStep->GetTotalEnergyDeposit(); // total energy deposit in MeV
	step_length = aStep->GetStepLength(); // step length to get dE/dx
	if(step_length>0.){ // to avoid NaNs
	    Edep_Birks = (Edep)/(1+0.0954*(Edep/step_length));//Birks Law with KB = 0.0954
	    depth = aStep->GetPreStepPoint()->GetPosition() + G4UniformRand()*(aStep->GetPostStepPoint()->GetPosition() - aStep->GetPreStepPoint()->GetPosition()); //randomize the point of energy dep.
	    rel_depth = depth[0] - theDetector->CalorimeterEntrance; // depth relative to calorimeter entrance 
	    calHist->Fill(rel_depth, Edep);//write out both a quenched and unquenched curve for debugging (This was to check the Birks constant)
	    calHistBirks->Fill(rel_depth,Edep_Birks);
	}
}


void Analysis::DoseAnalysis(G4double x, G4double y, G4double z, G4double dose){ // Score the dose to a voxel
	doseHist->Fill(x,y,z,dose);
}

void Analysis::Save(){
  f1->cd();
  t->Write("",TObject::kOverwrite);
  calHist->Write("spotHist",TObject::kOverwrite);
  calHistBirks->Write("spotHist_birks",TObject::kOverwrite);
  doseHist->Write();
  f1->Close();
}
