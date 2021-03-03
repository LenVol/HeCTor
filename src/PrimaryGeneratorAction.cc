#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4Alpha.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4IonTable.hh"
#include <math.h>
#include "G4NavigationHistory.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "TFile.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Analysis.hh"


#include "fstream"
#include "iostream"
#include "string"
#include "vector"
#include "time.h"

using namespace std;

PrimaryGeneratorAction* PrimaryGeneratorAction::theGenerator = NULL;

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String RifiFile, int nbPrimaries, int A, double y, double z):ntotal(nbPrimaries),Anumber(A),y0TPS(y),z0TPS(z)//:ENER(energy)
{
  G4cout << "The beam spot is " << y0TPS <<","<<z0TPS << G4endl;
  theGenerator = this;
  theDetector = DetectorConstruction::GetInstance();
  nProtonsGenerated = 0;
  particleGun = new G4ParticleGun();

  InitializeRippleFilter(RifiFile);
  EnergyRand = new G4RandGeneral(mPDF,250,0);
  fieldSizeY = 2*theDetector->NbinsX*theDetector->halfX;   
  fieldSizeZ = 2*theDetector->NbinsZ*theDetector->halfZ;
  for(int i = 0; i<250; i++ ) cout << EnergyRifi[i] << " " << mPDF[i] << endl; // for checking everything from the log files  
/*  if(Anumber == 1){
        particle = G4Proton::Definition();
  }*/
  if(Anumber == 4){
        //particle = G4Alpha::Definition();
	if(RifiFile.find("IES001") != std::string::npos){
          FWHM = 8.5*mm;
	}
	else if (RifiFile.find("IES002") != std::string::npos){
	  FWHM = 8.1*mm;
	}
	else {
	  FWHM = 7*mm;
	}
  }
  else{
        //particle = G4IonTable::GetIonTable()->GetIon(Anumber/2,Anumber,0);
        if(RifiFile.find("IES001") != std::string::npos){
          FWHM = 8.5*mm;
        }
        else if (RifiFile.find("IES002") != std::string::npos){
          FWHM = 9.*mm;
        }
        else {
          FWHM = 8*mm;
        }
  }
  G4cout << "The beam particle has A number of " << Anumber << G4endl;
  G4cout << "The Focus of the beam is " << FWHM << "mm FWHM" << G4endl;

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  theGenerator = NULL;
  delete particleGun;
}


void PrimaryGeneratorAction::InitializeRippleFilter(G4String RifiFile){
  G4cout << "Initializing ripple filter..." << RifiFile.data() << endl;
  std::string line;
  std::ifstream RIFI(RifiFile);
  double data[2];
  linecnt = 0;
  if (RIFI.is_open()){
    while(getline(RIFI, line)) {
      stringstream ss(line);
      for(int i=0;i<2;i++) ss >> data[i];
      EnergyRifi[linecnt] = data[0];
      mPDF[linecnt] = data[1];
      linecnt++;
    }
  }
  else{
    cout << "Could not open Rifi spectrum: Aborting simulation!" << endl;
    exit(0);
  }
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(Anumber==1) particle = G4Proton::Definition();
  else if(Anumber==4) particle = G4Alpha::Definition();
  else particle = G4IonTable::GetIonTable()->GetIon(Anumber/2,Anumber,0);
  particleGun->SetParticleDefinition(particle);
  int Eid = 250*EnergyRand->fire();
  Einit = Anumber*EnergyRifi[Eid]*MeV;

  // Ideal pencil beam, we neglect the beam divergence here 
  px0 = 1; 
  py0 = 0.;
  pz0 = 0.;

  x0 = -102.*cm; // start at BAMS

//  y0TPS = 12.*mm; //relative to isocenter //TODO: for some reason y -> -y and z->-z in the Geant4 representation of the CT scan
//  z0TPS = 0.*mm; 
   
  y0 = G4RandGauss::shoot(y0TPS,FWHM/2.355);//3.44 for ADAM C, 2.97 for He//Pencil beam with 8mm FWHM around pencil beam position
  z0 = G4RandGauss::shoot(z0TPS,FWHM/2.355);//3.44 for ADAM C, 2.97 for He//assuming equal spread in x and y (slight difference between the two fro real HIT beam)

//  y0 = fieldSizeY*G4UniformRand()-fieldSizeY/2 + theDetector->shift.y();
//  z0 = fieldSizeZ*G4UniformRand()-fieldSizeZ/2 + theDetector->shift.z();

  Position = G4ThreeVector(x0,y0,z0);
  Momentum = G4ThreeVector(px0,py0,pz0);
  Momentum.setMag(1);

  particleGun->SetParticleEnergy(Einit);
  particleGun->SetParticleMomentumDirection(Momentum);
  particleGun->SetParticlePosition(Position);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  nProtonsGenerated++;
  if(nProtonsGenerated%10000==0){
    G4cout << "********************************************************************" << G4endl;
    G4cout << "Total Progress: " << float(nProtonsGenerated)/float(ntotal)*100 << "%" << G4endl;
    G4cout << "********************************************************************" << G4endl;
  }

}

vector<G4double> PrimaryGeneratorAction::linspace(double a, double b, double step) {
  vector<double> array;
  while(a <= b) {
    array.push_back(a);
    a += step;
  }
  return array;
}

vector<string> PrimaryGeneratorAction::split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
 
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
 
  return internal;
}





