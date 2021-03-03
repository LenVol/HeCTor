#ifndef Analysis_hh
#define Analysis_hh
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TTree.h"
#include "TString.h"
#include "G4Step.hh"
#include "TH1.h"
#include "TProfile.h"
#include "TH3.h"

using namespace std;

class G4Step;
class TFile ;
class TTree ;
class TH3D  ;
class PrimaryGeneratorAction;
class DetectorConstruction;
class SteppingAction;
class Analysis
{
public:

  Analysis(G4String);
  ~Analysis();
  static inline Analysis* GetInstance() { return theAnalysis; }
  void analyseHit(G4Step*,G4String);
  void calorimeterHit(G4Step* );
  void DoseAnalysis(G4double, G4double, G4double, G4double); 
  
  TTree  *t;
  TFile *f1;
  TH1D* calHist;
  TH1D* calHistBirks;
  TH3D* doseHist;
  void Save();
  void EndOfSpotAction(int);

  //vector<double>* tracks_X = new vector<double>;
  //vector<double>* tracks_Y = new vector<double>;
  //vector<double>* tracks_Z = new vector<double>;
  //vector<double>* tracks_E = new vector<double>;
  //vector<double>* Radlen   = new vector<double>;
  //vector<TString>*  mat_name;

private:

  static Analysis* theAnalysis;
  PrimaryGeneratorAction *theGenerator;
  DetectorConstruction   *theDetector ;
  SteppingAction         *theSteppingAction;
  G4double cubicVolume;
  Double_t x0,y0,z0,y0TPS,z0TPS,px0,py0,pz0;
  Double_t x1,y1,z1,px1,py1,pz1;
  Double_t Einit,Estop;
  G4int Id;
  TString pName;
  G4double Edep;
  G4double rel_depth;
  G4ThreeVector depth;
  G4double step_length;
  G4double Edep_Birks;
  G4double maxRangeSpot;

  ofstream spotRangeFile;
  
};
#endif
