#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4ParticleDefinition.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "TTree.h"
#include "G4ThreeVector.hh"
#include "Analysis.hh"

#include "Randomize.hh"

using namespace std;
class G4ParticleGun;
class G4Event;
class DetectorConstruction;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

  PrimaryGeneratorAction( G4String , int, int, double, double);
  ~PrimaryGeneratorAction();
  //G4double ENER; //ESPR, ANGU_X, ANGU_Y, CORR_X, CORR_Y, SPOT_CX, SPOT_CY, SPOT_CZ, SPOT_X, SPOT_Y, SPOT_Z, RAD;
  
  int ntotal;
  int nbSpots;
  double maxEnergy;
  int maxPart;
  int Anumber; 
  int linecnt;

  vector<G4double> linspace(double , double , double );
  vector<string> split(string, char);

  void GeneratePrimaries(G4Event* );
  void InitializePlan(string );
  void InitializeRippleFilter(G4String);
  static G4String GetPrimaryName() ;                
  static inline PrimaryGeneratorAction* GetInstance() { return theGenerator; }

  Double_t x0,y0,z0,y0TPS,z0TPS,px0,py0,pz0;
  Double_t x1,y1,z1,px1,py1,pz1;
  Double_t Einit,Estop;
  Int_t   Id;
  

  G4RandGeneral * mDistriGeneral;
  G4ThreeVector Position;
  G4ThreeVector Momentum;
  G4ThreeVector isocenter;

  G4int nProtonsGenerated;   
  G4int spot;

  G4ParticleDefinition* particle;
  G4double FWHM;
  
  vector<G4double> spotPosZ;
  vector<G4double> spotPosY;
  vector<G4double> spotEnergy;
  vector<G4double> spotWeight;
  vector<G4int> nbIonsToGenerate;

  G4double mPDF[250];
  G4double EnergyRifi[250];
  G4RandGeneral *EnergyRand; 
 
  private:

  static PrimaryGeneratorAction* theGenerator;
  G4ParticleGun*  	       particleGun;  //pointer a to G4 service class
  DetectorConstruction* theDetector;  
  Analysis* theAnalysis;
  G4double Z_Position;
  G4String PSD_Path;
  G4String PSD_Name;
  G4String PARTICLE;
  G4double fieldSizeZ,fieldSizeY;
};

#endif



