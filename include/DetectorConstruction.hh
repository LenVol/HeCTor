#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "OrganicMaterial.hh"
#include "Analysis.hh"
#include <map>

class OrganicMaterial;
class Analysis;
class G4Box;
class G4VSolid;
class G4NistManager;
class PrimaryGeneratorAction;
class G4PhantomParameterisation;
class G4PVParameterised;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVParameterised;
class G4SubtractionSolid;
class G4Material;
class TH3S;
class TH3D;
class TFile;


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(G4String,G4double);
  ~DetectorConstruction();
  virtual G4VPhysicalVolume* Construct();  
  static inline DetectorConstruction* GetInstance() { return theDetector; }

  G4PhantomParameterisation* GetParam() { return param; }
  G4double theAngle;
  G4String thePhantom;
  G4int NbinsX,NbinsY,NbinsZ;
  G4double Xmax,Xmin,halfX;
  G4double Ymax,Ymin,halfY;
  G4double Zmin,Zmax,halfZ;
  G4double midX,midY,midZ ;
  G4double CalorimeterEntrance;
  G4double degrader_thick;
  G4double thickX;
  G4double phantom_degrader_distance;
  G4ThreeVector shift;
  G4ThreeVector isocenter;

  bool withCT;

  TH3D* hu;
  vector<G4Material*> theMaterialList;
  G4Material* water;
  G4Material* air;
  map<G4int,G4int> hu2id;
  vector<G4int> huList;
private:
  static DetectorConstruction* theDetector;
  OrganicMaterial* theOrganicMaterial;

  G4PhantomParameterisation* param;
  G4PVParameterised* phantomPhys;

  map<G4int,G4double> hu2density;
  TFile* f;
  G4Material* mat;

  TH3D* rho;
};
#endif

