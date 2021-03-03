#ifndef OrganicMaterial_hh
#define OrganicMaterial_hh
#include "G4ElementVector.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <set>

using namespace std;
class G4Material;
class G4NistManager;
class OrganicMaterial {
public :
  OrganicMaterial();
  ~OrganicMaterial();
  
  G4Material* ConstructMaterial(G4String Name,G4double density);

  G4Material* ConstructCompositeMaterial(G4String, G4String, G4double,G4double);
  G4double ConvertToMassDensity( std::map<G4Element*, G4double> , G4double);
  static inline OrganicMaterial* GetInstance() { 
    if(!theMaterial) theMaterial = new OrganicMaterial();
    return theMaterial; }

  G4Material* mat;
  G4Material* water;
  G4Element* elH, *elC, *elN, *elO, *elF, *elNa, *elMg, *elSi, *elP, *elS, *elCl, *elAr, *elK, *elCa, *elFe, *elZn, *elBa;
  std::map<G4String,G4Material*> theMaterialList;

private :

  G4NistManager* man;
  static OrganicMaterial* theMaterial;

  G4int nel;
};

#endif

