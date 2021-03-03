#include "MyPhantomParameterisation.hh"
#include "G4PhantomParameterisation.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

#include "globals.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VVolumeMaterialScanner.hh"
#include "G4GeometryTolerance.hh"
#include "G4Colour.hh"

G4Material* MyPhantomParameterisation::ComputeMaterial(const G4int repNo,
	                                                G4VPhysicalVolume *currentVol,
	                                          const G4VTouchable *parentTouch)  {

	G4Material* material = G4PhantomParameterisation::ComputeMaterial(repNo, currentVol, parentTouch);

	static G4VisAttributes* visAtts1 = new G4VisAttributes(G4VisAttributes::Invisible);
	static G4VisAttributes* visAtts2 = new G4VisAttributes(G4Colour::Grey()); //Lung Material
	static G4VisAttributes* visAtts3 = new G4VisAttributes(G4Colour::Magenta()); //low density tissue
	static G4VisAttributes* visAtts4 = new G4VisAttributes(G4Colour::Red()); // high density tissue
	static G4VisAttributes* visAtts5 = new G4VisAttributes(G4Colour::Yellow); //light bones
	static G4VisAttributes* visAtts6 = new G4VisAttributes(G4Colour::White()); //Bones
	visAtts2->SetForceSolid(true);
	visAtts3->SetForceSolid(true);
	visAtts4->SetForceSolid(true);
	visAtts5->SetForceSolid(true);
	visAtts6->SetForceSolid(true);

	G4LogicalVolume* log = currentVol->GetLogicalVolume();
	G4Material* mat = this->GetMaterial(repNo);

	if(mat->GetDensity()<=0.292*CLHEP::g/CLHEP::cm3) {
	   log->SetVisAttributes(visAtts1); // Air

	} else if (mat->GetDensity()<=0.483*CLHEP::g/CLHEP::cm3){
	   log->SetVisAttributes(visAtts2); // Lung
	} else if (mat->GetDensity()<=1.000*CLHEP::g/CLHEP::cm3){
	   log->SetVisAttributes(visAtts3); // light soft tissue
	} else if (mat->GetDensity()<=1.099*CLHEP::g/CLHEP::cm3){
	   log->SetVisAttributes(visAtts4); // dense soft tissue
	} else if (mat->GetDensity()<=1.285*CLHEP::g/CLHEP::cm3){
	   log->SetVisAttributes(visAtts5); // light bones
	} else{
	   log->SetVisAttributes(visAtts6); // dense bones
	}
	return material;

} 

