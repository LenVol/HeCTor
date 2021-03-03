#include <vector>

#include "G4Types.hh"
#include "G4VPVParameterisation.hh"
#include "G4AffineTransform.hh"
#include "G4PhantomParameterisation.hh"
#include "G4Material.hh"

class G4VPhysicalVolume;
class G4VTouchable; 
class G4VSolid;
class G4Material;

// Dummy forward declarations ...

class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Polycone;
class G4Polyhedra;

class MyPhantomParameterisation : public G4PhantomParameterisation
{
  public: 
  
    virtual G4Material* ComputeMaterial(const G4int repNo, 
                                              G4VPhysicalVolume *currentVol,
                                        const G4VTouchable *parentTouch=0);
};


