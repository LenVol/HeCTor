#include "HadrontherapyStepMax.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyStepMax::HadrontherapyStepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),MaxChargedStep(DBL_MAX)
{
}
 
/////////////////////////////////////////////////////////////////////////////
HadrontherapyStepMax::~HadrontherapyStepMax() { }

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyStepMax::IsApplicable(const G4ParticleDefinition& particle) 
{ 
  return (particle.GetPDGCharge() != 0.);
}

/////////////////////////////////////////////////////////////////////////////    
void HadrontherapyStepMax::SetMaxStep(G4double step) {MaxChargedStep = step;}

/////////////////////////////////////////////////////////////////////////////
G4double HadrontherapyStepMax::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                  G4double,
                                                  G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double ProposedStep = DBL_MAX;

  if((MaxChargedStep > 0.) &&
     (aTrack.GetVolume() != 0) &&
     (aTrack.GetVolume()->GetName() == "DetectorPhys"))
     ProposedStep = MaxChargedStep;

  return ProposedStep;
}

/////////////////////////////////////////////////////////////////////////////
G4VParticleChange* HadrontherapyStepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

