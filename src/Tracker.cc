#include "Tracker.hh"
#include "Analysis.hh"

Tracker::Tracker(G4String name):G4VSensitiveDetector(name),theName(name)
{
  theAnalysis = Analysis::GetInstance();  
}

G4bool Tracker::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary && aStep->GetTrack()->GetTrackID()==1){
    theAnalysis->analyseHit(aStep, theName);
    //aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }
        
  return true;
}
