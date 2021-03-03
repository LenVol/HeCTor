#include "Calorimeter.hh"
#include "Analysis.hh"

Calorimeter::Calorimeter(G4String name):G4VSensitiveDetector(name),theName(name)
{
  theAnalysis = Analysis::GetInstance();  
}

G4bool Calorimeter::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  theAnalysis->calorimeterHit(aStep);
  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){ // don't propagate the track after the calorimeter. 
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  } 
     
  return true;
}


