#include "SensitiveDetector.hh"
#include "Analysis.hh"

SensitiveDetector::SensitiveDetector(G4String name):G4VSensitiveDetector(name),theName(name)
{
  theAnalysis = Analysis::GetInstance();  
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary && aStep->GetTrack()->GetTrackID()==1){
    theAnalysis->analyseHit(aStep, theName);
  }
/*
  if (theName=="calorimeter"){
    theAnalysis->calorimeterHit(aStep);
    if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }
  }
*/      
  return true;
}
