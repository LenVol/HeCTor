#include "TH3S.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4RegularNavigationHelper.hh"
#include "SteppingAction.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Analysis.hh"
#include "DetectorConstruction.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"

SteppingAction *SteppingAction::theSteppingAction=NULL;
SteppingAction::~SteppingAction() 
{theSteppingAction=NULL;}


SteppingAction::SteppingAction()
{
  theSteppingAction=this;
  theAnalysis = Analysis::GetInstance();
  theDetector = DetectorConstruction::GetInstance();
  cubeVolume = (theDetector->halfX*mm*theDetector->halfY*mm*theDetector->halfZ*mm)*8;
  G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetPushVerbosity(0); // to avoid the waring about stuck tracks which is common for CT geometries (see SLAC)
}


void SteppingAction::UserSteppingAction(const G4Step* aStep) // Compute the dose to each voxel on the step level // TODO: Merge this with the analysis class
{

 G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
 G4VPVParameterisation* physParam = physVol->GetParameterisation(); //check if we are in a parameterized volume

 if(physParam){ //score dose only in parameterized volume
  edep = aStep->GetTotalEnergyDeposit(); 
  density = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetMaterial()->GetDensity();
  dose    = edep / ( density*cubeVolume); 
  //dose	  = dose*aStep->GetPreStepPoint()->GetWeight(); //Used in G4DoseDeposit3D primitive scorer (not sure what it does)
  pos = aStep->GetPreStepPoint()->GetPosition();
  theAnalysis->DoseAnalysis(pos[0],pos[1],pos[2],dose/gray);
 }

}

  



