#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "G4PhantomParameterisation.hh"
#include "SteppingAction.hh"
#include "TH3D.h"
#include "EventAction.hh"


EventAction *EventAction::theEventAction=NULL;
EventAction::EventAction() {
  theDetector = DetectorConstruction::GetInstance();
  theSteppingAction = SteppingAction::GetInstance();
  theAnalysis = Analysis::GetInstance();
}

EventAction::~EventAction() {
  theEventAction=NULL;
}

void EventAction::BeginOfEventAction(const G4Event*) {
  theSteppingAction->Edep.clear();
}

void EventAction::EndOfEventAction(const G4Event*) {

}
