#ifndef EventAction_hh
#define EventAction_hh

#include "G4UserEventAction.hh"
#include <map>

class DetectorConstruction;
class SteppingAction;

class EventAction : public G4UserEventAction {
public:
  EventAction();
  ~EventAction();
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  static inline EventAction* GetInstance() { return theEventAction; }

private:
  static EventAction* theEventAction;
  DetectorConstruction* theDetector;
  Analysis * theAnalysis;
  SteppingAction* theSteppingAction;
  map<G4int,G4double>::iterator it,it_end;

};

#endif

