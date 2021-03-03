#ifndef SENSITIVEDETECTOR_HH_
#define SENSITIVEDETECTOR_HH_


#include "G4VSensitiveDetector.hh"


#include "G4Step.hh"
#include "G4Track.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "Analysis.hh"
class Analysis;
class SensitiveDetector : public G4VSensitiveDetector
{

public:
  SensitiveDetector(G4String);
  virtual ~SensitiveDetector(){};
  G4String theName;
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent* ){};
  void clear(){};
  void DrawAll(){};
  void PrintAll(){};
private:
  Analysis* theAnalysis;


};




#endif 
