#include "EmPhysics_pCT.hh"
#include "HadrontherapyPhysicsList.hh"
//#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "OrganicMaterial.hh"
#include "Analysis.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"
//#include "ParallelWorldConstruction.hh"

#include <TROOT.h>


#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4Alpha.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"
#include "G4EmCalculator.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"

#ifdef VIS
#include "G4VisExecutive.hh"

#endif
#include <iostream>
#include <fstream>
#include "time.h"
using namespace std;


void calcRSP(OrganicMaterial*);
void calcX0(DetectorConstruction*, PrimaryGeneratorAction *);


int main(int argc,char** argv) {
  gROOT->ProcessLine("#include <vector>");
  time_t seed;
  time(&seed);
  CLHEP::RanecuEngine *theRanGenerator = new CLHEP::RanecuEngine;  
  theRanGenerator->setSeed(seed);
  CLHEP::HepRandom::setTheEngine(theRanGenerator);


  ifstream fileStream;
  fileStream.open(argv[1]);
  if(!fileStream.is_open()){ cout << "Abort!" << endl; exit(0);}
  typedef std::map<std::string, std::string> ConfigInfo;
  ConfigInfo configValues;
  std::string line;
  while (std::getline(fileStream, line))
  {
    std::istringstream is_line(line);
    std::string key;
    if (std::getline(is_line, key, '='))
    {
	std::string value;
	if (key[0] == '#')
	    continue;

	if (std::getline(is_line, value))
	{
	    configValues[key] = value;
	}

    }
  } 
 
  cout << configValues["nbPrimaries"] << endl;

  G4int nParticles = stoi(configValues["nbPrimaries"]);
  G4String phantom = configValues["Phantom"];
  G4String Rifi    = configValues["RIFIFile"]; // Relative path to rifi spectrum 
  G4double angle   = stof(configValues["Angle"]); // Gantry angle (0° here is -90° in treatment plan)
  G4String name    = configValues["OFILE"]; // Name of the output file
  G4int ANumber    = stoi(configValues["A"]); // Particle to be used if 4 Helium, else Carbon
  G4double yTPS    = stof(configValues["PosY"]);
  G4double zTPS    = stof(configValues["PosZ"]);
   

/*  //Reading in the parameters from command line
  G4int nParticles = atoi(argv[1]); // Number of primaries to be used
  G4String phantom = argv[2]; // In case of CT scans, relative path to phantom.root as outputted from config.py
  G4String Rifi	   = argv[3]; // Relative path to rifi spectrum 
  G4double angle   = atof(argv[4]); // Gantry angle (0° here is -90° in treatment plan)
  G4String name    = argv[5]; // Name of the output file
  G4int ANumber    = atoi(argv[6]); // Particle to be used if 4 Helium, else Carbon
*/


  G4RunManager* runManager = new G4RunManager;

  // Physics
  runManager->SetUserInitialization(new HadrontherapyPhysicsList()); // Custom Physics list (QGSP_BIC_HP slightly modified) following the extended Hadrontherapy example of Geant4
/**************  
  G4PhysListFactory *physListFactory = new G4PhysListFactory();  // Reference phyiscs list contains wrong paramterisation of helium cross sections
  G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("QGSP_BIC_HP_EMZ"); 
  runManager->SetUserInitialization(physicsList);
***************/
    
  // Constructing the detector
  DetectorConstruction* myDC = new DetectorConstruction(phantom,angle);
  OrganicMaterial *theMaterialList = OrganicMaterial::GetInstance();

  // Beam source
  PrimaryGeneratorAction* theGenerator = new PrimaryGeneratorAction(Rifi, nParticles, ANumber, yTPS,zTPS);
  runManager->SetUserAction( theGenerator );


  // User analysis
  Analysis* theAnalysis    = new Analysis(name);
  runManager->SetUserAction( new SteppingAction() );

  runManager->SetUserInitialization( myDC );
  runManager->Initialize();
 
/**************
//Additional debug command lines
  //G4UImanager * UImanager = G4UImanager::GetUIpointer();
  //UImanager->ApplyCommand("/process/eLoss/fluct false");
  //UImanager->ApplyCommand("/run/verbose 1");
  //UImanager->ApplyCommand("/event/verbose 1");
  //UImanager->ApplyCommand("/tracking/verbose 1");
  //G4UIExecutive * ui = new G4UIExecutive(argc,argv);
  //ui->SessionStart();
***************/

//Visualisation -> Not Recommended with full CT simulations! Use one or two test-slices instead
#ifdef VIS
  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4cout << " UI session starts ..." << G4endl;
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  UImanager->ApplyCommand("/control/execute vis.mac");
  ui->SessionStart();
#endif

  runManager->BeamOn( nParticles );
  theAnalysis->Save();
  calcRSP(theMaterialList); //Calculate the RSP of the materials used in the sim
  
  delete runManager;
  return 0;
}

// To quickly calculate the stopping power of materials in the sim
void calcRSP(OrganicMaterial *theMaterial){
  G4ParticleDefinition* particle = G4Alpha::Definition(); //stoppin power defined with Alpha here
  G4EmCalculator* emCal = new G4EmCalculator;
  G4NistManager *man = G4NistManager::Instance();

  ofstream myfile;
  myfile.open ("RSPList.txt");
  for(auto itr=theMaterial->theMaterialList.begin(); itr!=theMaterial->theMaterialList.end(); itr++) {
    G4int I = 0;
    G4double tot =0;
    G4Material* water = man->FindOrBuildMaterial("G4_WATER");
    G4Material* mat = itr->second;
    for(int j=1000;j<33000;j++){
      G4double dedx_w = emCal->ComputeElectronicDEDX( double(j)/100*MeV,particle,water);
      G4double dedx_b = emCal->ComputeElectronicDEDX( double(j)/100*MeV,particle,mat);
      tot +=dedx_b/dedx_w;
      I+=1;
    }
    G4double RSP = tot/I;
    myfile<<" "<<mat->GetDensity()/(g/cm3)<<" "<<RSP<<" "<<mat->GetName()<<endl;
  }
  myfile.close();

}




