#include "TFile.h"
#include "TH3D.h"
#include "TH3S.h"
#include "G4Proton.hh"
#include "G4ThreeVector.hh"
#include "G4PhantomParameterisation.hh"
#include "MyPhantomParameterisation.hh"
#include "Calorimeter.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4EmCalculator.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Region.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4VSolid.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "OrganicMaterial.hh"
#include "Analysis.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "SensitiveDetector.hh"
#include "DetectorConstruction.hh"
#include <set>
#include <stdlib.h> 

DetectorConstruction* DetectorConstruction::theDetector=NULL;
DetectorConstruction::~DetectorConstruction(){theDetector = NULL; }

DetectorConstruction::DetectorConstruction(G4String theName, G4double angle):theAngle(angle),thePhantom(theName)
{  
  theDetector = this;   
  theOrganicMaterial = OrganicMaterial::GetInstance();
  air   = theOrganicMaterial->ConstructMaterial("Air",0.001155);
  water = theOrganicMaterial->water;

  G4String rootFile = ".root";
  withCT = (thePhantom.find(rootFile) != string::npos); //check if the phantom to be simulated is a root file
  if(withCT){
	  G4cout << "Loading the CT scan from " << thePhantom.data() << "!" << G4endl;
	  f      = new TFile(Form("%s",thePhantom.data()),"update");
	  hu     = (TH3D*)f->Get("hu"); // Hounsfield unit Histogram
	  rho    = (TH3D*)f->Get("rho"); // Electron densities

	  NbinsX =  hu->GetNbinsX();
	  NbinsY =  hu->GetNbinsY();
	  NbinsZ =  hu->GetNbinsZ();
	  Xmax   =  hu->GetXaxis()->GetXmax();
	  Ymax   =  hu->GetYaxis()->GetXmax();
	  Zmax   =  hu->GetZaxis()->GetXmax();
	  Xmin   =  hu->GetXaxis()->GetXmin();
	  Ymin   =  hu->GetYaxis()->GetXmin();
	  Zmin   =  hu->GetZaxis()->GetXmin();
	  midX   =  (Xmax+Xmin)/2.;
	  midY   =  (Ymax+Ymin)/2.;
	  midZ   =  (Zmax+Zmin)/2.;

	  halfX  =  hu->GetXaxis()->GetBinWidth(1)/2.;
	  halfY  =  hu->GetYaxis()->GetBinWidth(1)/2.;
	  halfZ  =  hu->GetZaxis()->GetBinWidth(1)/2.;

	  // Calculate the shift to place the phantom according to the CT coordinate (Geant4 otherwise places everything relative to 0)
	  if(thePhantom.find("THORAX")!=string::npos){
	      shift = G4ThreeVector(0,-midX,midZ);
	      G4cout<<"Centroid shift from (0,0,0) ="<<shift<<G4endl;
	      isocenter = G4ThreeVector(0.,140.74 - 165.5, 411. - 80.1); 
	  }
	  else if (thePhantom.find("PETer")!=string::npos){
	      shift = G4ThreeVector(midX,midY,midZ);
	      G4cout<<"Centroid shift from (0,0,0) ="<<shift<<G4endl;
	      isocenter = G4ThreeVector(Xmin+245,Ymin + 249.5*mm, Zmax - 160*mm); //ADAMPETER //154.8, 245.1 // ADAMPETER use bead markers, ot isocentre!
	  }
	  else{
	      shift  =  G4ThreeVector(0,midY,midZ); // Here I wanted the phantom in beam direction to be centered on 0 
              G4cout<<"Centroid shift from (0,0,0) ="<<shift<<G4endl;
	      isocenter = G4ThreeVector(0,-472.512 + 248.5*mm, -460.1 - 169.9*mm); // ADAM Isocenter position -> Isocenter relative to min of scan (minPos + IsocenterPos) 
	  }
	  // Kept relative here for easier debugging
	  shift  = shift - isocenter;


	  G4cout<<"Alignment with Isocenter hift from (0,0,0) ="<< shift <<G4endl;
	  cout<<"NbinsX - Xmax - Xmin :"<<NbinsX<<" "<<Xmax<<" "<<Xmin<<endl;
	  cout<<"NbinsY - Ymax - Ymin :"<<NbinsY<<" "<<Ymax<<" "<<Ymin<<endl;
	  cout<<"NbinsZ - Zmax - Zmin :"<<NbinsZ<<" "<<Zmax<<" "<<Zmin<<endl;

	  //Gammex Calibration for the Scanner used
	  //vector< G4String > MaterialName  = {"Air","LN450","AP6Adipose","Breast","Water","NHMuscle3","Brain","Liver","InnerBone","MineralBone","CB230","CB250","CorticalBone"};
	  //vector< G4double > Threshold     = {0.,0.4,0.9,0.96,1.0,1.02,1.05,1.07,1.09,1.1,1.27,1.47,1.695}; // Electron densities used as thresholds in this case

	  // Materials and intervals from Huenemohr et al. 2014
          //vector< G4String > MaterialName  = {"Air","NHLungDeflated","NHAdipose3","NHAdipose1","NHMuscle1","NHMuscle3","NHCartilage","NHFemurWholeSpecimen","NHCorticalBone"};
          //vector< G4double > Threshold     = {0.,-741.,-98.,-55.,40.,44.,102.,702.,1524.}; // HU thresholds


	  // Phantom Materials
	  vector< G4String > MaterialName;
	  vector< G4double > Threshold;
	  if(thePhantom.find("THORAX")!=string::npos){
	    MaterialName  = {"ADAMAgarose2","ADAMSilicone"}; // Lung materials (PMMA overlaps with Silicone)
            Threshold     = {-30.,80.}; //HU thresholds

	  } 
	  else {
	    MaterialName  = {"ADAMPeanutOil","ADAMAgarose1","ADAMAgarose2","ADAMInnerBone","ADAMGypsum"}; // ADAM materials
	    Threshold     = {-160.,-30.,35.,80.,285.}; //HU thresholds
	 }
	  map<int,double> ids; //A map of the hu and their density counterpart //TODO: This can be done directly from CT scan, instead of the workaround through TH3D and config.py
	  for (G4int k=0;k<NbinsZ;k++) {
	    for (G4int j=0;j<NbinsY;j++) {
	      for (G4int i=0;i<NbinsX;i++) ids.insert(pair<double,double>(hu->GetBinContent(i+1,j+1,k+1),rho->GetBinContent(i+1,j+1,k+1)) ) ;
	    }
	  }

	  for(auto it=ids.begin();it!=ids.end();++it){

	    G4double huUnit   = it->first;
	    G4double density  = it->second;

	    // Find the closest material based on the density threshold in the material list
	    int lower = lower_bound(Threshold.begin(), Threshold.end(), huUnit)-Threshold.begin()-1;
	    //Construct composite material with the relevant density
	    if(lower<0) mat = air;                                   // Smaller than the smallest point
	    else{
		mat = theOrganicMaterial->ConstructMaterial(MaterialName.at(lower),density);
	    }
	    hu2id[huUnit] = theMaterialList.size(); // Indices later read in in the detector construction to assign the material to each voxel
	    huList.push_back(huUnit);
	    theMaterialList.push_back(mat);
	  }

	  G4cout << ids.size() << " materials " << G4endl;
    }
    else G4cout << "Thee PMMA phantom will be used!" << G4endl;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  G4cout << "DetectorConstruction::Construct" << G4endl;

  G4Colour* colour = new G4Colour();

  //Cubic Air world
  G4double world_size          = 3*m;  

  G4Box* hallWorld             = new G4Box("World",world_size,world_size,world_size);
  G4LogicalVolume* logicWorld  = new G4LogicalVolume(hallWorld,air,"World");                                    
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"physWorld",0,false,0);
  G4VisAttributes* hall_att    = new G4VisAttributes(G4Colour(1,0,0));
  hall_att->SetVisibility(true);
  logicWorld->SetVisAttributes(hall_att);

  //Declaration of some general variables
  G4double MylarScintDist = 2*cm; //Estimated, but without the efficiency of light output detection as function of the lateral/vertical position, this is irrelevant anyways
  G4double PMMADetecDist = 2*cm; // Distance between last slab of PMMA and the detector entrance window; 
  G4double PMMA_end = 0.; // Rear end of PMMA
   
  //PMMA Material
  G4Material* pmma_49mm = theOrganicMaterial->ConstructMaterial("pmma",1.1742); //for the 49.6mm slab
  G4Material *pmma_2mm = theOrganicMaterial->ConstructMaterial("pmma",1.148); // For the 2mm 
  G4Material *pmma_1mm = theOrganicMaterial->ConstructMaterial("pmma",1.132); // for the 1 mm 
  G4Material *pmma_5mm = theOrganicMaterial->ConstructMaterial("pmma",1.133); //for te 5 mm 
  G4Material *pmma_9mm = theOrganicMaterial->ConstructMaterial("pmma",1.170); //for the 9.7mm
  G4Material *pmma_10mm = theOrganicMaterial->ConstructMaterial("pmma",1.112); //for the 10.7mm
  G4Material *pmma_20mm = pmma_49mm;                                               // for the 19.55mm, same RSP as 49.9mm slab
  G4Material *pmma_51mm = theOrganicMaterial->ConstructMaterial("pmma",1.167); // for the 51.1mm
  G4Material *pmma_50mm = theOrganicMaterial->ConstructMaterial("pmma",1.1762);//For the 50.7mm 
  G4VisAttributes *degrader_att = new G4VisAttributes(colour->Cyan());
  degrader_att->SetVisibility(true);
  degrader_att->SetForceSolid(true);

  if( withCT){
	  //-------------------------
	  // CT scan phantom
	  //-------------------------

	  size_t* materialIDs          = new size_t[NbinsX*NbinsY*NbinsZ];
	  unsigned int p               = 0;

	  for (G4int k=0;k<NbinsZ;k++) {
	  for (G4int j=0;j<NbinsY;j++) {
	  for (G4int i=0;i<NbinsX;i++) {
	    G4int huUnit = hu->GetBinContent(i+1,j+1,k+1);
	    materialIDs[p] = hu2id[huUnit];
	    p++;
	  }
	  }
	  }


	  // Container volume
	  G4RotationMatrix* rot1;
	  if(thePhantom.find("PETer")!=string::npos){
	    rot1 = new G4RotationMatrix(); // If the angle is to be varied (e.g. for ADAM PETEr)
	    rot1->rotateY(theAngle*pi/180.);
	  }
	  else if (thePhantom.find("THORAX")!=string::npos){
	    rot1 = new G4RotationMatrix(); // If the angle is to be varied (e.g. for THORAX)
            rot1->rotateZ((90+theAngle)*pi/180.);
	  }
	  else{
            rot1 = new G4RotationMatrix(); // For ADAM no rotation is used
            rot1->rotateZ(0*pi/180.);
	  }
	  G4Box* cont_vol                      = new G4Box("cont_vol",NbinsX*halfX,NbinsY*halfY,NbinsZ*halfZ);
	  G4LogicalVolume * cont_log           = new G4LogicalVolume(cont_vol,air,"cont_log",0,0,0);

	  G4VPhysicalVolume* cont_phys         = new G4PVPlacement(rot1,shift,cont_log,"cont_phys",logicWorld,false,0);
	  G4VisAttributes* cont_att = new G4VisAttributes(G4Colour(0,1,1));
	  cont_att->SetVisibility(true);
	  cont_log->SetVisAttributes(cont_att);

	  // Phantom voxel volume 

	  G4Box* phantomVoxel_vol              = new G4Box("phantomVoxel_vol",halfX,halfY,halfZ);
	  G4LogicalVolume *phantomVoxel_log    = new G4LogicalVolume(phantomVoxel_vol,air,"phantomVoxel_log",0,0,0); //Material overwritten anyways

	  // Parameterisation of the CT phantom

	  G4cout << "Make a parameterised phantom" << G4endl;

	  param = new MyPhantomParameterisation(); // "Standard" for creating a CT geometry in Geant4, more powerful but more complicated wys possible (G4NestedParameterisation eg.)
	  param->SetVoxelDimensions(halfX,halfY,halfZ); // Phantom voxel half widths
	  param->SetNoVoxel(NbinsX,NbinsY,NbinsZ); // Phantom number of voxels
	  param->SetMaterials(theMaterialList); // Materials to be contained in the phantom
	  param->SetMaterialIndices(materialIDs); // Materials to be assigned to each voxel number
	  param->SetSkipEqualMaterials(false); //Would speed up simulation, as voxel edges can be skipped, but makes dose calculation wrong
	  param->BuildContainerSolid(cont_phys);  //Build the container
	  param->CheckVoxelsFillContainer(cont_vol->GetXHalfLength(),cont_vol->GetYHalfLength(),cont_vol->GetZHalfLength()); // Check if voxels fill the mother completely
	  phantomPhys = new G4PVParameterised("phantomPhys",phantomVoxel_log,cont_log,kUndefined,param->GetNoVoxel(),param); //kUndefined means smart voxel optimization is used
	  //phantomPhys->CheckOverlaps(); //Takes forever to do for full phantom
	  phantomPhys->SetRegularStructureId(1); //Speeds up simulation because RegNavigator is used


	  // Degrader
	  //before phantom
	  G4double degrader_0_thick = 19.57*mm;
          G4Box* degrader_0_vol = new G4Box("degrader_0_vol", degrader_0_thick/2.,15*cm/2.,15*cm/2.);
          G4LogicalVolume *degrader_0_log = new G4LogicalVolume(degrader_0_vol,pmma_20mm,"degrader_0_log",0,0,0);
          degrader_0_log->SetVisAttributes(degrader_att);

	  //after phantom
	  G4double degrader_1_thick = 10.7*mm;
          G4Box* degrader_1_vol = new G4Box("degrader_1_vol", degrader_1_thick/2.,15*cm/2.,15*cm/2.);
          G4LogicalVolume *degrader_1_log = new G4LogicalVolume(degrader_1_vol,pmma_10mm,"degrader_1_log",0,0,0);
          degrader_1_log->SetVisAttributes(degrader_att);

	  G4double degrader_2_thick = 9.7*mm;
          G4Box* degrader_2_vol = new G4Box("degrader_2_vol", degrader_2_thick/2.,15*cm/2.,15*cm/2.);
          G4LogicalVolume *degrader_2_log = new G4LogicalVolume(degrader_2_vol,pmma_9mm,"degrader_2_log",0,0,0);
          degrader_2_log->SetVisAttributes(degrader_att);

	  G4double degrader_3_thick = 19.55*mm;
          G4Box* degrader_3_vol = new G4Box("degrader_3_vol", degrader_3_thick/2.,15*cm/2.,15*cm/2.);
          G4LogicalVolume *degrader_3_log = new G4LogicalVolume(degrader_3_vol,pmma_20mm,"degrader_3_log",0,0,0);
          degrader_3_log->SetVisAttributes(degrader_att);

	  G4double degrader_4_thick = 51.10*mm;
	  G4Box* degrader_4_vol = new G4Box("degrader_4_vol", degrader_4_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_4_log = new G4LogicalVolume(degrader_4_vol,pmma_51mm,"degrader_4_log",0,0,0);
	  degrader_4_log->SetVisAttributes(degrader_att);

	  G4double degrader_5_thick = 50.7*mm;
	  G4Box* degrader_5_vol = new G4Box("degrader_5_vol", degrader_5_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_5_log = new G4LogicalVolume(degrader_5_vol,pmma_50mm,"degrader_5_log",0,0,0);
	  degrader_5_log->SetVisAttributes(degrader_att);

          G4double degrader_6_thick = 49.6*mm;
          G4Box* degrader_6_vol = new G4Box("degrader_6_vol", degrader_6_thick/2.,15*cm/2.,15*cm/2.);
          G4LogicalVolume *degrader_6_log = new G4LogicalVolume(degrader_6_vol,pmma_49mm,"degrader_6_log",0,0,0);
          degrader_6_log->SetVisAttributes(degrader_att);

	  G4double end_thick = 1*cm;
	   G4ThreeVector pos;
	  if(thePhantom.find("PETer") != string::npos){
		cout << "ADAM PETER is assumed" << endl;
		//Degrader before phantom
		pos = G4ThreeVector(-NbinsX*halfX*mm -  degrader_0_thick/2. + shift.x() - 1*cm, 0.,0.);
                new G4PVPlacement(0, pos, "build_phys",degrader_0_log,physWorld,false,0);

		//Degrader after phantom
               /* pos = G4ThreeVector(NbinsX*halfX*mm+ end_thick + degrader_1_thick/2. + shift.x(), 0.,0.);
                new G4PVPlacement(0, pos, "end_phys_1",degrader_1_log,physWorld,false,0);
		end_thick+=degrader_1_thick;
                pos = G4ThreeVector(NbinsX*halfX*mm+ end_thick + degrader_2_thick/2. + shift.x(), 0.,0.);
                new G4PVPlacement(0, pos, "end_phys_2",degrader_2_log,physWorld,false,0);
                end_thick+=degrader_2_thick;
   		pos = G4ThreeVector(NbinsX*halfX*mm+ end_thick + degrader_3_thick/2. + shift.x(), 0.,0.);
                new G4PVPlacement(0, pos, "end_phys_3",degrader_3_log,physWorld,false,0);
		end_thick+=degrader_3_thick;
		*/
		pos = G4ThreeVector(NbinsX*halfX*mm+ end_thick + degrader_4_thick/2. + shift.x(), 0.,0.);
          	new G4PVPlacement(0, pos, "end_phys_4",degrader_4_log,physWorld,false,0);
          	end_thick+=degrader_4_thick;
          	pos = G4ThreeVector(NbinsX*halfX*mm+end_thick+degrader_5_thick/2. + shift.x(), 0.,0.);
          	new G4PVPlacement(0, pos, "end_phys_5",degrader_5_log,physWorld,false,0);
          	end_thick += degrader_5_thick;
             	pos = G4ThreeVector(NbinsX*halfX*mm+end_thick+degrader_6_thick/2. + shift.x(), 0.,0.);
                new G4PVPlacement(0, pos, "end_phys_6",degrader_6_log,physWorld,false,0);
                end_thick += degrader_6_thick;

          	PMMA_end = NbinsX*halfX*mm + end_thick+ shift.x(); //Position of the last PMMA slab relative to the isocentre


	  }
	  else if(thePhantom.find("THORAX")!=string::npos){
		//No degrader used!
	        PMMA_end = NbinsY*halfY*mm;
	  }
	  else{
	  	pos = G4ThreeVector(NbinsX*halfX*mm+ degrader_5_thick/2. + shift.x(), 0.,0.); 
	  	new G4PVPlacement(0, pos, "end_phys_1",degrader_5_log,physWorld,false,0);
	  	end_thick+=degrader_5_thick;
	  	pos = G4ThreeVector(NbinsX*halfX*mm+end_thick+degrader_6_thick/2. + shift.x(), 0.,0.);
	  	new G4PVPlacement(0, pos, "end_phys_2",degrader_6_log,physWorld,false,0);
	  	end_thick += degrader_6_thick;
	  	PMMA_end = NbinsX*halfX*mm + end_thick+ shift.x(); //Position of the last PMMA slab relative to the isocentre
	  }
  }
  //MAKE THE PMMA PHANTOM
  else {
  //For the PMMA the Detecor entrance window was located in the laser isocenter (here (0.,0.,0.))
	  //-------------------------
	  // phantom
	  //-------------------------

	  //Build up
	  G4double build_thick = 49.6*mm;
	  G4Box* build_vol = new G4Box("degrader_vol", build_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *build_log = new G4LogicalVolume(build_vol,pmma_49mm,"degrader_log",0,0,0);
	  build_log->SetVisAttributes(degrader_att);

	  //Gap 
	  G4double gap_thick = 1.07*mm; //TODO MAKE THE GAP A VARIABLE
	  G4VSolid * gap_sheet = new G4Box("gap_sheet", gap_thick/2., 15*cm/2.,15*cm/2);
	  G4VSolid * gap_gap = new G4Box("gap_gap", gap_thick/2., 15*cm/2., 2*mm/2.); //TODO MAKE THE GAP A VARIABLE
	  G4SubtractionSolid* gap_vol = new G4SubtractionSolid("gap_vol",gap_sheet, gap_gap,0,G4ThreeVector(0,0,0));
	  G4LogicalVolume* gap_log = new G4LogicalVolume(gap_vol,pmma_1mm,"gap_log");
	  gap_log->SetVisAttributes(degrader_att);

	  G4double degrader_1_thick = 2.2*mm; //TODO MAKE THE GAP A VARIABLE
	  G4Box* degrader_1_vol = new G4Box("degrader_1_vol", degrader_1_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_1_log = new G4LogicalVolume(degrader_1_vol,pmma_2mm,"degrader_1_log",0,0,0);
	  degrader_1_log->SetVisAttributes(degrader_att);

	  G4double degrader_2_thick = 5.2*mm; //TODO: MAKE The GAP A VARIABLE
	  G4Box* degrader_2_vol = new G4Box("degrader_2_vol", degrader_2_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_2_log = new G4LogicalVolume(degrader_2_vol,pmma_5mm,"degrader_2_log",0,0,0);
	  degrader_2_log->SetVisAttributes(degrader_att);

	  G4double degrader_3_thick = 10.7*mm;
	  G4Box* degrader_3_vol = new G4Box("degrader_3_vol", degrader_3_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_3_log = new G4LogicalVolume(degrader_3_vol,pmma_10mm,"degrader_3_log",0,0,0);
	  degrader_3_log->SetVisAttributes(degrader_att);

	  G4double degrader_4_thick = 19.55*mm;
	  G4Box* degrader_4_vol = new G4Box("degrader_4_vol", degrader_4_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_4_log = new G4LogicalVolume(degrader_4_vol,pmma_20mm,"degrader_4_log",0,0,0);
	  degrader_4_log->SetVisAttributes(degrader_att);

	  G4double degrader_5_thick = 51.10*mm;
	  G4Box* degrader_5_vol = new G4Box("degrader_2_vol", degrader_5_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_5_log = new G4LogicalVolume(degrader_5_vol,pmma_51mm,"degrader_5_log",0,0,0);
	  degrader_5_log->SetVisAttributes(degrader_att);

	  G4double degrader_6_thick = 50.7*mm;
	  G4Box* degrader_6_vol = new G4Box("degrader_6_vol", degrader_6_thick/2.,15*cm/2.,15*cm/2.);
	  G4LogicalVolume *degrader_6_log = new G4LogicalVolume(degrader_6_vol,pmma_50mm,"degrader_6_log",0,0,0);
	  degrader_6_log->SetVisAttributes(degrader_att);

	  G4double end_thick = 0;
	  new G4PVPlacement(0, G4ThreeVector(-degrader_6_thick/2. - PMMADetecDist, 0.,0.), "end_phys_6",degrader_6_log,physWorld,false,0);
	  end_thick+=degrader_6_thick;
	  new G4PVPlacement(0, G4ThreeVector(-degrader_5_thick/2. - end_thick - PMMADetecDist, 0.,0.), "end_phys_5",degrader_5_log,physWorld,false,0);
	  end_thick+=degrader_5_thick;

	/*If the Gap were to be placed in the rear of the phantom (distal to C peak)
	  // new G4PVPlacement(0, G4ThreeVector(-end_thick - 2*cm - gap_thick/2., 0.,0.), "gap_phys",gap_log,physWorld,false,0);
	  //  end_thick+=gap_thick;
	*/
	  new G4PVPlacement(0, G4ThreeVector(-degrader_4_thick/2. - end_thick - PMMADetecDist, 0.,0.), "end_phys_4",degrader_4_log,physWorld,false,0);
	  end_thick+=degrader_4_thick;
	  new G4PVPlacement(0, G4ThreeVector(-degrader_3_thick/2. - end_thick - PMMADetecDist, 0.,0.), "end_phys_3",degrader_3_log,physWorld,false,0);
	  end_thick+=degrader_3_thick;
	  new G4PVPlacement(0, G4ThreeVector(-degrader_2_thick/2. - end_thick - PMMADetecDist, 0.,0.), "end_phys_2",degrader_2_log,physWorld,false,0);
	  end_thick+=degrader_2_thick;
	  new G4PVPlacement(0, G4ThreeVector(-degrader_1_thick/2. - end_thick - PMMADetecDist, 0.,0.), "end_phys_1",degrader_1_log,physWorld,false,0);
	  end_thick+=degrader_1_thick;
	  //Gap just before the first 5 cm block is placed (after a 5cmm build up in beam direction)
	  new G4PVPlacement(0, G4ThreeVector(-end_thick - PMMADetecDist - gap_thick/2., 0.,0.), "gap_phys",gap_log,physWorld,false,0);
	  end_thick+=gap_thick;
	  //Build up beofre gap
	  new G4PVPlacement(0, G4ThreeVector(-end_thick - PMMADetecDist - build_thick/2., 0.,0.), "build_phys",build_log,physWorld,false,0);
	  PMMA_end = -PMMADetecDist; // By construction the PMMA slabs end at -2*cm
}

  //-------------------------
  // Ideal trackers for single event tests/particle imaging
  //-------------------------
  G4double HalfY  = 30*cm; //NbinsY*halfY;
  G4double HalfZ  = 30*cm; //NbinsZ*halfZ;
  G4double thickTracker = 1*mm; //will contain air anyways

  G4Box* sd_vol                      = new G4Box("rad_vol2",thickTracker/2.,HalfY,HalfZ);
  G4LogicalVolume * sd_log           = new G4LogicalVolume(sd_vol,air,"sd_log",0,0,0); // air tracker  

  G4VisAttributes* sd_att = new G4VisAttributes(G4Colour(1,1,1));
  sd_att->SetVisibility(false);
  sd_log->SetVisAttributes(sd_att);

  G4ThreeVector posSD = G4ThreeVector(PMMA_end + thickTracker/2.,0.,0.); //directly at end of the Phantom (for now) TODO: Make this a parallel geometry for simpler phantom rotation
  new G4PVPlacement(0,posSD,"sd_phys",sd_log,physWorld,false,0);
  SensitiveDetector* sd2               = new SensitiveDetector("sd2");
  sd_log->SetSensitiveDetector(sd2);

  //------------------------
  // Calorimeter
  //------------------------
  // Calorimeter entrance at isocenter
  G4double calThickX = 12*cm; // Approximate Thickness of UCL telescope (anything after 12 cm is not relevant anyways)
  G4double calWidth  = 10*cm; // Transversal Aperture

  
  CalorimeterEntrance = PMMA_end+PMMADetecDist+MylarScintDist;
  

  //Mylar window
  G4Material* mylar             = theOrganicMaterial->ConstructMaterial("mylar",1.3925); // Mylar density taken from pCT simualtion (Giacometti et al. 2017 "software platform.." Med. Phys.)
  G4Box* mylar_vol              = new G4Box("mylar_vol",0.01*mm/2.,calWidth/2.,calWidth/2.);
  G4LogicalVolume* mylar_log    = new G4LogicalVolume(mylar_vol,mylar,"mylar_log",0,0,0);
  G4VisAttributes* mylar_att    = new G4VisAttributes(colour->Black());
  mylar_att->SetVisibility(true);
  mylar_att->SetForceSolid(true);
  mylar_log->SetVisAttributes(mylar_att);
 
  G4ThreeVector posMylar = G4ThreeVector(PMMA_end+PMMADetecDist,0.,0.);
  new G4PVPlacement(0, posMylar, "mylar_phys",mylar_log,physWorld,false,0);

  //Sensitive detector
  G4Material* polystyrene             = theOrganicMaterial->ConstructMaterial("polystyrene",1.0325); //Density optmised to fit RSP of detector
  G4Box* cal_vol                      = new G4Box("cal_vol",calThickX/2.,calWidth/2.,calWidth/2.);
  G4LogicalVolume * cal_log           = new G4LogicalVolume(cal_vol,polystyrene,"cal_log",0,0,0); // water calorimeter   

  // Fixing the step
  G4double stepLimitValueRange = .05*mm;
  cal_log->SetUserLimits(new G4UserLimits(stepLimitValueRange));

  //Make this a region with finer cuts (sets 0.5mm production)
  G4Region* energy_calorimeter = new G4Region("energy_calorimeter");
  energy_calorimeter->AddRootLogicalVolume(cal_log);


  G4VisAttributes* cal_att = new G4VisAttributes(colour->Yellow());
  cal_att->SetVisibility(true);
  cal_att->SetForceSolid(true);
  cal_log->SetVisAttributes(cal_att);

  G4ThreeVector posCal  = G4ThreeVector(CalorimeterEntrance+calThickX/2.,0.,0.);
  new G4PVPlacement(0,posCal,"cal_phys",cal_log,physWorld,false,0);


  Calorimeter* calorimeter = new Calorimeter("calorimeter");
  cal_log->SetSensitiveDetector(calorimeter);

  


  return physWorld;

}




