//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B2aDetectorConstruction.cc 101658 2016-11-21 09:00:41Z gcosmo $
//
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class
 
#include "MuDDetectorConstruction.hh"
#include "MuDDetectorMessenger.hh"
#include "MuDTrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* B2aDetectorConstruction::fMagFieldMessenger = 0;

B2aDetectorConstruction::B2aDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfSiPM(0),
 fLogicTarget(NULL), fLogicStopper(NULL), fLogicSiPM(NULL),
 fTargetMaterial(NULL), fSiPMMaterial(NULL), fStopperMaterial(NULL),
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new B2aDetectorMessenger(this);

  fNbOfSiPM = 3;
  fLogicSiPM = new G4LogicalVolume*[fNbOfSiPM];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete [] fLogicSiPM; 
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2aDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::DefineMaterials()
{
 // Material definition 

 G4NistManager* nistManager = G4NistManager::Instance();

 // Air defined using NIST Manager
 nistManager->FindOrBuildMaterial("G4_AIR");

 // Led defined using NIST Manager  
 fStopperMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");

 // Aluminium defined using NIST Manager
 fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Al");

 // Silicon defined using NIST Manager
 fSiPMMaterial = nistManager->FindOrBuildMaterial("G4_Si");

 // Germanium defined using NIST Manager
 fGeDetMaterial = nistManager->FindOrBuildMaterial("G4_Ge");

 // vacuum and air defined using NIST Manager
 fVacuum  = new G4Material("vacuum", 
		1, 1.008*g/mole,1.e-25*g/cm3,kStateGas,300*kelvin,3.e-18*pascal); 
 fAir  = G4Material::GetMaterial("G4_AIR");

 // Print materials
 G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{
 //-------------------------visualization attributes-----------------------------------
 G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
 G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
 G4VisAttributes* geDetVisAtt = new G4VisAttributes(G4Colour(1,1,1));
 G4VisAttributes* deadlayerVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
 G4VisAttributes* alVisAtt = new G4VisAttributes(G4Colour(.8,1,1));
 G4VisAttributes* vacuumVisAtt = new G4VisAttributes(G4Colour(0,0,0));
 geDetVisAtt->SetForceSolid(true);
 vacuumVisAtt->SetForceSolid(true);
 deadlayerVisAtt->SetForceSolid(true);
 alVisAtt->SetForceSolid(true);


 //------------------------------------------------------------------------------------
 //--------------------------parameters of all solids----------------------------------
 // o = outer, i = inner
 //Silicon Photomultipliers (SiPM)
 G4double siPMXY = 4.5*cm;
 G4double siPMZ  = 1.2*cm;

 //for all tubs
 G4double startAngle = 0.*deg;
 G4double spanningAngle = 360.*deg;

 //Target
 G4double tarRadius = 3.*cm;
 G4double tarLength = 3.*cm;

 //Stopper
 G4double stopoRadius = 4.*cm;
 G4double stopLength = 9.*cm;

 //Germaniumdetector (GeDet)
 G4double geDetHoleRadius = 0.5 *cm;
 G4double geDetHoleLength = 4.65 *cm;
 G4double geDetCrystaloRadius = 3.05 *cm;
 G4double geDetCrystalLength = 6.45 *cm;
 
 //deadlayer of GeDet 
 G4double geDetIDThickness = .00003 *cm;
 G4double geDetIDoRadius = geDetHoleRadius + geDetIDThickness;
 G4double geDetIDLength = geDetHoleLength + geDetIDThickness;
 G4double geDetODThickness = .07 *cm;	
 
 //active crystal of GeDet
 G4double geDetActiveoRadius = geDetCrystaloRadius - geDetODThickness;
 G4double geDetActiveLength = geDetCrystalLength - geDetODThickness;
 
 //inner Al shield of GeDet
 G4double geDetISThickness = 0.1 *cm;
 G4double geDetISDistance = 0.02 *cm;
 G4double geDetISLength = geDetCrystalLength + 2*geDetISThickness + .005 *cm;
 G4double geDetISiRadius = geDetCrystaloRadius + geDetISDistance;
 G4double geDetISoRadius = geDetISiRadius + geDetISThickness;

 //outer Al shield of GeDet
 G4double geDetOSThickness = 0.15 *cm;
 G4double geDetOSDistance = 0.6 *cm;
 G4double geDetOSiRadius = geDetISoRadius + geDetOSDistance;
 G4double geDetOSoRadius = geDetOSiRadius + geDetOSThickness;
 G4double geDetOSLength = geDetISLength + 2*geDetISThickness + geDetOSDistance;

 // World
 G4double distanceWorldStopper = 25 *cm - stopLength;
 G4double worldXY = 22 *cm;
 G4double worldZ  = 1.1*(stopLength + distanceWorldStopper + tarLength + 3*siPMZ);

 //-------------------------positons of all solids-------------------------------------
 //Stopper
 G4ThreeVector posStop = G4ThreeVector(0, 0,
			   -0.45*worldZ + .5*stopLength + distanceWorldStopper);
 //Detectors
 G4ThreeVector posSiPM0 = posStop + G4ThreeVector(0, 0, .5*stopLength + .5*siPMZ);
 G4ThreeVector posSiPM1 = posSiPM0 + G4ThreeVector(0, 0, siPMZ);   
 G4ThreeVector posSiPM2 = posSiPM1 + G4ThreeVector(0, 0, siPMZ + tarLength); 
 //Target
 G4ThreeVector posTar = posSiPM1 + G4ThreeVector(0, 0, .5*(siPMZ + tarLength)); 
 //Ge-Det
 G4ThreeVector posGeDetOS = posTar + G4ThreeVector(0, 0.5*geDetOSLength + tarRadius, 0);
 //all following positions relative to the mother volume on line above 
 G4ThreeVector posGeDetOuterVacuum = G4ThreeVector(0, 0, -0.5*geDetOSThickness);
 G4ThreeVector posGeDetIS = G4ThreeVector(0, 0, -0.5*geDetOSDistance);
 G4ThreeVector posGeDetInnerVacuum = G4ThreeVector(0, 0, 0);
 G4ThreeVector posGeDetISHole = G4ThreeVector(0, 0,
				-0.5*(geDetISLength-geDetISThickness));
 G4ThreeVector posGeDetCrystal = G4ThreeVector(0, 0, 0);
 G4ThreeVector posGeDetActive = G4ThreeVector(0, 0, -0.5*geDetODThickness);
 G4ThreeVector posGeDetID = G4ThreeVector(0, 0, -0.5*(geDetActiveLength-geDetIDLength));
 G4ThreeVector posGeDetHole = G4ThreeVector(0, 0, -0.5*(geDetIDThickness));

 //---------------------------------debug option---------------------------------------
 SetCheckOverlaps(true);
 bool debug = false;
 G4double displacement = 0.0001 *mm;
 G4ThreeVector offset = G4ThreeVector(0,0,0);
 if (debug)
 {
	spanningAngle = 180. *deg;
	offset = G4ThreeVector(0,-0.1*cm,-0.*cm);
 }

 //------------------------------------------------------------------------------------
 // --------------Definitions of Solids, Logical Volumes, Physical Volumes-------------

 //----------------------------------World---------------------------------------------
 G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldZ);

 G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

 G4Box* worldS
    = new G4Box("World",                                    //its name
                0.5*worldXY, 0.5*worldXY, 0.5*worldZ); //its size


 G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 fAir,     //its material
                 "World"); //its name
 worldLV->SetVisAttributes(boxVisAtt);

 //  Must place the World Physical volume unrotated at (0,0,0).
 G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 


 //----------------------------------Target--------------------------------------------
 G4Tubs* targetS
    = new G4Tubs("Target",0.,tarRadius,0.5*tarLength,startAngle,spanningAngle);
 fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
 fLogicTarget ->SetVisAttributes(alVisAtt);

 new G4PVPlacement(0,posTar,fLogicTarget,"Target",worldLV,false,0,fCheckOverlaps);

 G4cout << "Target is " << tarLength/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;

/*
 //----------------------------------Stopper-------------------------------------------
 G4Tubs* stopS  
    = new G4Tubs("Stopper",
                 0, 
                 stopoRadius,
                 0.5*stopLength,
                 startAngle, 
                 spanningAngle);

 fLogicStopper                        
     = new G4LogicalVolume(stopS,fStopperMaterial,"Stopper");
            
 new G4PVPlacement(0,posStop,fLogicStopper,"Stopper",worldLV,false,0,fCheckOverlaps);          //overlaps checking
*/
  

 //----------------------------------SiPM----------------------------------------------
 G4cout << "There are " << fNbOfSiPM << " SiPM in the tracker region. " << G4endl;
  
 G4Box* detS = new G4Box("Detector_S",0.5*siPMXY, 0.5*siPMXY, 0.5*siPMZ);

 for (G4int copyNo=0; copyNo<fNbOfSiPM; copyNo++)
 {
      fLogicSiPM[copyNo] =
	new G4LogicalVolume(detS,fSiPMMaterial,"SensitiveDetector_LV",0,0,0);
      fLogicSiPM[copyNo]->SetVisAttributes(chamberVisAtt);
 }

 /*
 //Test
 G4Box* detS2 = 
  new G4Box("Detector_S",                       //its name
       siPMXY+0.5*cm, siPMXY+0.5*cm, 0.5*siPMZ);     //its size
 fLogicSiPM[2] =
	new G4LogicalVolume(detS2,fSiPMMaterial,"SensitiveDetector_LV",0,0,0);
 */

 // Detector 0
 new G4PVPlacement(0,posSiPM0,fLogicSiPM[0],"Detector",worldLV,false,0,
								fCheckOverlaps);
 // Detector 1
     new G4PVPlacement(0,posSiPM1,fLogicSiPM[1],"Detector",worldLV,false,1,
								fCheckOverlaps);          

  // Detector 2
    new G4PVPlacement(0,posSiPM2,fLogicSiPM[2],"Detector",worldLV,false,2,
								fCheckOverlaps);         
         

 //-----------------------------------GeDet--------------------------------------------
 //build layer by layer from the outmost to the inmost
 //each layer is the daughter volume of the previous one

 //outer shielding
 G4Tubs* GeDetOS = new G4Tubs("GeDetOuterAlCapS", 0.,geDetOSoRadius,
			 0.5*geDetOSLength,startAngle,spanningAngle);

 fLogicGeDetOS  = new G4LogicalVolume(GeDetOS,fTargetMaterial,"GeDetOS");
 fLogicGeDetOS->SetVisAttributes(alVisAtt); 
 
 //rotation of the outmost layer and therfore all following layers (daughter volumes)
 //init rotation matrix as 1I-matrix
 G4RotationMatrix rotm  = G4RotationMatrix();
 //rotate rotation matrix 90° around the x-axis  
 rotm.rotateX(90*deg); 
 //rotation and position in one variable
 G4Transform3D transformGeDet = G4Transform3D(rotm,posGeDetOS);

 new G4PVPlacement(transformGeDet,fLogicGeDetOS, "GeDetOS", worldLV,
						false,0,fCheckOverlaps); 

 //vacuum between outer and inner shielding
 G4Tubs* GeDetOuterVacuum = new G4Tubs("GeDetOuterVacuum", 0.,geDetOSiRadius,
				 0.5*(geDetOSLength-geDetOSThickness),
				 startAngle,spanningAngle);

 fLogicGeDetOuterVacuum  = 
		new G4LogicalVolume(GeDetOuterVacuum,fVacuum,"GeDetOuterVacuum");
 fLogicGeDetOuterVacuum->SetVisAttributes(vacuumVisAtt); 

 new G4PVPlacement(0,posGeDetOuterVacuum,fLogicGeDetOuterVacuum,
		"GeDetOuterVacuum",fLogicGeDetOS,false,0,fCheckOverlaps); 

 //inner shielding
 G4Tubs* GeDetIS = new G4Tubs("GeDetIS", 0.,geDetISoRadius,0.5*geDetISLength,
				 startAngle,spanningAngle);
 fLogicGeDetIS  = new G4LogicalVolume(GeDetIS,fTargetMaterial,"GeDetIS");
 fLogicGeDetIS->SetVisAttributes(alVisAtt); 

 new G4PVPlacement(0,offset+posGeDetIS,fLogicGeDetIS,
		"GeDetIS",fLogicGeDetOuterVacuum,false,0,fCheckOverlaps); 

 //vacuum between inner shielding and crystal
 G4Tubs* GeDetInnerVacuum = new G4Tubs("GeDetInnerVacuum", 0.,geDetISiRadius,
				 0.5*geDetISLength-geDetISThickness,
				 startAngle,spanningAngle);

 fLogicGeDetInnerVacuum  = 
		new G4LogicalVolume(GeDetInnerVacuum,fVacuum,"GeDetInnerVacuum");
 fLogicGeDetInnerVacuum->SetVisAttributes(vacuumVisAtt); 

 new G4PVPlacement(0,posGeDetInnerVacuum,fLogicGeDetInnerVacuum,
		"GeDetInnerVacuum",fLogicGeDetIS,false,0,fCheckOverlaps); 
 
 //hole for cooling finger in the inner shielding
 G4Tubs* GeDetISHole = new G4Tubs("GeDetISHole", 0.,geDetHoleRadius,
				 0.5*geDetISThickness,
				 startAngle,spanningAngle);

 fLogicGeDetISHole  = 
		new G4LogicalVolume(GeDetISHole,fVacuum,"GeDetISHole");
 fLogicGeDetISHole->SetVisAttributes(vacuumVisAtt); 

 new G4PVPlacement(0,posGeDetISHole,fLogicGeDetISHole,
		"GeDetISHole",fLogicGeDetIS,false,0,fCheckOverlaps); 

 //crystal = outer deadlayer
 G4Tubs* GeDetCrystal = new G4Tubs("GeDetCrystal", 0.,geDetCrystaloRadius,
				 0.5*geDetCrystalLength,
				 startAngle,spanningAngle);

 fLogicGeDetCrystal  = 
		new G4LogicalVolume(GeDetCrystal,fGeDetMaterial,"GeDetCrystal");
 fLogicGeDetCrystal->SetVisAttributes(deadlayerVisAtt); 

 new G4PVPlacement(0,offset+posGeDetCrystal,fLogicGeDetCrystal,
		"GeDetCrystal",fLogicGeDetInnerVacuum,false,0,fCheckOverlaps); 

 //active crystal
 G4Tubs* GeDetActive = new G4Tubs("GeDetActive", 0.,geDetActiveoRadius,
				 0.5*geDetActiveLength,
				 startAngle,spanningAngle);

 fLogicGeDetActive  = 
		new G4LogicalVolume(GeDetActive,fGeDetMaterial,"SensitiveDetector_LV");
 fLogicGeDetActive->SetVisAttributes(geDetVisAtt); 

 new G4PVPlacement(0,offset+posGeDetActive,fLogicGeDetActive,
		"GeDetActive",fLogicGeDetCrystal,false,3,fCheckOverlaps);

 //inner deadlayer
 G4Tubs* GeDetID = new G4Tubs("GeDetID", 0.,geDetIDoRadius,
				 0.5*geDetIDLength,
				 startAngle,spanningAngle);

 fLogicGeDetID  = new G4LogicalVolume(GeDetID,fGeDetMaterial,"GeDetID");
 fLogicGeDetID->SetVisAttributes(deadlayerVisAtt); 

 new G4PVPlacement(0,offset+posGeDetID,fLogicGeDetID,
		"GeDetID",fLogicGeDetActive,false,0,fCheckOverlaps);
 
 //hole for colling finger in the crystal
 G4Tubs* GeDetHole = new G4Tubs("GeDetHole", 0.,geDetHoleRadius,
				 0.5*geDetHoleLength,
				 startAngle,spanningAngle);

 fLogicGeDetHole  = new G4LogicalVolume(GeDetHole,fVacuum,"GeDetHole");
 fLogicGeDetHole->SetVisAttributes(vacuumVisAtt); 

 new G4PVPlacement(0,offset+posGeDetHole,fLogicGeDetHole,
		"GeDetHole",fLogicGeDetID,false,0,fCheckOverlaps); 



 //-----------------------Example of user step limits----------------------------------
 //
 // Below is an example of how to set tracking constraints in a given
 // logical volume
 //
 // Sets a max step length in the tracker region, with G4StepLimiter

 G4double maxStep = 0.5*siPMXY;
 fStepLimit = new G4UserLimits(maxStep);
 fLogicSiPM[0]->SetUserLimits(fStepLimit);
 fLogicSiPM[1]->SetUserLimits(fStepLimit);
 fLogicSiPM[2]->SetUserLimits(fStepLimit);
 fLogicGeDetActive->SetUserLimits(fStepLimit);
 /// Set additional contraints on the track, with G4UserSpecialCuts
 ///
 /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
 /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
 ///                                           maxLength,
 ///                                           maxTime,
 ///                                           minEkin));


 // Always return the physical world
 return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "SensitiveDetector_LV".
  SetSensitiveDetector("SensitiveDetector_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fSiPMMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fSiPMMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfSiPM; copyNo++) {
            if (fLogicSiPM[copyNo]) fLogicSiPM[copyNo]->
                                               SetMaterial(fSiPMMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
