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
 fNbOfChambers(0),
 fLogicTarget(NULL), fLogicStopper(NULL), fLogicChamber(NULL),
 fTargetMaterial(NULL), fChamberMaterial(NULL), fStopperMaterial(NULL),
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new B2aDetectorMessenger(this);

  fNbOfChambers = 3;
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete [] fLogicChamber; 
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
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Si");

  // Germanium defined using NIST Manager
  fGeDetMaterial = nistManager->FindOrBuildMaterial("G4_Ge");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{

  // Visualization attributes
  G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* geDetVisAtt = new G4VisAttributes(G4Colour(.7,0.7,.7));
  G4VisAttributes* deadlayerVisAtt = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes* alVisAtt = new G4VisAttributes(G4Colour(.8,1,1));
  geDetVisAtt->SetForceSolid(true);
  deadlayerVisAtt->SetForceSolid(true);
  alVisAtt->SetForceSolid(true);

  //parameters for all solids
  //Si-Detectors
  G4double det_sizeXY = 4.5*cm;
  G4double det_sizeZ  = 1.2*cm;

  //Target
  G4double tar_radius = 3.*cm;
  G4double tar_hz = 3.*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;

  //Stopper
  G4double stop_radius = 4.*cm;
  G4double stop_hz = 20.*cm;

  //Germaniumdetector (GeDet)
  G4double geDetHole_radius = 0.6 *cm;
  G4double geDetHole_hz = 4.7 *cm;
  G4double geDet_radius = 2.88 *cm;
  G4double geDet_hz = 6.2 *cm;

  //deadlayer of GeDet 
  G4double deadlayer_hz = .12 *cm;
  G4double deadlayer_radius = geDet_radius + deadlayer_hz;

  //AlCap of GeDet
  G4double alCap_dist = 0.1 *cm;
  G4double alCapF_hz = 0.1 *cm;
  G4double alCapS_hz = geDet_hz + deadlayer_hz + alCap_dist + .5 *cm;
  G4double alCap_innerRadius = deadlayer_radius + alCap_dist;
  G4double alCap_outerRadius = deadlayer_radius + alCap_dist + alCapF_hz;

   // World
  G4double dist_z = 25 *cm - stop_hz;
  G4double world_sizeXY = 22 *cm;
  G4double world_sizeZ  = 1.1*(stop_hz + dist_z + tar_hz + 3*det_sizeZ);

  //all positions
  //Stopper
  G4ThreeVector pos_stop = G4ThreeVector(0, 0,
			   -0.45*world_sizeZ + .5*stop_hz + dist_z);
  //Detectors
  G4ThreeVector pos0 = pos_stop + G4ThreeVector(0,0,
		        .5*stop_hz + .5*det_sizeZ);
  G4ThreeVector pos1 = pos0 + G4ThreeVector(0, 0, det_sizeZ);   
  G4ThreeVector pos2 = pos1 + G4ThreeVector(0, 0, det_sizeZ + tar_hz); 
  //Target
  G4ThreeVector pos_tar = pos1 + G4ThreeVector(0, 0, .5*(det_sizeZ + tar_hz)); 
  //Ge-Detector
  G4ThreeVector pos_alCap = pos_tar + 
		G4ThreeVector(0, .5*alCapF_hz + tar_radius, 0);
  G4ThreeVector pos_deadlayer = pos_alCap + 
                G4ThreeVector(0, alCap_dist + .5*(deadlayer_hz+alCapF_hz), 0);
  G4ThreeVector pos_geDet = pos_deadlayer + 
                G4ThreeVector(0, .5*(deadlayer_hz + geDet_hz), 0);

  // Definitions of Solids, Logical Volumes, Physical Volumes
  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(world_sizeZ);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  G4Box* worldS
    = new G4Box("World",                                    //its name
                0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size

 G4Material* air  = G4Material::GetMaterial("G4_AIR");
 G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 air,      //its material
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


// Target
  G4Tubs* targetS
    = new G4Tubs("Target",0.,tar_radius,0.5*tar_hz,startAngle,spanningAngle);
  fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
  fLogicTarget ->SetVisAttributes(alVisAtt);

  new G4PVPlacement(0,               // no rotation
                    pos_tar,         // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << tar_hz/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;
  
//Stopper
  G4Tubs* stopS  
    = new G4Tubs("Stopper",
                 0, 
                 stop_radius,
                 0.5*stop_hz,
                 startAngle, 
                 spanningAngle);

    fLogicStopper                        
     = new G4LogicalVolume(stopS,         //its solid
                        fStopperMaterial,          //its material
                        "Stopper");           //its name
            
  new G4PVPlacement(0,                       //no rotation
                    pos_stop,                    //at position
                    fLogicStopper,             //its logical volume
                    "Stopper",                //its name
                    worldLV,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking

  
  

 
 
 
 
  // Tracker segments

  G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
         << G4endl;
 
  
  //Si-Detectors
  //Solid
G4Box* detS =    
    new G4Box("Detector_S",                       //its name
       0.5*det_sizeXY, 0.5*det_sizeXY, 0.5*det_sizeZ);     //its size

 //Logical Volumes
for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++)
{
      fLogicChamber[copyNo] =
	new G4LogicalVolume(detS,fChamberMaterial,"Chamber_LV",0,0,0);
      fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);
}
/*
//Test
G4Box* detS2 = 
  new G4Box("Detector_S",                       //its name
       det_sizeXY+0.5*cm, det_sizeXY+0.5*cm, 0.5*det_sizeZ);     //its size
fLogicChamber[2] =
	new G4LogicalVolume(detS2,fChamberMaterial,"Chamber_LV",0,0,0);
*/
  //Placements
  // Detector 0
  new G4PVPlacement(0,                       //no rotation
                    pos0,                    //at position
                    fLogicChamber[0],             //its logical volume
                    "Detector",                //its name
                    worldLV,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking
 // Detector 1
     new G4PVPlacement(0,                     
                    pos1,                  
                    fLogicChamber[1],            
                    "Detector",              
                    worldLV,                
                    false,                  
                    1,                      
                    fCheckOverlaps);          

  // Detector 2
    new G4PVPlacement(0,                      
                    pos2,                   
                    fLogicChamber[2],         
                    "Detector",                
                    worldLV,               
                    false,                
                    2,                       
                    fCheckOverlaps);         
         
//Germaniumdetektor

 //init rotation matrix as 1I-matrix
  G4RotationMatrix rotm  = G4RotationMatrix();
  //rotate rotation matrix 90° around the x-axis  
  rotm.rotateX(90*deg); 
     
  //rotation and position in one variable
  G4Transform3D transformGeDet = G4Transform3D(rotm,pos_geDet);

  //massiv cylinder for the active part of the GeDet
  G4Tubs* GeDetCyl
    = new G4Tubs("GeDetTop", 0.,geDet_radius,
		 0.5*geDet_hz,startAngle,spanningAngle);

  //massiv cylinder for the hole in the GeDet
  G4Tubs* GeDetHole
    = new G4Tubs("GeDetBottom",0 ,geDetHole_radius,
		 0.5*geDetHole_hz,startAngle,spanningAngle);

  //position of the hole relative to the massiv cylinder
  //small offsets to avoid visualisation error
  G4ThreeVector geDetHoleTrans = G4ThreeVector(0.0001 *mm, 0,
			     .5*(geDetHole_hz-geDet_hz)-0.001 *cm);

  //subtraction of hole from the active part and init. of the resulting solid
  G4SubtractionSolid* GeDet = new G4SubtractionSolid("GeDet", GeDetCyl,
			      GeDetHole, 0 , geDetHoleTrans);
  fLogicGeDet
    = new G4LogicalVolume(GeDet, fGeDetMaterial,"Chamber_LV",0,0,0);

  fLogicGeDet ->SetVisAttributes(geDetVisAtt);
 
  new G4PVPlacement(transformGeDet,  // rotation and position
                    fLogicGeDet,   
                    "GeDet",     
                    worldLV,       
                    false,          
                    3,              
		    fCheckOverlaps); 
  

//Deadlayers
  //rotation of deadlayer just like GeDet and position of deadlayer
  G4Transform3D transformDeadlayer = G4Transform3D(rotm,pos_deadlayer);

  //Front Deadlayer
  G4Tubs* deadlayerFront
    = new G4Tubs("DeadlayerF", 0.,deadlayer_radius, 0.5*deadlayer_hz,
		 startAngle,spanningAngle);

  //Side Deadlayer
  G4Tubs* deadlayerSide
    = new G4Tubs("DeadlayerS", geDet_radius,deadlayer_radius, 0.5*geDet_hz,
		 startAngle,spanningAngle);

  //position of the side relative to the front
  G4ThreeVector dlSideTrans = G4ThreeVector(0, 0,-0.5*geDet_hz - 0.1 *cm);

  //union of the side and the front of the deadlayer to one solid
  G4UnionSolid* deadlayer = new G4UnionSolid("Deadlayer", deadlayerFront,
					     deadlayerSide, 0, dlSideTrans);
  fLogicDeadlayer
    = new G4LogicalVolume(deadlayer, fGeDetMaterial,"Deadlayer_LV",0,0,0);
  fLogicDeadlayer->SetVisAttributes(deadlayerVisAtt);

  new G4PVPlacement(transformDeadlayer,       // rotation and position
                    fLogicDeadlayer,    
                    "Deadlayer",        
                    worldLV,         
                    false,           
                    0,               
                    fCheckOverlaps); 
  
   
    
//Aluminium cap
  //rotation of aluminium cap just like GeDet and position of AlCap
  G4Transform3D transformAlCap = G4Transform3D(rotm,pos_alCap);

  //Front Cap
  G4Tubs* AlCapFront
    = new G4Tubs("AlCapF", 0.,alCap_outerRadius, 0.5*alCapF_hz,
		 startAngle,spanningAngle);

  //Side Cap
  G4Tubs* AlCapSide
    = new G4Tubs("AlCapS", alCap_innerRadius,alCap_outerRadius, 
		 0.5*alCapS_hz,startAngle,spanningAngle);

  //position of the side relative to the front
  G4ThreeVector alSideTrans = G4ThreeVector(0, 0.0001 *cm,-0.5*geDet_hz );

  //union of the side and the front of the deadlayer to one solid
  G4UnionSolid* alCap = new G4UnionSolid("AlCap", AlCapFront,
					     AlCapSide, 0, alSideTrans);
  fLogicAlCap
    = new G4LogicalVolume(alCap, fTargetMaterial,"AlCap_LV",0,0,0);
  fLogicAlCap ->SetVisAttributes(alVisAtt);
  new G4PVPlacement(transformAlCap,       // rotation and position
                    fLogicAlCap,    
                    "AlCap",        
                    worldLV,         
                    false,           
                    0,               
                    fCheckOverlaps); 
  


  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.5*det_sizeXY;
  fStepLimit = new G4UserLimits(maxStep);
  fLogicChamber[0]->SetUserLimits(fStepLimit);
  fLogicChamber[1]->SetUserLimits(fStepLimit);
  fLogicChamber[2]->SetUserLimits(fStepLimit);
  fLogicGeDet->SetUserLimits(fStepLimit);
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
  // of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

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

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
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
