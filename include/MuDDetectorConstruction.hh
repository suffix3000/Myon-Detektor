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
// $Id: B2aDetectorConstruction.hh 73722 2013-09-09 10:23:05Z gcosmo $
//
/// \file B2aDetectorConstruction.hh
/// \brief Definition of the B2aDetectorConstruction class

#ifndef B2aDetectorConstruction_h
#define B2aDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class B2aDetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class B2aDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B2aDetectorConstruction();
    virtual ~B2aDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Set methods
    void SetTargetMaterial (G4String );
    void SetChamberMaterial(G4String );
    void SetMaxStep (G4double );
    void SetCheckOverlaps(G4bool );

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    // data members
    G4int fNbOfSiPM;				//number of SiPM
    //poiters to logical volumes 
    G4LogicalVolume*   fLogicTarget;     	//target
    G4LogicalVolume*   fLogicStopper;   	//stopper
    G4LogicalVolume*   fLogicGeDetActive;       //active region of GeDet
    G4LogicalVolume*   fLogicGeDetCrystal;	//germanium crystal = outer deadlayer
    G4LogicalVolume*   fLogicGeDetID;  		//inner deadlayer
    G4LogicalVolume*   fLogicGeDetHole; 	//hole for cooling finger in crystal
    G4LogicalVolume*   fLogicGeDetOS; 		//outer shielding
    G4LogicalVolume*   fLogicGeDetIS;   	//inner shielding
    G4LogicalVolume*   fLogicGeDetOuterVacuum;  //outer vacuum between OS and IS
    G4LogicalVolume*   fLogicGeDetInnerVacuum;  //inner vacuum between IS and crystal
    G4LogicalVolume*   fLogicGeDetISHole;	//hole for cooling finger in IS	
    G4LogicalVolume**  fLogicSiPM;   		//SiPM

    // pointer to materials
    G4Material*        fTargetMaterial; 	//target
    G4Material*        fSiPMMaterial;   	//SiPM
    G4Material*        fStopperMaterial; 	//stopper
    G4Material*        fGeDetMaterial;   	//germanium detector
    G4Material*        fVacuum;			//vaccum	
    G4Material*	       fAir;			//air
    G4UserLimits* fStepLimit;         		// pointer to user step limits

    B2aDetectorMessenger*  fMessenger;   	// messenger

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                         // magnetic field messenger
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
