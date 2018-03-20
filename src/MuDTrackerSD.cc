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
// $Id: B2TrackerSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "MuDAnalysis.hh"
#include "MuDRunAction.hh"
#include "MuDTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::B2TrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::~B2TrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new B2TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*) 
{

  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
 
  if (edep==0.) return false;

  B2TrackerHit* newHit = new B2TrackerHit();

 
 
  newHit->SetChamberNb (aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());

  newHit->SetEdep(edep);
 
  fHitsCollection->insert(newHit);
 
  newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent*)
{ 
 G4int nofHits = fHitsCollection->entries();
 if ( verboseLevel>0 ) { 
    
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
     }
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  G4int HitsD0 = 0;
  G4int HitsD1 = 0;
  G4int HitsD2 = 0;
  G4int HitsD3 = 0;
  G4int DetNb = 0;
  G4double eMess = 0;
  G4double eDep = 0;
  for ( G4int i=0; i<nofHits; i++ ){
    eDep = (*fHitsCollection)[i]->GetEdep();  
    DetNb = (*fHitsCollection)[i]->GetChamberNb();
    if(DetNb == 0){HitsD0++;}
    if(DetNb == 1){HitsD1++;}
    if(DetNb == 2){HitsD2++;}
    if(DetNb == 3){
      HitsD3++;
      eMess += eDep;
    }
    
  }
  //Data output in root file
  if(HitsD0 > 0 and HitsD1 > 0 and HitsD2 == 0){
    analysisManager->FillH1(0, eMess);
    analysisManager->FillNtupleDColumn(0, eMess);
    analysisManager->AddNtupleRow();  
  }
  else{eMess = 0;}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
