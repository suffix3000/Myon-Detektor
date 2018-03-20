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
// $Id: B2PrimaryGeneratorAction.cc 97979 2016-06-30 09:36:20Z gcosmo $
//
/// \file B2PrimaryGeneratorAction.cc
/// \brief Implementation of the B2PrimaryGeneratorAction class

#include "MuDPrimaryGeneratorAction.hh"
#include "MuDAnalysis.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <fstream>
using namespace std;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::B2PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{ 
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::~B2PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double worldHalfZ = 0;
  G4double worldHalfXY = 0;
  G4LogicalVolume* worldLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = NULL;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox )
  {
      worldHalfZ = worldBox->GetZHalfLength();
      worldHalfXY = worldBox->GetXHalfLength();
  }

  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }
  
  //ParticleGun
  //Charge Ratio/primary particle definition
  G4float chargeratio = 1.27;
  G4String name;
  if (G4UniformRand() > 1/(1 + chargeratio)) name="mu+";
  else name = "mu-";
  G4ParticleDefinition* particleDefinition 
	 = G4ParticleTable::GetParticleTable()->FindParticle(name);
  fParticleGun->SetParticleDefinition(particleDefinition);

  //Position
  G4double size = 1.; 
  G4double x0 = 2 * size * worldHalfXY * (G4UniformRand() - 0.5);
  G4double y0 = 2 * size * worldHalfXY * (G4UniformRand() - 0.5);
  G4double z0 = -worldHalfZ;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  //Energy
  //Konstant integriertes Spektrum nach ExPacs
  G4int ene_bins = 79;
  G4double ene_dis[] = {4.307e-12,6.771e-12,1.037e-11,1.640e-11,2.504e-11,
			3.934e-11,6.049e-11,9.452e-11,1.459e-10,2.279e-10,
			3.495e-10,5.495e-10,8.418e-10,1.328e-09,2.022e-09,
			3.171e-09,4.860e-09,7.563e-09,1.160e-08,1.799e-08,
			2.741e-08,4.255e-08,6.440e-08,1.003e-07,1.497e-07,
			2.292e-07,3.416e-07,5.138e-07,7.590e-07,1.127e-06,
			1.630e-06,2.399e-06,3.422e-06,5.015e-06,7.066e-06,
			1.028e-05,1.462e-05,2.102e-05,2.972e-05,4.194e-05,
			5.699e-05,7.770e-05,1.012e-04,1.330e-04,1.657e-04,
			2.088e-04,2.523e-04,3.043e-04,3.560e-04,4.136e-04,
			4.640e-04,5.217e-04,5.566e-04,5.986e-04,6.030e-04,
			6.036e-04,5.707e-04,5.277e-04,4.637e-04,3.961e-04,
			3.191e-04,2.531e-04,1.888e-04,1.407e-04,9.800e-05,
			6.836e-05,4.565e-05,3.036e-05,1.957e-05,1.256e-05,
			7.845e-06,4.969e-06,3.036e-06,1.903e-06,1.115e-06,
			6.343e-07,3.492e-07,1.944e-07,1.064e-07};
  CLHEP::RandGeneral randG_ene(ene_dis, ene_bins,0);
  G4double randnumb =  randG_ene.fire();
  //From flat random to exp. random with two parameters
  G4double a = -4.48295;
  G4double b = 18.1893;
  G4double energy = exp(a + b*randnumb);
  fParticleGun->SetParticleEnergy(energy);

  //Direction
  G4int ang_bins = 100;
  G4double ang_dis[] = {1.57E-2,1.57E-2,1.57E-2,1.57E-2,1.56E-2,1.56E-2,1.55E-2,
			1.55E-2,1.54E-2,1.54E-2,1.53E-2,1.52E-2,1.51E-2,1.50E-2,
			1.49E-2,1.48E-2,1.47E-2,1.46E-2,1.44E-2,1.43E-2,1.41E-2,
			1.40E-2,1.38E-2,1.37E-2,1.35E-2,1.33E-2,1.31E-2,1.30E-2,
			1.28E-2,1.26E-2,1.24E-2,1.22E-2,1.20E-2,1.17E-2,1.15E-2,
			1.13E-2,1.11E-2,1.09E-2,1.06E-2,1.04E-2,1.02E-2,9.93E-3,
			9.69E-3,9.45E-3,9.20E-3,8.96E-3,8.72E-3,8.47E-3,8.22E-3,
			7.98E-3,7.73E-3,7.48E-3,7.24E-3,6.99E-3,6.75E-3,6.50E-3,
			6.26E-3,6.02E-3,5.78E-3,5.54E-3,5.31E-3,5.08E-3,4.85E-3,
			4.62E-3,4.40E-3,4.18E-3,3.96E-3,3.75E-3,3.54E-3,3.34E-3,
			3.14E-3,2.94E-3,2.75E-3,2.57E-3,2.39E-3,2.21E-3,2.05E-3,
			1.88E-3,1.72E-3,1.57E-3,1.43E-3,1.29E-3,1.16E-3,1.03E-3,
			9.13E-4,8.01E-4,6.96E-4,5.98E-4,5.07E-4,4.24E-4,3.48E-4,
			2.79E-4,2.17E-4,1.63E-4,1.17E-4,7.87E-5,4.78E-5,2.45E-5,
			9.04E-6,1.29E-6};
  CLHEP::RandGeneral randG_ang(ang_dis, ang_bins);
  G4double theta = randG_ang.fire()*CLHEP::pi/2;
  G4double phi =  G4UniformRand()*2*CLHEP::pi;
  G4double px = cos(phi)*sin(theta);
  G4double py = sin(phi)*sin(theta);
  G4double pz = cos(theta);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));

  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
