# Myon-Detektor
Geant4-Simulation eines Detektoraufbaus zur Vermessung von Einfangreaktionen atmosphärischer Myonen in unterschiedlichen Materialien.

Ziele der Simulation:
1. Untersuchung des Einflusses eines Bleiblocks (Stopper) auf die Anzahl gestoppter Myonen im Target
2. Untersuchung des durch Myon-Einfangreaktionen erzeugten Spektrums
3. Untersuchung der Kerzerfälle, welche auf einen Myon-Einfang folgen.

Die Geometrie der Simulation besteht aus einem Target aus Aluminium (oder anderen Materialien), drei Silizium-Photomultiplier (SiPM), einem Germanium-Detektor (GeDet) und einem Stopper aus Blei.
In dem Target sollen die Myonen gestoppt und eingefangen werden. Um zu überprüfen ob dies geschehen ist, befinden sich über dem Target zwei SiPMs und darunter ein weiterer SiPM. Geben die beiden oberen SiPM ein Signal und der untere nicht (Antikoinzidenz), ist das Myon mutmaßlich im Target gestoppt wurden.
Auf die Antikoinzidenz wird der GeDet aktiviert und vermisst die aus dem Einfangprozess resultierende Strahlung.
Der Stopper aus Blei soll die Myonen "abbremsen" und so die Anzahl der niederenergetischen Myonen erhöhen, da diese eine höhere Wahrscheinlichkeit haben im Target gestoppt zu werden. 
Die Energiesverteilung für die atmosphärischen Myonen wurde von ExPacs übernommen. Die Winkelverteilung ist die bekannte Cosinus-Quadrat-Verteilung.

Zur Erstellung der Simulation wird geant4/10.03.p01-gcc710_withQT verwendet. 

Wichtige Dateien:

.cc:
Muon-Detector.cc
MuDDetectorConstruction.cc     MuDRunAction.cc
MuDDetectorMessenger.cc        MuDTrackerHit.cc
MuDEventAction.cc              MuDTrackerSD.cc
MuDActionInitialization.cc     MuDPrimaryGeneratorAction.cc
MuDDetectorConstruction.cc     MuDRunAction.cc

.hh:
MuDDetectorConstruction.hh    MuDTrackerHit.hh
MuDDetectorMessenger.hh       MuDTrackerSD.hh
MuDEventAction.hh             MuDActionInitialization.hh 
MuDPrimaryGeneratorAction.hh  MuDAnalysis.hh         
MuDRunAction.hh
