

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"


using namespace Pythia8;
using namespace HepMC3;

//==========================================================================

int main() {

  // Interface for conversion from Pythia8::Event to HepMC event.
  // Specify file where HepMC events will be stored.

  
  Pythia8ToHepMC toHepMC("PYTHIA8_eic_275_10_50000Events_100GeV_decays_off.hepmc");
  //toHepMC.setNewFile("test.hepmc");

  //100GeV cross sec - 2.15541e-04, total 1.30538e+05
  //figure out cross sections!

  // Generator. Process selection. RHIC initialization. Histogram.
  Pythia pythia;

  
  pythia.readString("Beams:idA = -11");
pythia.readString("Beams:idB = 2212");
pythia.readString("Beams:frameType = 2");
pythia.readString("Beams:eA = 10");
pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
pythia.readString("SpaceShower:dipoleRecoil = on");
pythia.readString("Beams:eB = 275.");
pythia.readString("PhaseSpace:pTHatMax = 1");
  
  //Decays Off
  
    pythia.readString(("111:mayDecay = false")); //pi0
    pythia.readString(("211:mayDecay = false")); //pi+
    pythia.readString(("221:mayDecay = false")); //eta
    pythia.readString(("321:mayDecay = false")); //K+
    pythia.readString(("310:mayDecay = false")); //Kshort
    pythia.readString(("130:mayDecay = false")); //Klong
    pythia.readString(("3122:mayDecay = false"));//Lambda0
    pythia.readString(("3212:mayDecay = false"));//Sigma0
    pythia.readString(("3112:mayDecay = false"));//Sigma-
    pythia.readString(("3222:mayDecay = false"));//Sigma+
    pythia.readString(("3312:mayDecay = false"));//Xi-
    pythia.readString(("3322:mayDecay = false"));//Xi0
    pythia.readString(("3334:mayDecay = false"));//Omega-

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 50000; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

    // Construct new empty HepMC event, fill it and write it out.
    toHepMC.writeNextEvent( pythia );

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  
  double crossSection = pythia.info.sigmaGen(0);
  double weightsumm = pythia.info.weightSum();

  std::cout << std::setprecision(5)<< "Estimated cross section: " << (crossSection) << " and total: " << weightsumm << std::endl;

  // Done.
  return 0;
}
