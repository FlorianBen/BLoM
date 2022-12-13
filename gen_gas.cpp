#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewSignal.hh"

int main(int argc, char *argv[]) {
  using namespace Garfield;
  MediumMagboltz gas;
  const size_t nE = 100;
  const double emin = 100.;
  const double emax = 100000.;
  // Flag to request logarithmic spacing.
  constexpr bool useLog = true;
  gas.SetFieldGrid(emin, emax, nE, useLog);

  // Set composition
  gas.SetComposition("N2", 100);
  // Set the temperature [K].
  gas.SetTemperature(293.15);
  // Set the pressure [Torr].
  gas.SetPressure(1.1 * 760.);

  // Number of collision
  const int ncoll = 10;

  // Run Magboltz to generate the gas table.
  gas.GenerateGasTable(ncoll);

  // Write to file
  gas.WriteGasFile("data/n2.gas");

  ViewMedium view;
  view.SetMedium(&gas);
  view.PlotElectronVelocity('e');
}
