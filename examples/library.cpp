// Synthetic seismograms from C++: the example in README.md, which the test
// suite builds and runs.

#include <YSpec/Synthetics.hpp>
#include <cmath>
#include <iostream>

int main() {
  auto model = YSpec::Model::ReadDeck(YSPEC_EXAMPLES_DIR "/prem.200");
  auto settings = YSpec::Settings{};
  settings.lMax = 200;
  settings.fMin = 0.2;  // mHz
  settings.fMax = 15;
  settings.length = 60;            // minutes
  settings.timeStep = 20;          // seconds
  settings.sourceDepth = 29.5;     // km
  settings.receiverDepth = 3;      // km: the sea floor of this model
  settings.moment = {.rr = 1e18};  // N m
  settings.receivers = {{0, 30}, {0, 60}};
  auto seismograms = YSpec::ComputeSeismograms(model, settings);
  double z = seismograms(0, 0, 100);  // receiver 0, vertical, sample 100

  std::cout << "vertical displacement at 30 degrees, 2000 s: " << z << " m\n";
  return std::isfinite(z) && z != 0 ? 0 : 1;
}
