// Run-time configuration read from environment variables, plus the
// per-event truth record and output buffers shared between the primary
// generator (which knows the truth), the sensitive detectors, EventAction
// (which fills the ntuples) and RunAction (which books them).
//
// Every setting here defaults to the pre-existing behaviour, so a run with
// none of these variables set produces the same output as before:
//
//   OUTPUT_FILE        ROOT output file name (default "dragon_hits.root").
//   OUTPUT_MODE        "steps"   (default): legacy per-step Bgo/Dsssd trees.
//                      "compact": one row per event in BgoEvent (summed per
//                                 crystal, zero-edep steps dropped) plus
//                                 Truth and RunInfo trees; no step trees.
//                      "both":    everything.
//   EVENT_PRINT        1/0: per-event "HITS BGO ..." stdout line. Default 1
//                      in "steps"/"both" mode, 0 in "compact" mode.
//   GAS_PRESSURE_TORR  target gas pressure (default 5.0 Torr). Every gas
//                      volume (cell, collars, pumping chain) scales with it.
//   GAS_TEMPERATURE_K  target gas temperature. Unset (default): the NIST
//                      reference density is scaled by pressure only (the
//                      legacy behaviour). Set: density additionally scales
//                      by kNistReferenceTemperatureK / T (ideal gas).
//
// Gamma-library mode (--gamma-library), see PrimaryGeneratorAction:
//   LIBRARY_CELLS            text file, one cell per line: "E0_MeV z_cm beta"
//                            ('#' comments allowed).
//   LIBRARY_EVENTS_PER_CELL  events simulated per cell (default 1000).
//
// Recoil acceptance scan (--track-reaction with RECOIL_SCAN=1), see
// PrimaryGeneratorAction: the recoil alone (no gammas), with its momentum
// replaced by p_compound*(1+delta) along direction (theta_x, theta_y),
// sampled uniformly in |theta_x|,|theta_y| <= RECOIL_SCAN_THETA_MAX (rad,
// default 0.035) and |delta| <= RECOIL_SCAN_DELTA_MAX (default 0.05). The
// separator retune (if any) still uses unperturbed reaction events.
#pragma once

#include <string>
#include <vector>

#include "G4ThreeVector.hh"

namespace RunConfig {

enum class OutputMode { kSteps, kCompact, kBoth };

// G4_H / G4_He NIST densities are tabulated at 20 C, 1 atm.
constexpr double kNistReferenceTemperatureK = 293.15;

std::string OutputFile();
OutputMode Mode();
bool WriteStepTrees();    // kSteps or kBoth
bool WriteCompactTrees(); // kCompact or kBoth
bool EventPrint();
double GasPressureTorr();
double GasTemperatureK();  // 0 if unset (legacy: no temperature scaling)
// Multiplier applied to the NIST STP density for the full-pressure target
// gas: (P/760) * (T_ref/T, or 1 if GAS_TEMPERATURE_K is unset).
double TargetGasDensityFractionOfStp();
bool RecoilScan();
double RecoilScanThetaMax();
double RecoilScanDeltaMax();
long RandomSeed();  // RANDOM_SEED, or -1 if unset

// Filled by PrimaryGeneratorAction for every event, read by EventAction.
// Gamma quantities are ordered by emission in the cascade.
struct Truth {
  int mode = 0;        // 0 = reaction, 1 = gamma library
  int cellIndex = -1;  // gamma-library cell, -1 otherwise
  double vertexXCm = 0.0, vertexYCm = 0.0, vertexZCm = 0.0;
  double beamKeMeV = 0.0;    // beam kinetic energy at the reaction
  double ecmMeV = 0.0;       // centre-of-mass energy above threshold
  double emitterBeta = 0.0;  // speed (v/c) of the compound system / emitter
  G4ThreeVector recoilPMeV;  // final ground-state recoil momentum, lab
  double recoilEMeV = 0.0;   // final recoil total energy, lab
  std::vector<float> gammaELabMeV, gammaEStarMeV;
  std::vector<float> gammaCosThetaLab, gammaPhiLab;    // lab direction
  std::vector<float> gammaCosThetaStar, gammaPhiStar;  // in emitting parent's rest frame
  std::vector<int> levelFrom, levelTo;                 // cascade step, -1 = resonance
  // Where the recoil (track 1) was last seen: its last step's post-step
  // volume (nullptr = left the world) and z. Filled by SteppingAction.
  const void* recoilEndVolume = nullptr;  // G4VPhysicalVolume*
  double recoilEndZCm = 0.0;

  void Clear();
};
Truth& CurrentTruth();

// Bound by reference to vector columns at ntuple creation (RunAction) and
// refilled per event (EventAction); must outlive the run, hence globals.
struct CompactBuffers {
  std::vector<int> crystalID;
  std::vector<float> edepMeV;
  std::vector<float> tFirstNs;
};
CompactBuffers& Buffers();

// Ntuple IDs as returned by G4AnalysisManager::CreateNtuple (they depend on
// which trees the output mode books). -1 = not booked.
struct NtupleIds {
  int bgoSteps = -1, dsssdSteps = -1, bgoEvent = -1, truth = -1, libTruth = -1, runInfo = -1;
};
NtupleIds& Ids();

// Gamma-library bookkeeping, set by PrimaryGeneratorAction, recorded in RunInfo.
struct LibraryInfo {
  int nCells = 0;
  int eventsPerCell = 0;
};
LibraryInfo& Library();

}  // namespace RunConfig
