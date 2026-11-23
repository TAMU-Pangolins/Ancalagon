// Shared ntuple/column layout for the "Bgo" and "Dsssd" ROOT TTrees --
// RunAction.cc creates columns in exactly this order (G4AnalysisManager
// assigns ntuple/column IDs by creation order, not by name), EventAction.cc
// fills them; kept in one place so the two can't drift out of sync.
#pragma once

namespace HitNtuple {

// Ntuple IDs themselves depend on OUTPUT_MODE -- see RunConfig::Ids().

// Bgo columns.
constexpr int kBgoEventID = 0;
constexpr int kBgoCrystalID = 1;
constexpr int kBgoEdepMeV = 2;
constexpr int kBgoXCm = 3;
constexpr int kBgoYCm = 4;
constexpr int kBgoZCm = 5;
constexpr int kBgoTimeNs = 6;
constexpr int kBgoTrackID = 7;
constexpr int kBgoParticle = 8;

// Dsssd columns.
constexpr int kDsssdEventID = 0;
constexpr int kDsssdEdepMeV = 1;
constexpr int kDsssdXCm = 2;
constexpr int kDsssdYCm = 3;
constexpr int kDsssdZCm = 4;
constexpr int kDsssdTimeNs = 5;
constexpr int kDsssdTrackID = 6;
constexpr int kDsssdParticle = 7;

// BgoEvent scalar columns (vector columns sit between them but are
// filled through RunConfig::Buffers(), not by index).
constexpr int kEvEventID = 0;
constexpr int kEvNCrystals = 1;
constexpr int kEvDsssdNHits = 5;
constexpr int kEvDsssdEdep = 6;
constexpr int kEvDsssdRecoil = 7;

// Truth scalar columns.
constexpr int kTrEventID = 0;
constexpr int kTrMode = 1;
constexpr int kTrCellIndex = 2;
constexpr int kTrVx = 3;
constexpr int kTrVy = 4;
constexpr int kTrVz = 5;
constexpr int kTrBeamKe = 6;
constexpr int kTrEcm = 7;
constexpr int kTrBeta = 8;
constexpr int kTrRecoilPx = 9;
constexpr int kTrRecoilPy = 10;
constexpr int kTrRecoilPz = 11;
constexpr int kTrRecoilE = 12;

// LibTruth columns (gamma-library mode).
constexpr int kLtCosThetaStar = 0;
constexpr int kLtPhiStar = 1;

// RunInfo columns.
constexpr int kRiNEvents = 0;
constexpr int kRiSeed = 1;
constexpr int kRiPressure = 2;
constexpr int kRiTemperature = 3;
constexpr int kRiDensityFrac = 4;
constexpr int kRiNCells = 5;
constexpr int kRiEventsPerCell = 6;

}  // namespace HitNtuple
