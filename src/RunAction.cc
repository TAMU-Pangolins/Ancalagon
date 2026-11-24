#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "HitNtuple.hh"
#include "RunConfig.hh"

RunAction::RunAction() {
  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(0);
  RunConfig::NtupleIds& ids = RunConfig::Ids();

  if (RunConfig::WriteStepTrees()) {
    // "Bgo": one row per BGO hit (BgoHit, see BgoHit.hh) -- column order
    // here must match HitNtuple.hh's kBgo* indices exactly.
    ids.bgoSteps = analysisManager->CreateNtuple("Bgo", "BGO crystal (SCNT) hits");
    analysisManager->CreateNtupleIColumn("eventID");   // kBgoEventID
    analysisManager->CreateNtupleIColumn("crystalID");  // kBgoCrystalID, 1-30
    analysisManager->CreateNtupleDColumn("edepMeV");    // kBgoEdepMeV
    analysisManager->CreateNtupleDColumn("x_cm");        // kBgoXCm
    analysisManager->CreateNtupleDColumn("y_cm");        // kBgoYCm
    analysisManager->CreateNtupleDColumn("z_cm");        // kBgoZCm
    analysisManager->CreateNtupleDColumn("t_ns");        // kBgoTimeNs
    analysisManager->CreateNtupleIColumn("trackID");    // kBgoTrackID
    analysisManager->CreateNtupleSColumn("particle");   // kBgoParticle
    analysisManager->FinishNtuple();

    // "Dsssd": one row per DSSSD hit (DsssdHit, see DsssdHit.hh) -- no
    // per-strip channel yet (see README), just the whole slab's own hits.
    ids.dsssdSteps = analysisManager->CreateNtuple("Dsssd", "DSSSD hits");
    analysisManager->CreateNtupleIColumn("eventID");  // kDsssdEventID
    analysisManager->CreateNtupleDColumn("edepMeV");   // kDsssdEdepMeV
    analysisManager->CreateNtupleDColumn("x_cm");       // kDsssdXCm
    analysisManager->CreateNtupleDColumn("y_cm");       // kDsssdYCm
    analysisManager->CreateNtupleDColumn("z_cm");       // kDsssdZCm
    analysisManager->CreateNtupleDColumn("t_ns");       // kDsssdTimeNs
    analysisManager->CreateNtupleIColumn("trackID");   // kDsssdTrackID
    analysisManager->CreateNtupleSColumn("particle");  // kDsssdParticle
    analysisManager->FinishNtuple();
  }

  if (RunConfig::WriteCompactTrees()) {
    // "BgoEvent": one row per simulated event (including events with no
    // BGO energy), energy summed per crystal over edep>0 steps. Column
    // order must match HitNtuple.hh's kEv* indices (vector columns take
    // an index too).
    RunConfig::CompactBuffers& buf = RunConfig::Buffers();
    ids.bgoEvent = analysisManager->CreateNtuple("BgoEvent", "BGO energy per crystal, per event");
    analysisManager->CreateNtupleIColumn("eventID");                   // 0 kEvEventID
    analysisManager->CreateNtupleIColumn("nCrystals");                 // 1 kEvNCrystals
    analysisManager->CreateNtupleIColumn("crystalID", buf.crystalID);  // 2 vector, 1-30
    analysisManager->CreateNtupleFColumn("edepMeV", buf.edepMeV);      // 3 vector, MeV
    analysisManager->CreateNtupleFColumn("tFirstNs", buf.tFirstNs);    // 4 vector, ns
    analysisManager->CreateNtupleIColumn("dsssdNHits");                // 5 kEvDsssdNHits
    analysisManager->CreateNtupleFColumn("dsssdEdepMeV");              // 6 kEvDsssdEdep
    analysisManager->CreateNtupleIColumn("dsssdRecoil");               // 7 kEvDsssdRecoil
    analysisManager->FinishNtuple();

    if (RunConfig::Library().nCells > 0) {
      // Gamma-library mode: lean truth. Everything else is fixed by the
      // cell (cellIndex = eventID / eventsPerCell; E0, z, beta from the
      // cells file, which the packer reads alongside this file); only the
      // emission direction in the emitter frame is random. Lab direction
      // and energy follow from these and beta.
      ids.libTruth = analysisManager->CreateNtuple("LibTruth", "Gamma-library emission truth");
      analysisManager->CreateNtupleFColumn("cosThetaStar");  // 0 kLtCosThetaStar
      analysisManager->CreateNtupleFColumn("phiStar");       // 1 kLtPhiStar, [-pi, pi)
      analysisManager->FinishNtuple();
    } else {
    // "Truth": one row per simulated event (RunConfig::Truth).
    RunConfig::Truth& t = RunConfig::CurrentTruth();
    ids.truth = analysisManager->CreateNtuple("Truth", "Generated event truth");
    analysisManager->CreateNtupleIColumn("eventID");      // 0 kTrEventID
    analysisManager->CreateNtupleIColumn("mode");         // 1 kTrMode
    analysisManager->CreateNtupleIColumn("cellIndex");    // 2 kTrCellIndex
    analysisManager->CreateNtupleFColumn("vertexX_cm");   // 3 kTrVx
    analysisManager->CreateNtupleFColumn("vertexY_cm");   // 4 kTrVy
    analysisManager->CreateNtupleFColumn("vertexZ_cm");   // 5 kTrVz
    analysisManager->CreateNtupleDColumn("beamKeMeV");    // 6 kTrBeamKe
    analysisManager->CreateNtupleDColumn("ecmMeV");       // 7 kTrEcm
    analysisManager->CreateNtupleDColumn("emitterBeta");  // 8 kTrBeta
    analysisManager->CreateNtupleDColumn("recoilPx");     // 9 kTrRecoilPx, MeV/c
    analysisManager->CreateNtupleDColumn("recoilPy");     // 10 kTrRecoilPy
    analysisManager->CreateNtupleDColumn("recoilPz");     // 11 kTrRecoilPz
    analysisManager->CreateNtupleDColumn("recoilE");      // 12 kTrRecoilE, total, MeV
    analysisManager->CreateNtupleFColumn("gammaELab", t.gammaELabMeV);
    analysisManager->CreateNtupleFColumn("gammaEStar", t.gammaEStarMeV);
    analysisManager->CreateNtupleFColumn("gammaCosThetaLab", t.gammaCosThetaLab);
    analysisManager->CreateNtupleFColumn("gammaPhiLab", t.gammaPhiLab);
    analysisManager->CreateNtupleFColumn("gammaCosThetaStar", t.gammaCosThetaStar);
    analysisManager->CreateNtupleFColumn("gammaPhiStar", t.gammaPhiStar);
    analysisManager->CreateNtupleIColumn("levelFrom", t.levelFrom);
    analysisManager->CreateNtupleIColumn("levelTo", t.levelTo);
    analysisManager->CreateNtupleSColumn("recoilEndVolume");  // 21 kTrRecoilEndVolume
    analysisManager->CreateNtupleFColumn("recoilEndZ_cm");    // 22 kTrRecoilEndZ
    analysisManager->FinishNtuple();
    }

    // "RunInfo": one row per run.
    ids.runInfo = analysisManager->CreateNtuple("RunInfo", "Run configuration");
    analysisManager->CreateNtupleIColumn("nEvents");         // 0 kRiNEvents
    analysisManager->CreateNtupleDColumn("randomSeed");      // 1 kRiSeed (-1 = unset)
    analysisManager->CreateNtupleDColumn("pressureTorr");    // 2 kRiPressure
    analysisManager->CreateNtupleDColumn("temperatureK");    // 3 kRiTemperature (0 = unset)
    analysisManager->CreateNtupleDColumn("densityFracStp");  // 4 kRiDensityFrac
    analysisManager->CreateNtupleIColumn("nLibraryCells");   // 5 kRiNCells
    analysisManager->CreateNtupleIColumn("eventsPerCell");   // 6 kRiEventsPerCell
    analysisManager->FinishNtuple();
  }
}

void RunAction::BeginOfRunAction(const G4Run*) {
  G4AnalysisManager::Instance()->OpenFile(RunConfig::OutputFile());
}

void RunAction::EndOfRunAction(const G4Run* run) {
  auto* analysisManager = G4AnalysisManager::Instance();
  const int id = RunConfig::Ids().runInfo;
  if (id >= 0) {
    using namespace HitNtuple;
    analysisManager->FillNtupleIColumn(id, kRiNEvents, run->GetNumberOfEvent());
    analysisManager->FillNtupleDColumn(id, kRiSeed, static_cast<double>(RunConfig::RandomSeed()));
    analysisManager->FillNtupleDColumn(id, kRiPressure, RunConfig::GasPressureTorr());
    analysisManager->FillNtupleDColumn(id, kRiTemperature, RunConfig::GasTemperatureK());
    analysisManager->FillNtupleDColumn(id, kRiDensityFrac,
                                       RunConfig::TargetGasDensityFractionOfStp());
    analysisManager->FillNtupleIColumn(id, kRiNCells, RunConfig::Library().nCells);
    analysisManager->FillNtupleIColumn(id, kRiEventsPerCell, RunConfig::Library().eventsPerCell);
    analysisManager->AddNtupleRow(id);
  }
  analysisManager->Write();
  analysisManager->CloseFile();
}
