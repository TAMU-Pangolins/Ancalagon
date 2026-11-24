#include "EventAction.hh"

#include <cstdio>
#include <map>

#include "BgoHit.hh"
#include "DsssdHit.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "HitNtuple.hh"
#include "RunConfig.hh"

void EventAction::EndOfEventAction(const G4Event* event) {
  auto* analysisManager = G4AnalysisManager::Instance();
  const G4int eventID = event->GetEventID();
  const RunConfig::NtupleIds& ids = RunConfig::Ids();
  const bool print = RunConfig::EventPrint();

  // Per-crystal sums for the compact BgoEvent tree (edep>0 steps only).
  std::map<int, std::pair<double, double>> crystalSums;  // id -> (edep MeV, first t ns)
  int dsssdNHits = 0;
  double dsssdEdepMeV = 0.0;
  bool dsssdRecoil = false;

  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();

  const G4int bgoHcID = hce ? sdManager->GetCollectionID("BgoSD/BgoHitsCollection") : -1;
  if (bgoHcID >= 0) {
    auto* hits = static_cast<BgoHitsCollection*>(hce->GetHC(bgoHcID));
    if (hits) {
      double totalEdepMeV = 0.0;
      if (print) std::printf("HITS BGO   n=%zu crystals=[", hits->entries());
      for (std::size_t i = 0; i < hits->entries(); ++i) {
        const BgoHit* hit = (*hits)[i];
        totalEdepMeV += hit->edepMeV;
        if (print) std::printf("%s%d", i == 0 ? "" : ",", hit->crystalID);

        if (hit->edepMeV > 0.0) {
          auto [it, inserted] =
              crystalSums.try_emplace(hit->crystalID, std::make_pair(0.0, hit->timeNs));
          it->second.first += hit->edepMeV;
          if (hit->timeNs < it->second.second) it->second.second = hit->timeNs;
        }

        if (ids.bgoSteps >= 0) {
          using namespace HitNtuple;
          const int id = ids.bgoSteps;
          analysisManager->FillNtupleIColumn(id, kBgoEventID, eventID);
          analysisManager->FillNtupleIColumn(id, kBgoCrystalID, hit->crystalID);
          analysisManager->FillNtupleDColumn(id, kBgoEdepMeV, hit->edepMeV);
          analysisManager->FillNtupleDColumn(id, kBgoXCm, hit->position.x() / cm);
          analysisManager->FillNtupleDColumn(id, kBgoYCm, hit->position.y() / cm);
          analysisManager->FillNtupleDColumn(id, kBgoZCm, hit->position.z() / cm);
          analysisManager->FillNtupleDColumn(id, kBgoTimeNs, hit->timeNs);
          analysisManager->FillNtupleIColumn(id, kBgoTrackID, hit->trackID);
          analysisManager->FillNtupleSColumn(id, kBgoParticle, hit->particleName);
          analysisManager->AddNtupleRow(id);
        }
      }
      if (print) std::printf("] edep_total=%.6f MeV\n", totalEdepMeV);
    }
  }

  const G4int dsssdHcID = hce ? sdManager->GetCollectionID("DsssdSD/DsssdHitsCollection") : -1;
  if (dsssdHcID >= 0) {
    auto* hits = static_cast<DsssdHitsCollection*>(hce->GetHC(dsssdHcID));
    if (hits) {
      double totalEdepMeV = 0.0;
      for (std::size_t i = 0; i < hits->entries(); ++i) {
        const DsssdHit* hit = (*hits)[i];
        totalEdepMeV += hit->edepMeV;
        if (hit->trackID == 1) dsssdRecoil = true;  // the recoil is the first primary

        if (ids.dsssdSteps >= 0) {
          using namespace HitNtuple;
          const int id = ids.dsssdSteps;
          analysisManager->FillNtupleIColumn(id, kDsssdEventID, eventID);
          analysisManager->FillNtupleDColumn(id, kDsssdEdepMeV, hit->edepMeV);
          analysisManager->FillNtupleDColumn(id, kDsssdXCm, hit->position.x() / cm);
          analysisManager->FillNtupleDColumn(id, kDsssdYCm, hit->position.y() / cm);
          analysisManager->FillNtupleDColumn(id, kDsssdZCm, hit->position.z() / cm);
          analysisManager->FillNtupleDColumn(id, kDsssdTimeNs, hit->timeNs);
          analysisManager->FillNtupleIColumn(id, kDsssdTrackID, hit->trackID);
          analysisManager->FillNtupleSColumn(id, kDsssdParticle, hit->particleName);
          analysisManager->AddNtupleRow(id);
        }
      }
      dsssdNHits = static_cast<int>(hits->entries());
      dsssdEdepMeV = totalEdepMeV;
      if (print) std::printf("HITS DSSSD n=%zu edep_total=%.6f MeV\n", hits->entries(), totalEdepMeV);
    }
  }

  if (ids.bgoEvent >= 0) {
    using namespace HitNtuple;
    RunConfig::CompactBuffers& buf = RunConfig::Buffers();
    buf.crystalID.clear();
    buf.edepMeV.clear();
    buf.tFirstNs.clear();
    for (const auto& [crystal, sums] : crystalSums) {
      buf.crystalID.push_back(crystal);
      buf.edepMeV.push_back(static_cast<float>(sums.first));
      buf.tFirstNs.push_back(static_cast<float>(sums.second));
    }
    const int id = ids.bgoEvent;
    analysisManager->FillNtupleIColumn(id, kEvEventID, eventID);
    analysisManager->FillNtupleIColumn(id, kEvNCrystals, static_cast<int>(buf.crystalID.size()));
    analysisManager->FillNtupleIColumn(id, kEvDsssdNHits, dsssdNHits);
    analysisManager->FillNtupleFColumn(id, kEvDsssdEdep, static_cast<float>(dsssdEdepMeV));
    analysisManager->FillNtupleIColumn(id, kEvDsssdRecoil, dsssdRecoil ? 1 : 0);
    analysisManager->AddNtupleRow(id);
  }

  if (ids.libTruth >= 0) {
    using namespace HitNtuple;
    const RunConfig::Truth& t = RunConfig::CurrentTruth();
    analysisManager->FillNtupleFColumn(ids.libTruth, kLtCosThetaStar, t.gammaCosThetaStar.at(0));
    analysisManager->FillNtupleFColumn(ids.libTruth, kLtPhiStar, t.gammaPhiStar.at(0));
    analysisManager->AddNtupleRow(ids.libTruth);
  }

  if (ids.truth >= 0) {
    using namespace HitNtuple;
    const RunConfig::Truth& t = RunConfig::CurrentTruth();
    const int id = ids.truth;
    analysisManager->FillNtupleIColumn(id, kTrEventID, eventID);
    analysisManager->FillNtupleIColumn(id, kTrMode, t.mode);
    analysisManager->FillNtupleIColumn(id, kTrCellIndex, t.cellIndex);
    analysisManager->FillNtupleFColumn(id, kTrVx, static_cast<float>(t.vertexXCm));
    analysisManager->FillNtupleFColumn(id, kTrVy, static_cast<float>(t.vertexYCm));
    analysisManager->FillNtupleFColumn(id, kTrVz, static_cast<float>(t.vertexZCm));
    analysisManager->FillNtupleDColumn(id, kTrBeamKe, t.beamKeMeV);
    analysisManager->FillNtupleDColumn(id, kTrEcm, t.ecmMeV);
    analysisManager->FillNtupleDColumn(id, kTrBeta, t.emitterBeta);
    analysisManager->FillNtupleDColumn(id, kTrRecoilPx, t.recoilPMeV.x());
    analysisManager->FillNtupleDColumn(id, kTrRecoilPy, t.recoilPMeV.y());
    analysisManager->FillNtupleDColumn(id, kTrRecoilPz, t.recoilPMeV.z());
    analysisManager->FillNtupleDColumn(id, kTrRecoilE, t.recoilEMeV);
    const auto* endVol = static_cast<const G4VPhysicalVolume*>(t.recoilEndVolume);
    analysisManager->FillNtupleSColumn(id, kTrRecoilEndVolume,
                                       endVol ? endVol->GetName() : G4String("OutOfWorld"));
    analysisManager->FillNtupleFColumn(id, kTrRecoilEndZ, static_cast<float>(t.recoilEndZCm));
    analysisManager->AddNtupleRow(id);
  }
}
