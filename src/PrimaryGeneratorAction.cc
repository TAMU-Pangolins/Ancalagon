#include "PrimaryGeneratorAction.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "GasStoppingPower.hh"
#include "Randomize.hh"
#include "ReactionConfig.hh"
#include "RunConfig.hh"
#include "TargetChamber.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(double x0Cm, double y0Cm, double momentumMeV,
                                                double z0Cm) {
  fGun = new G4ParticleGun(1);
  G4ParticleDefinition* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGun->SetParticleDefinition(proton);
  fGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  fGun->SetParticleMomentum(momentumMeV * MeV);
  fGun->SetParticlePosition(G4ThreeVector(x0Cm * cm, y0Cm * cm, z0Cm * cm));
}

PrimaryGeneratorAction::PrimaryGeneratorAction(double x0Cm, double y0Cm, double momentumMeV,
                                                double z0Cm, int ionZ, int ionA,
                                                int ionChargeState, double angleXDeg)
    : fNeedsIonLookup(true), fIonZ(ionZ), fIonA(ionA), fIonChargeState(ionChargeState) {
  fGun = new G4ParticleGun(1);
  const double angleXRad = angleXDeg * CLHEP::pi / 180.0;
  fGun->SetParticleMomentumDirection(
      G4ThreeVector(std::sin(angleXRad), 0.0, std::cos(angleXRad)));
  fGun->SetParticleMomentum(momentumMeV * MeV);
  fGun->SetParticlePosition(G4ThreeVector(x0Cm * cm, y0Cm * cm, z0Cm * cm));
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const std::string& reactionFilePath, double x0Cm,
                                               double y0Cm, double z0Cm)
    : fIsReaction(true),
      fConfig(ReactionConfig::Load(reactionFilePath)),
      fReactionX0Cm(x0Cm),
      fReactionY0Cm(y0Cm),
      fReactionZ0Cm(z0Cm) {
  fReaction = new ReactionKinematics(fConfig);
  fRecoilZ = fConfig.recoil.Z;
  fRecoilA = fConfig.recoil.A;
  fRecoilChargeState = fConfig.recoilChargeState;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const std::vector<LibraryCell>& cells,
                                               int eventsPerCell, double x0Cm, double y0Cm)
    : fIsLibrary(true),
      fLibraryCells(cells),
      fEventsPerCell(eventsPerCell),
      fLibraryX0Cm(x0Cm),
      fLibraryY0Cm(y0Cm) {
  if (fLibraryCells.empty()) throw std::runtime_error("gamma library: no cells");
  if (fEventsPerCell <= 0) throw std::runtime_error("gamma library: eventsPerCell must be > 0");
  for (const auto& c : fLibraryCells) {
    if (!(c.e0MeV > 0.0) || c.beta < 0.0 || c.beta >= 1.0) {
      throw std::runtime_error("gamma library: bad cell (need E0>0, 0<=beta<1)");
    }
  }
  RunConfig::Library().nCells = static_cast<int>(fLibraryCells.size());
  RunConfig::Library().eventsPerCell = fEventsPerCell;
}

std::vector<PrimaryGeneratorAction::LibraryCell> PrimaryGeneratorAction::LoadLibraryCells(
    const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("gamma library: cannot open cells file '" + path + "'");
  std::vector<LibraryCell> cells;
  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream ss(line);
    LibraryCell c{};
    if (!(ss >> c.e0MeV)) continue;  // blank/comment line
    if (!(ss >> c.zCm >> c.beta)) {
      throw std::runtime_error("gamma library: line " + std::to_string(lineNo) + " of '" + path +
                               "' needs 'E0_MeV z_cm beta'");
    }
    cells.push_back(c);
  }
  return cells;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
  delete fReaction;
  delete fStoppingPower;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  RunConfig::CurrentTruth().Clear();
  if (fIsLibrary) {
    GenerateLibraryPrimaries(event);
    return;
  }
  if (fIsReaction) {
    GenerateReactionPrimaries(event);
    return;
  }
  if (fNeedsIonLookup) {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);
    fGun->SetParticleDefinition(ion);
    fGun->SetParticleCharge(fIonChargeState * eplus);
    fNeedsIonLookup = false;
  }
  fGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::GenerateReactionPrimaries(G4Event* event) {
  double vertexZCm = fReactionZ0Cm;
  double beamKineticEnergyMeV = fReaction->BeamKineticEnergyMeV();

  // BKIN given: track the beam's own real energy loss through the target
  // gas (see ReactionConfig.hh/GasStoppingPower.hh) instead of firing
  // every event from the same fixed point at the same fixed (energy-
  // loss-free) beam energy.
  if (fConfig.beamEntranceKineticEnergyMeV > 0.0) {
    if (!fStoppingPowerReady) {
      G4Material* targetGas = G4Material::GetMaterial(TargetChamber::TargetGasMaterialName());
      fStoppingPower =
          new GasStoppingPower(fConfig.beam.Z, fConfig.beam.A, targetGas,
                                fConfig.beamEntranceKineticEnergyMeV,
                                2.0 * TargetChamber::BeamPathHalfLengthCm());
      fStoppingPowerReady = true;
    }
    beamKineticEnergyMeV = fReaction->SampleBeamKineticEnergyMeV(
        fStoppingPower->ExitKineticEnergyMeV(), fStoppingPower->EntranceKineticEnergyMeV());
    const double meanDepthCm = fStoppingPower->DepthForKineticEnergyCm(beamKineticEnergyMeV);

    // Straggling (see GasStoppingPower's own header comment): the mean
    // dE/dx curve alone says every ion of the same entrance energy
    // reaches beamKineticEnergyMeV at exactly meanDepthCm -- real energy
    // loss is a statistical process, so a real ion actually reaches it
    // at a depth that fluctuates around that mean. Model this by
    // perturbing the *depth-crossing* energy (not the reaction's own
    // kinematics, which stay at the Breit-Wigner-sampled
    // beamKineticEnergyMeV above -- that's still the energy the reaction
    // itself happens at) by a Gaussian sample with Bohr's straggling
    // sigma at that depth, then re-finding the depth a real ion reaches
    // *that* energy at.
    const double sigmaMeV = fStoppingPower->StragglingSigmaMeV(meanDepthCm);
    double straggledKeMeV = beamKineticEnergyMeV + G4RandGauss::shoot(0.0, sigmaMeV);
    straggledKeMeV = std::min(std::max(straggledKeMeV, fStoppingPower->ExitKineticEnergyMeV()),
                               fStoppingPower->EntranceKineticEnergyMeV());
    const double depthCm = fStoppingPower->DepthForKineticEnergyCm(straggledKeMeV);
    vertexZCm = -TargetChamber::BeamPathHalfLengthCm() + depthCm;
  }

  const ReactionKinematics::Event ev = fReaction->GenerateEvent(beamKineticEnergyMeV);

  const G4ThreeVector vertexPos(fReactionX0Cm * cm, fReactionY0Cm * cm, vertexZCm * cm);
  auto* vertex = new G4PrimaryVertex(vertexPos, 0.0);

  G4ParticleDefinition* recoilIon = G4IonTable::GetIonTable()->GetIon(fRecoilZ, fRecoilA, 0.0);
  G4ThreeVector recoilP = ev.recoilMomentumMeV;
  if (fRecoilScan) {
    // Acceptance scan: the compound momentum (all along +z), rescaled by
    // (1+delta) and pointed along (theta_x, theta_y); no gammas.
    const double pCompound = std::sqrt(std::pow(beamKineticEnergyMeV + fReaction->BeamMassMeV(), 2) -
                                       std::pow(fReaction->BeamMassMeV(), 2));
    const double tx = fScanThetaMax * (2.0 * G4UniformRand() - 1.0);
    const double ty = fScanThetaMax * (2.0 * G4UniformRand() - 1.0);
    const double delta = fScanDeltaMax * (2.0 * G4UniformRand() - 1.0);
    recoilP = G4ThreeVector(tx, ty, 1.0).unit() * (pCompound * (1.0 + delta));
  }
  auto* recoil = new G4PrimaryParticle(recoilIon, recoilP.x() * MeV, recoilP.y() * MeV,
                                        recoilP.z() * MeV);
  recoil->SetCharge(fRecoilChargeState * eplus);
  vertex->SetPrimary(recoil);

  // One gamma per cascade step actually sampled this event (see
  // ReactionKinematics::GenerateEvent) -- usually 1, but 2 for this
  // pilot's own reaction whenever the resonance decays through its
  // level-1 intermediate state first (see reactions/o15ag_19ne.reaction).
  for (const G4ThreeVector& pGamma : ev.gammaMomentaMeV) {
    if (fRecoilScan) break;
    auto* gamma = new G4PrimaryParticle(G4Gamma::Gamma(), pGamma.x() * MeV, pGamma.y() * MeV,
                                         pGamma.z() * MeV);
    vertex->SetPrimary(gamma);
  }

  event->AddPrimaryVertex(vertex);

  RunConfig::Truth& truth = RunConfig::CurrentTruth();
  truth.mode = 0;
  truth.vertexXCm = fReactionX0Cm;
  truth.vertexYCm = fReactionY0Cm;
  truth.vertexZCm = vertexZCm;
  truth.beamKeMeV = beamKineticEnergyMeV;
  truth.ecmMeV = ev.ecmMeV;
  truth.emitterBeta = ev.compoundBeta;
  truth.recoilPMeV = recoilP;
  truth.recoilEMeV = std::sqrt(recoilP.mag2() + std::pow(fReaction->RecoilGroundMassMeV(), 2));
  for (std::size_t i = 0; i < (fRecoilScan ? 0 : ev.gammaMomentaMeV.size()); ++i) {
    const G4ThreeVector& p = ev.gammaMomentaMeV[i];
    truth.gammaELabMeV.push_back(static_cast<float>(p.mag()));
    truth.gammaCosThetaLab.push_back(static_cast<float>(p.cosTheta()));
    truth.gammaPhiLab.push_back(static_cast<float>(p.phi()));
    truth.gammaEStarMeV.push_back(static_cast<float>(ev.gammaEnergyStarMeV[i]));
    truth.gammaCosThetaStar.push_back(static_cast<float>(ev.gammaDirectionStar[i].cosTheta()));
    truth.gammaPhiStar.push_back(static_cast<float>(ev.gammaDirectionStar[i].phi()));
    truth.levelFrom.push_back(ev.levelFrom[i]);
    truth.levelTo.push_back(ev.levelTo[i]);
  }
}

void PrimaryGeneratorAction::GenerateLibraryPrimaries(G4Event* event) {
  const int cellIndex = event->GetEventID() / fEventsPerCell;
  if (cellIndex >= static_cast<int>(fLibraryCells.size())) {
    throw std::runtime_error("gamma library: more events requested than nCells*eventsPerCell");
  }
  const LibraryCell& cell = fLibraryCells[cellIndex];

  // Isotropic in the emitter rest frame, boosted along +z.
  const double cosThetaStar = 2.0 * G4UniformRand() - 1.0;
  const double sinThetaStar = std::sqrt(std::max(0.0, 1.0 - cosThetaStar * cosThetaStar));
  const double phiStar = 2.0 * CLHEP::pi * G4UniformRand() - CLHEP::pi;  // [-pi, pi), like G4ThreeVector::phi
  const double beta = cell.beta;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  const double eStar = cell.e0MeV;
  const double px = eStar * sinThetaStar * std::cos(phiStar);
  const double py = eStar * sinThetaStar * std::sin(phiStar);
  const double pz = gamma * eStar * (cosThetaStar + beta);
  const G4ThreeVector pLab(px, py, pz);

  const G4ThreeVector vertexPos(fLibraryX0Cm * cm, fLibraryY0Cm * cm, cell.zCm * cm);
  auto* vertex = new G4PrimaryVertex(vertexPos, 0.0);
  vertex->SetPrimary(new G4PrimaryParticle(G4Gamma::Gamma(), px * MeV, py * MeV, pz * MeV));
  event->AddPrimaryVertex(vertex);

  RunConfig::Truth& truth = RunConfig::CurrentTruth();
  truth.mode = 1;
  truth.cellIndex = cellIndex;
  truth.vertexXCm = fLibraryX0Cm;
  truth.vertexYCm = fLibraryY0Cm;
  truth.vertexZCm = cell.zCm;
  truth.emitterBeta = beta;
  truth.gammaELabMeV.push_back(static_cast<float>(pLab.mag()));
  truth.gammaCosThetaLab.push_back(static_cast<float>(pLab.cosTheta()));
  truth.gammaPhiLab.push_back(static_cast<float>(pLab.phi()));
  truth.gammaEStarMeV.push_back(static_cast<float>(eStar));
  truth.gammaCosThetaStar.push_back(static_cast<float>(cosThetaStar));
  truth.gammaPhiStar.push_back(static_cast<float>(phiStar));
  truth.levelFrom.push_back(-1);
  truth.levelTo.push_back(0);
}
