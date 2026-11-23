#include "RunConfig.hh"

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace RunConfig {

namespace {
const char* Env(const char* name) {
  const char* v = std::getenv(name);
  return (v && *v) ? v : nullptr;
}
}  // namespace

std::string OutputFile() {
  const char* v = Env("OUTPUT_FILE");
  return v ? std::string(v) : std::string("dragon_hits.root");
}

OutputMode Mode() {
  static const OutputMode mode = [] {
    const char* v = Env("OUTPUT_MODE");
    if (!v) return OutputMode::kSteps;
    const std::string s(v);
    if (s == "steps") return OutputMode::kSteps;
    if (s == "compact") return OutputMode::kCompact;
    if (s == "both") return OutputMode::kBoth;
    throw std::runtime_error("OUTPUT_MODE must be steps, compact or both (got '" + s + "')");
  }();
  return mode;
}

bool WriteStepTrees() { return Mode() != OutputMode::kCompact; }
bool WriteCompactTrees() { return Mode() != OutputMode::kSteps; }

bool EventPrint() {
  static const bool print = [] {
    const char* v = Env("EVENT_PRINT");
    if (v) return std::atoi(v) != 0;
    return Mode() != OutputMode::kCompact;
  }();
  return print;
}

double GasPressureTorr() {
  static const double p = [] {
    const char* v = Env("GAS_PRESSURE_TORR");
    const double val = v ? std::atof(v) : 5.0;
    if (!(val > 0.0)) throw std::runtime_error("GAS_PRESSURE_TORR must be > 0");
    return val;
  }();
  return p;
}

double GasTemperatureK() {
  static const double t = [] {
    const char* v = Env("GAS_TEMPERATURE_K");
    if (!v) return 0.0;
    const double val = std::atof(v);
    if (!(val > 0.0)) throw std::runtime_error("GAS_TEMPERATURE_K must be > 0");
    return val;
  }();
  return t;
}

double TargetGasDensityFractionOfStp() {
  const double t = GasTemperatureK();
  const double tempFactor = t > 0.0 ? kNistReferenceTemperatureK / t : 1.0;
  return (GasPressureTorr() / 760.0) * tempFactor;
}

long RandomSeed() {
  const char* v = Env("RANDOM_SEED");
  return v ? std::atol(v) : -1;
}

void Truth::Clear() {
  *this = Truth();
}

Truth& CurrentTruth() {
  static Truth truth;
  return truth;
}

CompactBuffers& Buffers() {
  static CompactBuffers buffers;
  return buffers;
}

NtupleIds& Ids() {
  static NtupleIds ids;
  return ids;
}

LibraryInfo& Library() {
  static LibraryInfo info;
  return info;
}

}  // namespace RunConfig
