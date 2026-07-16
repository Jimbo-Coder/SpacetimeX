/* puncture.hxx */
/* (c) Liwei Ji 06/2024 */

#ifndef PUNCTURETRACKER_PUNCTURE_HXX
#define PUNCTURETRACKER_PUNCTURE_HXX

#include <loop_device.hxx>

#include <cctk.h>

#include <vector>

namespace PunctureTracker {

class PunctureContainer {
public:
  PunctureContainer() {}

  virtual ~PunctureContainer() = default; // Use default for trivial destructors

  void setNumPunctures() { numPunctures_ = location_[0].size(); }

  CCTK_INT getNumPunctures() { return numPunctures_; }

  std::vector<CCTK_REAL> &getTime() { return time_; }

  std::vector<CCTK_REAL> &getPreviousTime() { return previousTime_; }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getLocation() {
    return location_;
  }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getVelocity() {
    return velocity_;
  }

  std::vector<CCTK_REAL> &getMass() { return mass_; }

  std::vector<CCTK_REAL> &getEtaWeight() { return etaWeight_; }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getGroupLocation() {
    return groupLocation_;
  }

  std::vector<CCTK_REAL> &getGroupMass() { return groupMass_; }

  std::vector<CCTK_REAL> &getGroupEtaWeight() { return groupEtaWeight_; }

  std::vector<CCTK_INT> &getGroupMembership() { return groupMembership_; }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getBeta() { return beta_; }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getPreviousBeta() {
    return previousBeta_;
  }

  void updatePreviousTime(CCTK_ARGUMENTS);

  void interpolate(CCTK_ARGUMENTS);

  void evolve(CCTK_ARGUMENTS);

  void broadcast(CCTK_ARGUMENTS);

  void updateGroups(const bool trackMergers,
                    const CCTK_REAL mergerDistanceCoefficient);

private:
  CCTK_INT numPunctures_;
  std::vector<CCTK_REAL> time_, previousTime_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> location_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> velocity_;
  std::vector<CCTK_REAL> mass_, etaWeight_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> groupLocation_;
  std::vector<CCTK_REAL> groupMass_, groupEtaWeight_;
  std::vector<CCTK_INT> groupMembership_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> beta_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> previousBeta_;
};

} // namespace PunctureTracker

#endif // #ifndef PUNCTURETRACKER_PUNCTURE_HXX
