/* puncture.cxx */
/* (c) Liwei Ji 06/2024 */

#include "puncture.hxx"

#include <cctk_Parameters.h>
#include <util_Table.h>

#include <cmath>
#include <vector>

namespace PunctureTracker {

static int findRoot(std::vector<int> &parent, int n) {
  while (parent[n] != n) {
    parent[n] = parent[parent[n]];
    n = parent[n];
  }
  return n;
}

static void joinRoots(std::vector<int> &parent, const int a, const int b) {
  const int rootA = findRoot(parent, a);
  const int rootB = findRoot(parent, b);
  if (rootA != rootB) {
    parent[rootB] = rootA;
  }
}

void PunctureContainer::updatePreviousTime(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  for (int n = 0; n < numPunctures_; ++n) {
    previousTime_[n] = time_[n];
    time_[n] = cctk_time;
  }
}

void PunctureContainer::interpolate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  // Only processor 0 interpolates
  const CCTK_INT nPoints = (CCTK_MyProc(cctkGH) == 0) ? numPunctures_ : 0;

  // Interpolation coordinates
  const void *interpCoords[Loop::dim] = {
      location_[0].data(), location_[1].data(), location_[2].data()};

  // Interpolated variables
  const CCTK_INT nInputArrays = 3;
  const CCTK_INT inputArrayIndices[nInputArrays] = {
      CCTK_VarIndex("ADMBaseX::betax"), CCTK_VarIndex("ADMBaseX::betay"),
      CCTK_VarIndex("ADMBaseX::betaz")};

  CCTK_POINTER outputArrays[nInputArrays] = {beta_[0].data(), beta_[1].data(),
                                             beta_[2].data()};

  // DriverInterpolate arguments that aren't currently used
  const int coordSystemHandle = 0;
  const CCTK_INT interpCoordsTypeCode = 0;
  const CCTK_INT outputArrayTypes[1] = {0};

  const int interpHandle = CCTK_InterpHandle("CarpetX");
  if (interpHandle < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Can't get interpolation handle");
    return;
  }

  // Create parameter table for interpolation
  const int paramTableHandle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (paramTableHandle < 0) {
    CCTK_VERROR("Can't create parameter table: %d", paramTableHandle);
  }

  // Set interpolation order in the parameter table
  int ierr = Util_TableSetInt(paramTableHandle, interp_order, "order");
  if (ierr < 0) {
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  }

  if (CCTK_EQUALS(interp_level_mode, "max_level")) {
    if (interp_max_level < 0) {
      CCTK_ERROR("interp_level_mode='max_level' requires interp_max_level >= 0");
    }

    ierr = Util_TableSetInt(paramTableHandle, interp_max_level, "max_level");
    if (ierr < 0) {
      CCTK_VERROR("Can't set max_level in parameter table: %d", ierr);
    }
  }

  // Perform the interpolation
  ierr = DriverInterpolate(cctkGH, Loop::dim, interpHandle, paramTableHandle,
                           coordSystemHandle, nPoints, interpCoordsTypeCode,
                           interpCoords, nInputArrays, inputArrayIndices,
                           nInputArrays, outputArrayTypes, outputArrays);

  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Interpolation error");
  }

  // Destroy the parameter table
  Util_TableDestroy(paramTableHandle);
}

void PunctureContainer::evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0) {
    for (int n = 0; n < numPunctures_; ++n) {
      const CCTK_REAL dt = time_[n] - previousTime_[n];
      for (int i = 0; i < Loop::dim; ++i) {
        const CCTK_REAL beta_eff =
            CCTK_EQUALS(time_integrator, "trapezoidal")
                ? 0.5 * (beta_[i][n] + previousBeta_[i][n])
                : beta_[i][n];
        location_[i][n] += dt * (-beta_eff);
        velocity_[i][n] = -beta_[i][n];
        previousBeta_[i][n] = beta_[i][n];
      }
    }
  }
}

void PunctureContainer::broadcast(CCTK_ARGUMENTS) {
  const CCTK_INT numComponents = 6;
  // 3 components for location, 3 components for velocity
  std::vector<CCTK_REAL> buffer(numComponents * numPunctures_);

  if (CCTK_MyProc(cctkGH) == 0) {
    for (int i = 0; i < Loop::dim; ++i) {
      for (int n = 0; n < numPunctures_; ++n) {
        buffer[i * numPunctures_ + n] = location_[i][n];
        buffer[(i + Loop::dim) * numPunctures_ + n] = velocity_[i][n];
      }
    }
  }

  int mpiError =
      MPI_Bcast(buffer.data(), buffer.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (mpiError != MPI_SUCCESS) {
    CCTK_VINFO("MPI_Bcast failed with error code %d", mpiError);
    MPI_Abort(MPI_COMM_WORLD, mpiError);
  }

  for (int i = 0; i < Loop::dim; ++i) {
    for (int n = 0; n < numPunctures_; ++n) {
      location_[i][n] = buffer[i * numPunctures_ + n];
      velocity_[i][n] = buffer[(i + Loop::dim) * numPunctures_ + n];
    }
  }
}

void PunctureContainer::updateGroups(
    const bool trackMergers, const CCTK_REAL mergerDistanceCoefficient) {
  groupLocation_[0].clear();
  groupLocation_[1].clear();
  groupLocation_[2].clear();
  groupMass_.clear();
  groupEtaWeight_.clear();
  groupMembership_.assign(numPunctures_, -1);

  if (numPunctures_ <= 0) {
    return;
  }

  std::vector<int> parent(numPunctures_);
  for (int n = 0; n < numPunctures_; ++n) {
    parent[n] = n;
  }

  if (trackMergers) {
    for (int i = 0; i < numPunctures_; ++i) {
      for (int j = i + 1; j < numPunctures_; ++j) {
        CCTK_REAL distance2 = 0.0;
        for (int d = 0; d < Loop::dim; ++d) {
          const CCTK_REAL dx = location_[d][i] - location_[d][j];
          distance2 += dx * dx;
        }
        const CCTK_REAL distance = std::sqrt(distance2);
        const CCTK_REAL mergerDistance =
            mergerDistanceCoefficient * (mass_[i] + mass_[j]);
        if (distance < mergerDistance) {
          joinRoots(parent, i, j);
        }
      }
    }
  }

  std::vector<int> groupRoot;
  for (int n = 0; n < numPunctures_; ++n) {
    const int root = findRoot(parent, n);
    int group = -1;
    for (int g = 0; g < int(groupRoot.size()); ++g) {
      if (groupRoot[g] == root) {
        group = g;
        break;
      }
    }

    if (group < 0) {
      group = int(groupRoot.size());
      groupRoot.push_back(root);
      for (int d = 0; d < Loop::dim; ++d) {
        groupLocation_[d].push_back(0.0);
      }
      groupMass_.push_back(0.0);
      groupEtaWeight_.push_back(0.0);
    }

    const CCTK_REAL oldMass = groupMass_[group];
    const CCTK_REAL newMass = oldMass + mass_[n];
    for (int d = 0; d < Loop::dim; ++d) {
      groupLocation_[d][group] =
          (groupLocation_[d][group] * oldMass + location_[d][n] * mass_[n]) /
          newMass;
    }
    groupEtaWeight_[group] =
        (groupEtaWeight_[group] * oldMass + etaWeight_[n] * mass_[n]) /
        newMass;
    groupMass_[group] = newMass;
    groupMembership_[n] = group;
  }
}

} // namespace PunctureTracker
