#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "config.hpp"
#include "mfem.hpp"
#include "source.hpp"

#include <string>
#include <vector>

static const char* DEFAULT_FILE_NAME = "no-file";

class ReceiversSet;



class Parameters
{
public:
  Parameters();
  ~Parameters();

  void init(int argc, char **argv);

  double sx, sy, sz; ///< size of the computational domain
  int nx, ny, nz; ///< number of cells in x- and y-directions
  double T; ///< simulation time
  double dt; ///< time step
  int order; ///< finite element order

  double rho, vp, vs; ///< homogeneous media properties

  const char *rhofile; ///< heterogeneous media properties
  const char *vpfile;
  const char *vsfile;

  double *rho_array, *vp_array, *vs_array; ///< arrays of values describing
                                           ///< media properties

  double damp_layer; ///< thickness of a damping layer
  double damp_power; ///< power in damping coefficient functions
  int topsurf; ///< top surface: 0 absorbing, 1 free

  Source source; ///< source of the wave

  int step_snap; ///< time step for outputting snapshots (every *th time step)
  int snapshot_format; ///< 0 - binary, 1 - VTS

  const char* method; ///< method to be used - FEM or SEM

  const char *extra_string; ///< for naming output results

  const char *receivers_file; ///< file describing sets of receivers
  std::vector<ReceiversSet*> sets_of_receivers; ///< receivers


private:
  Parameters(const Parameters&); // no copies
  Parameters& operator=(const Parameters&); // no copies
};

#endif // PARAMETERS_HPP
