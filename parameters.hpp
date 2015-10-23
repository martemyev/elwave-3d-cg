#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <string>
#include <vector>

static const char* DEFAULT_FILE_NAME = "no-file";

class ReceiversSet;


class Source
{
public:
  Source();
  ~Source();

  mfem::Vertex location;
  double frequency;
  /**
   * The direction of the point force component of the source. The values 1-3
   * mean the direction along one of the axis: 1 OX, 2 OY, 3 OZ.
   */
  int direction;
  double ricker_scale;
  double point_force_scale;
  double gauss_support;
  double Mxx, Mxy, Mxz, Myy, Myz, Mzz; // components of a moment tensor
  int type; // 0 - Delta function, 1 - Gaussian function

  void AddOptions(mfem::OptionsParser& args);

  double Ricker(double t) const;
  void PointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void MomentTensorSource(const mfem::Vector& x, mfem::Vector& f) const;

private:
  void compute_direction_vector();
  void DeltaPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void GaussPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void DivDeltaMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
  void DivGaussMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
};



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

  /**
   * arrays of values describing media properties
   */
  double *rho_array, *vp_array, *vs_array;

  double damp_layer; ///< thickness of a damping layer
  double damp_power; ///< power in damping coefficient functions
  int topsurf; ///< top surface: 0 absorbing, 1 free

  Source source; ///< source of the wave

  int step_snap; ///< time step for outputting snapshots (every *th time step)

  int method; ///< 0 - FEM, 1 - SEM

  const char *extra_string; ///< for naming output results

  const char *receivers_file; ///< file describing sets of receivers
  std::vector<ReceiversSet*> sets_of_receivers; ///< receivers


private:
  Parameters(const Parameters&); // no copies
  Parameters& operator=(const Parameters&); // no copies
};

#endif // PARAMETERS_HPP
