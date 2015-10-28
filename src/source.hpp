#ifndef SOURCE_HPP
#define SOURCE_HPP

#include "config.hpp"
#include "mfem.hpp"

class Source
{
public:
  Source();
  ~Source() { }

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
  double GaussFirstDerivative(double t) const;
  void PointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void MomentTensorSource(const mfem::Vector& x, mfem::Vector& f) const;

private:
  void compute_direction_vector();
  void DeltaPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void GaussPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void DivDeltaMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
  void DivGaussMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
};

#endif // SOURCE_HPP
