#ifndef SOURCE_HPP
#define SOURCE_HPP

#include "config.hpp"
#include "mfem.hpp"

class Parameters;



/**
 * Parameters describing the source and related functions.
 */
class Source
{
public:
  Source();
  ~Source() { }

  enum SourceType { POINT_FORCE, MOMENT_TENSOR, PLANE_WAVE };

  mfem::Vertex location;
  double frequency;
  int direction; ///< The direction of the point force component of the source.
                 ///< The values 1-3 mean the direction along one of the axis:
                 ///< 1 OX, 2 OY, 3 OZ.
  double scale;
  double Mxx, Mxy, Mxz, Myy, Myz, Mzz; ///< components of a moment tensor
  const char *type_string; ///< "pointforce", "momenttensor", "planewave"
  const char *spatial_function; ///< "delta", "gauss"
  double gauss_support; ///< size of the support for the "gauss" spatial function

  SourceType type;

  void AddOptions(mfem::OptionsParser& args);
  void check_and_update_parameters();

  double Ricker(double t) const;
  double GaussFirstDerivative(double t) const;
  void PointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void MomentTensorSource(const mfem::Vector& x, mfem::Vector& f) const;
  void PlaneWaveSource(const Parameters& param, const mfem::Vector& x,
                       mfem::Vector& f) const;

private:
  void DeltaPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void GaussPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void DivDeltaMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
  void DivGaussMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
};



/**
 * Implementation of a vector point force type of source.
 */
class VectorPointForce: public mfem::VectorCoefficient
{
public:
  VectorPointForce(int dim, const Source& s);
  ~VectorPointForce() { }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip);

private:
  const Source& source;
};



/**
 * Implementation of a moment tensor type of source.
 */
class MomentTensorSource: public mfem::VectorCoefficient
{
public:
  MomentTensorSource(int dim, const Source& s);
  ~MomentTensorSource() { }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip);

private:
  const Source& source;
};



/**
 * Implementation of a plane wave type of source.
 */
class PlaneWaveSource: public mfem::VectorCoefficient
{
public:
  PlaneWaveSource(int dim, const Parameters& p);
  ~PlaneWaveSource() { }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip);

private:
  const Parameters& param;
};

#endif // SOURCE_HPP
