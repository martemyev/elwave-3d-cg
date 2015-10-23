#ifndef ELASTIC_WAVE3D_HPP
#define ELASTIC_WAVE3D_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <vector>

class Parameters;
class Source;



/**
 * Cell-wise constant coefficient.
 */
class CWConstCoefficient : public mfem::Coefficient
{
public:
  CWConstCoefficient(double *array, bool own = 1)
    : val_array(array), own_array(own)
  { }

  virtual ~CWConstCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    return val_array[T.ElementNo];
  }

protected:
  double *val_array;
  bool own_array;
};



/**
 * A coefficient obtained with multiplication of a cell-wise constant
 * coefficient and a function.
 */
class CWFunctionCoefficient : public CWConstCoefficient
{
public:
  CWFunctionCoefficient(double(*F)(const mfem::Vector&, const Parameters&),
                        const Parameters& _param,
                        double *array, bool own = 1)
    : CWConstCoefficient(array, own)
    , Function(F)
    , param(_param)
  { }
  virtual ~CWFunctionCoefficient() { }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    const double cw_coef = val_array[T.ElementNo];
    double x[3];
    mfem::Vector transip(x, 3);
    T.Transform(ip, transip);
    const double func_val = (*Function)(transip, param);
    return cw_coef * func_val;
  }

protected:
  double(*Function)(const mfem::Vector&, const Parameters&);
  const Parameters& param;
};



/**
 * Implementation of a vector point force as a component of a source term.
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
 * Implementation of a moment tensor density as a component of a source term.
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




class ElasticWave2D
{
public:
  ElasticWave2D(const Parameters& _param);
  ~ElasticWave2D();

  void run();

private:
  const Parameters& param;

  /**
   * Finite Element Method (FEM) (non-diagonal mass matrix) with Absorbing
   * Layers by Increasing Damping (ALID) for implementation of absorbing
   * boundary condition.
   */
  void run_FEM_ALID();

  /**
   * Spectral Element Method (SEM) (diagonal mass matrix) with Stiffness
   * Reduction Method (SRM) for implementation of absorbing boundary condition.
   */
  void run_SEM_SRM();
};



mfem::Vector compute_solution_at_points(const std::vector<mfem::Vertex>& points,
                                        const std::vector<int>& cells_containing_points,
                                        const mfem::GridFunction& U);

void show_SRM_damp_weights(const Parameters& param);

void cells_containing_vertices(const mfem::Mesh& mesh, int nx, int ny,
                               double sx, double sy, std::vector<int>& cells);

mfem::Vector get_nodal_values(const std::vector<int>& cells,
                              const mfem::Mesh& mesh,
                              const mfem::GridFunction& U, int vdim);

#endif // ELASTIC_WAVE3D_HPP
