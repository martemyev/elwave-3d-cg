#include "elastic_wave3D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
#include "utilities.hpp"

#include <fstream>

using namespace std;
using namespace mfem;



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
VectorPointForce::VectorPointForce(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{ }

void VectorPointForce::Eval(Vector &V, ElementTransformation &T,
                            const IntegrationPoint &ip)
{
  double x[3];
  Vector transip(x, 3);
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.PointForce(transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
MomentTensorSource::MomentTensorSource(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{ }

void MomentTensorSource::Eval(Vector &V, ElementTransformation &T,
                              const IntegrationPoint &ip)
{
  double x[3];
  Vector transip(x, 3);
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.MomentTensorSource(transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
ElasticWave2D::ElasticWave2D(const Parameters &_param)
  : param(_param)
{ }

ElasticWave2D::~ElasticWave2D() { }

void ElasticWave2D::run()
{
  if (param.method == 0) // FEM
    run_FEM_ALID();
  else if (param.method == 1) // SEM
    run_SEM_SRM();
  else MFEM_ABORT("Unknown method to be used");
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Vector compute_solution_at_points(const vector<Vertex>& points,
                                  const vector<int>& cells_containing_points,
                                  const GridFunction& U)
{
  MFEM_ASSERT(points.size() == cells_containing_points.size(), "Sizes mismatch");
  const int vdim = U.VectorDim();
  MFEM_ASSERT(vdim == N_ELAST_COMPONENTS, "Dimensions mismatch");
  Vector U_at_points(vdim*points.size());
  Vector values(vdim);
  IntegrationPoint ip;
  for (size_t p = 0; p < points.size(); ++p)
  {
    ip.x = points[p](0);
    ip.y = points[p](1);
    ip.z = points[p](2);
    U.GetVectorValue(cells_containing_points[p], ip, values);
    MFEM_ASSERT(values.Size() == vdim, "Unexpected vector size");
    for (int d = 0; d < vdim; ++d)
      U_at_points(vdim*p+d) = values(d);
  }
  return U_at_points;
}

void cells_containing_vertices(const Mesh& mesh, int nx, int ny, double sx,
                               double sy, vector<int>& cells)
{
  cells.resize(mesh.GetNV());

  for (int v = 0; v < mesh.GetNV(); ++v)
  {
    const double *vert_coord = mesh.GetVertex(v);
    const Vertex vert(vert_coord[0], vert_coord[1], vert_coord[2]);
    cells[v] = find_element(nx, ny, sx, sy, vert, true);
  }
}

Vector get_nodal_values(const vector<int>& cells, const Mesh& mesh,
                        const GridFunction& U, int vdim_select)
{
  const int nv = mesh.GetNV();
  Vector nodal_values(nv);
  const int vdim = U.VectorDim();
  MFEM_ASSERT(vdim == N_ELAST_COMPONENTS, "Dimensions mismatch");
  Vector values(vdim);
  IntegrationPoint ip;

  for (int v = 0; v < nv; ++v)
  {
    const double *vert = mesh.GetVertex(v);
    ip.x = vert[0];
    ip.y = vert[1];
    ip.z = vert[2];
    U.GetVectorValue(cells[v], ip, values);
    nodal_values(v) = values(vdim_select-1);
  }

  return nodal_values;
}



