#include "elastic_wave3D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
#include "utilities.hpp"

#include <fstream>

using namespace std;
using namespace mfem;



ElasticWave2D::ElasticWave2D(const Parameters &_param)
  : param(_param)
{ }



void ElasticWave2D::run()
{
  if (!strcmp(param.method, "fem") || !strcmp(param.method, "FEM"))
    run_FEM_ALID();
  else if (!strcmp(param.method, "sem") || !strcmp(param.method, "SEM"))
    run_SEM_SRM();
  else
    MFEM_ABORT("Unknown method to be used: " + string(param.method));
}



//------------------------------------------------------------------------------
//
// Auxiliary useful functions
//
//------------------------------------------------------------------------------
Vector compute_function_at_point(double sx, double sy, double sz,
                                 int nx, int ny, int nz, const Mesh& mesh,
                                 const Vertex& point, int cell,
                                 const GridFunction& U)
{
  const double hx = sx / nx;
  const double hy = sy / ny;
  const double hz = sz / nz;

  const Element *element = mesh.GetElement(cell);
  Array<int> vert_indices;
  element->GetVertices(vert_indices);
  const double *vert0 = mesh.GetVertex(vert_indices[0]); // min coords

  IntegrationPoint ip;
  ip.x = (point(0) - vert0[0]) / hx; // transfer to the reference space [0,1]^3
  ip.y = (point(1) - vert0[1]) / hy;
  ip.z = (point(2) - vert0[2]) / hz;

  Vector values;
  U.GetVectorValue(cell, ip, values);

  return values;
}



Vector compute_function_at_points(double sx, double sy, double sz,
                                  int nx, int ny, int nz, const Mesh& mesh,
                                  const vector<Vertex>& points,
                                  const vector<int>& cells_containing_points,
                                  const GridFunction& U)
{
  MFEM_ASSERT(points.size() == cells_containing_points.size(), "Sizes mismatch");
  Vector U_at_points(N_ELAST_COMPONENTS*points.size());

  for (size_t p = 0; p < points.size(); ++p)
  {
    Vector values = compute_function_at_point(sx, sy, sz,
                                              nx, ny, nz, mesh, points[p],
                                              cells_containing_points[p], U);
    MFEM_ASSERT(values.Size() == N_ELAST_COMPONENTS, "Unexpected vector size");
    for (int c = 0; c < N_ELAST_COMPONENTS; ++c)
      U_at_points(p*N_ELAST_COMPONENTS+c) = values(c);
  }
  return U_at_points;
}



void output_snapshots(int time_step, const string& snapshot_filebase,
                      const Parameters& param, const GridFunction& U,
                      const GridFunction& V)
{
  Vector u_x, u_y, u_z, v_x, v_y, v_z;
  U.GetNodalValues(u_x, 1); // components of displacement
  U.GetNodalValues(u_y, 2);
  U.GetNodalValues(u_z, 3);
  V.GetNodalValues(v_x, 1); // components of velocity
  V.GetNodalValues(v_y, 2);
  V.GetNodalValues(v_z, 3);

  string tstep = d2s(time_step,0,0,0,6), fname;
  if (param.snapshot_format == 0) // binary format
  {
    const int n_values = u_x.Size(); // the same for all components
    fname = snapshot_filebase + "_Ux_t" + tstep + ".bin";
    write_binary(fname.c_str(), n_values, u_x);
    fname = snapshot_filebase + "_Uy_t" + tstep + ".bin";
    write_binary(fname.c_str(), n_values, u_y);
    fname = snapshot_filebase + "_Uz_t" + tstep + ".bin";
    write_binary(fname.c_str(), n_values, u_z);
    fname = snapshot_filebase + "_Vx_t" + tstep + ".bin";
    write_binary(fname.c_str(), n_values, v_x);
    fname = snapshot_filebase + "_Vy_t" + tstep + ".bin";
    write_binary(fname.c_str(), n_values, v_y);
    fname = snapshot_filebase + "_Vz_t" + tstep + ".bin";
    write_binary(fname.c_str(), n_values, v_z);
  }
  else // VTS format
  {
    fname = snapshot_filebase + "_U_t" + tstep + ".vts";
    write_vts_vector(fname, "U", param.sx, param.sy, param.sz, param.nx,
                     param.ny, param.nz, u_x, u_y, u_z);
    fname = snapshot_filebase + "_V_t" + tstep + ".vts";
    write_vts_vector(fname, "V", param.sx, param.sy, param.sz, param.nx,
                     param.ny, param.nz, v_x, v_y, v_z);
  }
}



void output_seismograms(const Parameters& param, const Mesh& mesh,
                        const GridFunction &U, const GridFunction &V,
                        ofstream *seisU, ofstream *seisV)
{
  // for each set of receivers
  for (size_t rec = 0; rec < param.sets_of_receivers.size(); ++rec)
  {
    for (int c = 0; c < N_ELAST_COMPONENTS; ++c)
    {
      MFEM_ASSERT(seisU[rec*N_ELAST_COMPONENTS+c].is_open(), "The stream for "
                  "writing displacement seismograms is not open");
      MFEM_ASSERT(seisV[rec*N_ELAST_COMPONENTS+c].is_open(), "The stream for "
                  "writing velocity seismograms is not open");
    }

    const ReceiversSet *rec_set = param.sets_of_receivers[rec];
    // displacement at the receivers
    const Vector u = compute_function_at_points(param.sx, param.sy, param.sz,
                                                param.nx, param.ny, param.nz, mesh,
                                                rec_set->get_receivers(),
                                                rec_set->get_cells_containing_receivers(),
                                                U);
    // velocity at the receivers
    const Vector v = compute_function_at_points(param.sx, param.sy, param.sz,
                                                param.nx, param.ny, param.nz, mesh,
                                                rec_set->get_receivers(),
                                                rec_set->get_cells_containing_receivers(),
                                                V);

    MFEM_ASSERT(u.Size() == N_ELAST_COMPONENTS*rec_set->n_receivers() &&
                u.Size() == v.Size(), "Sizes mismatch");

    float val;
    for (int i = 0; i < u.Size(); i += N_ELAST_COMPONENTS)
    {
      for (int j = 0; j < N_ELAST_COMPONENTS; ++j)
      {
        val = u(i+j); // displacement
        seisU[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val),
                                                sizeof(val));

        val = v(i+j); // velocity
        seisV[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val),
                                                sizeof(val));
      }
    }
  }
}



