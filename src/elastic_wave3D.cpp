#include "elastic_wave3D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
#include "snapshots.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;



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
Vector compute_function_at_point(const Mesh& mesh, const Vertex& point, int cell,
                                 const GridFunction& U)
{
  const Element* element = mesh.GetElement(cell);
  const Hexahedron *hex = dynamic_cast<const Hexahedron*>(element);
  MFEM_VERIFY(hex, "The mesh element has to be a hexahedron");

  Array<int> vert_indices;
  hex->GetVertices(vert_indices);

  const double *vert0 = mesh.GetVertex(vert_indices[0]); // min coords
  double x0 = vert0[0], x1 = vert0[0];
  double y0 = vert0[1], y1 = vert0[1];
  double z0 = vert0[2], z1 = vert0[2];
  for (int i = 1; i < hex->GetNVertices(); ++i)
  {
    const double *vert = mesh.GetVertex(vert_indices[i]);
    x0 = std::min(x0, vert[0]);
    x1 = std::max(x1, vert[0]);
    y0 = std::min(y0, vert[1]);
    y1 = std::max(y1, vert[1]);
    z0 = std::min(z0, vert[2]);
    z1 = std::max(z1, vert[2]);
  }

  const double hx = x1 - x0;
  const double hy = y1 - y0;
  const double hz = z1 - z0;

  MFEM_VERIFY(hx > 0 && hy > 0 && hz > 0, "Size of the hex is wrong");

  IntegrationPoint ip;
  ip.x = (point(0) - x0) / hx; // transfer to the reference space [0,1]^3
  ip.y = (point(1) - y0) / hy;
  ip.z = (point(2) - z0) / hz;

  Vector values;
  U.GetVectorValue(cell, ip, values);

  return values;
}



Vector compute_function_at_points(const Mesh& mesh, int n_points,
                                  const Vertex *points,
                                  const int *cells_containing_points,
                                  const GridFunction& U)
{
//  MFEM_ASSERT(points.size() == cells_containing_points.size(), "Sizes mismatch");
  Vector U_at_points(N_ELAST_COMPONENTS*n_points);

  for (int p = 0; p < n_points; ++p)
  {
    Vector values = compute_function_at_point(mesh, points[p],
                                              cells_containing_points[p], U);
    MFEM_ASSERT(values.Size() == N_ELAST_COMPONENTS, "Unexpected vector size");
    for (int c = 0; c < N_ELAST_COMPONENTS; ++c)
      U_at_points(p*N_ELAST_COMPONENTS+c) = values(c);
  }
  return U_at_points;
}



void open_seismo_outs(ofstream* &seisU, ofstream* &seisV,
                      const Parameters &param, const string &method_name)
{
  const int n_rec_sets = param.sets_of_receivers.size();

  seisU = new ofstream[N_ELAST_COMPONENTS*n_rec_sets];
  seisV = new ofstream[N_ELAST_COMPONENTS*n_rec_sets];

  for (int r = 0; r < n_rec_sets; ++r)
  {
    const ReceiversSet *rec_set = param.sets_of_receivers[r];
    const string desc = rec_set->description();
    const string variable = rec_set->get_variable();

    if (variable.find("U") != string::npos) {
      for (int c = 0; c < N_ELAST_COMPONENTS; ++c) {
        string seismofile = (string)param.output_dir + "/" + SEISMOGRAMS_DIR +
                            method_name + param.extra_string + desc + "_u" +
                            d2s(c) + ".bin";
        seisU[r*N_ELAST_COMPONENTS + c].open(seismofile.c_str(), ios::binary);
        MFEM_VERIFY(seisU[r*N_ELAST_COMPONENTS + c], "File '" + seismofile +
                    "' can't be opened");
      }
    }

    if (variable.find("V") != string::npos) {
      for (int c = 0; c < N_ELAST_COMPONENTS; ++c) {
        string seismofile = (string)param.output_dir + "/" + SEISMOGRAMS_DIR +
                            method_name + param.extra_string + desc + "_v" +
                            d2s(c) + ".bin";
        seisV[r*N_ELAST_COMPONENTS + c].open(seismofile.c_str(), ios::binary);
        MFEM_VERIFY(seisV[r*N_ELAST_COMPONENTS + c], "File '" + seismofile +
                    "' can't be opened");
      }
    }
  } // loop for sets of receivers
}



void output_snapshots(int time_step, const string& snapshot_filebase,
                      const Parameters& param, const GridFunction& U,
                      const GridFunction& V, const Mesh& mesh)
{
  for (size_t sn = 0; sn < param.sets_of_snapshots.size(); ++sn)
  {
    const SnapshotsSet *snap_set = param.sets_of_snapshots[sn];
    const string desc = snap_set->description();
    const string variable = snap_set->get_variable();

    if (variable.find("U") != string::npos) {
      // displacement at the snapshot points
      const Vector u =
        compute_function_at_points(mesh, snap_set->n_snapshot_points(),
                                   snap_set->get_snapshot_points(),
                                   snap_set->get_cells_containing_snapshot_points(),
                                   U);
      snap_set->save_snapshot_vector(u, time_step, "U_CG", N_ELAST_COMPONENTS,
                                     snapshot_filebase + desc);
    }

    if (variable.find("V") != string::npos) {
      // particle velocity at the snapshot points
      const Vector v =
        compute_function_at_points(mesh, snap_set->n_snapshot_points(),
                                   snap_set->get_snapshot_points(),
                                   snap_set->get_cells_containing_snapshot_points(),
                                   V);
      snap_set->save_snapshot_vector(v, time_step, "V_CG", N_ELAST_COMPONENTS,
                                     snapshot_filebase + desc);
    }
  }

//  Vector u_x, u_y, u_z, v_x, v_y, v_z;
//  U.GetNodalValues(u_x, 1); // components of displacement
//  U.GetNodalValues(u_y, 2);
//  U.GetNodalValues(u_z, 3);
//  V.GetNodalValues(v_x, 1); // components of velocity
//  V.GetNodalValues(v_y, 2);
//  V.GetNodalValues(v_z, 3);

//  string tstep = d2s(time_step,0,0,0,6), fname;
//  if (param.snapshot_format == 0) // binary format
//  {
//    const int n_values = u_x.Size(); // the same for all components
//    fname = snapshot_filebase + "_Ux_t" + tstep + ".bin";
//    write_binary(fname.c_str(), n_values, u_x);
//    fname = snapshot_filebase + "_Uy_t" + tstep + ".bin";
//    write_binary(fname.c_str(), n_values, u_y);
//    fname = snapshot_filebase + "_Uz_t" + tstep + ".bin";
//    write_binary(fname.c_str(), n_values, u_z);
//    fname = snapshot_filebase + "_Vx_t" + tstep + ".bin";
//    write_binary(fname.c_str(), n_values, v_x);
//    fname = snapshot_filebase + "_Vy_t" + tstep + ".bin";
//    write_binary(fname.c_str(), n_values, v_y);
//    fname = snapshot_filebase + "_Vz_t" + tstep + ".bin";
//    write_binary(fname.c_str(), n_values, v_z);
//  }
//  else // VTS format
//  {
//    fname = snapshot_filebase + "_U_t" + tstep + ".vts";
//    write_vts_vector(fname, "U", param.grid.sx, param.grid.sy, param.grid.sz,
//                     param.grid.nx, param.grid.ny, param.grid.nz, u_x, u_y, u_z);
//    fname = snapshot_filebase + "_V_t" + tstep + ".vts";
//    write_vts_vector(fname, "V", param.grid.sx, param.grid.sy, param.grid.sz,
//                     param.grid.nx, param.grid.ny, param.grid.nz, v_x, v_y, v_z);
//  }
}



void output_seismograms(const Parameters& param, const Mesh& mesh,
                        const GridFunction &U, const GridFunction &V,
                        ofstream* &seisU, ofstream* &seisV)
{
  // for each set of receivers
  for (size_t rec = 0; rec < param.sets_of_receivers.size(); ++rec)
  {
    const ReceiversSet *rec_set = param.sets_of_receivers[rec];
    const string variable = rec_set->get_variable();

    // Displacement
    if (variable.find("U") != string::npos) {
      for (int c = 0; c < N_ELAST_COMPONENTS; ++c) {
        MFEM_VERIFY(seisU[rec*N_ELAST_COMPONENTS+c].is_open(), "The stream for "
                    "writing displacement seismograms is not open");
      }
      // displacement at the receivers
      const Vector u =
        compute_function_at_points(mesh, rec_set->n_receivers(),
                                   &(rec_set->get_receivers()[0]),
                                   &(rec_set->get_cells_containing_receivers()[0]), U);
      MFEM_ASSERT(u.Size() == N_ELAST_COMPONENTS*rec_set->n_receivers(),
                  "Sizes mismatch");
      for (int i = 0; i < u.Size(); i += N_ELAST_COMPONENTS) {
        for (int j = 0; j < N_ELAST_COMPONENTS; ++j) {
          float val = u(i+j); // displacement
          seisU[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val),
                                                  sizeof(val));
        }
      }
    }

    // Particle velocity
    if (variable.find("V") != string::npos) {
      for (int c = 0; c < N_ELAST_COMPONENTS; ++c) {
        MFEM_VERIFY(seisV[rec*N_ELAST_COMPONENTS+c].is_open(), "The stream for "
                    "writing velocity seismograms is not open");
      }
      // velocity at the receivers
      const Vector v =
        compute_function_at_points(mesh, rec_set->n_receivers(),
                                   &(rec_set->get_receivers()[0]),
                                   &(rec_set->get_cells_containing_receivers()[0]), V);
      MFEM_ASSERT(v.Size() == N_ELAST_COMPONENTS*rec_set->n_receivers(),
                  "Sizes mismatch");
      for (int i = 0; i < v.Size(); i += N_ELAST_COMPONENTS) {
        for (int j = 0; j < N_ELAST_COMPONENTS; ++j) {
          float val = v(i+j); // velocity
          seisV[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val),
                                                  sizeof(val));
        }
      }
    }
  } // loop over receiver sets
}



