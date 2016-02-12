#include "elastic_wave3D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"

#include <fstream>
#include <vector>

using namespace std;
using namespace mfem;

//#define OUTPUT_MASS_MATRIX

double mass_damp_weight(const mfem::Vector& point, const Parameters& param);
double stif_damp_weight(const mfem::Vector& point, const Parameters& param);



void ElasticWave2D::run_SEM_SRM()
{
  StopWatch chrono;

  chrono.Start();
  cout << "Mesh and FE space generation..." << flush;
  const bool generate_edges = true;
  Mesh mesh(param.grid.nx, param.grid.ny, param.grid.nz, Element::HEXAHEDRON,
            generate_edges, param.grid.sx, param.grid.sy, param.grid.sz);
  const int dim = mesh.Dimension();
  MFEM_VERIFY(dim == SPACE_DIM, "Unexpected mesh dimension");
  const int n_elements = param.grid.nx * param.grid.ny * param.grid.nz;
  MFEM_VERIFY(n_elements == mesh.GetNE(), "Unexpected number of mesh elements");

  FiniteElementCollection *fec = new H1_FECollection(param.order, dim);
  FiniteElementSpace fespace(&mesh, fec, dim); //, Ordering::byVDIM);
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "Number of unknowns: " << fespace.GetVSize() << endl;

  double *lambda_array = new double[n_elements];
  double *mu_array     = new double[n_elements];

  for (int i = 0; i < n_elements; ++i)
  {
    const double rho = param.media.rho_array[i];
    const double vp  = param.media.vp_array[i];
    const double vs  = param.media.vs_array[i];

    MFEM_VERIFY(rho > 1.0 && vp > 1.0 && vs > 1.0, "Incorrect media properties "
                "arrays");

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;
  }

  const bool own_array = false;
  CWConstCoefficient rho_coef(param.media.rho_array, own_array);
  CWFunctionCoefficient lambda_coef  (stif_damp_weight, param, lambda_array);
  CWFunctionCoefficient mu_coef      (stif_damp_weight, param, mu_array);
  CWFunctionCoefficient rho_damp_coef(mass_damp_weight, param,
                                      param.media.rho_array, own_array);

  IntegrationRule segment_GLL;
  create_segment_GLL_rule(param.order, segment_GLL);
  IntegrationRule hex_GLL(segment_GLL, segment_GLL, segment_GLL);

  cout << "Stif matrix..." << flush;
  ElasticityIntegrator *elast_int = new ElasticityIntegrator(lambda_coef, mu_coef);
  elast_int->SetIntRule(&hex_GLL);
  BilinearForm stif(&fespace);
  stif.AddDomainIntegrator(elast_int);
  stif.Assemble();
  stif.Finalize();
  const SparseMatrix& S = stif.SpMat();
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "Mass matrix..." << flush;
  VectorMassIntegrator *mass_int = new VectorMassIntegrator(rho_coef);
  mass_int->SetIntRule(&hex_GLL);
  BilinearForm mass(&fespace);
  mass.AddDomainIntegrator(mass_int);
  mass.Assemble();
  mass.Finalize();
  const SparseMatrix& M = mass.SpMat();
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

#if defined(OUTPUT_MASS_MATRIX)
  {
    cout << "Output mass matrix..." << flush;
    ofstream mout("mass_mat.dat");
    mass.PrintMatlab(mout);
    cout << "M.nnz = " << M.NumNonZeroElems() << endl;
    cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
    chrono.Clear();
  }
#endif

  cout << "Damp matrix..." << flush;
  VectorMassIntegrator *damp_int = new VectorMassIntegrator(rho_damp_coef);
  damp_int->SetIntRule(&hex_GLL);
  BilinearForm dampM(&fespace);
  dampM.AddDomainIntegrator(damp_int);
  dampM.Assemble();
  dampM.Finalize();
  SparseMatrix& D = dampM.SpMat();
  double omega = 2.0*M_PI*param.source.frequency; // angular frequency
  D *= 0.5*param.dt*omega;
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "RHS vector... " << flush;
  LinearForm b(&fespace);
  if (param.source.plane_wave)
  {
    PlaneWaveSource plane_wave_source(dim, param);
    VectorDomainLFIntegrator *plane_wave_int =
        new VectorDomainLFIntegrator(plane_wave_source);
    plane_wave_int->SetIntRule(&hex_GLL);
    b.AddDomainIntegrator(plane_wave_int);
    b.Assemble();
  }
  else
  {
    if (!strcmp(param.source.type, "pointforce"))
    {
      VectorPointForce vector_point_force(dim, param);
      VectorDomainLFIntegrator *point_force_int =
          new VectorDomainLFIntegrator(vector_point_force);
      point_force_int->SetIntRule(&hex_GLL);
      b.AddDomainIntegrator(point_force_int);
      b.Assemble();
    }
    else if (!strcmp(param.source.type, "momenttensor"))
    {
      MomentTensorSource momemt_tensor_source(dim, param);
      VectorDomainLFIntegrator *moment_tensor_int =
          new VectorDomainLFIntegrator(momemt_tensor_source);
      moment_tensor_int->SetIntRule(&hex_GLL);
      b.AddDomainIntegrator(moment_tensor_int);
      b.Assemble();
    }
    else MFEM_ABORT("Unknown source type: " + string(param.source.type));
  }
  cout << "||b||_L2 = " << b.Norml2() << endl;
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  Vector diagM; M.GetDiag(diagM); // mass matrix is diagonal
  Vector diagD; D.GetDiag(diagD); // damping matrix is diagonal

  for (int i = 0; i < diagM.Size(); ++i)
    MFEM_VERIFY(fabs(diagM[i]) > FLOAT_NUMBERS_EQUALITY_TOLERANCE,
                "There is a small (" + d2s(diagM[i]) + ") number (row "
                + d2s(i) + ") on the mass matrix diagonal");

  const string method_name = "SEM_";

  cout << "Open seismograms files..." << flush;
  ofstream *seisU; // for displacement
  ofstream *seisV; // for velocity
  open_seismo_outs(seisU, seisV, param, method_name);
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  GridFunction u_0(&fespace); // displacement
  GridFunction u_1(&fespace);
  GridFunction u_2(&fespace);
  GridFunction v_1(&fespace); // velocity
  u_0 = 0.0;
  u_1 = 0.0;
  u_2 = 0.0;

  const int n_time_steps = ceil(param.T / param.dt);
  const int tenth = 0.1 * n_time_steps;

  const string snapshot_filebase = method_name + param.extra_string;
  const int N = u_0.Size();

  cout << "N time steps = " << n_time_steps
       << "\nTime loop..." << endl;

  // the values of the time-dependent part of the source
  vector<double> time_values(n_time_steps);
  if (!strcmp(param.source.type, "pointforce")) {
    for (int time_step = 1; time_step <= n_time_steps; ++time_step) {
      const double cur_time = time_step * param.dt;
      time_values[time_step-1] = RickerWavelet(param.source,
                                               cur_time - param.dt);
    }
  } else if (!strcmp(param.source.type, "momenttensor")) {
    for (int time_step = 1; time_step <= n_time_steps; ++time_step) {
      const double cur_time = time_step * param.dt;
      time_values[time_step-1] = GaussFirstDerivative(param.source,
                                                      cur_time - param.dt);
    }
  } else MFEM_ABORT("Unknown source type: " + string(param.source.type));

  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    Vector y = u_1; y *= 2.0; y -= u_2;        // y = 2*u_1 - u_2

    Vector z0; z0.SetSize(N);                  // z0 = M * (2*u_1 - u_2)
    for (int i = 0; i < N; ++i) z0[i] = diagM[i] * y[i];

    Vector z1; z1.SetSize(N); S.Mult(u_1, z1);     // z1 = S * u_1
    Vector z2 = b; z2 *= time_values[time_step-1]; // z2 = timeval*source

    // y = dt^2 * (S*u_1 - timeval*source), where it can be
    // y = dt^2 * (S*u_1 - ricker*pointforce) OR
    // y = dt^2 * (S*u_1 - gaussfirstderivative*momenttensor)
    y = z1; y -= z2; y *= param.dt*param.dt;

    // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-timeval*source)
    Vector RHS = z0; RHS -= y;

    for (int i = 0; i < N; ++i) y[i] = diagD[i] * u_2[i]; // y = D * u_2

    // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-timeval*source) + D*u_2
    RHS += y;

    // (M+D)*x_0 = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b) + D*x_2
    for (int i = 0; i < N; ++i) u_0[i] = RHS[i] / (diagM[i]+diagD[i]);

    // velocity: v = du/dt, we use the central difference here
    v_1  = u_0;
    v_1 -= u_2;
    v_1 /= 2.0*param.dt;

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << u_0.Norml2() << endl;

    if (time_step % param.step_snap == 0)
      output_snapshots(time_step, snapshot_filebase, param, u_0, v_1);

    if (time_step % param.step_seis == 0)
      output_seismograms(param, mesh, u_0, v_1, seisU, seisV);

    u_2 = u_1;
    u_1 = u_0;
  }

  delete[] seisU;
  delete[] seisV;

  cout << "Time loop is over" << endl;

  delete fec;
}



double mass_damp_weight(const Vector& point, const Parameters& param)
{
  const double x = point(0);
  const double y = point(1);
  const double z = point(2);
  const bool left   = (!strcmp(param.bc.left,   "abs") ? true : false);
  const bool right  = (!strcmp(param.bc.right,  "abs") ? true : false);
  const bool bottom = (!strcmp(param.bc.bottom, "abs") ? true : false);
  const bool top    = (!strcmp(param.bc.top,    "abs") ? true : false);
  const bool front  = (!strcmp(param.bc.front,  "abs") ? true : false);
  const bool back   = (!strcmp(param.bc.back,   "abs") ? true : false);

  const double X0 = 0.0;
  const double Y0 = 0.0;
  const double Z0 = 0.0;
  const double X1 = param.grid.sx;
  const double Y1 = param.grid.sy;
  const double Z1 = param.grid.sz;
  const double layer = param.bc.damp_layer;
  const double power = param.bc.damp_power;

  // coef for the mass matrix in a damping region is computed
  // C_M = C_Mmax * x^p, where
  // p is typically 3,
  // x changes from 0 at the interface between damping and non-damping regions
  // to 1 at the boundary - the farthest damping layer
  // C_M in the non-damping region is 0

  double weight = 0.0;
  if (left && x - layer <= X0)
    weight += pow((X0-x+layer)/layer, power);
  else if (right && x + layer >= X1)
    weight += pow((x+layer-X1)/layer, power);

  if (bottom && y - layer <= Y0)
    weight += pow((Y0-y+layer)/layer, power);
  else if (top && y + layer >= Y1)
    weight += pow((y+layer-Y1)/layer, power);

  if (front && z - layer <= Z0)
    weight += pow((Z0-z+layer)/layer, power);
  else if (back && z + layer >= Z1)
    weight += pow((z+layer-Z1)/layer, power);

  return weight;
}



double stif_damp_weight(const Vector& point, const Parameters& param)
{
  const double x = point(0);
  const double y = point(1);
  const double z = point(2);
  const bool left   = (!strcmp(param.bc.left,   "abs") ? true : false);
  const bool right  = (!strcmp(param.bc.right,  "abs") ? true : false);
  const bool bottom = (!strcmp(param.bc.bottom, "abs") ? true : false);
  const bool top    = (!strcmp(param.bc.top,    "abs") ? true : false);
  const bool front  = (!strcmp(param.bc.front,  "abs") ? true : false);
  const bool back   = (!strcmp(param.bc.back,   "abs") ? true : false);

  const double X0 = 0.0;
  const double Y0 = 0.0;
  const double Z0 = 0.0;
  const double X1 = param.grid.sx;
  const double Y1 = param.grid.sy;
  const double Z1 = param.grid.sz;
  const double layer = param.bc.damp_layer;
  const double power = param.bc.damp_power+1;
  const double C0 = log(100.0);

  // coef for the stif matrix in a damping region is computed
  // C_K = exp(-C0*alpha(x)*k_inc*x), where
  // C0 = ln(100)
  // alpha(x) = a_Max * x^p
  // p is typically 3,
  // x changes from 0 to 1 (1 at the boundary - the farthest damping layer)
  // C_K in the non-damping region is 1

  double weight = 1.0;
  if (left && x - layer <= X0)
    weight *= exp(-C0*pow((X0-x+layer)/layer, power));
  else if (right && x + layer >= X1)
    weight *= exp(-C0*pow((x+layer-X1)/layer, power));

  if (bottom && y - layer <= Y0)
    weight *= exp(-C0*pow((Y0-y+layer)/layer, power));
  else if (top && y + layer >= Y1)
    weight *= exp(-C0*pow((y+layer-Y1)/layer, power));

  if (front && z - layer <= Z0)
    weight *= exp(-C0*pow((Z0-z+layer)/layer, power));
  else if (back && z + layer >= Z1)
    weight *= exp(-C0*pow((z+layer-Z1)/layer, power));

  return weight;
}



void show_SRM_damp_weights(const Parameters& param)
{
  const int nx = param.grid.nx;
  const int ny = param.grid.ny;
  const int nz = param.grid.nz;

  Vector mass_damp((nx+1)*(ny+1)*(nz+1));
  Vector stif_damp((nx+1)*(ny+1)*(nz+1));

  const double hx = param.grid.sx / nx;
  const double hy = param.grid.sy / ny;
  const double hz = param.grid.sz / nz;

  for (int iz = 0; iz < nz+1; ++iz)
  {
    const double z = (iz == nz ? param.grid.sz : iz*hz);
    for (int iy = 0; iy < ny+1; ++iy)
    {
      const double y = (iy == ny ? param.grid.sy : iy*hy);
      for (int ix = 0; ix < nx+1; ++ix)
      {
        const double x = (ix == nx ? param.grid.sx : ix*hx);
        Vector point(SPACE_DIM);
        point(0) = x;
        point(1) = y;
        point(2) = z;
        const double md = mass_damp_weight(point, param);
        const double sd = stif_damp_weight(point, param);
        const int index = iz*(nx+1)*(ny+1) + iy*(nx+1) + ix;
        mass_damp(index) = md;
        stif_damp(index) = sd;
      }
    }
  }

  string fname = "mass_damping_weights.vts";
  write_vts_scalar(fname, "mass_weights", param.grid.sx, param.grid.sy,
                   param.grid.sz, nx, ny, nz, mass_damp);

  fname = "stif_damping_weights.vts";
  write_vts_scalar(fname, "stif_weights", param.grid.sx, param.grid.sy,
                   param.grid.sz, nx, ny, nz, stif_damp);
}

