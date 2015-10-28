#include "elastic_wave3D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"

#include <fstream>

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
  bool generate_edges = 1;
  Mesh mesh(param.nx, param.ny, param.nz, Element::HEXAHEDRON, generate_edges,
            param.sx, param.sy, param.sz);
  const int dim = mesh.Dimension();
  MFEM_VERIFY(dim == SPACE_DIM, "Unexpected mesh dimension");
  const int n_elements = param.nx*param.ny*param.nz;
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
    const double rho = param.rho_array[i];
    const double vp  = param.vp_array[i];
    const double vs  = param.vs_array[i];

    MFEM_VERIFY(rho > 1.0 && vp > 1.0 && vs > 1.0, "Incorrect media properties "
                "arrays");

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;
  }

  CWConstCoefficient rho_coef(param.rho_array, 0);
  CWFunctionCoefficient lambda_coef  (stif_damp_weight, param, lambda_array);
  CWFunctionCoefficient mu_coef      (stif_damp_weight, param, mu_array);
  CWFunctionCoefficient rho_damp_coef(mass_damp_weight, param, param.rho_array, 0);

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

  VectorPointForce vector_point_force(dim, param.source);
  VectorDomainLFIntegrator *point_force_int = new VectorDomainLFIntegrator(vector_point_force);
  point_force_int->SetIntRule(&hex_GLL);

  MomentTensorSource momemt_tensor_source(dim, param.source);
  VectorDomainLFIntegrator *moment_tensor_int = new VectorDomainLFIntegrator(momemt_tensor_source);
  moment_tensor_int->SetIntRule(&hex_GLL);

  cout << "RHS vector..." << flush;
  LinearForm b(&fespace);
  b.AddDomainIntegrator(point_force_int);
  b.AddDomainIntegrator(moment_tensor_int);
  b.Assemble();
  cout << "||b||_L2 = " << b.Norml2() << endl;
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  Vector diagM; M.GetDiag(diagM); // mass matrix is diagonal
  Vector diagD; D.GetDiag(diagD); // damping matrix is diagonal

  for (int i = 0; i < diagM.Size(); ++i)
    MFEM_VERIFY(fabs(diagM[i]) > FLOAT_NUMBERS_EQUALITY_TOLERANCE,
                string("There is a small (") + d2s(diagM[i]) + ") number (row "
                + d2s(i) + ") on the mass matrix diagonal");

  const string method_name = "SEM_";

  cout << "Open seismograms files..." << flush;
  int n_rec_sets = param.sets_of_receivers.size();
  ofstream *seisU = new ofstream[N_ELAST_COMPONENTS*n_rec_sets]; // for displacement
  ofstream *seisV = new ofstream[N_ELAST_COMPONENTS*n_rec_sets]; // for velocity
  for (int r = 0; r < n_rec_sets; ++r)
  {
    const string desc = param.sets_of_receivers[r]->description();
#if defined(DEBUG_WAVE)
    cout << desc << "\n";
    param.sets_of_receivers[r]->print_receivers(mesh);
#endif
    for (int c = 0; c < N_ELAST_COMPONENTS; ++c)
    {
      string seismofile = method_name + param.extra_string + desc + "_u" + d2s(c) + ".bin";
      seisU[r*N_ELAST_COMPONENTS + c].open(seismofile.c_str(), ios::binary);
      MFEM_VERIFY(seisU[r*N_ELAST_COMPONENTS + c], "File '" + seismofile +
                  "' can't be opened");

      seismofile = method_name + param.extra_string + desc + "_v" + d2s(c) + ".bin";
      seisV[r*N_ELAST_COMPONENTS + c].open(seismofile.c_str(), ios::binary);
      MFEM_VERIFY(seisV[r*N_ELAST_COMPONENTS + c], "File '" + seismofile +
                  "' can't be opened");
    } // loop for components
  } // loop for sets of receivers
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

  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    const double cur_time = time_step * param.dt;

    const double r = param.source.Ricker(cur_time - param.dt);

    Vector y = u_1; y *= 2.0; y -= u_2;        // y = 2*u_1 - u_2

    Vector z0; z0.SetSize(N);                  // z0 = M * (2*u_1 - u_2)
    for (int i = 0; i < N; ++i) z0[i] = diagM[i] * y[i];

    Vector z1; z1.SetSize(N); S.Mult(u_1, z1); // z1 = S * u_1

    Vector z2 = b; z2 *= r;                    // z2 = r * b

    y = z1; y -= z2; y *= param.dt*param.dt;   // y = dt^2 * (S*u_1 - r*b)

    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b)

    for (int i = 0; i < N; ++i) y[i] = diagD[i] * u_2[i]; // y = D * u_2
    RHS += y;                                             // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
    // (M+D)*x_0 = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b) + D*x_2
    for (int i = 0; i < N; ++i) u_0[i] = RHS[i] / (diagM[i]+diagD[i]);

    // velocity
    v_1  = u_0;
    v_1 -= u_2;
    v_1 /= 2.0*param.dt;

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << u_0.Norml2() << endl;

    if (time_step % param.step_snap == 0)
      output_snapshots(time_step, snapshot_filebase, param, u_0, v_1);

    output_seismograms(param, mesh, u_0, v_1, seisU, seisV);

    u_2 = u_1;
    u_1 = u_0;
  }

  delete[] seisV;
  delete[] seisU;

  cout << "Time loop is over" << endl;

  delete fec;
}



double mass_damp_weight(const Vector& point, const Parameters& param)
{
  const double x = point(0);
  const double y = point(1);
  const double z = point(2);
  bool left = true, right = true, bottom = true, front = true, back = true;
  bool top = (param.topsurf == 0 ? true : false);

  const double X0 = 0.0;
  const double X1 = param.sx;
  const double Y0 = 0.0;
  const double Y1 = param.sy;
  const double Z0 = 0.0;
  const double Z1 = param.sz;
  const double layer = param.damp_layer;
  const double power = param.damp_power;

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
  bool left = true, right = true, bottom = true, front = true, back = true;
  bool top = (param.topsurf == 0 ? true : false);

  const double X0 = 0.0;
  const double X1 = param.sx;
  const double Y0 = 0.0;
  const double Y1 = param.sy;
  const double Z0 = 0.0;
  const double Z1 = param.sz;
  const double layer = param.damp_layer;
  const double power = param.damp_power; //+1;
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
  Vector mass_damp((param.nx+1)*(param.ny+1)*(param.nz+1));
  Vector stif_damp((param.nx+1)*(param.ny+1)*(param.nz+1));

  const double hx = param.sx / param.nx;
  const double hy = param.sy / param.ny;
  const double hz = param.sz / param.nz;

  for (int iz = 0; iz < param.nz+1; ++iz)
  {
    const double z = (iz == param.nz ? param.sz : iz*hz);
    for (int iy = 0; iy < param.ny+1; ++iy)
    {
      const double y = (iy == param.ny ? param.sy : iy*hy);
      for (int ix = 0; ix < param.nx+1; ++ix)
      {
        const double x = (ix == param.nx ? param.sx : ix*hx);
        Vector point(SPACE_DIM);
        point(0) = x;
        point(1) = y;
        point(2) = z;
        const double md = mass_damp_weight(point, param);
        const double sd = stif_damp_weight(point, param);
        const int index = iz*(param.nx+1)*(param.ny+1) + iy*(param.nx+1) + ix;
        mass_damp(index) = md;
        stif_damp(index) = sd;
      }
    }
  }

  string fname = "mass_damping_weights.vts";
  write_vts_scalar(fname, "mass_weights", param.sx, param.sy, param.sz,
                   param.nx, param.ny, param.nz, mass_damp);

  fname = "stif_damping_weights.vts";
  write_vts_scalar(fname, "stif_weights", param.sx, param.sy, param.sz,
                   param.nx, param.ny, param.nz, stif_damp);
}

