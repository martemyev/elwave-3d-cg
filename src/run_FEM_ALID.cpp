#include "elastic_wave3D.hpp"
#include "parameters.hpp"
#include "receivers.hpp"

#include <float.h>
#include <fstream>
#include <vector>

using namespace std;
using namespace mfem;

void compute_damping_weights(const Parameters& param, double *damping_weights);
void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2);



void ElasticWave2D::run_FEM_ALID()
{
  StopWatch chrono;

  chrono.Start();
  cout << "FE space generation..." << flush;

  const int dim = param.mesh->Dimension();
  const int n_elements = param.mesh->GetNE();

  FiniteElementCollection *fec = new H1_FECollection(param.order, dim);
  FiniteElementSpace fespace(param.mesh, fec, dim); //, Ordering::byVDIM);
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "Number of unknowns: " << fespace.GetVSize() << endl;

  double *lambda_array = new double[n_elements];
  double *mu_array     = new double[n_elements];

  double Rho[] = { DBL_MAX, DBL_MIN };
  double Vp[]  = { DBL_MAX, DBL_MIN };
  double Vs[]  = { DBL_MAX, DBL_MIN };
  double Lam[] = { DBL_MAX, DBL_MIN };
  double Mu[]  = { DBL_MAX, DBL_MIN };

  double *damp_weights = new double[n_elements];
  compute_damping_weights(param, damp_weights);

  double *rho_w_array   = new double[n_elements];
  double *lambda_w_array= new double[n_elements];
  double *mu_w_array    = new double[n_elements];

  for (int i = 0; i < n_elements; ++i)
  {
    const double rho = param.media.rho_array[i];
    const double vp  = param.media.vp_array[i];
    const double vs  = param.media.vs_array[i];
    const double w   = damp_weights[i];

    MFEM_VERIFY(rho > 1.0 && vp > 1.0 && vs > 1.0, "Incorrect media properties "
                "arrays");

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;

    rho_w_array[i]   = rho*w;
    lambda_w_array[i]= lambda_array[i]*w;
    mu_w_array[i]    = mu_array[i]*w;

    Rho[0] = std::min(Rho[0], rho);
    Rho[1] = std::max(Rho[1], rho);
    Vp[0]  = std::min(Vp[0], vp);
    Vp[1]  = std::max(Vp[1], vp);
    Vs[0]  = std::min(Vs[0], vs);
    Vs[1]  = std::max(Vs[1], vs);
    Lam[0] = std::min(Lam[0], lambda_array[i]);
    Lam[1] = std::max(Lam[1], lambda_array[i]);
    Mu[0]  = std::min(Mu[0], mu_array[i]);
    Mu[1]  = std::max(Mu[1], mu_array[i]);
  }

  std::cout << "Rho: min " << Rho[0] << " max " << Rho[1] << "\n";
  std::cout << "Vp: min "  << Vp[0]  << " max " << Vp[1] << "\n";
  std::cout << "Vs: min "  << Vs[0]  << " max " << Vs[1] << "\n";
  std::cout << "Lam: min " << Lam[0] << " max " << Lam[1] << "\n";
  std::cout << "Mu: min "  << Mu[0]  << " max " << Mu[1] << "\n";

  const bool own_array = false;
  CWConstCoefficient rho_coef(param.media.rho_array, own_array);
  CWConstCoefficient lambda_coef(lambda_array);
  CWConstCoefficient mu_coef(mu_array);

  CWConstCoefficient rho_w_coef(rho_w_array);
  CWConstCoefficient lambda_w_coef(lambda_w_array);
  CWConstCoefficient mu_w_coef(mu_w_array);

  cout << "Stif matrix..." << flush;
  BilinearForm stif(&fespace);
  stif.AddDomainIntegrator(new ElasticityIntegrator(lambda_coef, mu_coef));
  stif.Assemble();
  stif.Finalize();
  const SparseMatrix& S = stif.SpMat();
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "Mass matrix..." << flush;
  BilinearForm mass(&fespace);
  mass.AddDomainIntegrator(new VectorMassIntegrator(rho_coef));
  mass.Assemble();
  mass.Finalize();
  SparseMatrix& M = mass.SpMat();
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

  double alpha1, alpha2;
  get_damp_alphas(param.source.frequency, alpha1, alpha2);

  cout << "Damp matrix..." << flush;
  BilinearForm dampS(&fespace);
  dampS.AddDomainIntegrator(new ElasticityIntegrator(lambda_w_coef, mu_w_coef));
  dampS.Assemble();
  dampS.Finalize();
  SparseMatrix& D = dampS.SpMat();
  D *= alpha2;

  BilinearForm dampM(&fespace);
  dampM.AddDomainIntegrator(new VectorMassIntegrator(rho_w_coef));
  dampM.Assemble();
  dampM.Finalize();
  SparseMatrix& DM = dampM.SpMat();
  DM *= alpha1;

  D += DM;
  D *= 0.5*param.dt;
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "Sys matrix..." << flush;
  const SparseMatrix& CopyFrom = D;
  const int nnz = CopyFrom.NumNonZeroElems();
  const bool ownij  = false;
  const bool ownval = true;
  SparseMatrix Sys(CopyFrom.GetI(), CopyFrom.GetJ(), new double[nnz],
                   CopyFrom.Height(), CopyFrom.Width(), ownij, ownval,
                   CopyFrom.areColumnsSorted());
  Sys = 0.0;
  Sys += D;
  Sys += M;
  GSSmoother prec(Sys);
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  cout << "RHS vector... " << flush;
  LinearForm b(&fespace);
  if (param.source.plane_wave)
  {
    PlaneWaveSource plane_wave_source(dim, param);
    b.AddDomainIntegrator(new VectorDomainLFIntegrator(plane_wave_source));
    b.Assemble();
  }
  else
  {
    if (!strcmp(param.source.type, "pointforce"))
    {
      VectorPointForce vector_point_force(dim, param);
      b.AddDomainIntegrator(new VectorDomainLFIntegrator(vector_point_force));
      b.Assemble();
    }
    else if (!strcmp(param.source.type, "momenttensor"))
    {
      MomentTensorSource momemt_tensor_source(dim, param);
      b.AddDomainIntegrator(new VectorDomainLFIntegrator(momemt_tensor_source));
      b.Assemble();
    }
    else MFEM_ABORT("Unknown source type: " + string(param.source.type));
  }
  cout << "||b||_L2 = " << b.Norml2() << endl;
  cout << "done. Time = " << chrono.RealTime() << " sec" << endl;
  chrono.Clear();

  const string method_name = "FEM_";

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

  const string snapshot_filebase = (string)param.output_dir + "/" +
                                   SNAPSHOTS_DIR + method_name +
                                   param.extra_string;
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
    Vector y = u_1; y *= 2.0; y -= u_2;            // y = 2*u_1 - u_2
    Vector z0; z0.SetSize(N); M.Mult(y, z0);       // z0 = M * (2*u_1 - u_2)
    Vector z1; z1.SetSize(N); S.Mult(u_1, z1);     // z1 = S * u_1
    Vector z2 = b; z2 *= time_values[time_step-1]; // z2 = timeval*source

    // y = dt^2 * (S*u_1 - timeval*source), where it can be
    // y = dt^2 * (S*u_1 - ricker*pointforce) OR
    // y = dt^2 * (S*u_1 - gaussfirstderivative*momenttensor)
    y = z1; y -= z2; y *= param.dt*param.dt;

    // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-timeval*source)
    Vector RHS = z0; RHS -= y;

    D.Mult(u_2, y); // y = D * u_2

    // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-timeval*source) + D*u_2
    RHS += y;

    // (M+D)*u_0 = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
    PCG(Sys, prec, RHS, u_0, 0, 200, 1e-12, 0.0);

    // velocity: v = du/dt, we use the central difference here
    v_1  = u_0;
    v_1 -= u_2;
    v_1 /= 2.0*param.dt;

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << u_0.Norml2() << endl;

    if (time_step % param.step_snap == 0)
      output_snapshots(time_step, snapshot_filebase, param, u_0, v_1, *param.mesh);

    if (time_step % param.step_seis == 0)
      output_seismograms(param, *param.mesh, u_0, v_1, seisU, seisV);

    u_2 = u_1;
    u_1 = u_0;
  }

  delete[] seisU;
  delete[] seisV;

  cout << "Time loop is over" << endl;

  delete fec;
}



void compute_damping_weights(const Parameters& param, double *damping_weights)
{
  const bool left   = (!strcmp(param.bc.left,   "abs") ? true : false);
  const bool right  = (!strcmp(param.bc.right,  "abs") ? true : false);
  const bool bottom = (!strcmp(param.bc.bottom, "abs") ? true : false);
  const bool top    = (!strcmp(param.bc.top,    "abs") ? true : false);
  const bool front  = (!strcmp(param.bc.front,  "abs") ? true : false);
  const bool back   = (!strcmp(param.bc.back,   "abs") ? true : false);

  const double X0 = 0.;
  const double Y0 = 0.;
  const double Z0 = 0.;
  const double X1 = param.grid.sx;
  const double Y1 = param.grid.sy;
  const double Z1 = param.grid.sz;
  const double layer = param.bc.damp_layer;
  const double power = param.bc.damp_power;

  for (int el = 0; el < param.mesh->GetNE(); ++el)
  {
    const Element *element = param.mesh->GetElement(el);

    Array<int> vert_indices;
    element->GetVertices(vert_indices);

    double x = 0.; // coordinates of the center of the element
    double y = 0.;
    double z = 0.;
    for (int i = 0; i < vert_indices.Size(); ++i)
    {
      const double *v = param.mesh->GetVertex(vert_indices[i]);
      x += v[0];
      y += v[1];
      z += v[2];
    }
    x /= vert_indices.Size();
    y /= vert_indices.Size();
    z /= vert_indices.Size();

    double weight = 1e-12;

    if (left && x - layer < X0)
      weight += pow((X0-x+layer)/layer, power);
    else if (right && x + layer > X1)
      weight += pow((x+layer-X1)/layer, power);

    if (bottom && y - layer < Y0)
      weight += pow((Y0-y+layer)/layer, power);
    else if (top && y + layer > Y1)
      weight += pow((y+layer-Y1)/layer, power);

    if (front && z - layer < Z0)
      weight += pow((Z0-z+layer)/layer, power);
    else if (back && z + layer > Z1)
      weight += pow((z+layer-Z1)/layer, power);

    damping_weights[el] = weight;
  }
}



void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2)
{
  // These numbers have been obtained by Shubin Fu and Kai Gao, while they've
  // been PhD students at Texas A&M.
  const double w1 = 0.7 * source_frequency;
  const double w2 = 20. * source_frequency;
  const double xi1 = 4.5;
  const double xi2 = 0.6;

  alpha1 = 2.*w1*w2*(xi2*w1-xi1*w2)/(w1*w1-w2*w2);
  alpha2 = 2.*(xi1*w1-xi2*w2)/(w1*w1-w2*w2);
}
