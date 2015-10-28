#include "source.hpp"

using namespace std;
using namespace mfem;



Source::Source()
  : location(500.0, 500.0, 500.0)
  , frequency(10.0)
  , direction(2) // OY
  , ricker_scale(1.0)
  , point_force_scale(1.0)
  , gauss_support(1.0)
  , Mxx(1.0), Mxy(0.0), Mxz(0.0), Myy(1.0), Myz(0.0), Mzz(1.0) // explosive source
  , type(1)
{ }



void Source::AddOptions(OptionsParser& args)
{
  args.AddOption(&location(0), "-srcx", "--source-x", "x-coord of a source location");
  args.AddOption(&location(1), "-srcy", "--source-y", "y-coord of a source location");
  args.AddOption(&location(2), "-srcz", "--source-z", "z-coord of a source location");
  args.AddOption(&frequency, "-f0", "--frequency", "Central frequency of a source");
  args.AddOption(&direction, "-dir", "--direction", "Direction of point force component (1 OX, 2 OY, 3 OZ, 4 with auxpoint)");
  args.AddOption(&ricker_scale, "-rs", "--ricker-scale", "Factor for the Ricker wavelet");
  args.AddOption(&point_force_scale, "-pfs", "--point-force-scale", "Factor for the point force term of a source");
  args.AddOption(&gauss_support, "-gs", "--gauss-support", "Gauss support for Gaussian space functions of a source");
  args.AddOption(&Mxx, "-mxx", "--moment-tensor-xx", "xx-component of a moment tensor source");
  args.AddOption(&Mxy, "-mxy", "--moment-tensor-xy", "xy-component of a moment tensor source");
  args.AddOption(&Mxz, "-mxz", "--moment-tensor-xz", "xz-component of a moment tensor source");
  args.AddOption(&Myy, "-myy", "--moment-tensor-yy", "yy-component of a moment tensor source");
  args.AddOption(&Myz, "-myz", "--moment-tensor-yz", "yz-component of a moment tensor source");
  args.AddOption(&Mzz, "-mzz", "--moment-tensor-zz", "zz-component of a moment tensor source");
  args.AddOption(&type, "-st", "--source-type", "Type of spatial source distribution (0 delta, 1 gauss, 2 plane)");
}



double Source::Ricker(double t) const
{
  const double a  = M_PI*frequency*(t-1./frequency);
  return ricker_scale * (1.-2.*a*a)*exp(-a*a);
}



double Source::GaussFirstDerivative(double t) const
{
  const double a = M_PI*frequency*(t-1./frequency);
  return ricker_scale * (t-1./frequency)*exp(-a*a);
}



void Source::PointForce(const Vector& x, Vector& f) const
{
  if (type == 0) // Delta
    DeltaPointForce(x, f);
  else if (type == 1) // Gauss
    GaussPointForce(x, f);
  else
    mfem_error("Unknown source type");
}



void Source::MomentTensorSource(const Vector &x, Vector &f) const
{
  if (type == 0) // Delta
    DivDeltaMomentTensor(x, f);
  else if (type == 1) // Gauss
    DivGaussMomentTensor(x, f);
  else
    mfem_error("Unknown source type");
}



void Source::DeltaPointForce(const Vector& x, Vector& f) const
{
  const double tol = 1e-2;
  const double loc[] = { location(0), location(1), location(2) };
  double scale = 0.0;
  if (x.DistanceTo(loc) < tol)
    scale = point_force_scale;

  Vector direction_vector(3);
  direction_vector = 0.0;
  direction_vector(direction-1) = 1.0;

  f(0) = scale * direction_vector(0);
  f(1) = scale * direction_vector(1);
  f(2) = scale * direction_vector(2);
}



void Source::GaussPointForce(const Vector& x, Vector& f) const
{
  const double xdiff  = x(0)-location(0);
  const double ydiff  = x(1)-location(1);
  const double zdiff  = x(2)-location(2);
  const double xdiff2 = xdiff*xdiff;
  const double ydiff2 = ydiff*ydiff;
  const double zdiff2 = zdiff*zdiff;
  const double h2 = gauss_support*gauss_support;
  const double G = exp(-(xdiff2 + ydiff2 + zdiff2) / h2);

  Vector direction_vector(3);
  direction_vector = 0.0;
  direction_vector(direction-1) = 1.0;

  f(0) = point_force_scale * G * direction_vector(0);
  f(1) = point_force_scale * G * direction_vector(1);
  f(2) = point_force_scale * G * direction_vector(2);
}



void Source::DivDeltaMomentTensor(const Vector& x, Vector& f) const
{
  mfem_error("NOT implemented");
}



void Source::DivGaussMomentTensor(const Vector& x, Vector& f) const
{
  const double xdiff  = x(0)-location(0);
  const double ydiff  = x(1)-location(1);
  const double zdiff  = x(2)-location(2);
  const double xdiff2 = xdiff*xdiff;
  const double ydiff2 = ydiff*ydiff;
  const double zdiff2 = zdiff*zdiff;
  const double h2 = gauss_support*gauss_support;
  const double exp_val = exp(-(xdiff2 + ydiff2 + zdiff2) / h2);
  const double Gx = -2.*xdiff/h2 * exp_val;
  const double Gy = -2.*ydiff/h2 * exp_val;
  const double Gz = -2.*zdiff/h2 * exp_val;

  f(0) = Mxx*Gx + Mxy*Gy + Mxz*Gz;
  f(1) = Mxy*Gx + Myy*Gy + Myz*Gz;
  f(2) = Mxz*Gx + Myz*Gy + Mzz*Gz;
}

