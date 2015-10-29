#include "source.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;



Source::Source()
  : location(500.0, 500.0, 500.0)
  , frequency(10.0)
  , direction(2) // OY
  , scale(1.0)
  , Mxx(1.0), Mxy(0.0), Mxz(0.0), Myy(1.0), Myz(0.0), Mzz(1.0) // explosive source
  , type_string("pointforce")
  , spatial_function("gauss")
  , gauss_support(1.0)
  , type(POINT_FORCE)
{ }



void Source::AddOptions(OptionsParser& args)
{
  args.AddOption(&location(0), "-srcx", "--source-x", "x-coord of a source location");
  args.AddOption(&location(1), "-srcy", "--source-y", "y-coord of a source location");
  args.AddOption(&location(2), "-srcz", "--source-z", "z-coord of a source location");
  args.AddOption(&frequency, "-f0", "--source-frequency", "Central frequency of a source");
  args.AddOption(&direction, "-dir", "--source-direction", "Direction of point force or plane wave source (1 OX, 2 OY, 3 OZ, 4 with auxpoint)");
  args.AddOption(&scale, "-scale", "--source-scale", "Scaling factor for the source function");
  args.AddOption(&Mxx, "-mxx", "--moment-tensor-xx", "xx-component of a moment tensor source");
  args.AddOption(&Mxy, "-mxy", "--moment-tensor-xy", "xy-component of a moment tensor source");
  args.AddOption(&Mxz, "-mxz", "--moment-tensor-xz", "xz-component of a moment tensor source");
  args.AddOption(&Myy, "-myy", "--moment-tensor-yy", "yy-component of a moment tensor source");
  args.AddOption(&Myz, "-myz", "--moment-tensor-yz", "yz-component of a moment tensor source");
  args.AddOption(&Mzz, "-mzz", "--moment-tensor-zz", "zz-component of a moment tensor source");
  args.AddOption(&type_string, "-type", "--source-type", "Type of the source (pointforce, momenttensor, planewave)");
  args.AddOption(&spatial_function, "-spatial", "--source-spatial", "Spatial function of the source (delta, gauss)");
  args.AddOption(&gauss_support, "-gs", "--gauss-support", "Gauss support for 'gauss' spatial function of the source");
}



void Source::check_and_update_parameters()
{
  MFEM_VERIFY(frequency > 0, "Frequency (" + d2s(frequency) + ") must be >0");
  MFEM_VERIFY(direction == 1 || direction == 2 || direction == 3, "Unsupported "
              "direction of the source: " + d2s(direction));

  if (!strcmp(type_string, "pointforce"))
    type = POINT_FORCE;
  else if (!strcmp(type_string, "momenttensor"))
    type = MOMENT_TENSOR;
  else if (!strcmp(type_string, "planewave"))
    type = PLANE_WAVE;
  else
    MFEM_ABORT("Unknown source type: " + string(type_string));
}



double Source::Ricker(double t) const
{
  const double a = M_PI*frequency*(t-1./frequency);
  return scale * (1.-2.*a*a)*exp(-a*a);
}



double Source::GaussFirstDerivative(double t) const
{
  const double a = M_PI*frequency*(t-1./frequency);
  return scale * (t-1./frequency)*exp(-a*a);
}



void Source::PointForce(const Vector& x, Vector& f) const
{
  if (!strcmp(spatial_function, "delta"))
    DeltaPointForce(x, f);
  else if (!strcmp(spatial_function, "gauss"))
    GaussPointForce(x, f);
  else
    MFEM_ABORT("Unknown source spatial function: " + string(spatial_function));
}



void Source::MomentTensorSource(const Vector &x, Vector &f) const
{
  if (!strcmp(spatial_function, "delta"))
    DivDeltaMomentTensor(x, f);
  else if (!strcmp(spatial_function, "gauss"))
    DivGaussMomentTensor(x, f);
  else
    MFEM_ABORT("Unknown source spatial function: " + string(spatial_function));
}



void Source::PlaneWaveSource(const Vector &x, Vector &f) const
{
  const double py = x(1); // y-coordinate of the point

  f = 0.0;
  if (fabs(py - location(1)) < FLOAT_NUMBERS_EQUALITY_TOLERANCE)
    f(1) = 1.0;
}



void Source::DeltaPointForce(const Vector& x, Vector& f) const
{
  const double tol = 1e-2;
  const double loc[] = { location(0), location(1), location(2) };
  double value = 0.0;
  if (x.DistanceTo(loc) < tol)
    value = 1.0;

  f = 0.0;
  f(direction-1) = value;
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

  f = 0.0;
  f(direction-1) = G;
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



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
VectorPointForce::VectorPointForce(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{ }



void VectorPointForce::Eval(Vector &V, ElementTransformation &T,
                            const IntegrationPoint &ip)
{
  Vector transip;
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
  Vector transip;
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.MomentTensorSource(transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
PlaneWaveSource::PlaneWaveSource(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
  ,
{ }



void PlaneWaveSource::Eval(Vector &V, ElementTransformation &T,
                           const IntegrationPoint &ip)
{
  Vector transip;
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.PlaneWaveSource(transip, V);
}
