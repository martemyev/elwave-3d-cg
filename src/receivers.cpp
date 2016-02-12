#include "mfem.hpp"
#include "receivers.hpp"

using namespace mfem;

//==============================================================================
//
// ReceiversSet
//
//==============================================================================
ReceiversSet::ReceiversSet()
  : _variable(""),
    _n_receivers(0),
    _receivers(),
    _cells_containing_receivers()
{ }

void ReceiversSet::
find_cells_containing_receivers(int nx, int ny, int nz, double sx, double sy,
                                double sz)
{
  MFEM_VERIFY(!_receivers.empty(), "The receivers hasn't been distributed yet");

  _cells_containing_receivers.clear();
  _cells_containing_receivers.resize(_n_receivers);

  // we throw an exception if we don't find a cell containing a receiver
  const bool throw_exception = true;

  for (int p = 0; p < _n_receivers; ++p)
    _cells_containing_receivers[p] = find_element(sx, sy, sz, nx, ny, nz,
                                                  _receivers[p],
                                                  throw_exception);
}




//==============================================================================
//
// ReceiversLine
//
//==============================================================================
ReceiversLine::ReceiversLine()
  : _start(),
    _end()
{ }

void ReceiversLine::init(std::ifstream &in)
{
  MFEM_VERIFY(in.is_open(), "The stream for reading receivers is not open");

  std::string tmp;

  in >> _variable;
  in >> _n_receivers; getline(in, tmp);
  in >> _start(0) >> _start(1) >> _start(2); getline(in, tmp);
  in >> _end(0)   >> _end(1)   >> _end(2);   getline(in, tmp);

  MFEM_VERIFY(_n_receivers > 0, "The number of receivers (" + d2s(_n_receivers)+
              ") must be >0");
}

void ReceiversLine::distribute_receivers()
{
  const double x0 = _start(0);
  const double y0 = _start(1);
  const double z0 = _start(2);
  const double x1 = _end(0);
  const double y1 = _end(1);
  const double z1 = _end(2);

  _receivers.resize(_n_receivers);

  const double dx = (x1 - x0) / (_n_receivers-1);
  const double dy = (y1 - y0) / (_n_receivers-1);
  const double dz = (z1 - z0) / (_n_receivers-1);

  for (int i = 0; i < _n_receivers; ++i)
  {
    const double x = (i == _n_receivers-1 ? x1 : x0 + i*dx);
    const double y = (i == _n_receivers-1 ? y1 : y0 + i*dy);
    const double z = (i == _n_receivers-1 ? z1 : z0 + i*dz);
    _receivers[i] = Vertex(x, y, z);
  }
}

std::string ReceiversLine::description() const
{
  return "_rec_line_x" + d2s(_start(0)) + "_" + d2s(_end(0)) +
         "_y" + d2s(_start(1)) + "_" + d2s(_end(1)) +
         "_z" + d2s(_start(2)) + "_" + d2s(_end(2));
}




//==============================================================================
//
// Auxiliary
//
//==============================================================================
int find_element(double sx, double sy, double sz, int nx, int ny, int nz,
                 const Vertex &point, bool throw_exception)
{
  const double px = point(0); // coordinates of the point of interest
  const double py = point(1);
  const double pz = point(2);

  const double x0 = 0.0; // limits of the rectangular mesh
  const double x1 = sx;
  const double y0 = 0.0;
  const double y1 = sy;
  const double z0 = 0.0;
  const double z1 = sz;

  // check that the point is within the mesh
  const double tol = FIND_CELL_TOLERANCE;
  if (px < x0 - tol || px > x1 + tol ||
      py < y0 - tol || py > y1 + tol)
  {
    if (throw_exception)
      MFEM_ABORT("The given point [" + d2s(px) + "," + d2s(py) + "] doesn't "
                 "belong to the rectangular mesh");

    return -1; // to show that the point in not here
  }

  // since the elements of the cubic mesh are numerated in the following
  // way:
  // ^ Y          / Z
  // | -----------
  // |/ 6  / 7  /|
  // ----------- |
  // | 2  | 3  |/|
  // ----------- |
  // | 0  | 1  |/
  // -------------> X
  // we can simplify the search of the element containing the given point:

  int cellx = (px-x0) * nx / (x1 - x0);
  int celly = (py-y0) * ny / (y1 - y0);
  int cellz = (pz-z0) * nz / (z1 - z0);

  if (cellx) --cellx;
  if (celly) --celly;
  if (cellz) --cellz;

  return (cellz*nx*ny + celly*nx + cellx);
}
