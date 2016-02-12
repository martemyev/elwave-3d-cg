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
  MFEM_VERIFY(!_receivers.empty(), "The receivers haven't been distributed yet");

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

