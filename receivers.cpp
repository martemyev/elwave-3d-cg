#include "mfem.hpp"
#include "receivers.hpp"

using namespace mfem;

//==============================================================================
//
// ReceiversSet
//
//==============================================================================
ReceiversSet::ReceiversSet()
  : _n_receivers(0),
    _receivers(),
    _cells_containing_receivers()
{ }

ReceiversSet::ReceiversSet(const ReceiversSet& rec)
  : _n_receivers(rec._n_receivers),
    _receivers(rec._receivers),
    _cells_containing_receivers(rec._cells_containing_receivers)
{ }

ReceiversSet& ReceiversSet::operator =(const ReceiversSet& rec)
{
  _n_receivers = rec._n_receivers;
  _receivers = rec._receivers;
  _cells_containing_receivers = rec._cells_containing_receivers;
  return *this;
}

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

//void ReceiversSet::print_receivers(const Mesh& mesh, std::ostream& out) const
//{
//  out << "n_receivers = " << _n_receivers << "\n# x y ";
//  if (!_cells_containing_receivers.empty())
//    out << "cell\tx0\ty0\tx1\ty1";
//  out << "\n";
//  for (int r = 0; r < _n_receivers; ++r)
//  {
//    out << r+1 << "\t" << _receivers[r](0) << "\t" << _receivers[r](1) << "\t";
//    if (!_cells_containing_receivers.empty())
//    {
//      out << _cells_containing_receivers[r] << "\t";
//      const Element *cell = mesh.GetElement(_cells_containing_receivers[r]);
//      Array<int> vert_indices;
//      cell->GetVertices(vert_indices);
//      out << "vertices:\t";
//      for (int v = 0; v < vert_indices.Size(); ++v)
//      {
//        const double *vert = mesh.GetVertex(vert_indices[v]);
//        out << vert[0] << " " << vert[1] << "\t";
//      }
//    }
//    out << "\n";
//  }
//}



//==============================================================================
//
// ReceiversLine
//
//==============================================================================
ReceiversLine::ReceiversLine()
  : _start(),
    _end()
{ }

ReceiversLine::ReceiversLine(const ReceiversLine& rec)
  : ReceiversSet(rec),
    _start(rec._start),
    _end(rec._end)
{ }

ReceiversLine& ReceiversLine::operator =(const ReceiversLine& rec)
{
  ReceiversSet::operator =(rec);
  _start = rec._start;
  _end   = rec._end;
  return *this;
}

void ReceiversLine::init(std::ifstream &in)
{
  MFEM_VERIFY(in.is_open(), "The stream for reading receivers is not open");

  std::string tmp;

  in >> _n_receivers; getline(in, tmp);
  in >> _start(0); getline(in, tmp);
  in >> _end(0);   getline(in, tmp);
  in >> _start(1); getline(in, tmp);
  in >> _end(1);   getline(in, tmp);
  in >> _start(2); getline(in, tmp);
  in >> _end(2);   getline(in, tmp);

  MFEM_VERIFY(_n_receivers > 0, "The number of receivers (" + d2s(_n_receivers)+
              ") must be >0");
}

void ReceiversLine::distribute_receivers()
{
  const double x0 = _start(0);
  const double x1 = _end(0);
  const double y0 = _start(1);
  const double y1 = _end(1);
  const double z0 = _start(2);
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
// ReceiversPlane
//
//==============================================================================
ReceiversPlane::ReceiversPlane()
  : _n_receivers_1(0),
    _n_receivers_2(0),
    _plane(""),
    _plane_coord(0.0),
    _start(),
    _end()
{ }

ReceiversPlane::ReceiversPlane(const ReceiversPlane& rec)
  : ReceiversSet(rec),
    _n_receivers_1(rec._n_receivers_1),
    _n_receivers_2(rec._n_receivers_2),
    _plane(rec._plane),
    _plane_coord(rec._plane_coord),
    _start(rec._start),
    _end(rec._end)
{ }

ReceiversPlane& ReceiversPlane::operator =(const ReceiversPlane& rec)
{
  ReceiversSet::operator =(rec);
  _n_receivers_1 = rec._n_receivers_1;
  _n_receivers_2 = rec._n_receivers_2;
  _plane = rec._plane;
  _plane_coord = rec._plane_coord;
  _start = rec._start;
  _end   = rec._end;
  return *this;
}

void ReceiversPlane::init(std::ifstream &in)
{
  MFEM_VERIFY(in.is_open(), "The stream for reading receivers is not open");

  std::string tmp;

  in >> _n_receivers_1; getline(in, tmp);
  in >> _n_receivers_2; getline(in, tmp);
  in >> _plane; getline(in, tmp);
  in >> _plane_coord; getline(in, tmp);
  in >> _start(0); getline(in, tmp);
  in >> _end(0);   getline(in, tmp);
  in >> _start(1); getline(in, tmp);
  in >> _end(1);   getline(in, tmp);

  MFEM_VERIFY(_n_receivers_1 > 0 && _n_receivers_2 > 0, "The number of "
              "receivers (" + d2s(_n_receivers_1) + " and " +
              d2s(_n_receivers_2) + ") must be >0");
  MFEM_VERIFY(_plane == "XY" || _plane == "XZ" || _plane == "YZ", "Unknown "
              " plane of receivers (" + _plane + "). Supported: XY, XZ, YZ");

  _n_receivers = _n_receivers_1 * _n_receivers_2;
}

void ReceiversPlane::distribute_receivers()
{
  _receivers.resize(_n_receivers);

  const double c1_beg = _start(0); // coordinate 1
  const double c1_end = _end(0);
  const double c2_beg = _start(1); // coordinate 2
  const double c2_end = _end(1);
  const double dc1 = (c1_end - c1_beg) / (_n_receivers_1-1);
  const double dc2 = (c2_end - c2_beg) / (_n_receivers_2-1);
  for (int i = 0; i < _n_receivers_2; ++i)
  {
    const double c2 = (i == _n_receivers_2-1 ? c2_end : c2_beg + i*dc2);
    for (int j = 0; j < _n_receivers_1; ++j)
    {
      const double c1 = (j == _n_receivers_1-1 ? c1_end : c1_beg + j*dc1);
      if (_plane == "XY")
        _receivers[i*_n_receivers_1+j] = Vertex(c1, c2, _plane_coord);
      else if (_plane == "XZ")
        _receivers[i*_n_receivers_1+j] = Vertex(c1, _plane_coord, c2);
      else if (_plane == "YZ")
        _receivers[i*_n_receivers_1+j] = Vertex(_plane_coord, c1, c2);
      else MFEM_ABORT("Unknown receivers plane");
    }
  }
}

std::string ReceiversPlane::description() const
{
  return "_rec_plane_" + _plane + "_coord" + d2s(_plane_coord);
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

  int cellx = px * nx / (x1 - x0);
  int celly = py * ny / (y1 - y0);
  int cellz = pz * nz / (z1 - z0);

  if (cellx) --cellx;
  if (celly) --celly;
  if (cellz) --cellz;

  return (cellz*nx*ny + celly*nx + cellx);
}
