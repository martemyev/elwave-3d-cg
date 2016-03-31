#include "snapshots.hpp"
#include "utilities.hpp"

#include <cmath>
//#include <fstream>

//#include "mfem.hpp"

using namespace std;
using namespace mfem;

//==============================================================================
//
// SnapshotsSet
//
//==============================================================================
SnapshotsSet::SnapshotsSet()
  : _variable("")
  , _format("bin")
  , _n_snapshot_points(0)
  , _snapshot_points(nullptr)
  , _cells_contain_snapshot_points()
{ }

void SnapshotsSet::
find_cells_containing_snapshot_points(int nx, int ny, int nz, double sx,
                                      double sy, double sz)
{
  MFEM_VERIFY(_snapshot_points != nullptr, "The snapshot points haven't been "
              "distributed yet");

  delete[] _cells_contain_snapshot_points;
  _cells_contain_snapshot_points = new int[_n_snapshot_points];

  // we throw an exception if we don't find a cell containing a receiver
  const bool throw_exception = true;

  for (int p = 0; p < _n_snapshot_points; ++p)
    _cells_contain_snapshot_points[p] = find_element(sx, sy, sz, nx, ny, nz,
                                                     _snapshot_points[p],
                                                     throw_exception);
}

void SnapshotsSet::
save_snapshot_vector(const Vector &snapshot, int time_step,
                     const std::string &var_name, int n_components,
                     const std::string &output_file_base) const
{
  if (_format.find("bin") != std::string::npos)
    save_snapshot_vector_bin(snapshot, time_step, var_name, n_components,
                             output_file_base);
  if (_format.find("vts") != std::string::npos)
    save_snapshot_vector_vts(snapshot, time_step, var_name, n_components,
                             output_file_base);
}




//==============================================================================
//
// SnapshotsPlane
//
//==============================================================================
SnapshotsPlane::SnapshotsPlane()
  : _n_points_1(0), _n_points_2(0), _plane("unknown")
{ }

void SnapshotsPlane::init(std::ifstream &in)
{
  MFEM_VERIFY(in.is_open(), "The stream for reading snapshot plane is not open");

  std::string tmp;

  in >> _variable; getline(in, tmp);
  in >> _format;   getline(in, tmp);
  in >> _n_points_1 >> _n_points_2; getline(in, tmp);
  for (int v = 0; v < N_VERTICES; ++v) {
    in >> _vertices[v](0) >> _vertices[v](1) >> _vertices[v](2);
    getline(in, tmp);
  }

  _n_snapshot_points = _n_points_1*_n_points_2;
  MFEM_VERIFY(_n_snapshot_points > 0, "The number of snapshot points (" +
              d2s(_n_snapshot_points) + ") must be >0");
}

void SnapshotsPlane::distribute_snapshot_points()
{
  _snapshot_points = new mfem::Vertex[_n_snapshot_points];

  double x0 = _vertices[0](0);
  double y0 = _vertices[0](1);
  double z0 = _vertices[0](2);
  double x1 = _vertices[0](0);
  double y1 = _vertices[0](1);
  double z1 = _vertices[0](2);
  for (int v = 1; v < N_VERTICES; ++v)
  {
    x0 = std::min(x0, _vertices[v](0));
    y0 = std::min(y0, _vertices[v](1));
    z0 = std::min(z0, _vertices[v](2));
    x1 = std::max(x1, _vertices[v](0));
    y1 = std::max(y1, _vertices[v](1));
    z1 = std::max(z1, _vertices[v](2));
  }

  if (fabs(x1 - x0) < FLOAT_NUMBERS_EQUALITY_TOLERANCE)
    distribute_snapshot_points_YZ_plane(x0, y0, y1, z0, z1);
  else if (fabs(y1 - y0) < FLOAT_NUMBERS_EQUALITY_TOLERANCE)
    distribute_snapshot_points_XZ_plane(x0, x1, y0, z0, z1);
  else if (fabs(z1 - z0) < FLOAT_NUMBERS_EQUALITY_TOLERANCE)
    distribute_snapshot_points_XY_plane(x0, x1, y0, y1, z0);
  else MFEM_ABORT("Snapshot plane is not parallel to XY, XZ or YZ");
}

void SnapshotsPlane::distribute_snapshot_points_YZ_plane(double x,
                                                         double y0, double y1,
                                                         double z0, double z1)
{
  _plane = "YZ";

  const double dy = (y1 - y0) / (_n_points_1-1);
  const double dz = (z1 - z0) / (_n_points_2-1);

  int p = 0;
  for (int i = 0; i < _n_points_2; ++i)
  {
    double z = (i == _n_points_2-1 ? z1 : z0 + i*dz);
    for (int j = 0; j < _n_points_1; ++j)
    {
      double y = (j == _n_points_1-1 ? y1 : y0 + j*dy);
      _snapshot_points[p++] = Vertex(x, y, z);
    }
  }
}

void SnapshotsPlane::distribute_snapshot_points_XZ_plane(double x0, double x1,
                                                         double y,
                                                         double z0, double z1)
{
  _plane = "XZ";

  const double dx = (x1 - x0) / (_n_points_1-1);
  const double dz = (z1 - z0) / (_n_points_2-1);

  int p = 0;
  for (int i = 0; i < _n_points_2; ++i)
  {
    double z = (i == _n_points_2-1 ? z1 : z0 + i*dz);
    for (int j = 0; j < _n_points_1; ++j)
    {
      double x = (j == _n_points_1-1 ? x1 : x0 + j*dx);
      _snapshot_points[p++] = Vertex(x, y, z);
    }
  }
}

void SnapshotsPlane::distribute_snapshot_points_XY_plane(double x0, double x1,
                                                         double y0, double y1,
                                                         double z)
{
  _plane = "XY";

  const double dx = (x1 - x0) / (_n_points_1-1);
  const double dy = (y1 - y0) / (_n_points_2-1);

  int p = 0;
  for (int i = 0; i < _n_points_2; ++i)
  {
    double y = (i == _n_points_2-1 ? y1 : y0 + i*dy);
    for (int j = 0; j < _n_points_1; ++j)
    {
      double x = (j == _n_points_1-1 ? x1 : x0 + j*dx);
      _snapshot_points[p++] = Vertex(x, y, z);
    }
  }
}

std::string SnapshotsPlane::description() const
{
  const double x = _vertices[0](0);
  const double y = _vertices[0](1);
  const double z = _vertices[0](2);
  std::string coord = (_plane == "YZ" ? "_x" + d2s(x) :
                                        (_plane == "XZ" ? "_y" + d2s(y) :
                                                          "_z" + d2s(z)));
  return "_snap_plane" + _plane + coord;
}

void SnapshotsPlane::
save_snapshot_vector_vts(const Vector &snapshot, int time_step,
                         const std::string &var_name, int n_components,
                         const std::string &output_file_base) const
{
  // string representing a time step: 6 digits (first filled with 0)
  const std::string tstep_str = d2s(time_step, 0, 0, 0, 6);

  std::string fname = file_path(output_file_base) + "/" +
                      file_stem(output_file_base) + "_" + var_name + "_t" +
                      tstep_str + ".vts";

  std::ofstream out(fname.c_str(), std::ios::binary);
  MFEM_VERIFY(out, "File '" + fname + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"" << endianness() << "\">\n";
  if (_plane == "YZ") {
    out << "  <StructuredGrid WholeExtent=\"1 1 1 " << _n_points_1 << " 1 " << _n_points_2 << "\">\n";
    out << "    <Piece Extent=\"1 1 1 " << _n_points_1 << " 1 " << _n_points_2 << "\">\n";
  } else if (_plane == "XZ") {
    out << "  <StructuredGrid WholeExtent=\"1 " << _n_points_1 << " 1 1 1 " << _n_points_2 << "\">\n";
    out << "    <Piece Extent=\"1 " << _n_points_1 << " 1 1 1 " << _n_points_2 << "\">\n";
  } else if (_plane == "XY") {
    out << "  <StructuredGrid WholeExtent=\"1 " << _n_points_1 << " 1 " << _n_points_2 << " 1 1\">\n";
    out << "    <Piece Extent=\"1 " << _n_points_1 << " 1 " << _n_points_2 << " 1 1\">\n";
  } else MFEM_ABORT("Unknown plane: " + _plane);

  out << "      <PointData>\n";
  for (int c = 0; c < n_components; ++c) {
    out << "        <DataArray type=\"Float64\" Name=\"" << var_name << "_" << c << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";
    for (int p = c; p < snapshot.Size(); p += n_components) {
      out << snapshot(p) << " ";
    }
    out << "\n";
    out << "        </DataArray>\n";
  }
  out << "        <DataArray type=\"Float64\" Name=\"" << var_name << "_mag\" format=\"ascii\" NumberOfComponents=\"1\">\n";
  for (int p = 0; p < snapshot.Size(); p += n_components) {
    double mag = 0.;
    for (int c = 0; c < n_components; ++c)
      mag += snapshot(p+c)*snapshot(p+c);
    out << sqrt(mag) << " ";
  }
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";

  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int p = 0; p < _n_snapshot_points; ++p) {
    const Vertex &node = _snapshot_points[p];
    out << node(0) << " " << node(1) << " " << node(2) << " ";
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";
}






//==============================================================================
//
// SnapshotsVolume
//
//==============================================================================
SnapshotsVolume::SnapshotsVolume()
  : _min_coord(), _max_coord(), _n_points_x(0), _n_points_y(0), _n_points_z(0)
{ }

void SnapshotsVolume::init(std::ifstream &in)
{
  MFEM_VERIFY(in.is_open(), "The stream for reading snapshot plane is not open");

  std::string tmp;

  in >> _variable; getline(in, tmp);
  in >> _format;   getline(in, tmp);
  in >> _n_points_x >> _n_points_y >> _n_points_z; getline(in, tmp);
  in >> _min_coord(0) >> _min_coord(1) >> _min_coord(2); getline(in, tmp);
  in >> _max_coord(0) >> _max_coord(1) >> _max_coord(2); getline(in, tmp);

  _n_snapshot_points = _n_points_x*_n_points_y*_n_points_z;
  MFEM_VERIFY(_n_snapshot_points > 0, "The number of snapshot points (" +
              d2s(_n_snapshot_points) + ") must be >0");
}

void SnapshotsVolume::distribute_snapshot_points()
{
  _snapshot_points = new mfem::Vertex[_n_snapshot_points];

  double x0 = _min_coord(0);
  double y0 = _min_coord(1);
  double z0 = _min_coord(2);
  double x1 = _max_coord(0);
  double y1 = _max_coord(1);
  double z1 = _max_coord(2);

  const double dx = (x1 - x0) / (_n_points_x-1);
  const double dy = (y1 - y0) / (_n_points_y-1);
  const double dz = (z1 - z0) / (_n_points_z-1);

  int p = 0;
  for (int iz = 0; iz < _n_points_z; ++iz)
  {
    double z = (iz == _n_points_z-1 ? z1 : z0 + iz*dz);
    for (int iy = 0; iy < _n_points_y; ++iy)
    {
      double y = (iy == _n_points_y-1 ? y1 : y0 + iy*dy);
      for (int ix = 0; ix < _n_points_x; ++ix)
      {
        double x = (ix == _n_points_x-1 ? x1 : x0 + ix*dx);
        _snapshot_points[p++] = Vertex(x, y, z);
      }
    }
  }
}

std::string SnapshotsVolume::description() const
{
  return "_snap_vol";
}

void SnapshotsVolume::
save_snapshot_vector_vts(const Vector &snapshot, int time_step,
                         const std::string &var_name, int n_components,
                         const std::string &output_file_base) const
{
  // string representing a time step: 6 digits (first filled with 0)
  const std::string tstep_str = d2s(time_step, 0, 0, 0, 6);

  std::string fname = file_path(output_file_base) + "/" +
                      file_stem(output_file_base) + "_" + var_name + "_t" +
                      tstep_str + ".vts";

  std::ofstream out(fname.c_str(), std::ios::binary);
  MFEM_VERIFY(out, "File '" + fname + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"" << endianness() << "\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << _n_points_x << " 1 " << _n_points_y << " 1 " << _n_points_z << "\">\n";
  out << "    <Piece Extent=\"1 " << _n_points_x << " 1 " << _n_points_y << " 1 " << _n_points_z << "\">\n";

  out << "      <PointData>\n";
  for (int c = 0; c < n_components; ++c) {
    out << "        <DataArray type=\"Float64\" Name=\"" << var_name << "_" << c << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";
    for (int p = c; p < snapshot.Size(); p += n_components) {
      out << snapshot(p) << " ";
    }
    out << "\n";
    out << "        </DataArray>\n";
  }
  out << "        <DataArray type=\"Float64\" Name=\"" << var_name << "_mag\" format=\"ascii\" NumberOfComponents=\"1\">\n";
  for (int p = 0; p < snapshot.Size(); p += n_components) {
    double mag = 0.;
    for (int c = 0; c < n_components; ++c)
      mag += snapshot(p+c)*snapshot(p+c);
    out << sqrt(mag) << " ";
  }
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";

  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int p = 0; p < _n_snapshot_points; ++p) {
    const Vertex &node = _snapshot_points[p];
    out << node(0) << " " << node(1) << " " << node(2) << " ";
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";
}




void save_snapshot_vector_bin(const Vector &snapshot, int time_step,
                              const std::string &var_name, int n_components,
                              const std::string &output_file_base)
{
  // string representing a time step: 6 digits (first filled with 0)
  const std::string tstep_str = d2s(time_step, 0, 0, 0, 6);

  std::vector<std::string> fnames(n_components);
  for (int c = 0; c < n_components; ++c) {
    fnames[c] = file_path(output_file_base) + "/" +
                file_stem(output_file_base) + "_" + var_name + d2s(c) + "_t" +
                tstep_str + ".bin";
  }

  std::vector<std::ofstream> outs(n_components);
  for (int c = 0; c < n_components; ++c) {
    outs[c].open(fnames[c].c_str(), std::ios::binary);
    MFEM_VERIFY(outs[c], "File '" + fnames[c] + "' can't be opened");
  }

  for (int p = 0; p < snapshot.Size(); p += n_components) {
    for (int c = 0; c < n_components; ++c) {
      float val = snapshot(p+c);
      outs[c].write(reinterpret_cast<char*>(&val), sizeof(val));
    }
  }
}



