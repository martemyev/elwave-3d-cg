#ifndef SNAPSHOTS_HPP
#define SNAPSHOTS_HPP

#include "config_elwave.hpp"

#include <fstream>
//#include <vector>

#include "mfem.hpp"

//namespace mfem {
//  class Vector;
//  class Vertex;
//}


/**
 * Abstract class representing a set (a plane or a volume) where the snapshot of
 * a solution is taken and saved in a file in a specified format.
 */
class SnapshotsSet
{
public:

  virtual ~SnapshotsSet() { delete[] _snapshot_points; }

  /**
   * Find and save the numbers of coarse and fine cells containing the points
   * where a solution is computed to be saved.
   */
  void find_cells_containing_snapshot_points(int nx, int ny, int nz, double sx,
                                             double sy, double sz);

  std::string get_variable() const { return _variable; }

  int n_snapshot_points() const { return _n_snapshot_points; }

  const mfem::Vertex* get_snapshot_points() const
  { return _snapshot_points; }

  const int* get_cells_containing_snapshot_points() const
  { return _cells_contain_snapshot_points; }

  /**
   * Initialize the parameters of the snapshots set reading them from a given
   * and already open input stream (likely connected to a file).
   */
  virtual void init(std::ifstream &in) = 0;

  /**
   * Distribute the snapshot points in a corresponding set.
   */
  virtual void distribute_snapshot_points() = 0;

  /**
   * Description of the set of snapshots (to distinguish the derived sets).
   */
  virtual std::string description() const = 0;

  /**
   * Save a snapshot of a vector solution in a file.
   */
  void save_snapshot_vector(const mfem::Vector &snapshot, int time_step,
                            const std::string &var_name, int n_components,
                            const std::string &output_file_base) const;

protected:

  /**
   * Which variable to save: displacement (U), velocity (V) or both (UV)
   */
  std::string _variable;

  /**
   * Format in which a snapshot of a solution will be saved in a file.
   */
  std::string _format;

  /**
   * Total number of snapshot points in the set, where a solution is computed to
   * be saved.
   */
  int _n_snapshot_points;

  /**
   * The locations of the snapshot points.
   */
//  std::vector<mfem::Vertex> _snapshot_points;
  mfem::Vertex *_snapshot_points;

  /**
   * Numbers (indices) of grid cells containing the receivers.
   */
  int *_cells_contain_snapshot_points;

  /**
   * Save a snapshot of a vector solution in a file in VTS format
   */
  virtual void save_snapshot_vector_vts(
      const mfem::Vector &snapshot, int time_step,
      const std::string &var_name, int n_components,
      const std::string &output_file_base) const = 0;

  SnapshotsSet();
  SnapshotsSet(const SnapshotsSet&);
  SnapshotsSet& operator =(const SnapshotsSet&);
};




/**
 * Class representing a plane where snapshots of a solution will be taken.
 */
class SnapshotsPlane: public SnapshotsSet
{
public:
  SnapshotsPlane();
  ~SnapshotsPlane() { }
  void init(std::ifstream &in);
  void distribute_snapshot_points();
  std::string description() const;

protected:
  static const int N_VERTICES = 4; ///< Number of points to describe a plane
  mfem::Vertex _vertices[N_VERTICES];
  int _n_points; ///< Number of points in one direction, i.e. total number of
                 ///< snapshot points will be _n_points^2
  std::string _plane; ///< Description of the plane orientation

  void distribute_snapshot_points_YZ_plane(double x, double y0, double y1,
                                           double z0, double z1);
  void distribute_snapshot_points_XZ_plane(double x0, double x1, double y,
                                           double z0, double z1);
  void distribute_snapshot_points_XY_plane(double x0, double x1, double y0,
                                           double y1, double z);
  void save_snapshot_vector_vts(const mfem::Vector &snapshot, int time_step,
                                const std::string &var_name, int n_components,
                                const std::string &output_file_base) const;

  SnapshotsPlane(const SnapshotsPlane&);
  SnapshotsPlane& operator =(const SnapshotsPlane&);
};




/**
 * Class representing a volume where snapshots of a solution will be taken.
 */
class SnapshotsVolume: public SnapshotsSet
{
public:
  SnapshotsVolume();
  ~SnapshotsVolume() { }
  void init(std::ifstream &in);
  void distribute_snapshot_points();
  std::string description() const;

protected:
  mfem::Vertex _min_coord;
  mfem::Vertex _max_coord;
  int _n_points_x; ///< Number of points in x-direction
  int _n_points_y; ///< Number of points in y-direction
  int _n_points_z; ///< Number of points in z-direction

  void save_snapshot_vector_vts(const mfem::Vector &snapshot, int time_step,
                                const std::string &var_name, int n_components,
                                const std::string &output_file_base) const;

  SnapshotsVolume(const SnapshotsVolume&);
  SnapshotsVolume& operator =(const SnapshotsVolume&);
};



void save_snapshot_vector_bin(const mfem::Vector &snapshot, int time_step,
                              const std::string &var_name, int n_components,
                              const std::string &output_file_base);


#endif // SNAPSHOTS_HPP
