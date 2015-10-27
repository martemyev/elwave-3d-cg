#ifndef RECEIVERS_HPP
#define RECEIVERS_HPP

#include "config.hpp"
#include "utilities.hpp"

#include <fstream>
#include <vector>

namespace mfem { class Vertex; }

/**
 * Abstract class representing a set (a straight line, a circle, or other line)
 * of receivers (stations).
 */
class ReceiversSet
{
public:

  ReceiversSet();
  ReceiversSet(const ReceiversSet& rec);
  ReceiversSet& operator =(const ReceiversSet& rec);
  virtual ~ReceiversSet() { }

  /**
   * Find and save the numbers of cells containing the receivers.
   * @param nx - number of cells in x-direction
   * @param ny - number of cells in y-direction
   * @param nz - number of cells in z-direction
   * @param sx - size of domain in x-direction
   * @param sy - size of domain in y-direction
   * @param sz - size of domain in z-direction
   */
  void find_cells_containing_receivers(int nx, int ny, int nz, double sx,
                                       double sy, double sz);

  int n_receivers() const { return _n_receivers; }

  const std::vector<mfem::Vertex>& get_receivers() const
  { return _receivers; }

  const std::vector<int>& get_cells_containing_receivers() const
  { return _cells_containing_receivers; }

//  void print_receivers(const mfem::Mesh& mesh,
//                       std::ostream& out = std::cout) const;

  /**
   * Initialize the parameters of the receivers set reading them from a given
   * and already open input stream (likely connected to a file).
   */
  virtual void init(std::ifstream &in) = 0;

  /**
   * Distribute the receivers along the corresponding line.
   */
  virtual void distribute_receivers() = 0;

  /**
   * Description of the set of receivers (to distinguish the derived sets).
   */
  virtual std::string description() const = 0;


protected:

  /**
   * Number of receivers (stations) in the set.
   */
  int _n_receivers;

  /**
   * The locations of the receivers (stations).
   */
  std::vector<mfem::Vertex> _receivers;

  /**
   * Numbers of grid cells containing the receivers.
   */
  std::vector<int> _cells_containing_receivers;
};




/**
 * A class representing a straight line of receivers.
 */
class ReceiversLine: public ReceiversSet
{
public:
  ReceiversLine();
  ReceiversLine(const ReceiversLine& rec);
  ReceiversLine& operator =(const ReceiversLine& rec);
  ~ReceiversLine() { }
  void init(std::ifstream &in);
  void distribute_receivers();
  std::string description() const;
protected:
  mfem::Vertex _start; ///< beginning of line of recievers
  mfem::Vertex _end;   ///< end       of line of receivers
};




/**
 * A class representing a plane of receivers distributed at the vertices of the
 * Cartesian grid built over this plane.
 */
class ReceiversPlane: public ReceiversSet
{
public:
  ReceiversPlane();
  ReceiversPlane(const ReceiversPlane& rec);
  ReceiversPlane& operator =(const ReceiversPlane& rec);
  ~ReceiversPlane() { }
  void init(std::ifstream &in);
  void distribute_receivers();
  std::string description() const;
protected:
  int _n_receivers_1; ///< Number of receivers in one direction
  int _n_receivers_2; ///< Number of receivers in another direction
  std::string _plane; ///< To which plane this plane of receivers is parallel
  double _plane_coord;///< Coordinate of the axis orthogonal to the plane
  mfem::Vertex _start;///< Coordinates of the boundaries (begin) of the receivers plane
  mfem::Vertex _end;  ///< Coordinates of the boundaries (end) of the receivers plane
};



int find_element(double sx, double sy, double sz, int nx, int ny, int nz,
                 const mfem::Vertex &point, bool throw_exception);


#endif // RECEIVERS_HPP
