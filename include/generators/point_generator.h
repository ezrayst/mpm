#ifndef MPM_POINT_GENERATOR_H_
#define MPM_POINT_GENERATOR_H_

#include <ctime>

#include "io.h"
#include "mesh.h"

// JSON
using Json = nlohmann::json;

namespace mpm {
//! PointGenerator class
//! \brief Base class that defines generation of material points
//! \tparam Tdim Dimension
template <unsigned Tdim>
class PointGenerator {
 public:
  //! Constructor with mesh pointer and generator properties
  PointGenerator(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh) : mesh_{mesh} {}

  //! Destructor
  virtual ~PointGenerator(){};

  //! Delete copy constructor
  PointGenerator(const PointGenerator<Tdim>&) = delete;

  //! Delete assignement operator
  PointGenerator& operator=(const PointGenerator<Tdim>&) = delete;

  //! Generate material points
  virtual bool generate_points(const std::shared_ptr<mpm::IO>& io,
                               const Json& generator_props) = 0;

 protected:
  //! Mesh
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;

};  // PointGenerator class
}  // namespace mpm
#endif  // MPM_POINTGENERATOR_H_
