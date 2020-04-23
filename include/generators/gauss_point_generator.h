#ifndef MPM_GAUSS_POINT_GENERATOR_H_
#define MPM_GAUSS_POINT_GENERATOR_H_

#include "point_generator.h"

namespace mpm {

//! FilePointGenerator class
//! \brief Base class that defines generation of material points at cell gauss
//! locations \tparam Tdim Dimension
template <unsigned Tdim>
class GaussPointGenerator : public PointGenerator<Tdim> {

 public:
  //! Constructor with mesh pointer and generator properties
  GaussPointGenerator(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh);

  //! Virtual destructor
  ~GaussPointGenerator() override{};

  //! Delete copy constructor
  GaussPointGenerator(const GaussPointGenerator<Tdim>&) = delete;

  //! Delete assignment operator
  GaussPointGenerator& operator=(const GaussPointGenerator<Tdim>&) = delete;

  //! Generate material points
  bool generate_points(const std::shared_ptr<mpm::IO>& io,
                       const Json& generator_props) override;

 private:
  //! Mesh
  using PointGenerator<Tdim>::mesh_;
  //! Logger
  using PointGenerator<Tdim>::console_;

};  // FilePointGenerator class
}  // namespace mpm

#include "gauss_point_generator.tcc"

#endif  // MPM_GAUSSPOINTGENERATOR_H_
