// Constructor
template <unsigned Tdim>
mpm::GaussPointGenerator<Tdim>::GaussPointGenerator(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh)
    : mpm::PointGenerator<Tdim>(mesh) {

  //! Logger
  std::string logger = "gauss-point-generator" + std::to_string(Tdim) + "d";
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Generate points
template <unsigned Tdim>
bool mpm::GaussPointGenerator<Tdim>::generate_points(
    const std::shared_ptr<mpm::IO>& io, const Json& generator_props) {

  // Number of particles per dir
  unsigned nparticles_dir =
      generator_props["nparticles_per_dir"].template get<unsigned>();
  // Particle type
  auto particle_type =
      generator_props["particle_type"].template get<std::string>();
  // Material id
  unsigned material_id =
      generator_props["material_id"].template get<unsigned>();
  // Cell set id
  int cset_id = generator_props["cset_id"].template get<int>();

  bool status = mesh_->generate_material_points(nparticles_dir, particle_type,
                                                material_id, cset_id);
  if (!status)
    throw std::runtime_error(
        "Particle generation at cell gauss locations failed");

  return status;
}
