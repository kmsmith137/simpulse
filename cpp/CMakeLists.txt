add_library(simpulse SHARED
  constant_acceleration_phase_model.cpp
  keplerian_binary_phase_model.cpp
  phase_model_base.cpp
  single_pulse.cpp
  von_mises_profile.cpp
)

target_link_libraries(simpulse ${SIMPULSE_LIBS})

install(TARGETS simpulse DESTINATION lib)
