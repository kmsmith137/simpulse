add_library(libsimpulse SHARED
  constant_acceleration_phase_model.cpp
  phase_model_base.cpp
  single_pulse.cpp
  von_mises_profile.cpp
)

target_link_libraries(libsimpulse ${LIBSIMPULSE_LIBS})

install(TARGETS libsimpulse DESTINATION lib)
