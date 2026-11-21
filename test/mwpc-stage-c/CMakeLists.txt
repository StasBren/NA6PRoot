cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(NA6PMWPCStageC LANGUAGES CXX)

find_package(Garfield REQUIRED)

add_executable(mwpc_stage_c1_full_muon c1_full_muon_response.cxx)
target_link_libraries(mwpc_stage_c1_full_muon PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c1_full_muon PRIVATE cxx_std_17)

add_executable(mwpc_stage_c2_weighting_probe c2_numerical_weighting_probe.cxx)
target_link_libraries(mwpc_stage_c2_weighting_probe PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c2_weighting_probe PRIVATE cxx_std_17)

add_executable(mwpc_stage_c2_full_muon_numerical c2_full_muon_numerical_response.cxx)
target_link_libraries(mwpc_stage_c2_full_muon_numerical PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c2_full_muon_numerical PRIVATE cxx_std_17)

add_executable(mwpc_stage_c2_event_3d_dump c2_event_3d_dump.cxx)
target_link_libraries(mwpc_stage_c2_event_3d_dump PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c2_event_3d_dump PRIVATE cxx_std_17)


add_executable(mwpc_stage_c_wire_response c_wire_coordinate_response.cxx)
target_link_libraries(mwpc_stage_c_wire_response PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c_wire_response PRIVATE cxx_std_17)


add_executable(mwpc_stage_c_transport_box_map c_transport_box_map.cxx)
target_link_libraries(mwpc_stage_c_transport_box_map PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c_transport_box_map PRIVATE cxx_std_17)

add_executable(pi_garfield_drift_methods_benchmark pi_garfield_drift_methods_benchmark.cxx)
target_link_libraries(pi_garfield_drift_methods_benchmark PRIVATE Garfield::Garfield)
target_compile_features(pi_garfield_drift_methods_benchmark PRIVATE cxx_std_17)

add_executable(pi_garfield_attachment_field_probe pi_garfield_attachment_field_probe.cxx)
target_link_libraries(pi_garfield_attachment_field_probe PRIVATE Garfield::Garfield)
target_compile_features(pi_garfield_attachment_field_probe PRIVATE cxx_std_17)

add_executable(pi_fast_lut_existing_starts_mc pi_fast_lut_existing_starts_mc.cxx)
target_link_libraries(pi_fast_lut_existing_starts_mc PRIVATE Garfield::Garfield)
target_compile_features(pi_fast_lut_existing_starts_mc PRIVATE cxx_std_17)

add_executable(pi_fast_lut_attachment_hazards pi_fast_lut_attachment_hazards.cxx)
target_link_libraries(pi_fast_lut_attachment_hazards PRIVATE Garfield::Garfield)
target_compile_features(pi_fast_lut_attachment_hazards PRIVATE cxx_std_17)


# Stage C candidate P3 2100V: read-only 3D Elmer/Gmsh individual U/V
# stereo-strip near-anode FIELD potential+vector-gradient diagnostics.
# Never runs avalanches, ion transport, induced current or electronics.
add_executable(mwpc_stage_c_fem_stereo_near_wire_probe
  fem_stereo_strip_near_wire_probe.cxx)
target_link_libraries(mwpc_stage_c_fem_stereo_near_wire_probe PRIVATE
  Garfield::Garfield)
target_compile_features(mwpc_stage_c_fem_stereo_near_wire_probe PRIVATE cxx_std_17)


# Stage-C P3 2100V: dependency-free verified 10ns Ramo bin arithmetic.
# This test does not read any FEM or microscopic archive.
add_executable(mwpc_stage_c_time_binned_ramo_charge_test
  cathode-ion-reference/test_time_binned_ramo_charge.cc)
target_compile_features(mwpc_stage_c_time_binned_ramo_charge_test PRIVATE cxx_std_17)

# Independent K256 archive -> P3 ideal physical ion drift -> FEM U/V
# central selected-strip potential differences. PROVISIONAL weighting geometry
# until the curved-wire map and strip sharing are separately validated.
# This target never regenerates microscopic avalanches.
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
add_executable(mwpc_stage_c_stereo_ion_response_2100
  cathode-ion-reference/stereo_ion_response_from_ions_root.cc)
target_link_libraries(mwpc_stage_c_stereo_ion_response_2100 PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(mwpc_stage_c_stereo_ion_response_2100 PRIVATE cxx_std_17)

# Three-wire early cumulative-charge sensitivity on fixed analytic trajectories.
add_executable(mwpc_stage_c_p3_compare_early_ion_paths
  cathode-ion-reference/compare_p3_fixed_early_ion_paths.cc)
target_link_libraries(mwpc_stage_c_p3_compare_early_ion_paths PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(mwpc_stage_c_p3_compare_early_ion_paths PRIVATE cxx_std_17)

# Independent zero-mode analytical whole-cathode comparison grid.
add_executable(mwpc_stage_c_p3_whole_plane_fourier_reference
  cathode-ion-reference/p3_whole_plane_fourier_reference.cc)
target_link_libraries(mwpc_stage_c_p3_whole_plane_fourier_reference PRIVATE Garfield::Garfield)
target_compile_features(mwpc_stage_c_p3_whole_plane_fourier_reference PRIVATE cxx_std_17)

add_executable(mwpc_stage_c_p3_export_spectral_archive_paths
  cathode-ion-reference/p3_export_spectral_archive_paths.cc)
target_link_libraries(mwpc_stage_c_p3_export_spectral_archive_paths PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(mwpc_stage_c_p3_export_spectral_archive_paths PRIVATE cxx_std_17)

add_executable(mwpc_stage_c_p3_export_full_spectral_paths
  cathode-ion-reference/p3_export_full_spectral_paths.cc)
target_link_libraries(mwpc_stage_c_p3_export_full_spectral_paths PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(mwpc_stage_c_p3_export_full_spectral_paths PRIVATE cxx_std_17)
