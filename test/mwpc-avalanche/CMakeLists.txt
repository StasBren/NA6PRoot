cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(NA6PMWPCAvalanchePhaseA LANGUAGES CXX)

find_package(Garfield REQUIRED)

add_executable(mwpc_phase_a phaseA_single_electron.cxx)
target_link_libraries(mwpc_phase_a PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a PRIVATE cxx_std_17)

add_executable(mwpc_phase_a_collection phaseA_collection_scan.cxx)
target_link_libraries(mwpc_phase_a_collection PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_collection PRIVATE cxx_std_17)

add_executable(mwpc_phase_a_coordinates phaseA_coordinate_test.cxx)
target_compile_features(mwpc_phase_a_coordinates PRIVATE cxx_std_17)

add_executable(mwpc_phase_a_gain phaseA_gain_fluctuations.cxx)
target_link_libraries(mwpc_phase_a_gain PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_gain PRIVATE cxx_std_17)

add_executable(mwpc_phase_a_magnetic phaseA_magnetic_drift.cxx)
target_link_libraries(mwpc_phase_a_magnetic PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_magnetic PRIVATE cxx_std_17)

add_executable(mwpc_phase_a_heed phaseA_heed_muon.cxx)
target_link_libraries(mwpc_phase_a_heed PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_heed PRIVATE cxx_std_17)

add_executable(mwpc_phase_a_heed_avalanche phaseA_heed_avalanche.cxx)
target_link_libraries(mwpc_phase_a_heed_avalanche PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_heed_avalanche PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_weighting phaseB_strip_weighting.cxx)
target_link_libraries(mwpc_phase_b_weighting PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_weighting PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_single_ion phaseB_single_ion_signal.cxx)
target_link_libraries(mwpc_phase_b_single_ion PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_single_ion PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_electron_ion_pair phaseB_electron_ion_pair.cxx)
target_link_libraries(mwpc_phase_b_electron_ion_pair PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_electron_ion_pair PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_ion_cloud phaseB_ion_cloud_signal.cxx)
target_link_libraries(mwpc_phase_b_ion_cloud PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_ion_cloud PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_electron_ion_cloud phaseB_electron_ion_cloud_signal.cxx)
target_link_libraries(mwpc_phase_b_electron_ion_cloud PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_electron_ion_cloud PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_real_avalanche phaseB_real_avalanche_signal.cxx)
target_link_libraries(mwpc_phase_b_real_avalanche PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_real_avalanche PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_resolution_ensemble phaseB_resolution_ensemble.cxx)
target_link_libraries(mwpc_phase_b_resolution_ensemble PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_resolution_ensemble PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_strip_geometry_scan phaseB_strip_geometry_scan.cxx)
target_link_libraries(mwpc_phase_b_strip_geometry_scan PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_strip_geometry_scan PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_tilt_scan phaseB_tilt_scan.cxx)
target_link_libraries(mwpc_phase_b_tilt_scan PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_tilt_scan PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_opposite_cathode_plane_coupling phaseB_opposite_cathode_plane_coupling.cxx)
target_link_libraries(mwpc_phase_b_opposite_cathode_plane_coupling PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_opposite_cathode_plane_coupling PRIVATE cxx_std_17)

add_executable(mwpc_heed_visualization heed_muon_visualization.cxx)
target_link_libraries(mwpc_heed_visualization PRIVATE Garfield::Garfield)
target_compile_features(mwpc_heed_visualization PRIVATE cxx_std_17)


add_executable(mwpc_phase_a_avalanche_visualization phaseA_avalanche_visualization.cxx)
target_link_libraries(mwpc_phase_a_avalanche_visualization PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_avalanche_visualization PRIVATE cxx_std_17)


add_executable(mwpc_phase_a_heed_stats phaseA_heed_statistics.cxx)
target_link_libraries(mwpc_phase_a_heed_stats PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a_heed_stats PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_opposite_cathode_stereo phaseB_opposite_cathode_stereo.cxx)
target_link_libraries(mwpc_phase_b_opposite_cathode_stereo PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_opposite_cathode_stereo PRIVATE cxx_std_17)

add_executable(mwpc_phase_b_avalanche_angle_scan phaseB_avalanche_angle_scan.cxx)
target_link_libraries(mwpc_phase_b_avalanche_angle_scan PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_avalanche_angle_scan PRIVATE cxx_std_17)

# PI muon-driven B1: original immutable Geant4 gas chords + true TrackHeed
# clusters and every seed's microscopic ionization. New ±3mm 4-geometry code.
add_executable(mwpc_phase_b_native_muon_3p3 phaseB_native_muon_3p3_avalanche.cxx)
target_link_libraries(mwpc_phase_b_native_muon_3p3 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_native_muon_3p3 PRIVATE cxx_std_17)

# Research-only exact microscopic muon reference with enhanced independent
# 128-point per-primary ion-birth spatial/time reservoir; preserves original target.
add_executable(mwpc_phase_b_native_muon_3p3_spacetime_v1 phaseB_native_muon_3p3_spacetime_v1.cxx)
target_link_libraries(mwpc_phase_b_native_muon_3p3_spacetime_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_native_muon_3p3_spacetime_v1 PRIVATE cxx_std_17)

# Isolated source of the joint per-electron ion-cloud LUT.  The new modes
# save 128 sampled or ALL exact ion birth coordinates for heldout compression
# checks.  Never edit the pinned native-muon physics sources.
add_executable(mwpc_phase_b_joint_cloud_lut_v1 phaseB_native_muon_joint_cloud_lut_v1.cxx)
target_link_libraries(mwpc_phase_b_joint_cloud_lut_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_joint_cloud_lut_v1 PRIVATE cxx_std_17)

# Separate phase-C research-only macroscopic RKF/Polya baseline. The archived
# full microscopic native B1 target/physics and source are left untouched.
add_executable(mwpc_phase_c_rkf_polya_pilot phaseC_rkf_polya_pilot.cxx)
target_link_libraries(mwpc_phase_c_rkf_polya_pilot PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_c_rkf_polya_pilot PRIVATE cxx_std_17)

# Separate read-only physical comparison V2, preserves SHA-sealed V1 executable/source.
add_executable(mwpc_phase_c_rkf_ion_census_pilot_v2 phaseC_rkf_ion_census_pilot_v2.cxx)
target_link_libraries(mwpc_phase_c_rkf_ion_census_pilot_v2 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_c_rkf_ion_census_pilot_v2 PRIVATE cxx_std_17)

# Separate full-16k archived A100-start RKF-only macroscopic estimate benchmark.
add_executable(mwpc_phase_c_rkf_ion_census_full16k_v3 phaseC_rkf_ion_census_full16k_v3.cxx)
target_link_libraries(mwpc_phase_c_rkf_ion_census_full16k_v3 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_c_rkf_ion_census_full16k_v3 PRIVATE cxx_std_17)

# Isolated native Geant4 J/psi TrackHeed -> FULL Garfield benchmark for
# p3_noFW ONLY; per-primary and per-muon native timings. Does NOT replace
# the frozen four-geometry native/spacetime reference or their binaries.
add_executable(mwpc_phase_b_native_p3_noFW_jpsi_bench_v1 phaseB_native_muon_p3_noFW_jpsi_runtime_v1.cxx)
target_link_libraries(mwpc_phase_b_native_p3_noFW_jpsi_bench_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_native_p3_noFW_jpsi_bench_v1 PRIVATE cxx_std_17)

# Emergency high-Heed-multiplicity recovery ONLY. Physics identical to the
# frozen p3_noFW native J/psi executable except accepted-primary safety
# budget of 50000 (rather than 1000). Build this independent target alone.
add_executable(mwpc_phase_b_native_p3_noFW_jpsi_highheed_recovery_v1 phaseB_native_muon_p3_noFW_jpsi_highheed_recovery_v1.cxx)
target_link_libraries(mwpc_phase_b_native_p3_noFW_jpsi_highheed_recovery_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_b_native_p3_noFW_jpsi_highheed_recovery_v1 PRIVATE cxx_std_17)

# p3_noFW-only Stage C: two complete cathode planes (NO strips yet).
# Real 9-wire analytic drift and weighting fields, RKF positive-ion signal.
# Independent target. Does NOT alter any previous Garfield production binary.
add_executable(mwpc_phase_c_p3_noFW_two_cathode_ion_v1 phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx)
target_link_libraries(mwpc_phase_c_p3_noFW_two_cathode_ion_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_c_p3_noFW_two_cathode_ion_v1 PRIVATE cxx_std_17)

# SCIENCE-ONLY comparison: original frozen p3_noFW (no Penning) vs
# independent identical geometry/gas/TrackHeed/Microscopic clone with
# Ar(70)/CO2(30) Penning r=0.57; no production source/binary modification.
add_executable(mwpc_p3_noFW_jpsi_penning057_pilot_v1 phaseB_native_muon_p3_noFW_jpsi_penning057_pilot_v1.cxx)
target_link_libraries(mwpc_p3_noFW_jpsi_penning057_pilot_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_p3_noFW_jpsi_penning057_pilot_v1 PRIVATE cxx_std_17)

# Independent Penning-ON 10k PBS campaign executable. Only accepted TrackHeed
# primary safety bound is increased to 50000; physical avalanche cap unchanged.
add_executable(mwpc_p3_noFW_jpsi_penning057_highheed_v1 phaseB_native_muon_p3_noFW_jpsi_penning057_highheed_v1.cxx)
target_link_libraries(mwpc_p3_noFW_jpsi_penning057_highheed_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_p3_noFW_jpsi_penning057_highheed_v1 PRIVATE cxx_std_17)

# Independent 1000-real-muon truth pilot (same Penning-on Garfield physics).
# Output includes actual physical collecting wires and exact per-wire ion birth
# moment census. The old sealed 30k executables are never modified or relinked.
add_executable(mwpc_p3_noFW_penning057_wiretruth_v1
               phaseB_native_muon_p3_noFW_penning057_wiretruth_v1.cxx)
target_link_libraries(mwpc_p3_noFW_penning057_wiretruth_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_p3_noFW_penning057_wiretruth_v1 PRIVATE cxx_std_17)
# Isolated runtime-HV version; original 1800 V wire-truth target stays frozen.
add_executable(mwpc_p3_noFW_penning057_hvscan_v1
               phaseB_native_muon_p3_noFW_penning057_hvscan_v1.cxx)
target_link_libraries(mwpc_p3_noFW_penning057_hvscan_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_p3_noFW_penning057_hvscan_v1 PRIVATE cxx_std_17)

# Isolated high-HV microscopic wiretruth pilot; original voltage scan unchanged.
add_executable(mwpc_p3_noFW_penning057_hvhigh_v1 phaseB_native_muon_p3_noFW_penning057_hvhigh_v1.cxx)
target_link_libraries(mwpc_p3_noFW_penning057_hvhigh_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_p3_noFW_penning057_hvhigh_v1 PRIVATE cxx_std_17)

# Isolated uncapped 2100-V full-microscopic wire-truth run (cluster pilot / 1000 muons).
# Historical 1800-V and 1900-2250-V executables are NOT altered.
add_executable(mwpc_p3_noFW_penning057_2100_uncapped_v1
               phaseB_native_muon_p3_noFW_penning057_2100_uncapped_v1.cxx)
target_link_libraries(mwpc_p3_noFW_penning057_2100_uncapped_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_p3_noFW_penning057_2100_uncapped_v1 PRIVATE cxx_std_17)

# Isolated 2100 V uncapped microscopic ion-birth reference for future physical
# cathode-strip weighting-current calculations (NOT the strip-current solver).
# Writes per-callback buffered ROOT trees; original HV2100 target is untouched.
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist Graf Gpad Matrix)
add_executable(mwpc_p3_noFW_penning057_cathode_ion_reference_v1
               phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx)
target_link_libraries(mwpc_p3_noFW_penning057_cathode_ion_reference_v1
                      PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(mwpc_p3_noFW_penning057_cathode_ion_reference_v1 PRIVATE cxx_std_17)

add_executable(verify_cathode_ion_root
               ../mwpc-stage-c/cathode-ion-reference/verify_cathode_ion_root.cc)
target_link_libraries(verify_cathode_ion_root
                      PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(verify_cathode_ion_root PRIVATE cxx_std_17)

add_executable(audit_geant4_1000_root
               ../mwpc-stage-c/cathode-ion-reference/audit_geant4_1000_root.cc)
target_link_libraries(audit_geant4_1000_root
                      PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad)
target_compile_features(audit_geant4_1000_root PRIVATE cxx_std_17)

add_executable(compare_sampling_to_full_root
               ../mwpc-stage-c/cathode-ion-reference/compare_sampling_to_full_root.cc)
target_link_libraries(compare_sampling_to_full_root PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(compare_sampling_to_full_root PRIVATE cxx_std_17)

add_executable(select_detailed_1000_root
               ../mwpc-stage-c/cathode-ion-reference/select_detailed_1000_root.cc)
target_link_libraries(select_detailed_1000_root PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(select_detailed_1000_root PRIVATE cxx_std_17)

# Standalone M2a at 2100 V: real P3 9-wire drift + WHOLE-cathode Ramo
# control, independent of old 1800-V prototype; not a stereo-strip model.
add_executable(mwpc_phase_c_p3_noFW_2100_single_ion_control_v1
               phaseC_p3_noFW_2100_single_ion_control_v1.cxx)
target_link_libraries(mwpc_phase_c_p3_noFW_2100_single_ion_control_v1 PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_c_p3_noFW_2100_single_ion_control_v1 PRIVATE cxx_std_17)

# Stage C3 M2b: read immutable microscopic ROOT births, calculate whole-plane
# induced currents for weighted macroions and charge-only FULL archive baseline.
add_executable(mwpc_plane_response_from_ions_root
    ../mwpc-stage-c/cathode-ion-reference/plane_response_from_ions_root.cc)
target_include_directories(mwpc_plane_response_from_ions_root PRIVATE
    ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_plane_response_from_ions_root
    PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)
target_compile_features(mwpc_plane_response_from_ions_root PRIVATE cxx_std_17)


# Experimental fixed-field compact whole-plane kernel:
# train RKF ion responses from archived ion-birth positions, then compare
# raw IDW, low-rank PCA+IDW, and the non-usable projection oracle on
# disjoint held-out ion starts. NO stereo strips, electronics or Fast4 yet.
add_executable(mwpc_compact_plane_kernel_pilot
               ../mwpc-stage-c/cathode-ion-reference/compact_plane_kernel_pilot.cc)
target_include_directories(mwpc_compact_plane_kernel_pilot PRIVATE
               ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_compact_plane_kernel_pilot
               PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree
                       ROOT::Hist ROOT::Matrix)
target_compile_features(mwpc_compact_plane_kernel_pilot PRIVATE cxx_std_17)

# Stage C3 muon response: two WHOLE cathode planes, ion-only.
# Read one sealed microscopic ROOT ion-birth archive and a precomputed
# single-ion response template bank; compare to direct weighted Garfield RKF
# waveforms for the SAME sampled ion cloud. Saves ROOT/CSV/PNG.
add_executable(mwpc_plane_muon_kernel_compare
               ../mwpc-stage-c/cathode-ion-reference/muon_plane_kernel_compare.cc)
target_include_directories(mwpc_plane_muon_kernel_compare PRIVATE
               ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_plane_muon_kernel_compare
               PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree
                       ROOT::Hist ROOT::Graf ROOT::Gpad)
target_compile_features(mwpc_plane_muon_kernel_compare PRIVATE cxx_std_17)


# Electron contribution for SAME microscopic Geant4 muon (independent target).
# Copies unchanged uncapped 2100V Garfield physics selection, but enables
# electron-induced sensor currents for ALL Heed primaries and audits each
# avalanche seed, ion birth and collected electron census against frozen ROOT.
add_executable(mwpc_p3_2100_full_muon_electron_planes
  phaseC_p3_noFW_2100_muon_electron_response_v1.cxx)
target_link_libraries(mwpc_p3_2100_full_muon_electron_planes PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)
target_compile_features(mwpc_p3_2100_full_muon_electron_planes PRIVATE cxx_std_17)

# NO electronics filter: combine audited prompt electrons (0.05 ns) with
# already saved same-muon K512 direct/kernel ion currents at 100 ns,
# produce coarse e+ion ROOT and plots, preserve fine electron ROOT.
add_executable(mwpc_assemble_electron_ion_planes
  ../mwpc-stage-c/cathode-ion-reference/assemble_muon_electron_ion_planes.cc)
target_include_directories(mwpc_assemble_electron_ion_planes PRIVATE
  ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_assemble_electron_ion_planes PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad)
target_compile_features(mwpc_assemble_electron_ion_planes PRIVATE cxx_std_17)


# Stage C time-resolution physics study. The first target recomputes exactly
# the same 256 existing RKF training ion-start sites at finer (10ns) steps.
# The comparison tool REBINS fine signal by charge conservation and checks
# same K512 muon identity before comparing early 100ns/200ns totals.
add_executable(mwpc_rebuild_plane_templates_fine
  ../mwpc-stage-c/cathode-ion-reference/rebuild_plane_templates_fine.cc)
target_include_directories(mwpc_rebuild_plane_templates_fine PRIVATE
  ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_rebuild_plane_templates_fine PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)
target_compile_features(mwpc_rebuild_plane_templates_fine PRIVATE cxx_std_17)

add_executable(mwpc_compare_plane_timing_100ns_10ns
  ../mwpc-stage-c/cathode-ion-reference/compare_plane_timing_100ns_10ns.cc)
target_include_directories(mwpc_compare_plane_timing_100ns_10ns PRIVATE
  ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_compare_plane_timing_100ns_10ns PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad)
target_compile_features(mwpc_compare_plane_timing_100ns_10ns PRIVATE cxx_std_17)


# Read-only three independent K256 genuine muons (phi/jpsi/omega) run with
# unchanged 256 starting-position template RKF response at 10ns.
add_executable(mwpc_audit_independent_muon_kernel_10ns
  ../mwpc-stage-c/cathode-ion-reference/audit_independent_muon_kernel_10ns.cc)
target_include_directories(mwpc_audit_independent_muon_kernel_10ns PRIVATE
  ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_audit_independent_muon_kernel_10ns PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad)
target_compile_features(mwpc_audit_independent_muon_kernel_10ns PRIVATE cxx_std_17)

# Design baseline only: candidate two-sided 1.7mm stereo strips with analytic
# WIRE-FREE Laplace boundary solution. The actual wire-aware 3D strip field
# is still required before physical induced-current calculations.
add_executable(mwpc_stereo_strip_boundary_pilot
  ../mwpc-stage-c/cathode-ion-reference/plot_stereo_strip_boundary_pilot.cc)
target_include_directories(mwpc_stereo_strip_boundary_pilot PRIVATE
  ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_stereo_strip_boundary_pilot PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad)
target_compile_features(mwpc_stereo_strip_boundary_pilot PRIVATE cxx_std_17)


# Stage C proposal-facing 3D U/V stereo-cathode weighting-potential model:
# analytic two-plate strip base plus approximate grounded thin-anode-wire
# 3D line-element collocation, with Garfield analytic WHOLE-cathode
# weighting-potential reference and numerical boundary/convergence tests.
# THIS IS A FIELD STUDY, NOT an induced-current or electronics model.
add_executable(mwpc_wire_aware_stereo_strip_pilot
  ../mwpc-stage-c/cathode-ion-reference/wire_aware_stereo_strip_pilot.cc)
target_include_directories(mwpc_wire_aware_stereo_strip_pilot PRIVATE
  ../mwpc-stage-c/cathode-ion-reference)
target_link_libraries(mwpc_wire_aware_stereo_strip_pilot PRIVATE
  Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
  ROOT::Graf ROOT::Gpad)
target_compile_features(mwpc_wire_aware_stereo_strip_pilot PRIVATE cxx_std_17)
