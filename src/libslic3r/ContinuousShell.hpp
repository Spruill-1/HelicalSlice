#pragma once

#include "ExPolygon.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Flow.hpp"
#include "PrintConfig.hpp"
#include "Surface.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace Slic3r {

// Printable-width bounds for the adaptive bead width (C-lite) and the per-angle bead
// plan (C-full). Zero fields derive Arachne-consistent defaults from the flow's nozzle
// diameter (85% / 25% / 100% / 25%). The caller (LayerRegion) fills these from the
// object-scope min_bead_width / min_feature_size / wall_transition_* percent options so
// continuous shell honours the same limits as the Arachne generator.
struct ContinuousShellWidthLimits {
    double min_bead_width_mm   = 0.;   // narrowest printable bead
    double min_feature_size_mm = 0.;   // walls thinner than this are declined outright
    // C-full (continuous_shell_adaptive_beading) additions:
    double wall_transition_length_mm           = 0.;  // arc length of a ring-count transition ramp
    double wall_transition_filter_deviation_mm = 0.;  // window keep/drop hysteresis
};

class ContinuousShellPathBuilder
{
public:
    using WidthLimits = ContinuousShellWidthLimits;

    static bool build_cross_section(const Surfaces &surfaces, const Flow &flow, double layer_height, const PrintRegionConfig *config,
                                    ExtrusionEntityCollection &entities, std::string *reason = nullptr,
                                    bool reverse_path = false, size_t layer_index = 0, bool force_solid = false,
                                    const ContinuousShellWidthLimits &width_limits = {},
                                    // Set when the layer BUILT but took a weaker path than asked for
                                    // (adaptive beading skipped or declined, interior spiral declined).
                                    // Distinct from `reason`, which reports an outright failure.
                                    std::string *degrade_note = nullptr);

    // Collect the leaf ExtrusionPaths of a continuous-shell entity tree in traversal
    // order. Used by the whole-object Z planner (PrintObject::continuous_shell_z) at
    // the posContouring stage to reach the per-layer bead centreline.
    static void collect_paths(ExtrusionEntity *entity, std::vector<ExtrusionPath *> &out);

    // Assign a continuous per-point Z ramp across the concatenated leaf paths of one
    // layer's continuous-shell bead and mark them z_contoured, so the exporter emits a
    // genuine per-point non-planar Z instead of relying on the SpiralVase post-process.
    // z_diff(point) = -(1 - factor) * layer_height, factor = cumulative XY arc fraction.
    // min_z_rel (mm, relative to the layer's nominal Z) is a physical floor: no point,
    // including pre-loaded weave/interlock dips, is allowed below it.
    //
    // NOTE: a DEPTH-based profile (surface low, core high) was tried for uniform
    // per-column density but empirically breaks the reverse_path join continuity - it
    // puts a full-layer vertical jump at every layer seam - so we keep the arc-length
    // ramp (monotone-in-path).
    static void apply_continuous_z(ExtrusionEntity *entity, double layer_height, double min_z_rel = -1e9);
};

} // namespace Slic3r
