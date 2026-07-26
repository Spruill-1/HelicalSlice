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

class ContinuousShellPathBuilder
{
public:
    static bool build_cross_section(const Surfaces &surfaces, const Flow &flow, double layer_height, const PrintRegionConfig *config,
                                    ExtrusionEntityCollection &entities, std::string *reason = nullptr,
                                    bool reverse_path = false, size_t layer_index = 0, bool force_solid = false,
                                    // Interior islands of this layer that sit over AIR (a blind bore's
                                    // cap): reintroduced as holes for the ring march (the rings wrap
                                    // them like the bore below) and spanned by a serpentine BRIDGE
                                    // patch appended to the tour, so the next layer's rings close
                                    // over taut lines instead of sagging into the opening.
                                    const ExPolygons *bridge_over = nullptr);

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
