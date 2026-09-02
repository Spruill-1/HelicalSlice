#include "SpiralVase.hpp"
#include "GCode.hpp"
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>

namespace Slic3r {

namespace SpiralVaseHelpers {
/** == Smooth Spiral Helpers == */
/** Distance between a and b */
float distance(SpiralVase::SpiralPoint a, SpiralVase::SpiralPoint b) { return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2)); }

SpiralVase::SpiralPoint subtract(SpiralVase::SpiralPoint a, SpiralVase::SpiralPoint b)
{
    return SpiralVase::SpiralPoint(a.x - b.x, a.y - b.y);
}

SpiralVase::SpiralPoint add(SpiralVase::SpiralPoint a, SpiralVase::SpiralPoint b) { return SpiralVase::SpiralPoint(a.x + b.x, a.y + b.y); }

SpiralVase::SpiralPoint scale(SpiralVase::SpiralPoint a, float factor) { return SpiralVase::SpiralPoint(a.x * factor, a.y * factor); }

/** dot product */
float dot(SpiralVase::SpiralPoint a, SpiralVase::SpiralPoint b) { return a.x * b.x + a.y * b.y; }

/** Find the point on line ab closes to point c */
SpiralVase::SpiralPoint nearest_point_on_line(SpiralVase::SpiralPoint c, SpiralVase::SpiralPoint a, SpiralVase::SpiralPoint b, float& dist)
{
    SpiralVase::SpiralPoint ab      = subtract(b, a);
    SpiralVase::SpiralPoint ca      = subtract(c, a);
    float                   t       = dot(ca, ab) / dot(ab, ab);
    t                               = t > 1 ? 1 : t;
    t                               = t < 0 ? 0 : t;
    SpiralVase::SpiralPoint closest = SpiralVase::SpiralPoint(add(a, scale(ab, t)));
    dist                            = distance(c, closest);
    return closest;
}

/** Given a set of lines defined by points such as line[n] is the line from points[n] to points[n+1],
 *  find the closest point to p that falls on any of the lines */
SpiralVase::SpiralPoint nearest_point_on_lines(SpiralVase::SpiralPoint               p,
                                               std::vector<SpiralVase::SpiralPoint>* points,
                                               bool&                                 found,
                                               float&                                dist)
{
    if (points->size() < 2) {
        found = false;
        return SpiralVase::SpiralPoint(0, 0);
    }
    float                   min = std::numeric_limits<float>::max();
    SpiralVase::SpiralPoint closest(0, 0);
    for (unsigned long i = 0; i < points->size() - 1; i++) {
        float                   currentDist = 0;
        SpiralVase::SpiralPoint current     = nearest_point_on_line(p, points->at(i), points->at(i + 1), currentDist);
        if (currentDist < min) {
            min     = currentDist;
            closest = current;
            found   = true;
        }
    }
    dist = min;
    return closest;
}

// Uniform spatial grid over a layer's points so smooth-spiral's per-point nearest-segment
// query is O(1) local instead of O(N) linear. For continuous-shell each "layer" fed here is
// the whole wall fill (hundreds of points), so the original O(N^2) per-layer scan is a serial
// bottleneck. Cell size >= the smoothing radius, so a small fixed neighbourhood inspects only
// the segments that can be within range - behaviourally identical to the linear scan for the
// smoothing, which only uses matches within max_xy_smoothing. Each segment is rasterised into
// every cell it crosses, so long connectors are still found exactly.
struct SpiralGrid {
    const std::vector<SpiralVase::SpiralPoint>* P = nullptr;
    float cell = 1.f, inv = 1.f, minx = 0.f, miny = 0.f;
    int   nx = 1, ny = 1;
    std::vector<std::vector<int>> buckets;
    int cell_index(int gx, int gy) const { return gy * nx + gx; }
    void build(const std::vector<SpiralVase::SpiralPoint>* pts, float cell_size) {
        P = pts;
        cell = std::max(cell_size, 1e-3f); inv = 1.f / cell;
        minx = miny = std::numeric_limits<float>::max();
        float maxx = -std::numeric_limits<float>::max(), maxy = -std::numeric_limits<float>::max();
        for (const auto& p : *pts) { minx = std::min(minx, p.x); miny = std::min(miny, p.y); maxx = std::max(maxx, p.x); maxy = std::max(maxy, p.y); }
        nx = std::max(1, int((maxx - minx) * inv) + 1);
        ny = std::max(1, int((maxy - miny) * inv) + 1);
        buckets.assign(size_t(nx) * size_t(ny), {});
        for (int i = 0; i + 1 < int(pts->size()); ++i) {
            const SpiralVase::SpiralPoint a = (*pts)[i], b = (*pts)[i + 1];
            const int steps = std::max(1, int(distance(a, b) * inv) + 1);
            int last = -1;
            for (int s = 0; s <= steps; ++s) {
                const float t  = float(s) / float(steps);
                const int   gx = std::clamp(int((a.x + (b.x - a.x) * t - minx) * inv), 0, nx - 1);
                const int   gy = std::clamp(int((a.y + (b.y - a.y) * t - miny) * inv), 0, ny - 1);
                const int   ci = cell_index(gx, gy);
                if (ci != last) { buckets[ci].push_back(i); last = ci; }
            }
        }
    }
    SpiralVase::SpiralPoint nearest(SpiralVase::SpiralPoint q, bool& found, float& dist) const {
        found = false; float best = std::numeric_limits<float>::max();
        SpiralVase::SpiralPoint bp(0, 0);
        if (P == nullptr || P->size() < 2) { dist = best; return bp; }
        const int gx = std::clamp(int((q.x - minx) * inv), 0, nx - 1);
        const int gy = std::clamp(int((q.y - miny) * inv), 0, ny - 1);
        constexpr int R = 2; // cell >= smoothing radius, so a radius-2 neighbourhood always covers it
        for (int dy = -R; dy <= R; ++dy) for (int dx = -R; dx <= R; ++dx) {
            const int cx = gx + dx, cy = gy + dy;
            if (cx < 0 || cy < 0 || cx >= nx || cy >= ny) continue;
            for (int i : buckets[cell_index(cx, cy)]) {
                float d; SpiralVase::SpiralPoint c = nearest_point_on_line(q, (*P)[i], (*P)[i + 1], d);
                if (d < best) { best = d; bp = c; found = true; }
            }
        }
        dist = best; return bp;
    }
};
} // namespace SpiralVaseHelpers

std::string SpiralVase::process_layer(const std::string &gcode, bool last_layer)
{
    /*  This post-processor relies on several assumptions:
        - all layers are processed through it, including those that are not supposed
          to be transformed, in order to update the reader with the XY positions
        - each call to this method includes a full layer, with a single Z move
          at the beginning
        - each layer is composed by suitable geometry (i.e. a single complete loop)
        - loops were not clipped before calling this method  */
    
    // If we're not going to modify G-code, just feed it to the reader
    // in order to update positions.
    if (! m_enabled) {
        m_reader.parse_buffer(gcode);
        return gcode;
    }
    
    // Get total XY length for this layer by summing all extrusion moves.
    float total_layer_length = 0;
    float layer_height = 0;
    float z = 0.f;
    
    {
        //FIXME Performance warning: This copies the GCodeConfig of the reader.
        GCodeReader r = m_reader;  // clone
        bool set_z = false;
        r.parse_buffer(gcode, [&total_layer_length, &layer_height, &z, &set_z]
            (GCodeReader &reader, const GCodeReader::GCodeLine &line) {
            if (line.cmd_is("G1")) {
                if (line.extruding(reader)) {
                    total_layer_length += line.dist_XY(reader);
                } else if (line.has(Z)) {
                    layer_height += line.dist_Z(reader);
                    if (!set_z) {
                        z = line.new_Z(reader);
                        set_z = true;
                    }
                }
            }
        });
    }

    // Remove layer height from initial Z.
    z -= layer_height;

    std::vector<SpiralVase::SpiralPoint>* current_layer = new std::vector<SpiralVase::SpiralPoint>();
    std::vector<SpiralVase::SpiralPoint>* previous_layer = m_previous_layer;

    bool smooth_spiral = m_smooth_spiral;
    std::string new_gcode;
    std::string transition_gcode;
    float max_xy_dist_for_smoothing = m_max_xy_smoothing;
    //FIXME Tapering of the transition layer only works reliably with relative extruder distances.
    // For absolute extruder distances it will be switched off.
    // Tapering the absolute extruder distances requires to process every extrusion value after the first transition
    // layer.
    bool  transition_in = !m_continuous_shell_mode && m_transition_layer && m_config.use_relative_e_distances.value;
    bool  transition_out = !m_continuous_shell_mode && last_layer && m_config.use_relative_e_distances.value;

    float starting_flowrate  = float(m_config.spiral_starting_flow_ratio.value);
    float finishing_flowrate = float(m_config.spiral_finishing_flow_ratio.value);
    const float min_segment_length = std::max(float(EPSILON), 2 * float(m_config.resolution.value));

    float len = 0.f;
    SpiralVase::SpiralPoint last_point = previous_layer != NULL && previous_layer->size() >0? previous_layer->at(previous_layer->size()-1): SpiralVase::SpiralPoint(0,0);
    const bool continuous_shell_mode = m_continuous_shell_mode;
    // Build a spatial grid over the previous layer once so smooth-spiral's nearest-segment
    // query is O(1) local instead of scanning every previous-layer segment per point (the
    // continuous-shell O(N^2)-per-layer serial bottleneck).
    SpiralVaseHelpers::SpiralGrid prev_grid;
    // Cell size must be >= the smoothing radius for the fixed neighbourhood search to be
    // exhaustive, and is floored at 1 mm so a tiny (user-settable down to 0) smoothing
    // radius cannot blow the bucket allocation up to (bbox/epsilon)^2. A cell larger
    // than the radius is always correct - buckets just hold a few more segments. With a
    // zero radius no match can ever satisfy dist < radius, so skip the build entirely.
    if (smooth_spiral && previous_layer != NULL && previous_layer->size() >= 2 && max_xy_dist_for_smoothing > EPSILON)
        prev_grid.build(previous_layer, std::max(max_xy_dist_for_smoothing, 1.f));
    m_reader.parse_buffer(gcode, [this, &new_gcode, &z, total_layer_length, layer_height, transition_in, &len, &current_layer, &previous_layer, &transition_gcode, transition_out, smooth_spiral, &max_xy_dist_for_smoothing, &last_point, starting_flowrate, finishing_flowrate, min_segment_length, continuous_shell_mode, &prev_grid]
        (GCodeReader &reader, GCodeReader::GCodeLine line) {
        if (continuous_shell_mode &&
            (line.cmd_is("G0") || line.cmd_is("G1") || line.cmd_is("G2") || line.cmd_is("G3")) &&
            !line.extruding(reader) &&
            (line.has(X) || line.has(Y) || line.has(Z) || line.has(E)) &&
            // ...except the approach to the FIRST extrusion of the FIRST spiral layer. There
            // is no shell yet to be continuous with: stripping that travel made the print's
            // opening move extrude one step's E over the whole approach from the preamble
            // position - a 9-20 mm segment at an effective width of 0.001-0.03, the last
            // truly-zero moves the ban audit found. One real travel onto the bed is the
            // normal way to arrive; continuity begins where the shell does.
            !(previous_layer == NULL && len == 0.f)) {
            // A continuous shell cannot contain a travel, Z-hop, wipe, or retraction
            // between adjacent contour segments - drop them. Bare feedrate setters
            // ("G1 F...", no axis words) MUST pass through: they carry all in-layer
            // speed control and the ;_EXTRUDE_SET_SPEED markers the cooling buffer
            // rewrites (the extruding moves themselves carry no F).
            return;
        }
        if (line.cmd_is("G1")) {
            // Orca: Filter out retractions at layer change
            if (line.retracting(reader) || (line.extruding(reader) && line.dist_XY(reader) < min_segment_length)) return;
            if (line.has_z() && !(line.has_x() || line.has_y())) {
                // If this is the initial Z move of the layer, replace it with a
                // (redundant) move to the last Z of previous layer.
                line.set(Z, z);
                new_gcode += line.raw() + '\n';
                return;
            } else {
                float dist_XY = line.dist_XY(reader);
                if (line.has_x() || line.has_y()) { // Sometimes lines have X/Y but the move is to the last position
                    if (dist_XY > 0 && line.extruding(reader)) { // Exclude wipe and retract
                        const bool first_continuous_segment = continuous_shell_mode && len == 0.f;
                        // Nominal flow tracker for the layer-join honesty check. Joins and
                        // transition-scaled lines are excluded so one dishonest segment cannot
                        // teach the tracker its own dishonesty.
                        if (!first_continuous_segment && !transition_in && !transition_out
                            && dist_XY > EPSILON && line.e() > 0.f)
                            m_e_per_mm = m_e_per_mm <= 0.f ? line.e() / dist_XY
                                                           : 0.9f * m_e_per_mm + 0.1f * line.e() / dist_XY;
                        len += dist_XY;
                        float factor = len / total_layer_length;
                        if (transition_in){
                            // Transition layer, interpolate the amount of extrusion starting from spiral_vase_starting_flow_rate to 100%.
                            float starting_e_factor = starting_flowrate + (factor * (1.f - starting_flowrate));
                            line.set(E, line.e() * starting_e_factor, 5 /*decimal_digits*/);
                        } else if (transition_out) {
                            // We want the last layer to ramp down extrusion, but without changing z height!
                            // So clone the line before we mess with its Z and duplicate it into a new layer that ramps down E
                            // We add this new layer at the very end
                            // As with transition_in, the amount is ramped down from 100% to spiral_vase_finishing_flow_rate
                            GCodeReader::GCodeLine transitionLine(line);
                            float finishing_e_factor = finishing_flowrate + ((1.f -factor) * (1.f - finishing_flowrate));
                            transitionLine.set(E, line.e() * finishing_e_factor, 5 /*decimal_digits*/);
                            transition_gcode += transitionLine.raw() + '\n';
                        }
                        // This line is the core of Spiral Vase mode, ramp up the Z smoothly.
                        // In continuous-shell mode the whole-object Z planner has already
                        // written a genuine per-point Z into every extruding move (z_contoured),
                        // so SpiralVase must NOT overwrite it here - it only strips travels and
                        // smooths XY for this mode.
                        if (!continuous_shell_mode)
                            line.set(Z, z + factor * layer_height);
                        if (smooth_spiral) {
                            // Now we also need to try to interpolate X and Y
                            SpiralVase::SpiralPoint p(line.x(), line.y()); // Get current x/y coordinates
                            current_layer->push_back(p);       // Store that point for later use on the next layer
                            if (continuous_shell_mode && first_continuous_segment && previous_layer == NULL) {
                                // The FIRST processed layer's entry. The layer below it never
                                // went through this filter, so the reader has no idea where
                                // the head physically is - last_point is a zero-initialized
                                // lie, and scaling by it would print a 240 mm blob. This is
                                // the boundary where continuity BEGINS, the same boundary the
                                // approach travel below layer 0 crosses: emit the entry as a
                                // real travel, not as one resample step's E smeared over the
                                // gap (9-12 mm at an effective width of 0.004 on every shape).
                                line.set(E, 0.f, 5 /*decimal_digits*/);
                                last_point = p;
                            } else if (continuous_shell_mode && first_continuous_segment && previous_layer != NULL) {
                                // THE LAYER JOIN, MADE HONEST. The upstream writer suppresses
                                // the approach travel, so this segment arrives with its REAL
                                // geometry (the whole gap from wherever the head is, wipe
                                // included) while its E was authored for one resample step -
                                // 9-13 mm at an effective width of 0.004, the ban's last
                                // mid-air class, surviving every earlier fix because no
                                // stripped travel exists in the stream to compare against.
                                // The tell is the E-per-mm: a fraction of what every
                                // surrounding bead carries. Re-issue the E at the running
                                // nominal flow for the segment's true length; the old special
                                // case (retarget to last_point, measure ~0, CLEAR the line)
                                // just handed the gap to the next segment unpaid.
                                // The join line is honest for its PLANNED geometry - the input
                                // stream still holds the approach travel, the reader tracked
                                // it, dist_XY is the tiny planned step. It is the STRIP that
                                // creates the smear: the printer's head stays at the previous
                                // layer's end, so the segment physically spans the whole gap.
                                // Scale the E to that physical length, from the one position
                                // this filter knows to be true - last_point, the previous
                                // layer's final emitted target.
                                const float actual = SpiralVaseHelpers::distance(last_point, p);
                                if (std::getenv("CS_JOINDBG") != nullptr)
                                    fprintf(stderr, "CS_JOINDBG join: e=%.5f dist=%.3f actual=%.3f last=(%.2f,%.2f) p=(%.2f,%.2f)\n",
                                            line.e(), dist_XY, actual, last_point.x, last_point.y, p.x, p.y);
                                if (dist_XY > EPSILON && actual > dist_XY + float(EPSILON))
                                    line.set(E, line.e() * actual / dist_XY, 5 /*decimal_digits*/);
                                last_point = p;
                            } else if (previous_layer != NULL) {
                                bool        found    = false;
                                float       dist     = 0;
                                SpiralVase::SpiralPoint nearestp = prev_grid.nearest(p, found, dist);
                                if (found && dist < max_xy_dist_for_smoothing) {
                                    // Interpolate between the point on this layer and the point on the previous layer
                                    SpiralVase::SpiralPoint target =
                                        SpiralVaseHelpers::add(SpiralVaseHelpers::scale(nearestp, 1 - factor), SpiralVaseHelpers::scale(p, factor));
                                    // Continuous shell: the nearest-point pairing has a whole
                                    // LATTICE to choose from, and where a weaving ring passes a
                                    // neighbour ring of the previous layer this blend dragged
                                    // beads a full gap sideways - measured on the cone as
                                    // full-flow notches up to 0.8 mm deep at the phase-crossing
                                    // azimuths, made entirely here (the assembled path upstream
                                    // is clean). Capping the SEARCH radius instead traded them
                                    // for on/off kinks at the engage boundary, so the cap is on
                                    // the PULL itself: the reference stays (taper-following and
                                    // seam blending keep working), the correction saturates at
                                    // about a fifth of a bead, and the weave's per-layer phase
                                    // precession survives the export instead of being blended
                                    // away. (The layer-join first segment never reaches here
                                    // in continuous-shell mode - it has its own branch above.)
                                    if (continuous_shell_mode) {
                                        constexpr float pull_cap = 0.1f; // mm
                                        const float px = target.x - p.x, py = target.y - p.y;
                                        const float plen = std::sqrt(px * px + py * py);
                                        if (plen > pull_cap)
                                            target = SpiralVase::SpiralPoint(p.x + px * (pull_cap / plen),
                                                                             p.y + py * (pull_cap / plen));
                                    }

                                    // Remove tiny movement
                                    // We need to figure out the distance of this new line!
                                    float modified_dist_XY = SpiralVaseHelpers::distance(last_point, target);
                                    if (modified_dist_XY < min_segment_length)
                                        line.clear();
                                    else {
                                        line.set(X, target.x);
                                        line.set(Y, target.y);
                                        // Scale the extrusion amount according to change in length
                                        line.set(E, line.e() * modified_dist_XY / dist_XY, 5 /*decimal_digits*/);
                                        last_point = target;
                                    }
                                } else {
                                    last_point = p;
                                }
                            }
                        }
                        new_gcode += line.raw() + '\n';
                    }
                    return;
                    /*  Skip travel moves: the move to first perimeter point will
                        cause a visible seam when loops are not aligned in XY; by skipping
                        it we blend the first loop move in the XY plane (although the smoothness
                        of such blend depend on how long the first segment is; maybe we should
                        enforce some minimum length?).
                        When smooth_spiral is enabled, we're gonna end up exactly where the next layer should
                        start anyway, so we don't need the travel move */
                }
            }
        }
        new_gcode += line.raw() + '\n';
        if(transition_out) {
            transition_gcode += line.raw() + '\n';
        }
    });

    delete m_previous_layer;
    m_previous_layer = current_layer;
    
    return new_gcode + transition_gcode;
}

}
