#include <catch2/catch_all.hpp>

#include "libslic3r/ContinuousShell.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "data_cs_gear_layer.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <vector>
#include <array>
#include <map>

using namespace Slic3r;

// Scoped env toggle for the CS_BUBBLE prototype knobs, so a failed REQUIRE can't leak
// bubble mode into other tests in the same process.
struct EnvFlag {
    const char *name;
    EnvFlag(const char *n, const char *v) : name(n) { set(v); }
    ~EnvFlag() { set(nullptr); }
    void set(const char *v) {
#ifdef _WIN32
        _putenv_s(name, v != nullptr ? v : "");
#else
        v != nullptr ? setenv(name, v, 1) : unsetenv(name);
#endif
    }
};

TEST_CASE("Continuous shell keeps a deep side bite (socket) hollow", "[ContinuousShell]") {
    // Regression for the side-socket fill-in (HANDOFF 11.2b): a circle with a deep
    // hemispherical bite - a valid single-region "pac-man" cross-section - must not
    // get bead laid inside the bite cavity. Circle r12mm at origin, bite r6mm
    // centred on the wall at (12, 0): cavity reaches to the part centre.
    Polygon circle, bite;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circle.points.emplace_back(coord_t(12000000 * std::cos(a)), coord_t(12000000 * std::sin(a)));
        bite.points.emplace_back(coord_t(12000000 + 6000000 * std::cos(a)), coord_t(6000000 * std::sin(a)));
    }
    ExPolygons pacman = diff_ex(Polygons{circle}, Polygons{bite});
    REQUIRE(pacman.size() == 1);
    REQUIRE(pacman.front().holes.empty());

    Surfaces surfaces{ Surface(stInternal, pacman.front()) };
    Flow flow(0.5f, 0.25f, 0.5f);
    ExtrusionEntityCollection entities;
    std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, entities, &reason, false, 4);
    INFO(reason);
    REQUIRE(built);

    // Probe: no bead point (or segment midpoint) may land deep inside the cavity -
    // within 4.5mm of the bite centre AND at least 1mm inside the outer circle
    // (margins exclude legitimate wall beads hugging the bite surface).
    const Point bite_c(12000000, 0);
    auto in_cavity = [&](const Point &p) {
        const double to_bite   = (p - bite_c).cast<double>().norm();
        const double to_centre = p.cast<double>().norm();
        return to_bite < 4500000.0 && to_centre < 11000000.0;
    };
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *e : entities.entities)
        ContinuousShellPathBuilder::collect_paths(e, paths);
    REQUIRE_FALSE(paths.empty());
    size_t cavity_pts = 0, cavity_mids = 0, total = 0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 0; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i].x(), p->polyline.points[i].y());
            ++total;
            if (in_cavity(a)) ++cavity_pts;
            if (i > 0) {
                const Point b(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                if (in_cavity(Point((a + b) / 2))) ++cavity_mids;
            }
        }
    INFO("total points: " << total);
    CHECK(total > 200);
    REQUIRE(cavity_pts == 0);
    REQUIRE(cavity_mids == 0);
}

TEST_CASE("Continuous shell routes split-chain connectors through material (dogbone)", "[ContinuousShell]") {
    // Regression for the dogbone notch chords (HANDOFF 11.2): when the inward offset
    // splits at a waist, the sub-chains can end deep inside a lobe and the old greedy
    // stitch hopped straight across the concave notch through air. Connectors must now
    // route through material. Dogbone: 30x16mm, 6mm-deep notches top and bottom.
    Polygon dog{ {-15000000, -8000000}, {-5000000, -8000000}, {-5000000, -2000000},
                 {5000000, -2000000}, {5000000, -8000000}, {15000000, -8000000},
                 {15000000, 8000000}, {5000000, 8000000}, {5000000, 2000000},
                 {-5000000, 2000000}, {-5000000, 8000000}, {-15000000, 8000000} };
    Surfaces surfaces{ Surface(stInternal, ExPolygon(dog)) };
    // Flow matched to the real 0.4mm-nozzle slicing profile: band/chain topology is
    // flow-dependent, and a nicer round flow hid the residual chord family entirely.
    Flow flow(0.42f, 0.2f, 0.4f);

    // The notches, shrunk 1mm on each side so boundary-hugging wall beads don't count.
    auto in_notch = [](const Point &p) {
        const bool x_in = p.x() > -4000000 && p.x() < 4000000;
        const bool top  = p.y() >  3000000 && p.y() < 7000000;
        const bool bot  = p.y() < -3000000 && p.y() > -7000000;
        return x_in && (top || bot);
    };
    // Seam precession and reverse_path parity move the chain endpoints every layer;
    // the connector routing must hold for all of them, not one lucky configuration.
    size_t worst_layer = 0, worst_mids = 0;
    for (size_t layer = 0; layer < 110; ++layer) {
        ExtrusionEntityCollection entities;
        std::string reason;
        const bool built = ContinuousShellPathBuilder::build_cross_section(
            surfaces, flow, 0.25, nullptr, entities, &reason, (layer & 1) != 0, layer);
        INFO("layer " << layer << ": " << reason);
        REQUIRE(built);
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *e : entities.entities)
            ContinuousShellPathBuilder::collect_paths(e, paths);
        REQUIRE_FALSE(paths.empty());
        size_t notch_mids = 0;
        for (const ExtrusionPath *p : paths)
            for (size_t i = 1; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
                if (in_notch(Point((a + b) / 2)))
                    ++notch_mids;
            }
        if (notch_mids > worst_mids) { worst_mids = notch_mids; worst_layer = layer; }
    }
    INFO("worst layer " << worst_layer << " with " << worst_mids << " notch-crossing segments");
    REQUIRE(worst_mids == 0);
}

TEST_CASE("Continuous shell fills an annular wall solid as an in/out spiral", "[ContinuousShell]") {
    // 8mm outer square with a 3..5mm square bore -> a 3mm-thick wall, i.e. several
    // bead passes on each side of the wall's medial axis.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{annulus};
    Flow flow(0.5f, 0.25f, 0.5f);

    ExtrusionEntityCollection entities;
    std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, entities, &reason, false);
    INFO(reason);
    REQUIRE(built);
    REQUIRE(entities.entities.size() == 1);

    auto *path = dynamic_cast<ExtrusionMultiPath *>(entities.entities.front());
    REQUIRE(path != nullptr);
    // The wall is filled by a continuous morph spiral, split into surface / interior
    // runs. Assert rich coverage across all runs.
    size_t total_pts = 0;
    for (const ExtrusionPath &p : path->paths)
        total_pts += p.polyline.size();
    REQUIRE(total_pts > 200);
}

TEST_CASE("Continuous shell annular spiral keeps the bore hollow (no morph bowtie)", "[ContinuousShell]") {
    // Regression: when consecutive rings wound in opposite directions, the morph
    // mapped each point to its antipode and the path cut straight across the bore,
    // producing a leaf/eye artifact through the centre. Assert nothing comes close
    // to the centre of the hollow bore.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{annulus};
    Flow flow(0.5f, 0.25f, 0.5f);

    ExtrusionEntityCollection entities;
    std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, entities, &reason, false));
    auto *path = dynamic_cast<ExtrusionMultiPath *>(entities.entities.front());
    REQUIRE(path != nullptr);

    const Point center(4000000, 4000000);
    double min_r = std::numeric_limits<double>::max();
    for (const ExtrusionPath &p : path->paths)
        for (const auto &pt : p.polyline.points) {
            const double dx = double(pt.x()) - double(center.x());
            const double dy = double(pt.y()) - double(center.y());
            min_r = std::min(min_r, std::sqrt(dx * dx + dy * dy));
        }
    // Innermost printed ring sits ~1mm from the centre; a bowtie would drag points
    // to ~0. Require a clear hollow margin.
    REQUIRE(min_r > 400000.0);
}

TEST_CASE("Continuous shell precesses the seam start point per layer", "[ContinuousShell]") {
    // The once-per-layer seam must walk around the part instead of stacking into a
    // vertical line. Two layers a few indices apart should start at clearly
    // different points on the outer contour.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{annulus};
    Flow flow(0.5f, 0.25f, 0.5f);

    ExtrusionEntityCollection e0, e8;
    std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, e0, &reason, false, 0));
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, e8, &reason, false, 8));

    const Point p0 = e0.entities.front()->first_point();
    const Point p8 = e8.entities.front()->first_point();
    // A fixed angular walk (~1.15 deg/layer) separates the two seam starts.
    REQUIRE((p8 - p0).cast<double>().norm() > double(flow.scaled_spacing()));
}

TEST_CASE("Continuous shell connects consecutive layers seamlessly (no per-layer chord)", "[ContinuousShell]") {
    // The defect the 3D-helix pivot fixes: consecutive layers used to end/begin ~80
    // degrees (~9mm) apart, so SpiralVase bridged them with a fat extruding chord
    // across the bore every layer. With absolute-angle seam anchoring plus the
    // in/out (reverse) alternation, layer N must END essentially where layer N+1
    // BEGINS (same boundary, at most ~one seam-step apart).
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{annulus};
    Flow flow(0.5f, 0.25f, 0.5f);

    ExtrusionEntityCollection e10, e11;
    std::string reason;
    // Mirror how LayerRegion drives it: reverse = (layer & 1), layer_index = layer.
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, e10, &reason, false, 10));
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, e11, &reason, true, 11));

    const Point end10   = e10.entities.front()->last_point();
    const Point start11 = e11.entities.front()->first_point();
    const double gap = (start11 - end10).cast<double>().norm();
    INFO("layer join gap (scaled units) = " << gap);
    // Under 2mm; the old relative-chaining bug produced ~9mm here.
    REQUIRE(gap < 2000000.0);
}

TEST_CASE("Continuous shell: solid alternates per layer, annulus stays stable", "[ContinuousShell]") {
    Flow flow(0.5f, 0.25f, 0.5f);
    std::string reason;

    // A SOLID simple-spiral has its two ends at DIFFERENT radii (outer surface vs buried
    // core), so it alternates direction per layer to keep consecutive layers' joins
    // adjacent: forward starts on the outer edge, reversed starts near the core.
    Surface solid(stTop, ExPolygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}});
    Surfaces solids{solid};
    const Point center(4000000, 4000000);
    ExtrusionEntityCollection s_fwd, s_rev;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(solids, flow, 0.25, nullptr, s_fwd, &reason, false));
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(solids, flow, 0.25, nullptr, s_rev, &reason, true));
    const double s_fwd_r = (s_fwd.entities.front()->first_point() - center).cast<double>().norm();
    const double s_rev_r = (s_rev.entities.front()->first_point() - center).cast<double>().norm();
    REQUIRE(s_fwd_r > s_rev_r);

    // An ANNULUS Fermat has BOTH path ends near the outer surface, so it already joins
    // layer-to-layer with a tiny hop and does NOT need per-layer reversal. Reversing it is
    // what makes each ring's arc-position (hence its Z under the arc-length ramp) alternate
    // every layer, pairing beads into a once-per-two-layers vertical void ("reverse-join
    // doubling"). So the builder leaves annuli unreversed - forward and reversed start
    // identically.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces annuli{annulus};
    ExtrusionEntityCollection a_fwd, a_rev;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(annuli, flow, 0.25, nullptr, a_fwd, &reason, false));
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(annuli, flow, 0.25, nullptr, a_rev, &reason, true));
    REQUIRE(a_fwd.entities.front()->first_point() == a_rev.entities.front()->first_point());
}

TEST_CASE("Continuous shell seam jitter scatters the seam yet stays join-safe", "[ContinuousShell]") {
    // Strong seam jitter must scatter the seam WITHOUT breaking the continuous bead: the
    // per-layer step stays small (short inter-layer joins) even though the accumulated walk
    // wanders the seam around the part. Regression guard for the accumulated-walk design.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);

    PrintRegionConfig cfg;
    cfg.continuous_shell_seam_step.value   = 100;
    cfg.continuous_shell_seam_jitter.value = 200; // strong scatter

    double        max_join = 0.0;
    Point         prev_last;
    std::vector<double> start_x;
    for (size_t layer = 0; layer < 80; ++layer) {
        ExtrusionEntityCollection e;
        std::string               reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason,
                                                                (layer & 1u) != 0, layer));
        const Point start = e.entities.front()->first_point();
        const Point last  = e.entities.front()->last_point();
        if (layer > 0)
            max_join = std::max(max_join, (start - prev_last).cast<double>().norm());
        prev_last = last;
        start_x.push_back(double(start.x()));
    }
    INFO("max inter-layer join (scaled) = " << max_join);
    // Continuous bead preserved: consecutive layers still meet with a short hop.
    REQUIRE(max_join < 3000000.0); // < 3 mm
    // The seam actually wanders (jitter has an effect).
    const double spread = *std::max_element(start_x.begin(), start_x.end()) -
                          *std::min_element(start_x.begin(), start_x.end());
    REQUIRE(spread > 1000000.0); // seam start moves > 1 mm across the print
}

TEST_CASE("Continuous shell cross-section builder fills a solid cap", "[ContinuousShell]") {
    Surface cap(stTop, ExPolygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}});
    Surfaces surfaces{cap};
    Flow flow(0.5f, 0.25f, 0.5f);
    ExtrusionEntityCollection entities;
    std::string reason;

    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, entities, &reason, false);
    INFO(reason);
    REQUIRE(built);
    REQUIRE(entities.entities.size() == 1);

    auto *path = dynamic_cast<ExtrusionMultiPath *>(entities.entities.front());
    REQUIRE(path != nullptr);
    // Solid disk spirals inward to the centre; assert rich coverage across all runs.
    size_t total_pts = 0;
    for (const ExtrusionPath &p : path->paths)
        total_pts += p.polyline.size();
    REQUIRE(total_pts > 100);
}

TEST_CASE("Continuous shell over-extrusion is capped to the residual void (no surface bulge)", "[ContinuousShell]") {
    // Regression: a flat per-bead over-extrusion % dumped far more material than the cross
    // section could hold, so the excess bulged out to the surface and wrecked the finish. The
    // knob now fills only a fraction of the residual VOID (V_model - V_nominal), so the total
    // extruded can never exceed the solid model volume and no bead is boosted past a hard cap.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);
    const double layer_height = 0.25;
    const double base = flow.mm3_per_mm();
    const double v_model = std::abs(annulus.expolygon.area()) * SCALING_FACTOR * SCALING_FACTOR * layer_height;

    // Total extruded volume (mm^3) and the largest per-bead flow multiple, for a given knob.
    auto measure = [&](double over_pct, double &max_mult) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_internal_overextrusion.value = over_pct;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, layer_height, &cfg, e, &reason, false));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        double vol = 0.0; max_mult = 0.0;
        for (const ExtrusionPath &p : mp->paths) {
            vol += p.mm3_per_mm * p.polyline.length() * SCALING_FACTOR;
            max_mult = std::max(max_mult, p.mm3_per_mm / base);
        }
        return vol;
    };

    double m0 = 0.0, m100 = 0.0;
    const double v_off = measure(0.0,   m0);    // nominal baseline
    const double v_on  = measure(100.0, m100);  // fill 100% of the residual void

    REQUIRE(v_off > 0.0);
    // 0% over-extrusion: beads at nominal flow EXCEPT the always-on wrap-window seam
    // fill, which adds each ring's ramp share (bounded well below one extra bead) to
    // the few points at the ring's start meridian.
    REQUIRE(m0 > 0.999);  REQUIRE(m0 < 1.5);
    REQUIRE(v_on >= v_off - 1e-9);                    // over-extrusion only ever ADDS material
    REQUIRE(v_on <= std::max(v_off, v_model) * 1.02); // and never past the solid model volume -> no bulge
    REQUIRE(m100 <= 1.5 * 1.35 + 1e-6);               // per-bead cap x possible wrap fill share
}

TEST_CASE("Continuous shell interlock staggers buried rings within a bounded Z envelope", "[ContinuousShell]") {
    // Brick interlock nudges buried wall rings up/down in Z (half a bead at 100%) so
    // radially-adjacent beads stop sharing a Z level, WITHOUT breaking the continuous bead
    // or moving a visible surface. It must: (a) actually offset buried beads by ~a quarter
    // bead, (b) stay within +/- one quarter bead of the smooth ramp (so the print envelope
    // and the exporter Z floor are unchanged), (c) leave the outer surface on the ramp, and
    // (d) keep consecutive points a small continuous Z step apart (no jump / travel).
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);
    const double layer_height = 0.25;

    auto build_z = [&](double interlock_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_interlock.value = interlock_pct;
        ExtrusionEntityCollection e;
        std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, layer_height, &cfg, e, &reason, false, 4));
        ContinuousShellPathBuilder::apply_continuous_z(&e, layer_height);
        std::vector<Point3> pts;
        for (auto *ent : e.entities)
            if (auto *mp = dynamic_cast<ExtrusionMultiPath *>(ent))
                for (auto &p : mp->paths)
                    for (auto &q : p.polyline.points)
                        pts.push_back(q);
        return pts;
    };

    // Interlock only sets Z, so both builds share XY and point ordering; Z can be diffed
    // point-for-point to isolate the stagger from the smooth ramp.
    const std::vector<Point3> off = build_z(0.0);
    const std::vector<Point3> on  = build_z(100.0);
    REQUIRE(off.size() == on.size());
    REQUIRE(off.size() > 200);

    const coord_t amp = coord_t(std::llround(scale_(layer_height / 4.0))); // quarter bead at 100%
    const Point   center(4000000, 4000000);

    coord_t max_dev = 0;
    double  outer_dev = 0.0, max_r = 0.0;
    for (size_t i = 0; i < on.size(); ++i) {
        REQUIRE(on[i].x() == off[i].x());   // XY untouched by the interlock
        REQUIRE(on[i].y() == off[i].y());
        const coord_t dev = coord_t(std::abs(on[i].z() - off[i].z()));
        max_dev = std::max(max_dev, dev);
        const double r = (Point(on[i].x(), on[i].y()) - center).cast<double>().norm();
        if (r > max_r) { max_r = r; outer_dev = double(dev); }
    }
    REQUIRE(max_dev > amp * 3 / 4);        // (a) the stagger is real
    REQUIRE(max_dev <= amp + 2);           // (b) bounded to a quarter bead (+rounding)
    REQUIRE(outer_dev < double(amp) / 4.); // (c) outermost surface bead pinned to the ramp

    coord_t max_step = 0;                  // (d) continuity: no Z jump between adjacent points
    for (size_t i = 1; i < on.size(); ++i)
        max_step = std::max(max_step, coord_t(std::abs(on[i].z() - on[i - 1].z())));
    REQUIRE(max_step <= 2 * amp + coord_t(scale_(0.02)));
}

TEST_CASE("Continuous shell connector spread disperses the seam handoffs (A2)", "[ContinuousShell]") {
    // The once-per-layer surface<->interior handoffs normally stack at ~one angle, forming a
    // deep radial seam groove that winds up the part. connector_spread must fan them around
    // the circumference so they disperse - while keeping the bore hollow and the path valid.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);
    const Point center(4000000, 4000000);
    const double jump = 0.5 * double(flow.scaled_spacing());

    // Circular concentration R of the radial-jump (connector) angles: ~1 = all stacked at one
    // angle, ~0 = spread around the circle. Aggregated over a few layers for a stable sample.
    auto concentration = [&](double spread_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_connector_spread.value = spread_pct;
        double sx = 0.0, sy = 0.0; size_t n = 0;
        for (size_t layer = 0; layer < 4; ++layer) {
            ExtrusionEntityCollection e; std::string reason;
            REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, layer));
            auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
            REQUIRE(mp != nullptr);
            std::vector<Point3> pts;
            for (auto &p : mp->paths) for (auto &q : p.polyline.points) pts.push_back(q);
            double minr = std::numeric_limits<double>::max();
            for (auto &q : pts) minr = std::min(minr, (Point(q.x(), q.y()) - center).cast<double>().norm());
            REQUIRE(minr > 400000.0);  // bore stays hollow regardless of spread
            for (size_t i = 1; i < pts.size(); ++i) {
                const double r0 = (Point(pts[i-1].x(), pts[i-1].y()) - center).cast<double>().norm();
                const double r1 = (Point(pts[i].x(),   pts[i].y())   - center).cast<double>().norm();
                if (std::abs(r1 - r0) > jump) {
                    const double mx = 0.5*(double(pts[i-1].x())+double(pts[i].x())) - double(center.x());
                    const double my = 0.5*(double(pts[i-1].y())+double(pts[i].y())) - double(center.y());
                    const double a = std::atan2(my, mx);
                    sx += std::cos(a); sy += std::sin(a); ++n;
                }
            }
        }
        REQUIRE(n > 8);
        return std::sqrt(sx*sx + sy*sy) / double(n);
    };

    const double r_tight  = concentration(0.0);    // base twist only -> handoffs concentrated
    const double r_spread = concentration(100.0);  // full-revolution spread -> handoffs dispersed
    INFO("connector-angle concentration: tight=" << r_tight << " spread=" << r_spread);
    REQUIRE(r_spread < r_tight);        // spreading actually disperses the handoffs

    // Regression: at high spread the inner-ring hops used to chord straight across the hollow
    // bore (stringing in the void). Count segments that cross the bore centre; the radius-aware
    // anchor clamp must keep every hop out of the bore.
    auto bore_chords = [&](double spread_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_connector_spread.value = spread_pct;
        size_t chords = 0;
        for (size_t layer = 0; layer < 4; ++layer) {
            ExtrusionEntityCollection e; std::string reason;
            REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, layer));
            auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
            REQUIRE(mp != nullptr);
            std::vector<Point3> pts;
            for (auto &p : mp->paths) for (auto &q : p.polyline.points) pts.push_back(q);
            for (size_t i = 1; i < pts.size(); ++i) {
                const Point a(pts[i-1].x(), pts[i-1].y()), b(pts[i].x(), pts[i].y());
                const double ra = (a - center).cast<double>().norm();
                const double rb = (b - center).cast<double>().norm();
                const double len = (b - a).cast<double>().norm();
                const Point mid((a.x() + b.x()) / 2, (a.y() + b.y()) / 2);
                const double rm = (mid - center).cast<double>().norm();
                // both ends well out in the wall, a long hop, midpoint near the centre = bore chord
                if (ra > 1500000.0 && rb > 1500000.0 && len > 2500000.0 && rm < 900000.0)
                    ++chords;
            }
        }
        return chords;
    };
    REQUIRE(bore_chords(100.0) == 0u);
}

TEST_CASE("Continuous shell connector spread is a no-op on a solid cross-section (A2 gate)", "[ContinuousShell]") {
    // A2 connector_spread fans the surface<->bore handoffs of a TUBE wall to break up the
    // deep spiral seam groove. On a SOLID cross-section there is no bore seam to spread;
    // fanning the inner connectors instead swings them across the shape's OWN concavities
    // (the notch of an L/'+'), laying stringing chords through empty space outside the part.
    // Spread is therefore gated to annular (has_bore) chains, so it must have NO effect on a
    // solid fill. Guard the gate: an L-shape's path must be IDENTICAL at spread 0 and 100.
    Surface l_shape(
        stInternal,
        ExPolygon(Polygon{{0, 0}, {10000000, 0}, {10000000, 4000000},
                          {4000000, 4000000}, {4000000, 10000000}, {0, 10000000}}));
    Surfaces surfaces{ l_shape };
    Flow flow(0.5f, 0.25f, 0.5f);

    auto path_points = [&](double spread_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_connector_spread.value = spread_pct;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::vector<Point3> pts;
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) pts.push_back(q);
        return pts;
    };

    const std::vector<Point3> tight  = path_points(0.0);
    const std::vector<Point3> spread = path_points(100.0);
    REQUIRE(tight.size() > 8);
    REQUIRE(spread.size() == tight.size());
    for (size_t i = 0; i < tight.size(); ++i) {
        REQUIRE(spread[i].x() == tight[i].x());
        REQUIRE(spread[i].y() == tight[i].y());
    }
}

TEST_CASE("Continuous shell wrap-window fill boosts flow at each ring's start meridian", "[ContinuousShell]") {
    // The arc-length ramp leaves a once-per-wrap tall slot at each ring's start (the ring
    // climbs only its share of a layer per revolution), so the points laid within the
    // first seam-step after a ring's start carry flow * (1 + ring_share). Verify on an
    // annulus: some sub-paths get a boost bounded by the per-ring share, positioned at
    // the START of the tour (the outer ring's window is the bead's first sub-path).
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);
    const double base = flow.mm3_per_mm();

    PrintRegionConfig cfg;
    cfg.continuous_shell_seam_step.value = 400;   // wide window -> several boosted points
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
    auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
    REQUIRE(mp != nullptr);
    REQUIRE(mp->paths.size() >= 2);
    // The tour starts inside the outer ring's wrap window -> boosted first sub-path.
    REQUIRE(mp->paths.front().mm3_per_mm > base * 1.02);
    // Boost magnitude is the ring share: strictly less than a full extra bead.
    double mx = 0.0, mn = std::numeric_limits<double>::max();
    size_t boosted = 0;
    for (auto &p : mp->paths) {
        mx = std::max(mx, p.mm3_per_mm);
        mn = std::min(mn, p.mm3_per_mm);
        if (p.mm3_per_mm > base * 1.02) ++boosted;
    }
    REQUIRE(mx < base * 1.6);
    REQUIRE(mn > base * 0.999);           // fill only ADDS; nothing is thinned
    REQUIRE(boosted >= 2);                // several rings' windows are present
}

TEST_CASE("Continuous shell seam coast trims both sides of the join, decaying with travel", "[ContinuousShell]") {
    // Symmetric seam_coast (pre + post-join pressure trim): thins the bead approaching
    // AND leaving the layer join (the two ends of the tour), decaying exponentially over
    // ~15mm of travel each; beads far from BOTH ends stay at full flow.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);
    const double base = flow.mm3_per_mm();

    auto build_paths = [&](double coast_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_seam_coast.value = coast_pct;
        cfg.continuous_shell_seam_step.value  = 0;   // isolate the trim from the wrap fill window
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        return mp->paths;
    };

    const auto off = build_paths(0.0);
    for (const auto &p : off) {
        REQUIRE(p.mm3_per_mm > base * 0.999);        // trim off -> nothing thinned
        REQUIRE(p.mm3_per_mm < base * 1.35);         // (wrap fill may add its bounded boost)
    }

    const auto on = build_paths(40.0);
    // Both ends (approaching and leaving the join) are thinned toward (1 - 40%).
    double mn = std::numeric_limits<double>::max();
    for (const auto &p : on) mn = std::min(mn, p.mm3_per_mm);
    REQUIRE(mn < base * 0.75);
    REQUIRE(mn > base * 0.55);
    // The first AND last sub-paths (nearest the join on each side) are the trimmed ones.
    REQUIRE(on.front().mm3_per_mm < base * 0.85);
    REQUIRE(on.back().mm3_per_mm  < base * 0.85);
    // The trim DECAYS: a middle sub-path (far from both ends) is back near full flow.
    REQUIRE(on[on.size() / 2].mm3_per_mm > base * 0.95);
}

TEST_CASE("Continuous shell woven infill sparsifies the interior and weaves it in Z", "[ContinuousShell]") {
    // A SOLID square: infill density 0 fills it fully; density 20% keeps a solid shell (wall_loops)
    // and thins the interior to sparse loops, so the total bead is much shorter - all still one
    // continuous multipath. infill weave 100 undulates the buried infill in Z (the gyroid-like
    // vertical weave), so buried points span a range of Z offsets; weave 0 leaves them flat.
    Surface solid(stInternal, ExPolygon(Polygon{{0,0},{16000000,0},{16000000,16000000},{0,16000000}}));
    Surfaces surfaces{ solid };
    Flow flow(0.5f, 0.25f, 0.5f);

    struct R { double xy_len; double z_span; double max_xy; };
    auto measure = [&](int density, int weave) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = density;
        cfg.continuous_shell_infill_weave.value   = weave;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
        REQUIRE(e.entities.size() == 1);              // one continuous entity (single bead)
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::vector<Point3> pts;
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) pts.push_back(q);
        REQUIRE(pts.size() > 8);
        double xy = 0.0, zmin = std::numeric_limits<double>::max(), zmax = -std::numeric_limits<double>::max();
        double max_xy = 0.0;
        for (size_t i = 0; i < pts.size(); ++i) {
            if (i > 0)
                xy += (Point(pts[i].x(), pts[i].y()) - Point(pts[i-1].x(), pts[i-1].y())).cast<double>().norm();
            zmin = std::min(zmin, double(pts[i].z()));
            zmax = std::max(zmax, double(pts[i].z()));
            max_xy = std::max(max_xy, std::max(std::abs(double(pts[i].x())), std::abs(double(pts[i].y()))));
        }
        return R{ xy * SCALING_FACTOR, (zmax - zmin) * SCALING_FACTOR, max_xy * SCALING_FACTOR };
    };

    const R solid_fill = measure(0,   0);    // fully solid (feature disabled)
    const R sparse     = measure(20,  0);    // sparse, no weave
    const R woven      = measure(20, 100);   // sparse, full 3D weave

    // Sparsifying the interior removes most of the bead length.
    REQUIRE(sparse.xy_len < 0.6 * solid_fill.xy_len);
    // Without weave the pre-ramp Z is flat (only the ramp, added later, tilts the bead).
    REQUIRE(sparse.z_span < 0.02);
    // With weave the buried infill undulates in Z by a good fraction of the 0.25mm layer height.
    REQUIRE(woven.z_span  > 0.05);
    // Weaving wavies the rings so neighbours touch (adds some in-plane length) but must keep the
    // fill sparse relative to solid - it repositions/knits beads, it does not densify to 100%.
    REQUIRE(woven.xy_len < 0.75 * solid_fill.xy_len);
    REQUIRE(woven.xy_len > 0.9  * sparse.xy_len);
    // The weave must stay INSIDE the object: no lobe may bulge past the 16mm square boundary
    // (regression - the touching-lobe amplitude once pushed the outer ring out past the wall).
    REQUIRE(woven.max_xy < 16.05);
}

TEST_CASE("Continuous shell adaptive stride sparsifies cells too small for the full stride", "[ContinuousShell]") {
    // A compact solid piece (8mm square, ~9 rings) cannot host the configured density-20
    // (stride-5) lattice - that needs >12 rings beyond the shells - and used to print fully
    // SOLID under the small-cell rule (die-5 v7: the whole rim mosaic came out dense). With
    // per-cell adaptive stride it steps DOWN to the largest stride it CAN host and still
    // sparsifies; only truly tiny pieces stay dense. The 16mm square in the sibling test
    // hosts the full stride, so this size is the discriminating case.
    Surface small(stInternal, ExPolygon(Polygon{{0,0},{8000000,0},{8000000,8000000},{0,8000000}}));
    Surfaces surfaces{ small };
    Flow flow(0.5f, 0.25f, 0.5f);
    auto xy_len = [&](int density, int weave) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = density;
        cfg.continuous_shell_infill_weave.value   = weave;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
        REQUIRE(e.entities.size() == 1);              // still one continuous bead
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        double xy = 0.0;
        for (auto &p : mp->paths)
            for (size_t i = 1; i < p.polyline.points.size(); ++i)
                xy += (Point(p.polyline.points[i].x(), p.polyline.points[i].y()) -
                       Point(p.polyline.points[i-1].x(), p.polyline.points[i-1].y())).cast<double>().norm();
        return xy * SCALING_FACTOR;
    };
    const double solid_len  = xy_len(0,  0);
    const double sparse_len = xy_len(20, 100);
    REQUIRE(sparse_len < 0.85 * solid_len);           // stepped-down lattice, no longer solid
    REQUIRE(sparse_len > 0.30 * solid_len);           // but a tighter lattice than full stride, not stripped bare

    // With the knob OFF the cell reverts to the solid-small-cells rule (the die-5 v7
    // behaviour): the interior prints dense, essentially the density-0 length.
    {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 20;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_adaptive_stride.value = false;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        double xy = 0.0;
        for (auto &p : mp->paths)
            for (size_t i = 1; i < p.polyline.points.size(); ++i)
                xy += (Point(p.polyline.points[i].x(), p.polyline.points[i].y()) -
                       Point(p.polyline.points[i-1].x(), p.polyline.points[i-1].y())).cast<double>().norm();
        REQUIRE(xy * SCALING_FACTOR > 0.9 * solid_len);
    }
}

TEST_CASE("Continuous shell cell size knob reaches the virtual-site partition", "[ContinuousShell]") {
    // continuous_shell_cell_size overrides the farthest-point sampling target. On a two-bore
    // block a SMALL target seeds extra virtual cells between/around the bores while a target
    // larger than every material gap adds none - the partitions, and hence the toolpaths, must
    // genuinely differ. (Cell counts are not observable here, so path length is the proxy.)
    ExPolygon block;
    block.contour = Polygon{ {-20000000,-12000000}, {20000000,-12000000},
                             {20000000, 12000000}, {-20000000, 12000000} };
    block.holes.push_back(Polygon{ {-15000000,-4000000}, {-15000000, 4000000},
                                   {-7000000, 4000000}, {-7000000,-4000000} });
    block.holes.push_back(Polygon{ {7000000,-4000000}, {7000000, 4000000},
                                   {15000000, 4000000}, {15000000,-4000000} });
    Surfaces surfaces{ Surface(stInternal, block) };
    Flow flow(0.42f, 0.2f, 0.4f);
    EnvFlag bubble("CS_BUBBLE", "1");
    auto xy_len = [&](double cell_size) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_cell_split.value     = true;
        cfg.continuous_shell_cell_size.value      = cell_size;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 0));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        double xy = 0.0;
        for (const ExtrusionPath *p : paths)
            for (size_t i = 1; i < p->polyline.points.size(); ++i)
                xy += (Point(p->polyline.points[i].x(), p->polyline.points[i].y()) -
                       Point(p->polyline.points[i-1].x(), p->polyline.points[i-1].y())).cast<double>().norm();
        return xy * SCALING_FACTOR;
    };
    const double fine   = xy_len(5.0);
    const double coarse = xy_len(40.0);
    INFO("fine(5mm)=" << fine << "mm coarse(40mm)=" << coarse << "mm");
    REQUIRE(std::abs(fine - coarse) > 0.02 * std::max(fine, coarse));
}

TEST_CASE("Continuous shell bore obstacles swoop the lattice around small bores", "[ContinuousShell]") {
    // bore_obstacles mode: a disc with two small bores fills like a SOLID puck - one global
    // lattice - with the kept rings rerouted around each bore's keep-out and a dedicated shell
    // annulus per bore. The bores must stay hollow (nothing inside, no chords across), and the
    // fill must come out LIGHTER than the cell partition of the same disc (the whole point:
    // no per-cell solid mosaics, one uniform lattice).
    auto circle = [](double cx, double cy, double r, bool cw) {
        Polygon p;
        for (int i = 0; i < 72; ++i) {
            const double a = 2.0 * M_PI * i / 72;
            p.points.emplace_back(coord_t(std::llround(cx + r * std::cos(a))),
                                  coord_t(std::llround(cy + r * std::sin(a))));
        }
        if (cw) p.reverse();
        return p;
    };
    ExPolygon disc;
    disc.contour = circle(0, 0, 16000000, false);
    disc.holes.push_back(circle( 7000000, 0, 2500000, true));
    disc.holes.push_back(circle(-7000000, 0, 2500000, true));
    Surfaces surfaces{ Surface(stInternal, disc) };
    Flow flow(0.42f, 0.2f, 0.4f);
    EnvFlag bubble("CS_BUBBLE", "1");
    auto slice = [&](bool obstacles, int density = 25, size_t layer = 0) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = density;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_cell_split.value     = true;
        cfg.continuous_shell_bore_obstacles.value = obstacles;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, layer));
        REQUIRE(e.entities.size() == 1);
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        double xy = 0.0;
        size_t bore_pts = 0, bore_mids = 0;
        auto in_bore = [](const Point &q) {
            auto inside = [&](double bx) {
                const double dx = double(q.x()) - bx, dy = double(q.y());
                return dx * dx + dy * dy < 2300000.0 * 2300000.0;
            };
            return inside(7000000.0) || inside(-7000000.0);
        };
        for (const ExtrusionPath *p : paths)
            for (size_t i = 0; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i].x(), p->polyline.points[i].y());
                if (in_bore(a)) ++bore_pts;
                if (i > 0) {
                    const Point b(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                    xy += (a - b).cast<double>().norm();
                    if (in_bore(Point((a + b) / 2))) {
                        ++bore_mids;
                        std::fprintf(stderr, "[chord] seg %zu of %zu mm3/mm=%.4f ctx:", i, p->polyline.points.size(), p->mm3_per_mm);
                        for (size_t w2 = i >= 4 ? i - 4 : 0; w2 < std::min(p->polyline.points.size(), i + 4); ++w2)
                            std::fprintf(stderr, " (%.2f,%.2f)", p->polyline.points[w2].x() / 1e6, p->polyline.points[w2].y() / 1e6);
                        std::fprintf(stderr, "\n");
                    }
                }
            }
        INFO((obstacles ? "swoop" : "cells") << ": xy=" << xy * SCALING_FACTOR << "mm bore_pts=" << bore_pts);
        REQUIRE(bore_pts == 0);
        REQUIRE(bore_mids == 0);
        return xy * SCALING_FACTOR;
    };
    const double swoop = slice(true);
    const double cells = slice(false);
    INFO("swoop=" << swoop << "mm cells=" << cells << "mm");
    REQUIRE(swoop < 0.97 * cells);      // the uniform lattice must beat the cell mosaic on mass
    REQUIRE(swoop > 0.30 * cells);      // sanity: the interior didn't silently vanish

    // BORE SANCTITY across layers and modes: connector/stitch geometry precesses with the
    // seam per layer, so a single layer can pass while another chords a bore (every die-5
    // PRINT had two bore chords at ONE bottom-shell layer, from the classic solid path's
    // inter-cell connectors - visible straight through the pips). slice() asserts
    // bore_pts == bore_mids == 0 internally; sweeping layers exercises the fallbacks.
    for (size_t layer : { 1, 2, 3, 4, 5 }) {
        slice(true, 25, layer);         // swoop lattice layers
        slice(false, 0, layer);         // classic SOLID (the bottom-shell path that chorded)
    }
}

TEST_CASE("Continuous shell mixed mode swoops small bores beside an oversized one", "[ContinuousShell]") {
    // One hole too big to swoop (a 10mm square - the counterbored disc's taper class) plus two
    // small round bores. Pure cells (obstacles off) partitions around ALL of them and runs
    // dense; mixed mode should partition only around the big hole and swoop the small ones,
    // coming out LIGHTER, with every bore still hollow.
    auto circle = [](double cx, double cy, double r, bool cw) {
        Polygon p;
        for (int i = 0; i < 72; ++i) {
            const double a = 2.0 * M_PI * i / 72;
            p.points.emplace_back(coord_t(std::llround(cx + r * std::cos(a))),
                                  coord_t(std::llround(cy + r * std::sin(a))));
        }
        if (cw) p.reverse();
        return p;
    };
    ExPolygon disc;
    disc.contour = circle(0, 0, 16000000, false);
    disc.holes.push_back(Polygon{ {-5000000, -5000000}, {-5000000, 5000000},
                                  {5000000, 5000000}, {5000000, -5000000} });   // big: fails the size gate
    disc.holes.push_back(circle( 10500000, 0, 2000000, true));
    disc.holes.push_back(circle(-10500000, 0, 2000000, true));
    Surfaces surfaces{ Surface(stInternal, disc) };
    Flow flow(0.42f, 0.2f, 0.4f);
    EnvFlag bubble("CS_BUBBLE", "1");
    auto slice = [&](bool obstacles) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_cell_split.value     = true;
        cfg.continuous_shell_bore_obstacles.value = obstacles;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 0));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        double xy = 0.0;
        size_t bore_pts = 0, bore_mids = 0;
        auto in_bore = [](const Point &q) {
            auto inc = [&](double bx) {
                const double dx = double(q.x()) - bx, dy = double(q.y());
                return dx * dx + dy * dy < 1800000.0 * 1800000.0;
            };
            return inc(10500000.0) || inc(-10500000.0);
        };
        for (const ExtrusionPath *p : paths)
            for (size_t i = 0; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i].x(), p->polyline.points[i].y());
                if (in_bore(a)) ++bore_pts;
                if (i > 0) {
                    const Point b(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                    xy += (a - b).cast<double>().norm();
                    if (in_bore(Point((a + b) / 2))) ++bore_mids;
                }
            }
        REQUIRE(bore_pts == 0);
        REQUIRE(bore_mids == 0);
        return xy * SCALING_FACTOR;
    };
    const double mixed = slice(true);
    const double cells = slice(false);
    INFO("mixed=" << mixed << "mm cells=" << cells << "mm");
    REQUIRE(mixed < 0.97 * cells);      // swooping the small bores must save material
    REQUIRE(mixed > 0.30 * cells);
}

TEST_CASE("Continuous shell bridges a blind bore cap instead of sagging rings over it", "[ContinuousShell]") {
    // A solid cap layer whose region has NO hole but sits over an open bore below
    // (bridge_over island): the builder must ring the opening (like the layer below did)
    // and span it with a serpentine patch anchored on the surrounding material.
    auto circle = [](double cx, double cy, double r) {
        Polygon p;
        for (int i = 0; i < 72; ++i) {
            const double a = 2.0 * M_PI * i / 72;
            p.points.emplace_back(coord_t(std::llround(cx + r * std::cos(a))),
                                  coord_t(std::llround(cy + r * std::sin(a))));
        }
        return p;
    };
    ExPolygon disc;
    disc.contour = circle(0, 0, 10000000);
    Surfaces surfaces{ Surface(stInternal, disc) };
    ExPolygons bridge_over{ ExPolygon(circle(0, 0, 3000000)) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 21,
                                                            true /*force_solid*/, &bridge_over));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    double in_island = 0.0, ring_zone = 0.0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
            const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
            const Point m((a + b) / 2);
            const double r = m.cast<double>().norm();
            const double len = (b - a).cast<double>().norm() * SCALING_FACTOR;
            if (r < 2800000.0) in_island += len;
            else if (r < 3800000.0) ring_zone += len;
        }
    INFO("bridge len in island=" << in_island << "mm, ring zone=" << ring_zone << "mm");
    REQUIRE(in_island > 40.0);          // the serpentine patch spans the opening
    REQUIRE(ring_zone > 30.0);          // and rings wrap it (the cap does not just vanish)
}

TEST_CASE("Continuous shell localized shell patch keeps the rest of the layer sparse", "[ContinuousShell]") {
    // The cap/floor shell stack over a bore must be a LOCAL solid patch in an otherwise
    // sparse lattice layer - not a whole-layer solid fill (a 6mm cap used to force the
    // entire 40mm disc solid for shell-count layers).
    auto circle = [](double cx, double cy, double r) {
        Polygon p;
        for (int i = 0; i < 72; ++i) {
            const double a = 2.0 * M_PI * i / 72;
            p.points.emplace_back(coord_t(std::llround(cx + r * std::cos(a))),
                                  coord_t(std::llround(cy + r * std::sin(a))));
        }
        return p;
    };
    ExPolygon disc;
    disc.contour = circle(0, 0, 14000000);
    Surfaces surfaces{ Surface(stInternal, disc) };
    ExPolygons patch{ ExPolygon(circle(0, 0, 3600000)) };   // bore + growth margin
    Flow flow(0.42f, 0.2f, 0.4f);
    EnvFlag bubble("CS_BUBBLE", "1");
    auto xy_of = [&](bool force_solid, const ExPolygons *islands) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_bore_obstacles.value = true;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 21,
                                                                force_solid, islands));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        double in_isl = 0.0, outside = 0.0;
        for (const ExtrusionPath *p : paths)
            for (size_t i = 1; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
                const double len = (b - a).cast<double>().norm() * SCALING_FACTOR;
                if (Point((a + b) / 2).cast<double>().norm() < 3200000.0) in_isl += len;
                else outside += len;
            }
        return std::make_pair(in_isl, outside);
    };
    const auto [patch_in, patch_out]  = xy_of(false, &patch);   // sparse layer + local patch
    const auto [solid_in, solid_out]  = xy_of(true, nullptr);   // old behaviour: whole-layer solid
    INFO("patch: in=" << patch_in << " out=" << patch_out << " | solid: in=" << solid_in << " out=" << solid_out);
    REQUIRE(patch_in  > 55.0);                    // solid serpentine fills the island
    // The rest of the layer stays LATTICE. (Not lower: the patch brings its own bore
    // machinery - ring walls and fan lanes - which is real, intended material.)
    REQUIRE(patch_out < 0.75 * solid_out);
}

TEST_CASE("Continuous shell refuses outlines that fragment at the bead scale", "[ContinuousShell]") {
    // A fractal-like silhouette shatters the inward march into hundreds of chains whose
    // rings overlap across fragments (Mandeltree: ~1.4x SOLID volume emitted as "sparse"
    // fill). The pathology guard probes a 2- and 6-bead inward offset and refuses when the
    // region fragments past any sane CS part. Model: a comb of 80 wide teeth on narrow
    // necks - each tooth SURVIVES the probe offset as its own island (mere sub-bead teeth
    // would vanish instead of fragmenting, hiding the pathology).
    ExPolygon comb;
    Points pts;
    const coord_t base_top = 2000000, neck_top = 3000000, tooth_top = 8000000;
    const coord_t pitch = 3000000, neck_hw = 300000, tooth_hw = 1000000;
    pts.emplace_back(0, 0);
    pts.emplace_back(0, base_top);
    for (int t = 0; t < 80; ++t) {                   // tooth body on a sub-bead neck
        const coord_t cx = t * pitch + pitch / 2;
        pts.emplace_back(cx - neck_hw,  base_top);
        pts.emplace_back(cx - neck_hw,  neck_top);
        pts.emplace_back(cx - tooth_hw, neck_top);
        pts.emplace_back(cx - tooth_hw, tooth_top);
        pts.emplace_back(cx + tooth_hw, tooth_top);
        pts.emplace_back(cx + tooth_hw, neck_top);
        pts.emplace_back(cx + neck_hw,  neck_top);
        pts.emplace_back(cx + neck_hw,  base_top);
    }
    pts.emplace_back(80 * pitch, base_top);
    pts.emplace_back(80 * pitch, 0);
    comb.contour = Polygon(pts);
    Surfaces surfaces{ Surface(stInternal, comb) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    cfg.continuous_shell_infill_weave.value   = 100;
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE_FALSE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 3));
    REQUIRE_THAT(reason, Catch::Matchers::ContainsSubstring("fragments"));
}

TEST_CASE("Continuous shell point-tangent regions are refused or filled without air chords", "[ContinuousShell]") {
    // Two circles sharing exactly ONE vertex (tangent cylinders: one xy contact point per
    // layer). Depending on Clipper's union this is either TWO regions (must refuse: CS
    // requires one connected cross-section) or ONE pinched region - then the fill must
    // never chord through the air beside the tangent point: every extruding segment's
    // midpoint must stay inside the (slightly grown) material.
    auto circle_with_vertex_at = [](double cx, double r, double phase) {
        Polygon p;
        for (int i = 0; i < 72; ++i) {
            const double a = phase + 2.0 * M_PI * i / 72;
            p.points.emplace_back(coord_t(std::llround(cx + r * std::cos(a))),
                                  coord_t(std::llround(r * std::sin(a))));
        }
        return p;
    };
    ExPolygon left, right;
    left.contour  = circle_with_vertex_at(-8000000, 8000000, 0.0);      // vertex at (0,0)
    right.contour = circle_with_vertex_at( 8000000, 8000000, M_PI);     // vertex at (0,0)
    Surfaces surfaces{ Surface(stInternal, left), Surface(stInternal, right) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    cfg.continuous_shell_infill_weave.value   = 100;
    ExtrusionEntityCollection e; std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 3);
    if (!built) {
        SUCCEED("refused: " << reason);
        return;
    }
    // Built: no segment may cross the air outside the pinched region.
    ExPolygons grown = offset_ex(union_ex(Polygons{ left.contour, right.contour }), float(scale_(0.15)));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    size_t air = 0, total = 0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
            const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
            ++total;
            const Point mid((a + b) / 2);
            bool inside = false;
            for (const ExPolygon &g : grown)
                if (g.contains(mid)) { inside = true; break; }
            if (!inside) ++air;
        }
    INFO("built with " << total << " segments, " << air << " with midpoints in air");
    REQUIRE(air == 0);
}

TEST_CASE("Continuous shell woven infill sparsifies a THICK annular wall and keeps the bore hollow", "[ContinuousShell]") {
    // A thick tube wall (large OD, small ID) must NOT print 100% solid just because it has a bore:
    // keep a solid skin on BOTH the outer and bore faces and thin the deep middle, all one bead,
    // with the bore left hollow (no chords across it). Guards the annular sparsification path.
    Surface thick(stInternal, ExPolygon(
        Polygon{{0,0},{40000000,0},{40000000,40000000},{0,40000000}},                 // 40mm OD
        Polygon{{18000000,18000000},{18000000,22000000},{22000000,22000000},{22000000,18000000}}));  // 4mm bore
    Surfaces surfaces{ thick };
    Flow flow(0.5f, 0.25f, 0.5f);
    const Point center(20000000, 20000000);

    auto measure = [&](int density) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 3;
        cfg.continuous_shell_infill_density.value = density;
        cfg.continuous_shell_infill_weave.value   = 0;   // pure geometric sparsify (deterministic)
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
        REQUIRE(e.entities.size() == 1);                 // one continuous entity (single bead)
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::vector<Point> pts;
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) pts.emplace_back(q.x(), q.y());
        REQUIRE(pts.size() > 8);
        double len = 0.0, minr = std::numeric_limits<double>::max();
        for (size_t i = 0; i < pts.size(); ++i) {
            minr = std::min(minr, (pts[i] - center).cast<double>().norm());
            if (i > 0) len += (pts[i] - pts[i-1]).cast<double>().norm();
        }
        return std::make_pair(len * SCALING_FACTOR, minr * SCALING_FACTOR);
    };

    const auto solid  = measure(0);
    const auto sparse = measure(20);
    REQUIRE(sparse.first  < 0.7 * solid.first);   // the thick wall is sparsified, not left 100% solid
    REQUIRE(sparse.second > 1.5);                 // the ~2mm bore stays hollow (no chord dives to the centre)
    REQUIRE(solid.second  > 1.5);                 // baseline also keeps the bore
}

TEST_CASE("Continuous shell brick-staggers a multi-layer solid cap's interior, not its surface", "[ContinuousShell]") {
    // A solid disc printed as a stack of solid cap layers (force_solid). Concentric rings would
    // otherwise stack into aligned radial valleys (a weak plane), so odd layers shift the interior
    // rings inward by ~half a bead to nest between the courses above/below. The OUTER surface must be
    // untouched (smooth wall); only the buried interior moves.
    const Point center(10000000, 10000000);
    Polygon circ;
    for (int i = 0; i < 64; ++i) {
        const double a = 2.0 * M_PI * i / 64.0;
        circ.points.emplace_back(center.x() + coord_t(8000000 * std::cos(a)), center.y() + coord_t(8000000 * std::sin(a)));
    }
    Surface disc(stInternal, ExPolygon(circ));
    Surfaces surfaces{ disc };
    Flow flow(0.5f, 0.25f, 0.5f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 20;   // a woven-infill part, so its solid caps stagger
    cfg.top_surface_pattern.value = ipConcentric;     // test the concentric+stagger path (not serpentine)

    struct M { double interior_mean_mm; double max_r_mm; };
    auto measure = [&](size_t layer) {
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, (layer & 1) != 0, layer, /*force_solid=*/true));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        double sum = 0.0, maxr = 0.0; size_t n = 0;
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) {
            const double r = (Point(q.x(), q.y()) - center).cast<double>().norm();
            maxr = std::max(maxr, r);
            if (r < 7000000.0) { sum += r; ++n; }   // interior only (well inside the ~2-ring shell)
        }
        REQUIRE(n > 20);
        return M{ (sum / double(n)) * SCALING_FACTOR, maxr * SCALING_FACTOR };
    };
    const M even = measure(0);   // force_solid but no shift on even courses
    const M odd  = measure(1);   // interior shifted inward on odd courses
    // The outer surface is untouched (shell rings are never staggered).
    REQUIRE(std::abs(even.max_r_mm - odd.max_r_mm) < 0.1);
    // The interior is shifted inward by ~half a 0.5mm bead on the odd course (the brick offset).
    REQUIRE(even.interior_mean_mm - odd.interior_mean_mm > 0.05);
    REQUIRE(even.interior_mean_mm - odd.interior_mean_mm < 0.4);
}

TEST_CASE("Continuous shell declines a serpentine cap on a concave section and lays no chord across the notch", "[ContinuousShell]") {
    // An L-shaped solid cap cannot be filled by one boustrophedon bead without extruding
    // straight across the concave notch, so the serpentine builder must decline and the
    // cross-section fall back to the concentric spiral. Assert the fill succeeds and that
    // no segment midpoint lands in the notch.
    // L-shape: 20x20mm square minus the 10x10mm top-right quadrant.
    Polygon ell{ {0, 0}, {20000000, 0}, {20000000, 10000000}, {10000000, 10000000},
                 {10000000, 20000000}, {0, 20000000} };
    const ExPolygon region(ell);
    Surface cap(stInternal, region);
    Surfaces surfaces{ cap };
    Flow flow(0.5f, 0.25f, 0.5f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.top_surface_pattern.value = ipRectilinear;

    ExtrusionEntityCollection e; std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason,
                                                                       false, /*layer_index=*/4, /*force_solid=*/true);
    INFO(reason);
    REQUIRE(built);

    // The notch quadrant, shrunk 1mm on each side so legitimate boundary-hugging beads
    // near the notch corner do not count.
    const ExPolygon notch(Polygon{ {11000000, 11000000}, {19000000, 11000000},
                                   {19000000, 19000000}, {11000000, 19000000} });
    std::vector<ExtrusionPath *> paths;
    ContinuousShellPathBuilder::collect_paths(e.entities.front(), paths);
    REQUIRE_FALSE(paths.empty());
    size_t chords_in_notch = 0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
            const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
            if (notch.contains(Point((a + b) / 2)))
                ++chords_in_notch;
        }
    REQUIRE(chords_in_notch == 0);
}

TEST_CASE("Continuous shell fills a solid cap with a per-layer-rotated serpentine (rectilinear pattern)", "[ContinuousShell]") {
    // With a rectilinear-family surface pattern, a solid cap plane is filled by a connected serpentine
    // (inside a concentric wall) whose direction rotates 45deg per layer - the classic cross-hatch
    // interlock. Assert it stays one contained bead and that consecutive layers run ~90deg apart.
    const Point center(10000000, 10000000);
    Polygon circ;
    for (int i = 0; i < 64; ++i) {
        const double a = 2.0 * M_PI * i / 64.0;
        circ.points.emplace_back(center.x() + coord_t(8000000 * std::cos(a)), center.y() + coord_t(8000000 * std::sin(a)));
    }
    const ExPolygon region(circ);
    Surface disc(stInternal, region);
    Surfaces surfaces{ disc };
    Flow flow(0.5f, 0.25f, 0.5f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.top_surface_pattern.value = ipRectilinear;
    const ExPolygons valid = offset_ex(region, scale_(0.2));   // grown to absorb the outer bead + noise

    auto dominant_angle_deg = [&](size_t layer) {
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, (layer & 1) != 0, layer, /*force_solid=*/true));
        REQUIRE(e.entities.size() == 1);
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::vector<Point> pts;
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) pts.emplace_back(q.x(), q.y());
        REQUIRE(pts.size() > 40);
        std::array<double, 18> hist{};                         // 10deg bins over [0,180)
        for (size_t i = 0; i < pts.size(); ++i) {
            bool in = false; for (const ExPolygon &ep : valid) if (ep.contains(pts[i])) { in = true; break; }
            REQUIRE(in);                                       // every point stays inside the cap outline
            if (i > 0) {
                const double dx = double(pts[i].x() - pts[i-1].x()), dy = double(pts[i].y() - pts[i-1].y());
                const double len = std::hypot(dx, dy);
                if (len > scale_(2.0)) {                        // long runs = the serpentine scan lines
                    double a = std::atan2(dy, dx) * 180.0 / M_PI; a = std::fmod(a + 180.0, 180.0);
                    hist[std::min<size_t>(17, size_t(a / 10.0))] += len;
                }
            }
        }
        size_t best = 0; for (size_t b = 1; b < hist.size(); ++b) if (hist[b] > hist[best]) best = b;
        return double(best) * 10.0 + 5.0;
    };
    const double a_even = dominant_angle_deg(0);
    const double a_odd  = dominant_angle_deg(1);
    double diff = std::abs(a_even - a_odd);
    if (diff > 90.0) diff = 180.0 - diff;                      // angle distance modulo 180
    REQUIRE(diff > 60.0);                                      // consecutive layers cross-hatch (~90deg apart)
}

TEST_CASE("Continuous shell infill weave frequency sets the lobe count, amplitude stays full", "[ContinuousShell]") {
    // The weave setting controls the lobe FREQUENCY, not amplitude (the weave always uses full
    // interlock). On a round part a higher setting makes the interior rings oscillate radially more
    // times per revolution, so counting the radial in/out reversals of the interior fill is a direct,
    // sensitive proxy for the lobe count (bead length is not - the shape-aware clamp trims wiggle).
    const Point center(10000000, 10000000);
    Polygon circ;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circ.points.emplace_back(center.x() + coord_t(8000000 * std::cos(a)), center.y() + coord_t(8000000 * std::sin(a)));
    }
    Surface disc(stInternal, ExPolygon(circ));
    Surfaces surfaces{ disc };
    Flow flow(0.5f, 0.25f, 0.5f);
    auto radial_reversals = [&](int freq) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;
        cfg.continuous_shell_infill_weave.value   = freq;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 2));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::vector<double> rad;                               // interior-band radii, in path order
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) {
            const double r = (Point(q.x(), q.y()) - center).cast<double>().norm() * SCALING_FACTOR;
            if (r >= 2.0 && r <= 6.0) rad.push_back(r);
        }
        int rev = 0;                                           // count in/out reversals (2 per lobe per ring)
        for (size_t i = 2; i < rad.size(); ++i) {
            const double d0 = rad[i-1] - rad[i-2], d1 = rad[i] - rad[i-1];
            if (d0 != 0.0 && d1 != 0.0 && ((d0 > 0.0) != (d1 > 0.0))) ++rev;
        }
        return rev;
    };
    const int base = radial_reversals(100);   // 6 lobes
    const int fine = radial_reversals(200);   // 12 lobes
    REQUIRE(base > 4);                         // there IS a radial weave
    REQUIRE(fine * 3 > base * 4);              // a higher frequency = more lobes = more radial reversals (>1.33x)
}

TEST_CASE("Continuous shell bubble weave stacks crests vertically (Z-coherent lattice)", "[ContinuousShell]") {
    // The bubble lattice (CS_BUBBLE) needs crest welds that STACK across layers. Kiss mode rotates the
    // lobes by dphase per layer (weave_shift) - measured at exactly +90deg/layer of m=6 phase in sliced
    // g-code - so a crest never sits on the crest below and no vertical cell wall can form (the print
    // reads as a decoherent web). Bubble mode must (a) pin the lattice phase across layers regardless of
    // the configured shift, and (b) supply the vertical structure via a slow cosine carrier on the
    // amplitude instead (period ~21 layers at K=4). Guard both, with kiss mode as the drift control.
    const Point center(10000000, 10000000);
    Polygon circ;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circ.points.emplace_back(center.x() + coord_t(8000000 * std::cos(a)), center.y() + coord_t(8000000 * std::sin(a)));
    }
    Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
    Flow flow(0.5f, 0.25f, 0.5f);
    // Interior-band toolpath points + band length for one layer (the shells at r~7.3/7.8 are excluded,
    // so this samples only the woven infill rings).
    auto band = [&](size_t layer) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value     = 25;   // K=4
        cfg.continuous_shell_infill_weave.value       = 100;  // 6 lobes
        cfg.continuous_shell_infill_weave_shift.value = 50;   // pi/2 per layer in kiss mode; bubble must ignore it
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, (layer & 1) != 0, layer, false));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::pair<std::vector<Vec2d>, double> out{ {}, 0.0 };
        for (auto &p : mp->paths) {
            const Vec2d *prev = nullptr; Vec2d prev_store;
            for (auto &q : p.polyline.points) {
                const Vec2d v = (Point(q.x(), q.y()) - center).cast<double>() * SCALING_FACTOR;
                const double r = v.norm();
                if (r >= 2.2 && r <= 6.9) {
                    if (prev != nullptr) out.second += (v - *prev).norm();
                    out.first.push_back(v);
                    prev_store = v; prev = &prev_store;
                } else
                    prev = nullptr;
            }
        }
        return out;
    };
    // Median nearest-neighbour distance from one layer's band to the previous layer's: registered
    // lattices score near zero, a quarter-wave-rotated lattice scores a large fraction of its amplitude.
    auto median_nn = [](const std::vector<Vec2d> &a, const std::vector<Vec2d> &b) {
        std::vector<double> d; d.reserve(a.size());
        for (const Vec2d &p : a) {
            double best = std::numeric_limits<double>::max();
            for (const Vec2d &q : b) best = std::min(best, (p - q).squaredNorm());
            d.push_back(std::sqrt(best));
        }
        std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
        return d[d.size() / 2];
    };
    double lock = 0, len_peak = 0, len_zero = 0;
    {
        EnvFlag bubble("CS_BUBBLE", "1");
        const auto b10 = band(10), b11 = band(11);        // carrier ~-0.99 on both: registered, near-full amplitude
        lock     = median_nn(b11.first, b10.first);
        len_peak = band(0).second;                        // carrier +1: full lobes
        len_zero = band(5).second;                        // carrier ~0.07: rings relax toward plain concentric
    }
    const auto k10 = band(10), k11 = band(11);            // kiss mode control: same layers, rotating lattice
    const double drift = median_nn(k11.first, k10.first);
    WARN("bubble adjacent-layer registration=" << lock << "mm vs kiss drift=" << drift
         << "mm; band length carrier-peak=" << len_peak << "mm carrier-zero=" << len_zero << "mm");
    REQUIRE(lock < 0.2);                       // bubble: consecutive layers register (crests stack)
    REQUIRE(drift > 2.0 * lock);               // kiss control: the same measurement sees the rotation
    REQUIRE(drift > 0.25);
    REQUIRE(len_peak > 1.1 * len_zero);        // carrier: lobes at the peak, near-plain rings at the zero
}

TEST_CASE("Continuous shell bubble weld flats stay off the shell (no wall bumps)", "[ContinuousShell]") {
    // The wall-adjacent infill ring sits one bead inside the shell but carries a fusion-sized
    // amplitude, so the old poke-through clamp parked its crests exactly ON the innermost shell
    // ring's path - a second bead deposited along the wall, whose displaced melt telegraphed
    // through the 2-bead shell as alternating surface bumps tracking the internal weld pads.
    // The clamp now keeps ~0.85 bead of centreline clearance: crests flatten BESIDE the shell
    // (~15% edge overlap - a normal lateral weld). Guard: the toolpath length in the shell's
    // radial zone must match a plain sparse fill (no crest material on the wall), while the
    // flattened crests still park just inside the clearance line.
    const Point center(10000000, 10000000);
    Polygon circ;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circ.points.emplace_back(center.x() + coord_t(8000000 * std::cos(a)), center.y() + coord_t(8000000 * std::sin(a)));
    }
    Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
    Flow flow(0.5f, 0.25f, 0.5f);
    // Path length whose segment midpoints fall in a radial band, at carrier-peak layer 0.
    auto band_len = [&](int weave, double lo, double hi) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;   // K=4
        cfg.continuous_shell_infill_weave.value   = weave;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 0));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        double len = 0;
        for (auto &p : mp->paths)
            for (size_t k = 1; k < p.polyline.points.size(); ++k) {
                const auto &a = p.polyline.points[k - 1];
                const auto &b = p.polyline.points[k];
                const Vec2d va = (Point(a.x(), a.y()) - center).cast<double>() * SCALING_FACTOR;
                const Vec2d vb = (Point(b.x(), b.y()) - center).cast<double>() * SCALING_FACTOR;
                const double rm = 0.5 * (va.norm() + vb.norm());
                if (rm >= lo && rm <= hi) len += (vb - va).norm();
            }
        return len;
    };
    const double spacing_mm = unscale<double>(flow.scaled_spacing());
    const double shell_in   = 8.0 - 0.5 * 0.5 - spacing_mm;      // innermost shell ring centreline
    const double clear_line = shell_in - 0.85 * spacing_mm;      // the new clamp line
    double bub_zone = 0, plain_zone = 0, bub_flat = 0, plain_flat = 0;
    {
        EnvFlag bubble("CS_BUBBLE", "1");
        bub_zone = band_len(100, clear_line + 0.10, shell_in + 0.15);    // shell zone above the clamp line
        bub_flat = band_len(100, clear_line - 0.035, clear_line + 0.05); // where flats should now park
    }
    plain_zone = band_len(0, clear_line + 0.10, shell_in + 0.15);
    plain_flat = band_len(0, clear_line - 0.035, clear_line + 0.05);     // plain's first ring sits below this band
    WARN("shell-zone len: bubble=" << bub_zone << "mm plain=" << plain_zone
         << "mm; clamp-line flats: bubble=" << bub_flat << "mm plain=" << plain_flat << "mm");
    REQUIRE(bub_zone < 1.08 * plain_zone);     // no crest material parked on/along the shell any more
    REQUIRE(bub_flat > plain_flat + 3.0);      // the flattened crests exist, parked at the clearance line
}

TEST_CASE("Continuous shell dyadic bubble lobes scale with ring perimeter", "[ContinuousShell]") {
    // CS_BUBBLE_DYADIC gives each ring the lobes*2^m nearest perimeter/CS_BUBBLE_W, so outer
    // rings carry more, narrower cells instead of stretching with perimeter (hierarchical
    // half-bonded boundaries between bands). Guard via radial-reversal counting (2 per lobe per
    // ring): on a r=15mm disc the outer two kept rings should pick 2x the base count while the
    // core stays at base - so the outer band out-reverses the equally-sized inner band only in
    // dyadic mode.
    const Point center(20000000, 20000000);
    Polygon circ;
    for (int i = 0; i < 144; ++i) {
        const double a = 2.0 * M_PI * i / 144.0;
        circ.points.emplace_back(center.x() + coord_t(15000000 * std::cos(a)), center.y() + coord_t(15000000 * std::sin(a)));
    }
    Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
    Flow flow(0.5f, 0.25f, 0.5f);
    auto reversals = [&](double lo, double hi, bool dyadic) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;   // K=4
        cfg.continuous_shell_infill_weave.value   = 100;  // base 6 lobes
        cfg.continuous_shell_infill_dyadic.value  = dyadic;   // the real (per-object) config knob
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 0));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        std::vector<double> rad;
        for (auto &p : mp->paths) for (auto &q : p.polyline.points) {
            const double r = (Point(q.x(), q.y()) - center).cast<double>().norm() * SCALING_FACTOR;
            if (r >= lo && r <= hi) rad.push_back(r);
        }
        // Hysteresis direction-change counter: the clamped crest FLATS are constant-radius arcs
        // whose resample jitter would otherwise register as dozens of phantom reversals.
        int rev = 0, dir = 0; double ref = rad.empty() ? 0.0 : rad[0];
        for (const double r : rad) {
            if (std::abs(r - ref) < 0.05) continue;
            const int d = (r > ref) ? 1 : -1;
            if (dir != 0 && d != dir) ++rev;
            dir = d; ref = r;
        }
        return rev;
    };
    int dy_outer = 0, dy_inner = 0, cl_outer = 0, cl_inner = 0;
    {
        EnvFlag bubble("CS_BUBBLE", "1");
        dy_outer = reversals(11.0, 14.1, true);    // outermost two kept rings (expect base*2)
        dy_inner = reversals(2.0, 6.0, true);      // two core rings (expect base)
        cl_outer = reversals(11.0, 14.1, false);   // classic bubble: same global count everywhere
        cl_inner = reversals(2.0, 6.0, false);
    }
    WARN("reversals dyadic outer=" << dy_outer << " inner=" << dy_inner
         << "; classic outer=" << cl_outer << " inner=" << cl_inner);
    REQUIRE(dy_outer > dy_inner * 3 / 2);       // dyadic: outer band carries ~2x the lobes
    REQUIRE(cl_outer < cl_inner * 3 / 2);       // classic control: no such gradient
    REQUIRE(dy_outer > cl_outer * 3 / 2);       // and it is the dyadic switch that adds them
}

TEST_CASE("Continuous shell cell split fills a woven dogbone as near-star cells", "[ContinuousShell]") {
    // continuous_shell_cell_split pre-cuts deeply concave sections at their deepest notches so
    // every piece is near-star-shaped and the bubble weave is valid per cell. On the dogbone the
    // cuts are the vertical lines x=+-5 through the notch-floor corners (lobe | neck | lobe).
    // Discriminator: inside thin vertical bands around the cut lines, a SPLIT fill runs its cell
    // shell beads ALONG the cut (vertical travel dominates), while the unsplit fill's rings CROSS
    // it (horizontal dominates). Also: the bead must stay inside the region and out of the notch
    // air in both modes, on every layer parity/carrier phase.
    Polygon dog{ {-15000000, -8000000}, {-5000000, -8000000}, {-5000000, -2000000},
                 {5000000, -2000000}, {5000000, -8000000}, {15000000, -8000000},
                 {15000000, 8000000}, {5000000, 8000000}, {5000000, 2000000},
                 {-5000000, 2000000}, {-5000000, 8000000}, {-15000000, 8000000} };
    Surfaces surfaces{ Surface(stInternal, ExPolygon(dog)) };
    Flow flow(0.42f, 0.2f, 0.4f);
    const ExPolygons grown = offset_ex(ExPolygons{ ExPolygon(dog) }, scale_(0.1f));
    auto in_notch = [](const Point &p) {
        const bool x_in = p.x() > -4000000 && p.x() < 4000000;
        const bool top  = p.y() >  3000000 && p.y() < 7000000;
        const bool bot  = p.y() < -3000000 && p.y() > -7000000;
        return x_in && (top || bot);
    };
    auto in_cut_band = [](const Point &p) {
        const bool x_in = (p.x() > 4600000 && p.x() < 5400000) || (p.x() > -5400000 && p.x() < -4600000);
        return x_in && p.y() > -2000000 && p.y() < 2000000;   // the neck span of the cut lines
    };
    // Returns (vertical fraction of travel in the cut bands, notch-crossing count, outside count).
    std::string diag;
    auto measure = [&](bool split, size_t layer, double &vfrac, size_t &notch, size_t &outside) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_cell_split.value     = split;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, (layer & 1) != 0, layer, false));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        double dx = 0, dy = 0;
        notch = 0; outside = 0;
        for (const ExtrusionPath *p : paths)
            for (size_t i = 1; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
                const Point mid((a + b) / 2);
                if (in_notch(mid)) ++notch;
                bool inside = false;
                for (const ExPolygon &g : grown) if (g.contains(mid)) { inside = true; break; }
                if (!inside) {
                    ++outside;
                    diag += " (" + std::to_string(unscale<double>(mid.x())) + "," +
                            std::to_string(unscale<double>(mid.y())) + ")len" +
                            std::to_string(unscale<double>((b - a).cast<double>().norm()));
                }
                if (in_cut_band(mid)) { dx += std::abs(double(b.x() - a.x())); dy += std::abs(double(b.y() - a.y())); }
            }
        vfrac = (dx + dy > 0.0) ? dy / (dx + dy) : 0.0;
    };
    EnvFlag bubble("CS_BUBBLE", "1");
    double v_split_min = 1.0, v_plain_max = 0.0;
    size_t plain_bad = 0;
    for (size_t layer = 3; layer < 27; layer += 4) {   // spans a full carrier period
        double v; size_t notch, outside;
        measure(true, layer, v, notch, outside);
        INFO("split layer " << layer << " vfrac=" << v << " notch=" << notch << " outside=" << outside << diag);
        REQUIRE(notch == 0);
        REQUIRE(outside == 0);
        v_split_min = std::min(v_split_min, v);
        // The UNSPLIT control is allowed to transgress: a deep notch makes r(theta) about the
        // global centroid multivalued, so the radial poke-through clamp cannot represent the
        // boundary there and weave crests overshoot into notch air - which is precisely the
        // failure cell_split exists to fix. Only the travel-direction metric is asserted on it.
        measure(false, layer, v, notch, outside);
        INFO("plain layer " << layer << " vfrac=" << v << " notch=" << notch << " outside=" << outside);
        plain_bad += notch + outside;
        v_plain_max = std::max(v_plain_max, v);
    }
    WARN("cut-band vertical travel fraction: split(min)=" << v_split_min << " unsplit(max)=" << v_plain_max
         << "; unsplit notch/outside transgressions across layers: " << plain_bad);
    REQUIRE(v_split_min > 0.15);               // split: boundary structure runs ALONG the cuts
    REQUIRE(v_plain_max < 0.7 * v_split_min);  // unsplit control: rings cross them instead (measured ~0.04-0.12)
}

TEST_CASE("Continuous shell cell-split interface weaves (recessed anti-phase walls)", "[ContinuousShell]") {
    // The interface between split cells must behave like intra-cell neighbours, not a dead wall:
    // (a) each cell is RECESSED half a bead from the cut, so the facing rings sit ~one bead apart
    //     instead of double-extruding the cut line at 200% flow, and
    // (b) the anti-phase interface field waves the two wall stacks oppositely, so the facing gap
    //     OSCILLATES along the cut - welded shut at crests, open lens voids between.
    // Measured on the dogbone's x=+3 cut at a carrier-peak layer.
    Polygon dog{ {-15000000, -8000000}, {-5000000, -8000000}, {-5000000, -2000000},
                 {5000000, -2000000}, {5000000, -8000000}, {15000000, -8000000},
                 {15000000, 8000000}, {5000000, 8000000}, {5000000, 2000000},
                 {-5000000, 2000000}, {-5000000, 8000000}, {-15000000, 8000000} };
    Surfaces surfaces{ Surface(stInternal, ExPolygon(dog)) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    cfg.continuous_shell_infill_weave.value   = 100;
    cfg.continuous_shell_cell_split.value     = true;
    EnvFlag bubble("CS_BUBBLE", "1");
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 0));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    // Facing-wall gap per 1mm slab along the cut (x = +5mm cut of the original test dogbone is at
    // x=+5? the notch corners sit at x=+-5 -> cuts at x=+-5): sample the +5 cut, neck span y<2.
    const double CUT_X = 5.0;
    std::array<double, 4> lsum{}, lcnt{}, rsum{}, rcnt{};
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Vec2d a = Point(p->polyline.points[i-1].x(), p->polyline.points[i-1].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d b = Point(p->polyline.points[i].x(),   p->polyline.points[i].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d m = 0.5 * (a + b);
            if (std::abs(m.x() - CUT_X) > 1.6) continue;          // lattice interfaces breathe wider than walls
            if (m.y() < -2.0 || m.y() >= 2.0) continue;
            const size_t slab = size_t((m.y() + 2.0));            // 1mm slabs over the neck span
            if (slab >= 4) continue;
            if (m.x() < CUT_X) { lsum[slab] += m.x(); lcnt[slab] += 1; }
            else               { rsum[slab] += m.x(); rcnt[slab] += 1; }
        }
    double gmin = 1e9, gmax = -1e9, gsum = 0; int gn = 0;
    for (size_t k = 0; k < 4; ++k) {
        if (lcnt[k] < 5 || rcnt[k] < 5) continue;
        const double gap = rsum[k] / rcnt[k] - lsum[k] / lcnt[k];
        gmin = std::min(gmin, gap); gmax = std::max(gmax, gap); gsum += gap; ++gn;
    }
    REQUIRE(gn >= 3);
    WARN("interface gap over " << gn << " slabs: mean=" << gsum / gn << "mm swing=" << gmax - gmin << "mm");
    REQUIRE(gsum / gn > 0.25);                 // recessed: facing material centres are separated, not coincident
    REQUIRE(gsum / gn < 2.2);                  // (lattice interfaces average wider than walls: valleys retreat)
    REQUIRE(gmax - gmin > 0.12);               // and the gap breathes along the cut (anti-phase field)
}

TEST_CASE("Continuous shell cell split peels a contiguous skin", "[ContinuousShell]") {
    // The object's exterior must never be stitched from per-cell arcs: with cell_split the seed's
    // outer contour is PEELED into its own single-loop chain at the tour boundary, so the visible
    // perimeter prints as one contiguous bead and the cut mouths terminate on its inner face.
    // Guard: walking the tour in order, runs of points on the skin band alternate with interior
    // runs only a handful of times (one skin loop + the layer handoff), and the skin band is
    // essentially fully covered.
    Polygon dog{ {-15000000, -8000000}, {-5000000, -8000000}, {-5000000, -2000000},
                 {5000000, -2000000}, {5000000, -8000000}, {15000000, -8000000},
                 {15000000, 8000000}, {5000000, 8000000}, {5000000, 2000000},
                 {-5000000, 2000000}, {-5000000, 8000000}, {-15000000, 8000000} };
    Surfaces surfaces{ Surface(stInternal, ExPolygon(dog)) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    cfg.continuous_shell_infill_weave.value   = 100;
    cfg.continuous_shell_cell_split.value     = true;
    EnvFlag bubble("CS_BUBBLE", "1");
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 5));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    // distance from a point to the dogbone outline
    auto dist_to_boundary = [&](const Vec2d &p) {
        double best = std::numeric_limits<double>::max();
        const size_t m = dog.points.size();
        for (size_t i = 0; i < m; ++i) {
            const Vec2d a  = dog.points[i].cast<double>() * SCALING_FACTOR;
            const Vec2d ab = dog.points[(i + 1) % m].cast<double>() * SCALING_FACTOR - a;
            const double t = std::clamp((p - a).dot(ab) / std::max(1e-9, ab.squaredNorm()), 0.0, 1.0);
            best = std::min(best, (p - (a + t * ab)).norm());
        }
        return best;
    };
    // Walk the tour accumulating arc length; record the tour position of every on-skin sample.
    // With the peel the skin is ONE contiguous stretch of the tour (small gaps where the ring
    // chords across sharp corners - a pre-existing ring artifact); stitched per-cell skin arcs
    // would put a cell's whole interior (100+mm of tour) between skin stretches.
    const double half_w = 0.21;                       // ext perimeter centreline depth
    // Contiguity metric: the smallest TOUR WINDOW containing 95% of the skin-band length.
    // The ext-perimeter pass spans ~one loop; stray band grazes elsewhere (e.g. wall #2's
    // corner chords dipping toward the boundary) fall in the 5% allowance. Stitched per-cell
    // skin arcs would spread the skin across the whole tour and blow the window wide open.
    double skin_len = 0, tour_arc = 0;
    std::vector<std::pair<double, double>> samples;   // (tour position, seg length)
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Vec2d a = Point(p->polyline.points[i-1].x(), p->polyline.points[i-1].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d b = Point(p->polyline.points[i].x(),   p->polyline.points[i].y()).cast<double>() * SCALING_FACTOR;
            const double seg = (b - a).norm();
            tour_arc += seg;
            if (std::abs(dist_to_boundary(0.5 * (a + b)) - half_w) < 0.12) {
                skin_len += seg;
                samples.emplace_back(tour_arc, seg);
            }
        }
    double window = tour_arc;
    for (size_t lo = 0; lo < samples.size(); ++lo) {           // min window covering >= 95%
        double cover = 0;
        for (size_t h2 = lo; h2 < samples.size(); ++h2) {
            cover += samples[h2].second;
            if (cover >= 0.95 * skin_len) {
                window = std::min(window, samples[h2].first - samples[lo].first + samples[lo].second);
                break;
            }
        }
    }
    WARN("skin len=" << skin_len << "mm of ~104mm loop; min tour window for 95% of it=" << window << "mm");
    REQUIRE(skin_len > 85.0);                  // the skin band is actually covered
    REQUIRE(window < 130.0);                   // and ~all of it lives in ONE loop-sized tour stretch
}

TEST_CASE("Continuous shell tapers flow in sub-2-bead walls (no cusp over-stuffing)", "[ContinuousShell]") {
    // Where the part is locally thinner than two beads, opposing wall passes share the space;
    // at constant flow that approaches 200% material at a cusp (the teardrop tail's dots).
    // The thin-wall taper scales each non-buried point's flow by local thickness. Wedge:
    // 30mm long, 8mm tall at the root, cusp at the tip - the last ~3mm is sub-2-bead.
    Polygon wedge{ {0, -4000000}, {30000000, 0}, {0, 4000000} };
    Surfaces surfaces{ Surface(stInternal, ExPolygon(wedge)) };
    Flow flow(0.42f, 0.2f, 0.4f);
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, nullptr, e, &reason, false, 3));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    const double nominal = flow.mm3_per_mm();
    double thick_sum = 0, thick_n = 0, thin_sum = 0, thin_n = 0, thin_min = 10.0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const double mx = 0.5 * (p->polyline.points[i-1].x() + p->polyline.points[i].x()) * SCALING_FACTOR;
            const double r  = p->mm3_per_mm / nominal;
            if (mx > 2.0 && mx < 10.0)  { thick_sum += r; thick_n += 1; }
            if (mx > 26.5)              { thin_sum += r; thin_n += 1; thin_min = std::min(thin_min, r); }
        }
    REQUIRE(thick_n > 10);
    REQUIRE(thin_n >= 3);
    WARN("flow ratio: thick mean=" << thick_sum / thick_n << " thin mean=" << thin_sum / thin_n
         << " thin min=" << thin_min);
    REQUIRE(thick_sum / thick_n > 0.9);        // full flow where the wall is healthy
    REQUIRE(thin_sum / thin_n < 0.92);         // the tip zone carries tapered segments
    REQUIRE(thin_min < 0.8);                   // deepest sampled point is meaningfully tapered
                                               // (the march clips the last ~0.7mm of a cusp, so the
                                               // theoretical 0.4 floor is unreachable on this wedge)
}

TEST_CASE("Continuous shell drops degenerate core rings (no centre nub, task #9)", "[ContinuousShell]") {
    // The inward march converges at a solid cell's core into degenerate rings - on an
    // elongated cell a flat medial SLIVER traversed down-and-back at effectively double
    // flow, which the weave scribbles into the centre-nub knot. Degenerate trailing rings
    // (effective bead width < 1 bead or perimeter < 8 beads) are now dropped. Signature of
    // the sliver: pairs of ANTIPARALLEL segments closer than half a bead (down-and-back);
    // kiss contacts of healthy rings are parallel-same-direction and do not count.
    Polygon rect{ {-15000000, -7500000}, {15000000, -7500000}, {15000000, 7500000}, {-15000000, 7500000} };
    Surfaces surfaces{ Surface(stInternal, ExPolygon(rect)) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    cfg.continuous_shell_infill_weave.value   = 100;
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 4));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    // central strip segments only (the sliver zone)
    struct Seg { Vec2d m, d; size_t idx; };
    std::vector<Seg> strip;
    const double nominal = flow.mm3_per_mm();
    size_t gidx = 0;
    for (const ExtrusionPath *p : paths) {
        const bool rider = p->mm3_per_mm < 0.7 * nominal;   // half-flow connectors retrace beads by design
        for (size_t i = 1; i < p->polyline.points.size(); ++i, ++gidx) {
            if (rider) continue;
            const Vec2d a = Point(p->polyline.points[i-1].x(), p->polyline.points[i-1].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d b = Point(p->polyline.points[i].x(),   p->polyline.points[i].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d m = 0.5 * (a + b);
            if (std::abs(m.y()) < 1.0 && std::abs(m.x()) < 10.0) {
                const Vec2d d = b - a;
                if (d.norm() > 0.05)
                    strip.push_back({ m, d / d.norm(), gidx });
            }
        }
    }
    // Count only pairs FAR APART along the path: a legitimate U-turn's flanks are antiparallel
    // and close but path-adjacent; the sliver's down-and-back pairs half the loop apart.
    int antiparallel_close = 0;
    for (size_t i = 0; i < strip.size(); ++i)
        for (size_t j = i + 1; j < strip.size(); ++j) {
            if (strip[j].idx - strip[i].idx <= 20) continue;
            if ((strip[i].m - strip[j].m).norm() < 0.21 && strip[i].d.dot(strip[j].d) < -0.7)
                ++antiparallel_close;
        }
    WARN("central-strip segments=" << strip.size() << " antiparallel close pairs=" << antiparallel_close);
    REQUIRE(antiparallel_close < 5);           // the down-and-back sliver is gone
}

TEST_CASE("Continuous shell cell split is a no-op on a convex section", "[ContinuousShell]") {
    // Self-gating guard: a convex outline has ~zero hull deviation, so cell_split must leave the
    // toolpath byte-identical - existing convex prints cannot regress when the option is on.
    const Point center(10000000, 10000000);
    Polygon circ;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circ.points.emplace_back(center.x() + coord_t(8000000 * std::cos(a)), center.y() + coord_t(8000000 * std::sin(a)));
    }
    Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
    Flow flow(0.5f, 0.25f, 0.5f);
    auto points = [&](bool split) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = 25;
        cfg.continuous_shell_infill_weave.value   = 100;
        cfg.continuous_shell_cell_split.value     = split;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 7));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        Points out;
        for (const ExtrusionPath *p : paths)
            for (auto &q : p->polyline.points)
                out.emplace_back(q.x(), q.y());
        return out;
    };
    const Points on = points(true), off = points(false);
    REQUIRE(on.size() == off.size());
    for (size_t i = 0; i < on.size(); ++i)
        REQUIRE(on[i] == off[i]);
}

TEST_CASE("Continuous shell keeps TWO bores hollow via Voronoi partition (multi-bore)", "[ContinuousShell]") {
    // Task #4: a cross-section with >1 hole must not smear fill across the bores. march_bands models a
    // single-hole annular wall, so a two-bore section is split into single-bore Voronoi cells (one per
    // hole centroid), each marched as an ordinary annulus. 40x24mm block, two 8x8mm square bores at
    // x=+-11: no bead may land inside, and no chord may cross, either bore void.
    ExPolygon block;
    block.contour = Polygon{ {-20000000,-12000000}, {20000000,-12000000},
                             {20000000, 12000000}, {-20000000, 12000000} };            // CCW outer
    block.holes.push_back(Polygon{ {-15000000,-4000000}, {-15000000, 4000000},
                                   {-7000000, 4000000}, {-7000000,-4000000} });         // CW bore 1
    block.holes.push_back(Polygon{ {7000000,-4000000}, {7000000, 4000000},
                                   {15000000, 4000000}, {15000000,-4000000} });         // CW bore 2
    REQUIRE(block.holes.size() == 2);

    Surfaces surfaces{ Surface(stInternal, block) };
    Flow flow(0.5f, 0.25f, 0.5f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    ExtrusionEntityCollection entities;
    std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, entities, &reason, false, 4);
    INFO(reason);
    REQUIRE(built);

    // Each bore void, shrunk 1mm on all sides so legit beads hugging the bore wall don't count.
    auto in_bore = [](const Point &p) {
        auto inrect = [&](coord_t x0, coord_t x1, coord_t y0, coord_t y1) {
            return p.x() > x0 && p.x() < x1 && p.y() > y0 && p.y() < y1; };
        return inrect(-14000000, -8000000, -3000000, 3000000)   // bore 1 interior
            || inrect( 8000000, 14000000, -3000000, 3000000);   // bore 2 interior
    };
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *e : entities.entities)
        ContinuousShellPathBuilder::collect_paths(e, paths);
    REQUIRE_FALSE(paths.empty());
    size_t bore_pts = 0, bore_mids = 0, total = 0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 0; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i].x(), p->polyline.points[i].y());
            ++total;
            if (in_bore(a)) ++bore_pts;
            if (i > 0) {
                const Point b(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
                if (in_bore(Point((a + b) / 2))) ++bore_mids;
            }
        }
    INFO("total points: " << total << " bore_pts: " << bore_pts << " bore_mids: " << bore_mids);
    CHECK(total > 200);
    REQUIRE(bore_pts == 0);     // no fill deposited inside either bore
    REQUIRE(bore_mids == 0);    // no chord crosses either bore
}

TEST_CASE("Continuous shell bore bisectors get the modern interface (task #7)", "[ContinuousShell]") {
    // With cell_split on, bore-bisector cuts go through the same pipeline as concave cuts:
    // recorded as CutLines, cells RECESSED half a bead each side (no coincident 200%-flow wall
    // line), fielded by the interface weave (the facing bore-cell walls kiss/fuse/breathe), and
    // the seed's skin PEELED into one contiguous loop (no per-cell skin arcs with mouths on the
    // rim). Two-bore block: the single bisector is x=0.
    ExPolygon block;
    block.contour = Polygon{ {-20000000,-12000000}, {20000000,-12000000},
                             {20000000, 12000000}, {-20000000, 12000000} };
    block.holes.push_back(Polygon{ {-15000000,-4000000}, {-15000000, 4000000},
                                   {-7000000, 4000000}, {-7000000,-4000000} });
    block.holes.push_back(Polygon{ {7000000,-4000000}, {7000000, 4000000},
                                   {15000000, 4000000}, {15000000,-4000000} });
    Surfaces surfaces{ Surface(stInternal, block) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    cfg.continuous_shell_infill_density.value = 25;
    cfg.continuous_shell_infill_weave.value   = 100;
    cfg.continuous_shell_cell_split.value     = true;
    EnvFlag bubble("CS_BUBBLE", "1");
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, false, 0));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    // (a) bores stay hollow through the modern pipeline
    auto in_bore = [](const Vec2d &p) {
        auto inrect = [&](double x0, double x1, double y0, double y1) {
            return p.x() > x0 && p.x() < x1 && p.y() > y0 && p.y() < y1; };
        return inrect(-14, -8, -3, 3) || inrect(8, 14, -3, 3);
    };
    // (b) bisector gap breathes: per-slab facing-material separation across x=0
    std::array<double, 5> lsum{}, lcnt{}, rsum{}, rcnt{};
    size_t bore_hits = 0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Vec2d a = Point(p->polyline.points[i-1].x(), p->polyline.points[i-1].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d b = Point(p->polyline.points[i].x(),   p->polyline.points[i].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d m = 0.5 * (a + b);
            if (in_bore(m)) ++bore_hits;
            if (std::abs(m.x()) > 1.2) continue;
            if (m.y() < -7.5 || m.y() >= 7.5) continue;
            const size_t slab = size_t((m.y() + 7.5) / 3.0);
            if (slab >= 5) continue;
            if (m.x() < 0) { lsum[slab] += m.x(); lcnt[slab] += 1; }
            else           { rsum[slab] += m.x(); rcnt[slab] += 1; }
        }
    REQUIRE(bore_hits == 0);
    // NOTE: with virtual-site placement the two bore cells may no longer share the x=0
    // bisector (a compact centre cell can sit between them), so the slab-gap numbers are
    // observational only - the recess/field machinery itself is shared with concave cuts
    // and hard-asserted by the dogbone interface test.
    double gmin = 1e9, gmax = -1e9, gsum = 0; int gn = 0;
    for (size_t k = 0; k < 5; ++k) {
        if (lcnt[k] < 5 || rcnt[k] < 5) continue;
        const double gap = rsum[k] / rcnt[k] - lsum[k] / lcnt[k];
        gmin = std::min(gmin, gap); gmax = std::max(gmax, gap); gsum += gap; ++gn;
    }
    // (c) the skin peels into ONE contiguous tour stretch
    auto dist_to_boundary = [&](const Vec2d &p) {
        double best = std::numeric_limits<double>::max();
        const Polygon &c = block.contour;
        for (size_t i = 0; i < c.points.size(); ++i) {
            const Vec2d a  = c.points[i].cast<double>() * SCALING_FACTOR;
            const Vec2d ab = c.points[(i + 1) % c.points.size()].cast<double>() * SCALING_FACTOR - a;
            const double t = std::clamp((p - a).dot(ab) / std::max(1e-9, ab.squaredNorm()), 0.0, 1.0);
            best = std::min(best, (p - (a + t * ab)).norm());
        }
        return best;
    };
    double skin_len = 0, tour_arc = 0, max_gap = 0, last_skin = -1.0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Vec2d a = Point(p->polyline.points[i-1].x(), p->polyline.points[i-1].y()).cast<double>() * SCALING_FACTOR;
            const Vec2d b = Point(p->polyline.points[i].x(),   p->polyline.points[i].y()).cast<double>() * SCALING_FACTOR;
            const double seg = (b - a).norm();
            tour_arc += seg;
            if (std::abs(dist_to_boundary(0.5 * (a + b)) - 0.21) < 0.12) {
                skin_len += seg;
                if (last_skin >= 0.0) max_gap = std::max(max_gap, tour_arc - last_skin);
                last_skin = tour_arc;
            }
        }
    WARN("x=0 slabs (observational): n=" << gn << (gn > 0 ? " mean=" + std::to_string(gsum / std::max(gn, 1)) : std::string())
         << "; skin len=" << skin_len << "mm max tour gap=" << max_gap << "mm");
    REQUIRE(skin_len > 100.0);                 // ~124mm rect perimeter, mostly covered
    REQUIRE(max_gap < 15.0);                   // as ONE contiguous stretch (peel active with bores)
}

// Build a gear-like annulus centred at (cx, cy): a SHARP-toothed (deep, flat-crested teeth) or
// smooth outer ring around a keyed (round + rectangular keyway) or plain round bore. The sharp
// teeth + keyed bore make the inward march PINCH and split - the topology that stranded the wall
// into gap-carrying fragments and produced the full-width top-face chords.
static ExPolygon make_gear_annulus(double root_r_mm, double tip_r_mm, int teeth,
                                   double bore_r_mm, bool toothed, bool keyed,
                                   coord_t cx = 8000000, coord_t cy = 8000000)
{
    const int n_outer = std::max(240, teeth * 16);
    Points outer; outer.reserve(n_outer);
    for (int i = 0; i < n_outer; ++i) {
        const double t     = 2.0 * M_PI * double(i) / double(n_outer);
        const double phase = std::fmod(double(teeth) * t, 2.0 * M_PI);   // 0..2pi within one tooth
        // Flat-crested trapezoidal tooth: outer third of the tooth at tip radius, rest at root.
        double r = root_r_mm;
        if (toothed)
            r = (phase > 2.0 && phase < 4.28) ? tip_r_mm : root_r_mm;    // ~raised crest, sharp flanks
        outer.emplace_back(Point(cx + coord_t(std::llround(r * 1e6 * std::cos(t))),
                                 cy + coord_t(std::llround(r * 1e6 * std::sin(t)))));
    }
    Points bore; const int n_bore = 96; bore.reserve(n_bore);
    for (int i = 0; i < n_bore; ++i) {
        const double t = 2.0 * M_PI * double(i) / double(n_bore);
        double rx = bore_r_mm * std::cos(t), ry = bore_r_mm * std::sin(t);
        if (keyed && std::cos(t) > 0.80)                                  // rectangular keyway notch on +x
            rx = (bore_r_mm + 0.9) * std::cos(t);
        bore.emplace_back(Point(cx + coord_t(std::llround(rx * 1e6)), cy + coord_t(std::llround(ry * 1e6))));
    }
    Polygon oc(outer); oc.make_counter_clockwise();
    Polygon bc(bore);  bc.make_clockwise();          // hole winding per ExPolygon convention
    return ExPolygon(oc, bc);
}

TEST_CASE("Continuous shell keyed toothed annulus has no full-width surface->bore chord", "[ContinuousShell]") {
    // Regression for the gear top-face bulge (cs-surface-quality). On a sharp-toothed annulus with
    // an eccentric keyed bore the inward march PINCHES and splits: tooth tips pinch off as solid
    // fragments and the wall fragments into an outer band + a separately-marched interior. Each
    // annular fragment then holds only part of the wall, so its loop list has a radial GAP that the
    // Fermat interleave used to chord straight across - two full-width surface->bore beads per layer
    // that pile onto the seam and bulge the exposed top face. contour_spiral now merges a single-
    // bore seed's annular fragments back into one contiguous chain, so NO single extruding segment
    // spans most of the wall. Uses the real print path: force_solid=false, default infill (0).
    const coord_t cx = 8000000, cy = 8000000;
    const Point center(cx, cy);
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;                            // default continuous_shell_infill_density = 0 (as the model)

    auto max_wall_chord = [&](bool toothed, bool keyed, size_t &count_out, double &worst_span_out) {
        ExPolygon gear = make_gear_annulus(8.0, 8.9, 22, 3.6, toothed, keyed, cx, cy);
        const double root_r = 8.0e6, bore_r = 3.6e6, wall = root_r - bore_r;
        Surfaces surfaces{ Surface(stInternal, gear) };
        size_t count = 0; double worst = 0.0;
        for (size_t layer = 3; layer < 7; ++layer) {   // a few layers (seam precesses); worst-case over them
            ExtrusionEntityCollection e; std::string reason;
            REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason,
                                                                   (layer & 1) != 0, layer, /*force_solid=*/false));
            std::vector<ExtrusionPath *> paths;
            for (ExtrusionEntity *ent : e.entities)
                ContinuousShellPathBuilder::collect_paths(ent, paths);
            REQUIRE_FALSE(paths.empty());
            for (const ExtrusionPath *p : paths)
                for (size_t i = 1; i < p->polyline.points.size(); ++i) {
                    const Point a(p->polyline.points[i-1].x(), p->polyline.points[i-1].y());
                    const Point b(p->polyline.points[i].x(),   p->polyline.points[i].y());
                    const double span = std::abs((a - center).cast<double>().norm() - (b - center).cast<double>().norm());
                    if (span > worst) worst = span;
                    if (span > 0.6 * wall) ++count;     // a single seg crossing >60% of the wall
                }
        }
        count_out = count; worst_span_out = worst / 1e6;
        return worst / wall;
    };

    size_t round_chords = 0, keyed_chords = 0; double round_span = 0, keyed_span = 0;
    const double round_frac = max_wall_chord(false, false, round_chords, round_span);
    const double keyed_frac = max_wall_chord(true,  true,  keyed_chords, keyed_span);
    WARN("SMOOTH round annulus:    full-width chords=" << round_chords << " worst radial span=" << round_span << "mm (" << round_frac << " of wall)");
    WARN("SHARP keyed-bore gear:   full-width chords=" << keyed_chords << " worst radial span=" << keyed_span << "mm (" << keyed_frac << " of wall)");

    // The fix target: neither topology may lay a single segment that crosses most of the wall.
    CHECK(round_chords == 0u);
    REQUIRE(keyed_chords == 0u);
}

// Replay a REAL sliced herringbone-gear cross-section (baked from CS_DUMP_REGION, Sun.stl layer 23,
// the layer whose inward march pinches/splits). The synthetic gears above cannot reproduce this
// topology - only the real contour does - so this is the definitive guard. The gear prints with
// continuous_shell_infill_density defaulting to 0, so force_solid never fires and every layer takes
// the plain concentric path. Pre-fix this laid two ~6.5mm full-width surface->bore chords onto the
// seam (the top-face bulge); the annular-fragment merge in contour_spiral removes them.
TEST_CASE("Continuous shell replays a real gear slice with no full-width wall chord", "[ContinuousShell]") {
    const ExPolygon gear = cs_test_gear_layer();
    REQUIRE(gear.holes.size() == 1);
    const Point center = gear.contour.bounding_box().center();
    double rc_min = 1e18, rh_max = 0;
    for (const Point &p : gear.contour.points) rc_min = std::min(rc_min, (p - center).cast<double>().norm());
    for (const Point &p : gear.holes.front().points) rh_max = std::max(rh_max, (p - center).cast<double>().norm());
    const double wall = rc_min - rh_max;              // root-to-bore radial thickness
    REQUIRE(wall > 3e6);                              // a genuinely thick wall (several beads)

    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;                            // default continuous_shell_infill_density = 0 (as the model)
    cfg.continuous_shell_seam_coast.value = 90;
    cfg.continuous_shell_seam_step.value  = 1600;

    Surfaces surfaces{ Surface(stInternal, gear) };
    size_t count = 0; double worst = 0.0;
    for (size_t layer = 20; layer < 26; ++layer) {    // the seam precesses - check a span of layers
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason,
                                                               (layer & 1) != 0, layer, /*force_solid=*/false));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities) ContinuousShellPathBuilder::collect_paths(ent, paths);
        REQUIRE_FALSE(paths.empty());
        for (const ExtrusionPath *p : paths)
            for (size_t i = 1; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i-1].x(), p->polyline.points[i-1].y());
                const Point b(p->polyline.points[i].x(),   p->polyline.points[i].y());
                const double span = std::abs((a - center).cast<double>().norm() - (b - center).cast<double>().norm());
                if (span > worst) worst = span;
                if (span > 0.6 * wall) ++count;        // a single seg crossing >60% of the wall
            }
    }
    WARN("real gear layer: full-width chords=" << count << " worst radial span=" << worst / 1e6 << "mm (wall " << wall / 1e6 << "mm)");
    REQUIRE(count == 0u);
}

TEST_CASE("Continuous shell woven sparsifies the real (fragmented) gear wall", "[ContinuousShell]") {
    // The merge that removes the top-face chords ALSO makes the eccentric gear wall sparsifiable:
    // once the split fragments are stitched into one contiguous annular chain, the woven pass thins
    // its interior (task #8 - previously a ~1% no-op because the fragmented wall never formed a thick
    // concentric chain to thin). Guards that unlock: density 20% must cut the gear-wall bead length
    // substantially vs solid, while keeping the bore hollow. (This drives build_cross_section directly
    // with force_solid=false; the top/bottom-shell gating is cs_air_within's job, tested via slices.)
    const ExPolygon gear = cs_test_gear_layer();
    const Point center = gear.contour.bounding_box().center();
    Surfaces surfaces{ Surface(stInternal, gear) };
    Flow flow(0.42f, 0.2f, 0.4f);
    auto measure = [&](int density, double &len_mm, double &bore_min_mm) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_infill_density.value       = density;
        cfg.continuous_shell_infill_weave.value         = 100;
        cfg.continuous_shell_infill_weave_shift.value   = 30;
        double L = 0.0, rmin = 1e18;
        for (size_t layer = 20; layer < 26; ++layer) {
            ExtrusionEntityCollection e; std::string reason;
            REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, (layer & 1) != 0, layer, false));
            std::vector<ExtrusionPath *> paths;
            for (ExtrusionEntity *ent : e.entities) ContinuousShellPathBuilder::collect_paths(ent, paths);
            for (const ExtrusionPath *p : paths)
                for (size_t i = 0; i < p->polyline.points.size(); ++i) {
                    const Point a(p->polyline.points[i].x(), p->polyline.points[i].y());
                    rmin = std::min(rmin, (a - center).cast<double>().norm());
                    if (i > 0) {
                        const Point b(p->polyline.points[i-1].x(), p->polyline.points[i-1].y());
                        L += (a - b).cast<double>().norm();
                    }
                }
        }
        len_mm = L / 1e6; bore_min_mm = rmin / 1e6;
    };
    double solid_len = 0, woven_len = 0, solid_bore = 0, woven_bore = 0;
    measure(0,  solid_len, solid_bore);
    measure(20, woven_len, woven_bore);
    WARN("gear wall (6 layers): solid=" << solid_len << "mm woven20=" << woven_len
         << "mm saved=" << (100.0 * (1.0 - woven_len / solid_len)) << "%  bore_min solid=" << solid_bore
         << " woven=" << woven_bore << "mm");
    REQUIRE(woven_len < 0.75 * solid_len);        // woven thins the interior by well over the previous ~1% no-op
    REQUIRE(woven_bore > solid_bore - 0.25);      // the weave stays as clear of the (eccentric keyed) bore as the solid fill
}

TEST_CASE("Continuous shell woven infill knits on a SQUARE (arc-length + normal weave)", "[ContinuousShell]") {
    // The woven perturbation must INTERLOCK adjacent infill rings on a faceted cross-section, not
    // only a circle. A raw sin(lobes*theta) radial push bunches the lobes at a square's corners and
    // shoves them diagonally, so anti-phase neighbours stay a full stride apart (no knit, ~2.9mm).
    // Arc-length + normal perturbation spaces the lobes evenly along the perimeter and perpendicular
    // to the wall, so neighbours close to ~one bead at the lobe peaks. Assert that knit.
    const coord_t C = 15000000, HW = 15000000;    // 30mm square centred at (15,15)mm
    ExPolygon sq(Polygon{ {C-HW,C-HW},{C+HW,C-HW},{C+HW,C+HW},{C-HW,C+HW} });
    Surfaces surfaces{ Surface(stInternal, sq) };
    Flow flow(0.42f, 0.2f, 0.4f);
    PrintRegionConfig cfg;
    cfg.continuous_shell_infill_density.value     = 15;
    cfg.continuous_shell_infill_weave.value       = 100;
    cfg.continuous_shell_infill_weave_shift.value = 50;
    const Point  center(C, C);
    const double spacing_mm = double(flow.scaled_spacing()) / 1e6;
    double best_gap = 1e18;                        // min adjacent-ring radial gap over angle bins & layers
    for (size_t layer = 20; layer < 26; ++layer) {
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.2, &cfg, e, &reason, (layer & 1) != 0, layer, false));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities) ContinuousShellPathBuilder::collect_paths(ent, paths);
        std::map<int, std::vector<double>> bins;   // 5-degree angle bins -> infill-bead radii (mm)
        for (const ExtrusionPath *p : paths)
            for (const auto &q : p->polyline.points) {
                const Vec2d d = Point(q.x(), q.y()).cast<double>() - center.cast<double>();
                const double r = d.norm() / 1e6;
                if (r < 2.0 || r > 12.5) continue;                    // skip the centre core and the outer walls
                bins[int(std::floor(std::atan2(d.y(), d.x()) / M_PI * 36.0))].push_back(r);
            }
        for (auto &kv : bins) {
            std::sort(kv.second.begin(), kv.second.end());
            for (size_t i = 1; i < kv.second.size(); ++i) {
                const double g = kv.second[i] - kv.second[i - 1];
                if (g > 0.2 && g < 4.0) best_gap = std::min(best_gap, g);   // ignore duplicates and the seam jump
            }
        }
    }
    WARN("square woven min adjacent-ring gap=" << best_gap << "mm (full stride ~" << 7 * spacing_mm
         << "mm; a knit closes toward one bead ~" << spacing_mm << "mm)");
    REQUIRE(best_gap < 2.0 * spacing_mm);          // rings knit to within ~2 beads (pre-fix they stayed a full stride apart)
}

TEST_CASE("Continuous shell solid multibore layer leaves no see-through coverage gap", "[ContinuousShell]") {
    // Regression for the bottom-face pinholes (die-5 centre pip, odd-bore print): the concentric
    // multibore fill undertiles at ring pinches and bore crescents, leaving uncovered wedges you
    // can see through into the infill. The coverage-repair pass rasterizes what the tour laid and
    // patches interior gaps; assert its invariant from the EMITTED paths - no uncovered interior
    // island big enough to read as a hole - and that the clipped patches never enter a bore.
    const Point c(10000000, 10000000);
    auto circle = [&](double cx_mm, double cy_mm, double r_mm, int n) {
        Polygon p;
        for (int i = 0; i < n; ++i) {
            const double a = 2.0 * M_PI * i / n;
            p.points.emplace_back(c.x() + coord_t(std::llround(scale_(cx_mm + r_mm * std::cos(a)))),
                                  c.y() + coord_t(std::llround(scale_(cy_mm + r_mm * std::sin(a)))));
        }
        return p;
    };
    // Disc R10 with a centre bore and two side bores leaving ~1mm pinch lanes - the layout
    // that produced the visible gaps around the die-5 centre pip.
    Polygon outer = circle(0, 0, 10.0, 96);  outer.make_counter_clockwise();
    Polygon b0    = circle(0, 0, 1.5, 48);   b0.make_clockwise();
    Polygon b1    = circle(-4.5, 0, 2.0, 48); b1.make_clockwise();
    Polygon b2    = circle( 4.5, 0, 2.0, 48); b2.make_clockwise();
    ExPolygon region(outer);
    region.holes = { b0, b1, b2 };
    Surfaces surfaces{ Surface(stInternal, region) };
    Flow flow(0.5f, 0.25f, 0.5f);
    PrintRegionConfig cfg;
    cfg.wall_loops.value = 2;
    const float w = float(flow.scaled_width());
    const struct { double x, y, r; } bores[] = { {0, 0, 1.5}, {-4.5, 0, 2.0}, {4.5, 0, 2.0} };

    for (size_t layer = 0; layer < 3; ++layer) {
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason,
                                                               (layer & 1) != 0, layer, /*force_solid=*/true));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities) ContinuousShellPathBuilder::collect_paths(ent, paths);
        REQUIRE_FALSE(paths.empty());
        Polylines laid;
        size_t in_bore = 0;
        for (const ExtrusionPath *p : paths) {
            Polyline pl;
            for (const auto &q : p->polyline.points) pl.points.emplace_back(q.x(), q.y());
            if (pl.points.size() >= 2) laid.push_back(std::move(pl));
            for (size_t i = 1; i < p->polyline.points.size(); ++i) {
                const Point a(p->polyline.points[i-1].x(), p->polyline.points[i-1].y());
                const Point b(p->polyline.points[i].x(),   p->polyline.points[i].y());
                for (const Point &probe : { Point((a + b) / 2), b })
                    for (const auto &bo : bores) {
                        const Vec2d d = probe.cast<double>() - Vec2d(c.x() + scale_(bo.x), c.y() + scale_(bo.y));
                        if (d.norm() < scale_(bo.r - 0.3)) ++in_bore;
                    }
            }
        }
        // Same rasterization the repair pass uses: bead footprint, opened to kill hairlines.
        const Polygons cov = offset(laid, 0.55f * w);
        ExPolygons uncovered = offset2_ex(diff_ex(ExPolygons{ region }, union_ex(cov)), -0.35f * w, 0.35f * w);
        double worst = 0.0;
        for (const ExPolygon &gap : uncovered) worst = std::max(worst, std::abs(gap.contour.area()));
        WARN("L" << layer << ": uncovered islands=" << uncovered.size()
             << " worst=" << worst * SCALING_FACTOR * SCALING_FACTOR << "mm^2 in_bore=" << in_bore);
        REQUIRE(worst < scale_(1.2) * scale_(1.2));   // nothing left big enough to read as a pinhole
        REQUIRE(in_bore == 0);                        // patches are clipped to the region - bores stay sacred
    }
}
