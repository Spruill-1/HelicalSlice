#include <catch2/catch_all.hpp>

#include "libslic3r/ContinuousShell.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/PrintConfig.hpp"

#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <array>

using namespace Slic3r;

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
    Flow flow(0.5f, 0.25f, 0.5f);

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
    const double v_model = std::abs(annulus.expolygon.area()) * SCALING_FACTOR * SCALING_FACTOR * layer_height;

    // Total extruded volume (mm^3) and the largest per-bead flow MULTIPLIER, for a given knob.
    // The multiplier is measured against each path's OWN width baseline (adaptive bead width
    // legitimately raises mm3_per_mm on widened medial beads at multiplier 1.0).
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
            max_mult = std::max(max_mult, p.mm3_per_mm / flow.with_width(p.width).mm3_per_mm());
        }
        return vol;
    };

    double m0 = 0.0, m100 = 0.0;
    const double v_off = measure(0.0,   m0);    // nominal baseline
    const double v_on  = measure(100.0, m100);  // fill 100% of the residual void

    REQUIRE(v_off > 0.0);
    REQUIRE(m0 > 0.999);  REQUIRE(m0 < 1.001);        // 0% -> every bead at nominal flow
    REQUIRE(v_on >= v_off - 1e-9);                    // over-extrusion only ever ADDS material
    REQUIRE(v_on <= std::max(v_off, v_model) * 1.02); // and never past the solid model volume -> no bulge
    REQUIRE(m100 <= 1.5 + 1e-6);                      // per-bead safety cap respected
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

TEST_CASE("Continuous shell seam coast thins flow at the surface handoff (A3)", "[ContinuousShell]") {
    // seam_coast thins the extrusion across the surface<->interior handoffs to soften the
    // seam groove, while leaving the surface and interior beads at full flow.
    Surface annulus(
        stInternal,
        ExPolygon(
            Polygon{{0, 0}, {8000000, 0}, {8000000, 8000000}, {0, 8000000}},
            Polygon{{3000000, 3000000}, {3000000, 5000000}, {5000000, 5000000}, {5000000, 3000000}}));
    Surfaces surfaces{ annulus };
    Flow flow(0.5f, 0.25f, 0.5f);

    // Flow MULTIPLIER extremes against each path's own width baseline (adaptive bead width
    // legitimately raises mm3_per_mm on widened medial beads at multiplier 1.0).
    auto flow_extremes = [&](double coast_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_seam_coast.value = coast_pct;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false));
        auto *mp = dynamic_cast<ExtrusionMultiPath *>(e.entities.front());
        REQUIRE(mp != nullptr);
        double mn = std::numeric_limits<double>::max(), mx = 0.0;
        for (auto &p : mp->paths) {
            const double mult = p.mm3_per_mm / flow.with_width(p.width).mm3_per_mm();
            mn = std::min(mn, mult); mx = std::max(mx, mult);
        }
        return std::make_pair(mn, mx);
    };

    const auto off = flow_extremes(0.0);
    REQUIRE(off.first  > 0.999);          // no coast -> no flow modifier anywhere
    REQUIRE(off.second < 1.001);

    const auto on = flow_extremes(50.0);  // 50% coast
    REQUIRE(on.second  > 0.999);          // surface + interior beads keep full flow
    REQUIRE(on.second  < 1.001);
    REQUIRE(on.first   < 0.6);            // the seam crossings are thinned to ~0.5x
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

namespace {
// A deep crescent: circle r8 at the origin minus a r6 bite centred at (5.5, 0). The bite
// reaches past the part centre, so the outline FOLDS relative to any interior point - the
// canonical non-star-shaped section (the angle about the centroid backtracks along the
// bite arc near the horns). Used by the star-shape-gate regressions below.
ExPolygon make_deep_crescent() {
    Polygon circle, bite;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circle.points.emplace_back(coord_t(8000000 * std::cos(a)), coord_t(8000000 * std::sin(a)));
        bite.points.emplace_back(coord_t(5500000 + 6000000 * std::cos(a)), coord_t(6000000 * std::sin(a)));
    }
    ExPolygons cres = diff_ex(Polygons{circle}, Polygons{bite});
    REQUIRE(cres.size() == 1);
    REQUIRE(cres.front().holes.empty());
    return cres.front();
}

std::vector<Point3> cs_path_points(const Surfaces &surfaces, const Flow &flow, PrintRegionConfig &cfg,
                                   size_t layer_index, bool force_solid) {
    ExtrusionEntityCollection e; std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason,
                                                                       false, layer_index, force_solid);
    INFO(reason);
    REQUIRE(built);
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    std::vector<Point3> pts;
    for (const ExtrusionPath *p : paths)
        for (const Point3 &q : p->polyline.points) pts.push_back(q);
    REQUIRE(pts.size() > 8);
    return pts;
}
} // namespace

TEST_CASE("Continuous shell closes the medial deficit by widening the innermost beads", "[ContinuousShell]") {
    // Adaptive bead width (C-lite). A 2.193mm annular wall fits two band pairs (0.5mm beads
    // at ~0.446 spacing) and leaves a ~0.3mm residual void between the innermost bead edges -
    // the old code patched it with heuristic medial strips or left it open. Now the deficit
    // is measured (vanish-depth bisection) and closed by widening the innermost BURIED ring
    // pair by v/2 each, centrelines shifted v/4 toward the medial. The free surfaces must
    // keep nominal width, every sub-path's flow must be consistent with its declared width
    // (the preview contract), and the deposited volume must come out near the model's own.
    Polygon outer, inner;
    for (int i = 0; i < 240; ++i) {
        const double a = 2.0 * M_PI * i / 240.0;
        outer.points.emplace_back(coord_t(10000000 * std::cos(a)), coord_t(10000000 * std::sin(a)));
        inner.points.emplace_back(coord_t( 7807000 * std::cos(a)), coord_t( 7807000 * std::sin(a)));
    }
    ExPolygons annulus = diff_ex(Polygons{outer}, Polygons{inner});
    REQUIRE(annulus.size() == 1);
    REQUIRE(annulus.front().holes.size() == 1);
    const double model_area_mm2 = annulus.front().area() * SCALING_FACTOR * SCALING_FACTOR;

    Surfaces surfaces{ Surface(stInternal, annulus.front()) };
    Flow flow(0.5f, 0.25f, 0.5f);
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, e, &reason, false, 3));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    REQUIRE_FALSE(paths.empty());

    float  max_w = 0.f, first_w = paths.front()->width;
    double volume = 0.0;
    for (const ExtrusionPath *p : paths) {
        max_w = std::max(max_w, p->width);
        // Preview/flow contract: mm3_per_mm must be exactly what a Flow of this width
        // gives (no flow modifiers are configured here).
        REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i-1].x(), p->polyline.points[i-1].y());
            const Point b(p->polyline.points[i].x(),   p->polyline.points[i].y());
            volume += (b - a).cast<double>().norm() * SCALING_FACTOR * p->mm3_per_mm;
        }
    }
    // The path starts on the outer surface ring: it must stay at nominal width.
    REQUIRE_THAT(double(first_w), Catch::Matchers::WithinAbs(0.5, 0.001));
    // The innermost pair was widened by ~v/2 = ~0.15mm.
    REQUIRE(max_w > 0.56f);
    REQUIRE(max_w < 0.75f);
    // Total deposited volume is close to the model's solid volume (the old gap lost ~20%).
    const double model_volume = model_area_mm2 * 0.25;
    REQUIRE(volume > 0.90 * model_volume);
    REQUIRE(volume < 1.08 * model_volume);
}

TEST_CASE("Continuous shell scarf tapers junctions, spares surfaces, bounds sub-paths", "[ContinuousShell]") {
    // Option A: the in-plane scarf tapers bead width toward a floor across every ring-to-ring
    // junction so the once-per-revolution radial step blends instead of ridging. Free-surface
    // rings are never tapered (their outer flank would recede into a relief groove). Tapered
    // widths are quantized to 0.05mm buckets so the taper stays a handful of sub-paths.
    Polygon circ;
    for (int i = 0; i < 96; ++i) {
        const double a = 2.0 * M_PI * i / 96.0;
        circ.points.emplace_back(coord_t(8000000 * std::cos(a)), coord_t(8000000 * std::sin(a)));
    }
    Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
    Flow flow(0.5f, 0.25f, 0.5f);

    struct R { double vol = 0.; float min_w = 1e9f, first_w = 0.f; double first_len = 0.; size_t n_paths = 0; };
    auto measure = [&](double scarf_mm) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_scarf_length.value = scarf_mm;
        ExtrusionEntityCollection e; std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason, false, 3));
        std::vector<ExtrusionPath *> paths;
        for (ExtrusionEntity *ent : e.entities)
            ContinuousShellPathBuilder::collect_paths(ent, paths);
        REQUIRE_FALSE(paths.empty());
        R r;
        r.n_paths   = paths.size();
        r.first_w   = paths.front()->width;
        r.first_len = paths.front()->polyline.length() * SCALING_FACTOR;
        for (const ExtrusionPath *p : paths) {
            r.min_w = std::min(r.min_w, p->width);
            r.vol  += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
            // Width/flow consistency holds through the scarf (no flow modifiers configured).
            REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
        }
        return r;
    };

    const R off = measure(0.0);
    const R on  = measure(8.0);

    REQUIRE(off.min_w > 0.49f);                        // scarf off: nothing below nominal width
    REQUIRE(on.min_w  > 0.20f);                        // scarf on: tapered down to the floor...
    REQUIRE(on.min_w  < 0.31f);                        // ...which is max(w/2, nozzle/4) = 0.25
    // The outer surface ring is spared: the path starts with the full ring 0 loop at
    // nominal width (circumference of the r = 8 - w/2 centreline).
    REQUIRE_THAT(double(on.first_w), Catch::Matchers::WithinAbs(0.5, 0.001));
    REQUIRE(on.first_len > 0.9 * 2.0 * M_PI * 7.75);
    // The taper deficit is local and bounded; the buckets keep the sub-path count sane.
    REQUIRE(on.vol < off.vol);
    REQUIRE(on.vol > 0.85 * off.vol);
    REQUIRE(on.n_paths < 500);
    REQUIRE(on.n_paths > off.n_paths);                 // it did split the taper into runs
}

TEST_CASE("Continuous shell fits the bead to a sub-bead wall and declines below min feature size", "[ContinuousShell]") {
    // Adaptive bead width on the sub-bead medial fallback: a 0.3mm wall used to print at
    // the nominal 0.5mm width (66% over-extrusion); now it prints at the narrowest
    // printable bead (default 85% of the 0.5 nozzle = 0.425). A 0.1mm wall is below the
    // minimum feature size (default 25% of nozzle = 0.125) and must decline loudly.
    Flow flow(0.5f, 0.25f, 0.5f);
    auto thin_annulus = [](coord_t r_out, coord_t r_in) {
        Polygon o, i2;
        for (int i = 0; i < 180; ++i) {
            const double a = 2.0 * M_PI * i / 180.0;
            o.points.emplace_back(coord_t(r_out * std::cos(a)), coord_t(r_out * std::sin(a)));
            i2.points.emplace_back(coord_t(r_in  * std::cos(a)), coord_t(r_in  * std::sin(a)));
        }
        ExPolygons ann = diff_ex(Polygons{o}, Polygons{i2});
        REQUIRE(ann.size() == 1);
        return ann.front();
    };

    Surfaces s03{ Surface(stInternal, thin_annulus(5000000, 4700000)) };   // 0.3mm wall
    ExtrusionEntityCollection e1; std::string r1;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(s03, flow, 0.2, nullptr, e1, &r1, false, 2));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e1.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    REQUIRE_FALSE(paths.empty());
    for (const ExtrusionPath *p : paths) {
        REQUIRE_THAT(double(p->width), Catch::Matchers::WithinAbs(0.425, 0.01));   // fitted, not nominal
        REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
    }

    Surfaces s01{ Surface(stInternal, thin_annulus(5000000, 4900000)) };   // 0.1mm wall
    ExtrusionEntityCollection e2; std::string r2;
    REQUIRE_FALSE(ContinuousShellPathBuilder::build_cross_section(s01, flow, 0.2, nullptr, e2, &r2, false, 2));
    REQUIRE_FALSE(r2.empty());
}

TEST_CASE("Continuous shell sub-bead medial fallback: closed ring works, open blade declines", "[ContinuousShell]") {
    // A region thinner than one bead falls back to a medial-axis centreline. That is only
    // honest when the medial is a genuinely CLOSED loop (a thin tube wall): force-closing an
    // OPEN centreline would fabricate an extruding chord from its end straight back to its
    // start, through air. The open case must decline loudly, not print a phantom chord.
    Flow flow(0.5f, 0.25f, 0.5f);

    // Closed control: a 0.3mm-thick full annulus (thinner than the 0.5mm bead), built via
    // diff_ex so the hole gets the proper (CW) winding.
    Polygon outer, inner;
    for (int i = 0; i < 180; ++i) {
        const double a = 2.0 * M_PI * i / 180.0;
        outer.points.emplace_back(coord_t(5000000 * std::cos(a)), coord_t(5000000 * std::sin(a)));
        inner.points.emplace_back(coord_t(4700000 * std::cos(a)), coord_t(4700000 * std::sin(a)));
    }
    ExPolygons annulus = diff_ex(Polygons{outer}, Polygons{inner});
    REQUIRE(annulus.size() == 1);
    REQUIRE(annulus.front().holes.size() == 1);
    Surfaces ring_s{ Surface(stInternal, annulus.front()) };
    ExtrusionEntityCollection e1; std::string reason1;
    const bool ring_built = ContinuousShellPathBuilder::build_cross_section(ring_s, flow, 0.2, nullptr, e1, &reason1, false, 2);
    INFO(reason1);
    REQUIRE(ring_built);
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e1.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    REQUIRE_FALSE(paths.empty());
    double min_r = std::numeric_limits<double>::max(), max_r = 0.0;
    for (const ExtrusionPath *p : paths)
        for (const Point3 &q : p->polyline.points) {
            const double r = Point(q.x(), q.y()).cast<double>().norm() * SCALING_FACTOR;
            min_r = std::min(min_r, r);
            max_r = std::max(max_r, r);
        }
    REQUIRE(min_r > 4.5);    // the bead rides the medial: never dives across the bore
    REQUIRE(max_r < 5.1);    // and never bulges outside the wall

    // Open case: the same thin wall with a 60-degree mouth cut out (a C-strip). Its medial
    // is an open arc; the old code closed it with a ~2.9mm chord across the mouth.
    Polygon cpoly;
    const double a0 = 0.5, a1 = 2.0 * M_PI - 0.5;
    for (int i = 0; i <= 150; ++i) {
        const double a = a0 + (a1 - a0) * i / 150.0;
        cpoly.points.emplace_back(coord_t(5000000 * std::cos(a)), coord_t(5000000 * std::sin(a)));
    }
    for (int i = 150; i >= 0; --i) {
        const double a = a0 + (a1 - a0) * i / 150.0;
        cpoly.points.emplace_back(coord_t(4700000 * std::cos(a)), coord_t(4700000 * std::sin(a)));
    }
    Surfaces cstrip_s{ Surface(stInternal, ExPolygon(cpoly)) };
    ExtrusionEntityCollection e2; std::string reason2;
    const bool blade_built = ContinuousShellPathBuilder::build_cross_section(cstrip_s, flow, 0.2, nullptr, e2, &reason2, false, 2);
    REQUIRE_FALSE(blade_built);
    REQUIRE_FALSE(reason2.empty());
}

TEST_CASE("Continuous shell samples the free surface finely on a large section", "[ContinuousShell]") {
    // The resampled ring points ARE the toolpath - z_contoured beads are exempt from
    // downstream simplification - so the old fixed 256-point cap quantized a large part's
    // outer surface to ~1mm chords (and put polygon corners up to ~0.8mm off). Ring
    // sampling is now per-ring, targeting a half-bead chord on the free surfaces.
    // A r=30mm disc (188mm circumference) must get outer-surface segments well under
    // 0.5mm; the fixed cap would have produced ~0.74mm.
    Polygon circ;
    for (int i = 0; i < 720; ++i) {                         // fine input: not the limiting factor
        const double a = 2.0 * M_PI * i / 720.0;
        circ.points.emplace_back(coord_t(30000000 * std::cos(a)), coord_t(30000000 * std::sin(a)));
    }
    Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
    Flow flow(0.5f, 0.25f, 0.5f);
    ExtrusionEntityCollection e; std::string reason;
    REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, nullptr, e, &reason, false, 3));
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    REQUIRE_FALSE(paths.empty());
    // Longest segment whose BOTH ends lie within one bead of the outer surface (the
    // radial band excludes the seam hop down to the second ring).
    const double band = 29500000.0;                        // r > 29.5mm: outer ring only
    double max_seg = 0.0;
    size_t surface_segs = 0;
    for (const ExtrusionPath *p : paths)
        for (size_t i = 1; i < p->polyline.points.size(); ++i) {
            const Point a(p->polyline.points[i - 1].x(), p->polyline.points[i - 1].y());
            const Point b(p->polyline.points[i].x(), p->polyline.points[i].y());
            if (a.cast<double>().norm() > band && b.cast<double>().norm() > band) {
                max_seg = std::max(max_seg, (b - a).cast<double>().norm());
                ++surface_segs;
            }
        }
    REQUIRE(surface_segs > 400);                            // the outer ring is densely sampled
    REQUIRE(max_seg * SCALING_FACTOR < 0.5);                // half-bead-chord target (was ~0.74mm)
}

TEST_CASE("Continuous shell declines the weave and Z-wave on a non-star (crescent) section", "[ContinuousShell]") {
    // The woven infill phases and clamps every ring by its angle about the section centroid,
    // which is only meaningful on star-shaped sections. On this crescent the 360-bin shell
    // profile would average the two boundary branches a ray crosses, and "clamping" to that
    // fiction can drag a bead half the section width off its ring. The star-shape gate must
    // therefore decline the radial weave AND the Z-wave here - while KEEPING the sparse
    // decimation, which needs no angular parametrization.
    Surfaces surfaces{ Surface(stInternal, make_deep_crescent()) };
    Flow flow(0.5f, 0.25f, 0.5f);
    auto measure = [&](int density, int weave) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = density;
        cfg.continuous_shell_infill_weave.value   = weave;
        return cs_path_points(surfaces, flow, cfg, /*layer_index=*/3, /*force_solid=*/false);
    };
    auto xy_len = [](const std::vector<Point3> &pts) {
        double len = 0.0;
        for (size_t i = 1; i < pts.size(); ++i)
            len += (Point(pts[i].x(), pts[i].y()) - Point(pts[i-1].x(), pts[i-1].y())).cast<double>().norm();
        return len * SCALING_FACTOR;
    };
    const std::vector<Point3> solid  = measure(0,  0);
    const std::vector<Point3> sparse = measure(20, 0);
    const std::vector<Point3> woven  = measure(20, 100);

    // Sparse decimation still applies on the folded section (it drops rings, no parametrization).
    REQUIRE(xy_len(sparse) < 0.85 * xy_len(solid));
    // The gate declines the weave outright: with it requested, the path is IDENTICAL to the
    // plain sparse path - no point moved, no Z perturbation added.
    REQUIRE(woven.size() == sparse.size());
    coord_t max_dev = 0, z_span = 0;
    for (size_t i = 0; i < woven.size(); ++i) {
        max_dev = std::max(max_dev, coord_t(std::abs(woven[i].x() - sparse[i].x())));
        max_dev = std::max(max_dev, coord_t(std::abs(woven[i].y() - sparse[i].y())));
        z_span  = std::max(z_span,  coord_t(std::abs(woven[i].z())));
    }
    REQUIRE(max_dev == 0);
    REQUIRE(z_span == 0);
}

TEST_CASE("Continuous shell declines the solid-cap stagger on a non-star (crescent) section", "[ContinuousShell]") {
    // The odd-layer brick stagger shifts interior rings radially toward the centroid, clamped
    // by the same angle-binned bore profile as the weave - equally fictional on a folded
    // section. The gate must decline it: an odd solid-cap layer of a woven-infill part prints
    // identically to the stagger-free baseline (density 0). The star-shaped control for the
    // stagger actually MOVING the interior is the multi-layer solid cap test above.
    Surfaces surfaces{ Surface(stInternal, make_deep_crescent()) };
    Flow flow(0.5f, 0.25f, 0.5f);
    auto measure = [&](int density) {
        PrintRegionConfig cfg;
        cfg.wall_loops.value = 2;
        cfg.continuous_shell_infill_density.value = density;
        cfg.top_surface_pattern.value = ipConcentric;    // concentric cap path (serpentine declines anyway)
        return cs_path_points(surfaces, flow, cfg, /*layer_index=*/1, /*force_solid=*/true);
    };
    const std::vector<Point3> plain     = measure(0);    // stagger disabled by density=0
    const std::vector<Point3> staggered = measure(20);   // stagger requested on this odd layer
    REQUIRE(staggered.size() == plain.size());
    coord_t max_dev = 0;
    for (size_t i = 0; i < staggered.size(); ++i) {
        max_dev = std::max(max_dev, coord_t(std::abs(staggered[i].x() - plain[i].x())));
        max_dev = std::max(max_dev, coord_t(std::abs(staggered[i].y() - plain[i].y())));
    }
    REQUIRE(max_dev == 0);
}

// ---- Adaptive beading (C-full) ------------------------------------------------------------------

namespace {
// A circular annulus centred at the origin: outer radius r_out_mm, wall t_geom_mm thick.
ExPolygon make_annulus(double r_out_mm, double t_geom_mm, int samples = 240) {
    Polygon outer, inner;
    for (int i = 0; i < samples; ++i) {
        const double a = 2.0 * M_PI * i / double(samples);
        outer.points.emplace_back(coord_t(std::llround(scale_(r_out_mm) * std::cos(a))),
                                  coord_t(std::llround(scale_(r_out_mm) * std::sin(a))));
        inner.points.emplace_back(coord_t(std::llround(scale_(r_out_mm - t_geom_mm) * std::cos(a))),
                                  coord_t(std::llround(scale_(r_out_mm - t_geom_mm) * std::sin(a))));
    }
    ExPolygons ann = diff_ex(Polygons{outer}, Polygons{inner});
    REQUIRE(ann.size() == 1);
    REQUIRE(ann.front().holes.size() == 1);
    return ann.front();
}

// Build a cross-section and return its leaf paths (owned by `e`, so `e` must outlive use).
std::vector<ExtrusionPath *> cs_paths(const Surfaces &surfaces, const Flow &flow, PrintRegionConfig &cfg,
                                      ExtrusionEntityCollection &e, size_t layer_index = 3, bool force_solid = false) {
    std::string reason;
    const bool built = ContinuousShellPathBuilder::build_cross_section(surfaces, flow, 0.25, &cfg, e, &reason,
                                                                       false, layer_index, force_solid);
    INFO(reason);
    REQUIRE(built);
    std::vector<ExtrusionPath *> paths;
    for (ExtrusionEntity *ent : e.entities)
        ContinuousShellPathBuilder::collect_paths(ent, paths);
    REQUIRE_FALSE(paths.empty());
    return paths;
}

// Radius clusters of the path points inside a narrow angular wedge about the origin -
// the ring count seen by a ray at `azimuth`, with 0.15mm cluster tolerance.
std::vector<double> ray_ring_radii(const std::vector<ExtrusionPath *> &paths, double azimuth) {
    std::vector<double> radii;
    for (const ExtrusionPath *p : paths)
        for (const Point3 &q : p->polyline.points) {
            const double th = std::atan2(double(q.y()), double(q.x()));
            double d = th - azimuth;
            d = std::atan2(std::sin(d), std::cos(d));
            if (std::abs(d) < 0.05)
                radii.push_back(std::hypot(double(q.x()), double(q.y())) * SCALING_FACTOR);
        }
    std::sort(radii.begin(), radii.end());
    std::vector<double> clusters;
    for (const double r : radii)
        if (clusters.empty() || r - clusters.back() > 0.15)
            clusters.push_back(r);
        else
            clusters.back() = r;                       // keep the cluster's running max
    return clusters;
}

// Full geometric signature of a built cross-section: per point (x, y, z, width, flow).
struct CsSig { std::vector<std::tuple<coord_t, coord_t, coord_t, float, double>> pts; };
CsSig cs_signature(const Surfaces &surfaces, const Flow &flow, PrintRegionConfig &cfg,
                   size_t layer_index = 3, bool force_solid = false) {
    ExtrusionEntityCollection e;
    CsSig sig;
    for (const ExtrusionPath *p : cs_paths(surfaces, flow, cfg, e, layer_index, force_solid))
        for (const Point3 &q : p->polyline.points)
            sig.pts.emplace_back(q.x(), q.y(), q.z(), p->width, double(p->mm3_per_mm));
    return sig;
}
} // namespace

TEST_CASE("Adaptive beading plans a concentric annulus at nominal widths", "[ContinuousShell]") {
    // A wall that tiles exactly (6 beads: t_geom = 6*spacing + (width - spacing)) must plan
    // 6 rings, all at nominal width, with the surface rings pinned and volume conserved -
    // the degenerate uniform case C-lite validated, now produced by the plan.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const double t_geom = 6.0 * sp + (0.5 - sp);
    const ExPolygon ann = make_annulus(10.0, t_geom);
    const double model_volume = ann.area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;

    Surfaces surfaces{ Surface(stInternal, ann) };
    PrintRegionConfig cfg;
    cfg.continuous_shell_adaptive_beading.value = true;
    ExtrusionEntityCollection e;
    const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);

    double volume = 0.;
    for (const ExtrusionPath *p : paths) {
        // Exact tiling: every sub-path prints at nominal width (one 0.05 bucket).
        REQUIRE_THAT(double(p->width), Catch::Matchers::WithinAbs(0.5, 0.051));
        REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
        volume += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
    }
    REQUIRE_THAT(volume, Catch::Matchers::WithinRel(model_volume, 0.05));
    for (const double az : {0.3, 1.7, 2.9, -2.1, -0.8}) {
        const std::vector<double> rings = ray_ring_radii(paths, az);
        REQUIRE(rings.size() == 6);
        REQUIRE(rings.front() > 10.0 - t_geom);        // bore stays hollow
    }
}

TEST_CASE("Adaptive beading squeeze bias adds a ring early and never stretches past 10%", "[ContinuousShell]") {
    // The count law from the printed T4 microscopy: squeezed beads consolidate, stretched
    // beads void along the medial - so a wall 4.7 spacings thick plans FIVE rings (interior
    // squeezed below nominal) instead of four stretched ones, while a wall just past an
    // integer tiling (4.05 spacings) keeps FOUR near-nominal rings (a 1% stretch beats a
    // gratuitous squeeze). Surface rings stay pinned at nominal width in both.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();

    auto measure = [&](double t_sp_units) {
        const double t_geom = t_sp_units * sp + (0.5 - sp);
        const ExPolygon ann = make_annulus(10.0, t_geom);
        Surfaces surfaces{ Surface(stInternal, ann) };
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        ExtrusionEntityCollection e;
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);
        struct M { size_t rings; float min_w = 1e9f, max_w = 0.f; double vol = 0., model = 0.; float outer_w = 0.f, bore_w = 0.f; } m;
        m.model = ann.area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;
        m.rings = ray_ring_radii(paths, 0.3).size();
        const double r_outer_c = 10.0 - 0.25;              // outer surface centreline
        const double r_bore_c  = 10.0 - t_geom + 0.25;     // bore surface centreline
        for (const ExtrusionPath *p : paths) {
            m.min_w = std::min(m.min_w, p->width);
            m.max_w = std::max(m.max_w, p->width);
            m.vol  += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
            REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
            // Classify surface sub-paths by MEAN radius (runs share endpoints, so a
            // surface ring's run may carry one junction point from its neighbour).
            if (p->polyline.points.size() >= 16) {
                double mean_r = 0.;
                for (const Point3 &q : p->polyline.points)
                    mean_r += std::hypot(double(q.x()), double(q.y())) * SCALING_FACTOR;
                mean_r /= double(p->polyline.points.size());
                if (std::abs(mean_r - r_outer_c) < 0.12) m.outer_w = std::max(m.outer_w, p->width);
                if (std::abs(mean_r - r_bore_c)  < 0.12) m.bore_w  = std::max(m.bore_w,  p->width);
            }
        }
        return m;
    };

    const auto squeezed = measure(4.7);
    REQUIRE(squeezed.rings == 5);                       // the extra ring came in EARLY...
    REQUIRE(squeezed.min_w < 0.5f - 0.04f);             // ...so the interior runs squeezed
    REQUIRE(squeezed.max_w < 0.5f + 0.051f);            // and nothing stretches
    REQUIRE_THAT(double(squeezed.outer_w), Catch::Matchers::WithinAbs(0.5, 0.051));
    REQUIRE_THAT(double(squeezed.bore_w),  Catch::Matchers::WithinAbs(0.5, 0.051));
    REQUIRE_THAT(squeezed.vol, Catch::Matchers::WithinRel(squeezed.model, 0.05));

    const auto near_tile = measure(4.05);
    REQUIRE(near_tile.rings == 4);                      // no gratuitous densification
    REQUIRE(near_tile.min_w > 0.5f - 0.051f);           // everything within a bucket of nominal
    REQUIRE(near_tile.max_w < 0.5f + 0.06f);            // the 1% stretch, bucketed
    REQUIRE_THAT(near_tile.vol, Catch::Matchers::WithinRel(near_tile.model, 0.05));
}

TEST_CASE("Adaptive beading transitions one ring through an angular window on an eccentric wall", "[ContinuousShell]") {
    // An eccentric bore makes t(theta) cross ONE count threshold: the thick side plans
    // five rings, the thin side four, and the fifth (medial) ring is a partial arc whose
    // tips carry the plan's own width ramp down to the termination floor. Coverage must
    // stay gap-free (no missing bead line where the ring retires) and volume conserved.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const double t_base = 4.4 * sp + (0.5 - sp);         // mid transition band
    const double ecc    = 0.34 * sp;                     // +-0.34 spacings of variation
    Polygon outer, inner;
    for (int i = 0; i < 360; ++i) {
        const double a = 2.0 * M_PI * i / 360.0;
        outer.points.emplace_back(coord_t(std::llround(scale_(10.0) * std::cos(a))),
                                  coord_t(std::llround(scale_(10.0) * std::sin(a))));
        inner.points.emplace_back(coord_t(std::llround(scale_(ecc) + scale_(10.0 - t_base) * std::cos(a))),
                                  coord_t(std::llround(scale_(10.0 - t_base) * std::sin(a))));
    }
    ExPolygons ann = diff_ex(Polygons{outer}, Polygons{inner});
    REQUIRE(ann.size() == 1);
    REQUIRE(ann.front().holes.size() == 1);
    const double model_volume = ann.front().area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;

    Surfaces surfaces{ Surface(stInternal, ann.front()) };
    PrintRegionConfig cfg;
    cfg.continuous_shell_adaptive_beading.value = true;
    ExtrusionEntityCollection e;
    const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);

    double volume = 0.;
    float  min_w = 1e9f, max_w = 0.f;
    for (const ExtrusionPath *p : paths) {
        min_w = std::min(min_w, p->width);
        max_w = std::max(max_w, p->width);
        volume += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
        REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
    }
    CAPTURE(min_w, max_w, paths.size());
    REQUIRE(min_w > 0.29f);                              // nothing below the termination floor
    REQUIRE(min_w < 0.36f);                              // ...and the floored thin medial exists
    REQUIRE(max_w < 0.56f);                              // no stretch past ~10%, no C-lite widened medial
    REQUIRE_THAT(volume, Catch::Matchers::WithinRel(model_volume, 0.06));

    // The medial bead never terminates (tips gap - the printed T5 verdict): through the
    // thin sector it thins to the termination floor while its neighbours squeeze to
    // compensate, so EVERY azimuth sees five gap-free rings and the bead count changes
    // nowhere within the layer.
    for (int i = 0; i < 72; ++i) {
        const double az = -M_PI + (i + 0.5) * (2.0 * M_PI / 72.0);
        const std::vector<double> rings = ray_ring_radii(paths, az);
        CAPTURE(az);
        REQUIRE(rings.size() == 5);
        for (size_t j = 1; j < rings.size(); ++j)
            REQUIRE(rings[j] - rings[j - 1] < 1.3 * sp); // no missing-bead gap on any ray
    }
}

TEST_CASE("Adaptive beading retires a ring gradually across a taper via a width ramp", "[ContinuousShell]") {
    // The printed-T4 scenario: an axisymmetric wall thinning with height. The old banding
    // flipped 6->4 rings inside one 0.25mm layer; the plan instead retires ONE ring by
    // thinning it across the flip band (0.25 spacings of thickness, here ~30 layers):
    // the medial ring's width ramps down to the termination floor, holds there, and the
    // topology changes only at the single flip layer - a one-floor-width step instead of
    // a two-full-rings step. Counts stay uniform within every layer and consecutive
    // layers keep short joins (no bead ever terminates mid-air, so no long chords).
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const int    n_layers = 100;
    auto t_sp_of = [&](int layer) { return 4.65 - 0.0075 * layer; };   // spacings; crosses the 4->5 flip at ~4.2

    PrintRegionConfig cfg;
    cfg.continuous_shell_adaptive_beading.value = true;

    std::vector<size_t> counts(n_layers, 0);
    std::vector<float>  layer_min_w(n_layers, 1e9f);
    Point3 prev_last(0, 0, 0);
    bool   have_prev = false;
    for (int layer = 0; layer < n_layers; ++layer) {
        const double t_geom = t_sp_of(layer) * sp + (0.5 - sp);
        const ExPolygon ann = make_annulus(10.0, t_geom);
        Surfaces surfaces{ Surface(stInternal, ann) };
        ExtrusionEntityCollection e;
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e, size_t(layer));
        for (const ExtrusionPath *p : paths)
            layer_min_w[layer] = std::min(layer_min_w[layer], p->width);

        // Ring count is uniform around the circumference on every single layer.
        counts[layer] = ray_ring_radii(paths, 0.3).size();
        for (const double az : {1.7, 2.9, -2.1, -0.8}) {
            CAPTURE(layer, az);
            REQUIRE(ray_ring_radii(paths, az).size() == counts[layer]);
        }

        // Inter-layer join: the next layer's first point must sit near this layer's last.
        const Point3 first = paths.front()->polyline.points.front();
        if (have_prev) {
            const double gap = std::hypot(double(first.x() - prev_last.x()), double(first.y() - prev_last.y())) * SCALING_FACTOR;
            REQUIRE(gap < 3.0);
        }
        prev_last = paths.back()->polyline.points.back();
        have_prev = true;
    }

    // Exactly one count change over the whole taper, and it is single-step 5 -> 4.
    REQUIRE(counts.front() == 5);
    REQUIRE(counts.back() == 4);
    int flip_layer = -1;
    for (int layer = 0; layer + 1 < n_layers; ++layer)
        if (counts[layer + 1] != counts[layer]) {
            REQUIRE(counts[layer] == 5);
            REQUIRE(counts[layer + 1] == 4);
            REQUIRE(flip_layer == -1);                   // only once
            flip_layer = layer + 1;
        }
    REQUIRE(flip_layer > 5);
    REQUIRE(flip_layer < n_layers - 5);

    // While five rings exist, the retiring ring's width ramps down through several
    // buckets to the termination floor and pins there before the flip; it never dips
    // below the floor, and the ramp spans many layers (the flip spreads over height).
    std::set<float> buckets;
    for (int layer = 0; layer < flip_layer; ++layer) {
        REQUIRE(layer_min_w[layer] > 0.29f);
        buckets.insert(std::round(layer_min_w[layer] / 0.05f) * 0.05f);
        if (layer > 0)
            REQUIRE(layer_min_w[layer] <= layer_min_w[layer - 1] + 0.051f);   // non-increasing (bucket slack)
    }
    REQUIRE(buckets.size() >= 3);
    REQUIRE(layer_min_w[flip_layer - 1] < 0.31f);        // floor-pinned just before the flip
    REQUIRE(layer_min_w[flip_layer] > 0.4f);             // settled four-ring layers are near nominal
}

namespace {
// Count of abrupt radial jumps along the whole path - ring-to-ring hops and junction
// diagonals. The interior spiral's entire purpose is to make this collapse.
size_t junction_census(const std::vector<ExtrusionPath *> &paths, double spacing_mm) {
    size_t jumps = 0;
    bool   have_prev = false;
    double prev_r = 0.;
    for (const ExtrusionPath *p : paths)
        for (const Point3 &q : p->polyline.points) {
            const double r = std::hypot(double(q.x()), double(q.y())) * SCALING_FACTOR;
            if (have_prev && std::abs(r - prev_r) > 0.6 * spacing_mm)
                ++jumps;
            prev_r = r;
            have_prev = true;
        }
    return jumps;
}
} // namespace

TEST_CASE("Interior spiral replaces buried rings with a continuous double-arm spiral", "[ContinuousShell]") {
    // Phase B: an 8-bead wall (2 surfaces + 6 buried) prints its interior as one
    // double-arm spiral. The junction census (abrupt radial jumps) collapses from one
    // per ring to almost none; coverage stays gap-free; the merge and wedge widths
    // stay within [w_term, ~1.5 pitch]; volume and the width/flow contract hold.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const ExPolygon ann = make_annulus(12.0, 8.0 * sp + (0.5 - sp));
    const double model_volume = ann.area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;
    Surfaces surfaces{ Surface(stInternal, ann) };

    auto build = [&](bool spiral_on) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        cfg.continuous_shell_interior_spiral.value  = spiral_on;
        ExtrusionEntityCollection e;
        struct R { size_t census; double vol = 0.; float min_w = 1e9f, max_w = 0.f; std::vector<size_t> ray_counts; double max_gap = 0.; size_t n_paths; };
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);
        R r;
        r.n_paths = paths.size();
        r.census = junction_census(paths, sp);
        for (const ExtrusionPath *p : paths) {
            r.min_w = std::min(r.min_w, p->width);
            r.max_w = std::max(r.max_w, p->width);
            r.vol  += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
            REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
        }
        for (int i = 0; i < 36; ++i) {
            const double az = -M_PI + (i + 0.5) * (2.0 * M_PI / 36.0);
            const std::vector<double> radii = ray_ring_radii(paths, az);
            r.ray_counts.push_back(radii.size());
            for (size_t j = 1; j < radii.size(); ++j)
                r.max_gap = std::max(r.max_gap, radii[j] - radii[j - 1]);
        }
        return r;
    };

    const auto ring = build(false);
    const auto spir = build(true);
    CAPTURE(ring.census, spir.census, spir.min_w, spir.max_w, spir.max_gap, spir.n_paths);

    REQUIRE(ring.census >= 6);                          // ring mode: ~one hop per ring
    REQUIRE(spir.census <= 3);                          // spiral: tangent merges, no hops
    REQUIRE(spir.min_w > 0.29f);                        // merge/wedge floors hold
    REQUIRE(spir.max_w < 0.75f);                        // wedge widening bounded ~1.5 pitch
    REQUIRE_THAT(spir.vol, Catch::Matchers::WithinRel(model_volume, 0.07));
    REQUIRE(spir.n_paths < 1200);                       // width buckets keep runs bounded
    for (const size_t c : spir.ray_counts) {
        REQUIRE(c >= 5);                                // gap-free comb everywhere
        REQUIRE(c <= 9);
    }
    REQUIRE(spir.max_gap < 2.4 * sp);                   // worst gap = the U-turn wedge
}

TEST_CASE("Interior spiral holds across a taper and keeps short layer joins", "[ContinuousShell]") {
    // Across a thinning wall the spiral's slot count follows the ring plan's count law
    // (single-step) and the junction census stays collapsed on every layer. The
    // traversal is fixed outer->bore (reversal causes reverse-join doubling), so the
    // inter-layer join is a radial in-material crossing of the wall: bounded by the
    // wall thickness plus the seam's angular drift.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    PrintRegionConfig cfg;
    cfg.continuous_shell_adaptive_beading.value = true;
    cfg.continuous_shell_interior_spiral.value  = true;

    Point3 prev_last(0, 0, 0);
    bool   have_prev = false;
    for (int layer = 0; layer < 40; ++layer) {
        const double t_units = 8.4 - 0.03 * layer;           // crosses one flip band
        const ExPolygon ann = make_annulus(12.0, t_units * sp + (0.5 - sp));
        Surfaces surfaces{ Surface(stInternal, ann) };
        ExtrusionEntityCollection e;
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e, size_t(layer));
        CAPTURE(layer);
        REQUIRE(junction_census(paths, sp) <= 3);
        const Point3 first = paths.front()->polyline.points.front();
        if (have_prev) {
            const double gap = std::hypot(double(first.x() - prev_last.x()), double(first.y() - prev_last.y())) * SCALING_FACTOR;
            REQUIRE(gap < t_units * sp + 1.0);       // one radial wall crossing, plus drift
        }
        prev_last = paths.back()->polyline.points.back();
        have_prev = true;
    }
}

TEST_CASE("Interior spiral composes with an eccentric wall and ignores the interlock", "[ContinuousShell]") {
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();

    SECTION("eccentric wall: gap-free comb, conserved volume") {
        const double t_base = 8.0 * sp + (0.5 - sp);
        Polygon outer, inner;
        for (int i = 0; i < 360; ++i) {
            const double a = 2.0 * M_PI * i / 360.0;
            outer.points.emplace_back(coord_t(std::llround(scale_(12.0) * std::cos(a))),
                                      coord_t(std::llround(scale_(12.0) * std::sin(a))));
            inner.points.emplace_back(coord_t(std::llround(scale_(0.25 * sp) + scale_(12.0 - t_base) * std::cos(a))),
                                      coord_t(std::llround(scale_(12.0 - t_base) * std::sin(a))));
        }
        ExPolygons ann = diff_ex(Polygons{outer}, Polygons{inner});
        REQUIRE(ann.size() == 1);
        const double model_volume = ann.front().area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;
        Surfaces surfaces{ Surface(stInternal, ann.front()) };
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        cfg.continuous_shell_interior_spiral.value  = true;
        ExtrusionEntityCollection e;
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);
        double volume = 0.;
        for (const ExtrusionPath *p : paths) {
            volume += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
            REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
        }
        REQUIRE_THAT(volume, Catch::Matchers::WithinRel(model_volume, 0.08));
        REQUIRE(junction_census(paths, sp) <= 3);
        for (int i = 0; i < 36; ++i) {
            const double az = -M_PI + (i + 0.5) * (2.0 * M_PI / 36.0);
            const std::vector<double> radii = ray_ring_radii(paths, az);
            CAPTURE(az);
            REQUIRE(radii.size() >= 5);
            for (size_t j = 1; j < radii.size(); ++j)
                REQUIRE(radii[j] - radii[j - 1] < 2.4 * sp);
        }
    }
    SECTION("interlock drives the ramped arm's stagger in Z only, and no-ops on a pure spiral") {
        // The interior arm's Z staircase rides the brick-interlock knob when the arm holds
        // ring-like columns: at 0 (the default) the columns keep the nesting the
        // arc-length ramp gives them for free, which measured 0.17 of a layer on the T5
        // tube. A PURE spiral still staggers unconditionally - it is one helical thread,
        // with nothing else resisting radial cleavage - so the knob is inert there.
        const ExPolygon ann = make_annulus(12.0, 8.0 * sp + (0.5 - sp));
        Surfaces surfaces{ Surface(stInternal, ann) };
        auto build = [&](double interlock_pct, int ramp_pct) {
            PrintRegionConfig cfg;
            cfg.continuous_shell_adaptive_beading.value = true;
            cfg.continuous_shell_interior_spiral.value  = true;
            cfg.continuous_shell_junction_ramp.value    = ramp_pct;
            cfg.continuous_shell_interlock.value        = interlock_pct;
            return cs_signature(surfaces, flow, cfg);
        };
        const CsSig flat = build(0, 5), locked = build(100, 5);
        REQUIRE(flat.pts.size() == locked.pts.size());
        bool any_z = false;
        for (size_t i = 0; i < flat.pts.size(); ++i) {
            REQUIRE(std::get<0>(flat.pts[i]) == std::get<0>(locked.pts[i]));   // X untouched
            REQUIRE(std::get<1>(flat.pts[i]) == std::get<1>(locked.pts[i]));   // Y untouched
            REQUIRE(std::get<3>(flat.pts[i]) == std::get<3>(locked.pts[i]));   // width untouched
            if (std::get<2>(flat.pts[i]) != std::get<2>(locked.pts[i]))
                any_z = true;
        }
        REQUIRE(any_z);                                        // it does move Z on a ramped arm

        const CsSig sp_flat = build(0, 100), sp_locked = build(100, 100);
        REQUIRE(sp_flat.pts.size() == sp_locked.pts.size());
        for (size_t i = 0; i < sp_flat.pts.size(); ++i)
            REQUIRE(sp_flat.pts[i] == sp_locked.pts[i]);       // inert on a pure spiral
    }
}

TEST_CASE("Interior spiral declines below four buried rings", "[ContinuousShell]") {
    // A 5-bead wall has too little interior for the double arm; the option must
    // degrade to ring mode byte-identically.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const ExPolygon ann = make_annulus(10.0, 5.0 * sp + (0.5 - sp));
    Surfaces surfaces{ Surface(stInternal, ann) };
    PrintRegionConfig on, off;
    on.continuous_shell_adaptive_beading.value  = true;
    off.continuous_shell_adaptive_beading.value = true;
    on.continuous_shell_interior_spiral.value   = true;
    off.continuous_shell_interior_spiral.value  = false;
    const CsSig a = cs_signature(surfaces, flow, on);
    const CsSig b = cs_signature(surfaces, flow, off);
    REQUIRE(a.pts.size() == b.pts.size());
    for (size_t i = 0; i < a.pts.size(); ++i)
        REQUIRE(a.pts[i] == b.pts[i]);
}

TEST_CASE("Interior spiral carries a per-revolution incommensurate Z-stagger", "[ContinuousShell]") {
    // v5: a monotone spiral offsets radial neighbours by only h/k (smooth laminate
    // interfaces - three printed revisions all cleaved along them), so each buried
    // revolution holds one constant Z level from the golden low-discrepancy sequence
    // in [-h/2, h/2], stepped over a short arc past the seam azimuth. (The v4 sine
    // printed peel-free - mechanism confirmed - but oscillating interdigitation depth
    // left organized void rows; ring mode's interfaces hold one constant depth each.)
    // Assert: (a) the stagger is real and bounded to half a layer, (b) surfaces and
    // the tangent-merge zones stay pinned to the ramp, (c) it is continuous point to
    // point, (d) it is suppressed on layer 0, (e) on rays clear of the step zone the
    // adjacent-interface depths are the golden gaps (~0.38 h / 0.62 h - never
    // coplanar) and (f) each interface's level is CONSTANT across those rays.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const double layer_height = 0.25;
    const ExPolygon ann = make_annulus(14.0, 12.0 * sp + (0.5 - sp));   // deep wall: many buried revolutions
    Surfaces surfaces{ Surface(stInternal, ann) };

    auto build_raw_z = [&](size_t layer_idx) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        cfg.continuous_shell_interior_spiral.value  = true;
        cfg.continuous_shell_junction_ramp.value    = 100;   // the staircase is unconditional here
        ExtrusionEntityCollection e;
        std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, layer_height, &cfg, e, &reason, false, layer_idx));
        // apply_continuous_z NOT called: the points' Z is the raw wave alone.
        std::vector<Point3> pts;
        for (auto *ent : e.entities)
            if (auto *mp = dynamic_cast<ExtrusionMultiPath *>(ent))
                for (auto &p : mp->paths)
                    for (auto &q : p.polyline.points)
                        pts.push_back(q);
        return pts;
    };

    const std::vector<Point3> pts = build_raw_z(4);
    REQUIRE(pts.size() > 500);
    const coord_t amp = coord_t(std::llround(scale_(layer_height / 4.0)));

    // (a) real, and bounded by the WRAP range rather than the amplitude: levels are
    // nest-equivalent modulo a layer, so ramp-tilt compensation may legitimately place a
    // column anywhere within half a layer of nominal. The amplitude governs the gaps
    // BETWEEN columns (checked below), which is what sets how much they overlap.
    const coord_t wrap = coord_t(std::llround(scale_(layer_height / 2.0)));
    coord_t z_reach = 0;
    for (const Point3 &q : pts)
        z_reach = std::max(z_reach, coord_t(std::abs(q.z())));
    REQUIRE(z_reach > amp * 7 / 10);
    REQUIRE(z_reach <= wrap + 2);

    // (b) surface rings and merge-adjacent radii pinned to the ramp
    double r_min = std::numeric_limits<double>::max(), r_max = 0.;
    for (const Point3 &q : pts) {
        const double r = Point(q.x(), q.y()).cast<double>().norm() * SCALING_FACTOR;
        r_min = std::min(r_min, r);
        r_max = std::max(r_max, r);
    }
    for (const Point3 &q : pts) {
        const double r = Point(q.x(), q.y()).cast<double>().norm() * SCALING_FACTOR;
        if (r > r_max - 0.15 * sp || r < r_min + 0.15 * sp)
            REQUIRE(std::abs(q.z()) <= 2);
    }

    // (c) continuity: no jump between consecutive points. A ramped handoff concentrates
    // the level change into a short arc by design, so the bound is a tenth of a layer per
    // segment (a few percent of slope) rather than the pure spiral's gentler drift.
    coord_t max_step = 0;
    for (size_t i = 1; i < pts.size(); ++i)
        max_step = std::max(max_step, coord_t(std::abs(pts[i].z() - pts[i - 1].z())));
    REQUIRE(max_step < coord_t(std::llround(scale_(0.1 * layer_height))));

    // (d) first layer flat
    coord_t z0 = 0;
    for (const Point3 &q : build_raw_z(0))
        z0 = std::max(z0, coord_t(std::abs(q.z())));
    REQUIRE(z0 == 0);

    // (e) THE quantity that matters: the PHYSICAL height offset between radially adjacent
    // columns, which is the stagger plus the arc-length ramp's own rise (neighbours are a
    // revolution apart in path), folded modulo a layer because offsets are nest-equivalent
    // mod h. That offset sets how much two columns overlap - too little and the interfaces
    // are coplanar cleavage planes, too much and the columns barely touch. The design puts
    // it on the golden gaps, so every adjacent pair should land in a band around them and
    // none should be near coplanar. Measured on the FINAL Z (ramp applied), unlike the raw
    // stagger checked above.
    auto build_final_z = [&](size_t layer_idx) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        cfg.continuous_shell_interior_spiral.value  = true;
        cfg.continuous_shell_junction_ramp.value    = 100;   // the staircase is unconditional here
        ExtrusionEntityCollection e;
        std::string reason;
        REQUIRE(ContinuousShellPathBuilder::build_cross_section(surfaces, flow, layer_height, &cfg, e, &reason, false, layer_idx));
        ContinuousShellPathBuilder::apply_continuous_z(&e, layer_height);
        std::vector<Point3> out;
        for (auto *ent : e.entities)
            if (auto *mp = dynamic_cast<ExtrusionMultiPath *>(ent))
                for (auto &p : mp->paths)
                    for (auto &q : p.polyline.points)
                        out.push_back(q);
        return out;
    };
    const std::vector<Point3> fpts = build_final_z(4);
    auto ray_offsets = [&](double az) {
        std::vector<std::pair<double, double>> ray;     // (r_mm, z_mm)
        for (const Point3 &q : fpts) {
            const Vec2d d = Point(q.x(), q.y()).cast<double>();
            double a = std::atan2(d.y(), d.x()) - az;
            a = std::atan2(std::sin(a), std::cos(a));
            if (std::abs(a) < 0.02)
                ray.emplace_back(d.norm() * SCALING_FACTOR, double(q.z()) * SCALING_FACTOR);
        }
        std::sort(ray.begin(), ray.end());
        std::vector<double> offs;
        double prev_z = 0.;
        bool   have   = false;
        for (size_t i = 0; i < ray.size(); ++i) {
            if (i > 0 && ray[i].first - ray[i - 1].first <= 0.5 * sp)
                continue;                               // same crossing
            if (ray[i].first <= r_min + 2.5 * sp || ray[i].first >= r_max - 2.5 * sp) {
                have = false;                           // skip the merge-tapered ends
                continue;
            }
            if (have) {
                const double d = std::abs(ray[i].second - prev_z);
                const double m = d - layer_height * std::floor(d / layer_height);
                offs.push_back(std::min(m, layer_height - m));
            }
            prev_z = ray[i].second;
            have   = true;
        }
        return offs;
    };
    int good_rays = 0, rays = 0;
    for (int k = 0; k < 8; ++k) {
        const std::vector<double> offs = ray_offsets(-M_PI + (k + 0.5) * (2.0 * M_PI / 8.0));
        if (offs.size() < 3)
            continue;
        ++rays;
        bool ok = true;
        for (double o : offs)
            if (o < 0.10 * layer_height || o > 0.42 * layer_height)
                ok = false;
        if (ok)
            ++good_rays;
    }
    REQUIRE(rays >= 6);
    REQUIRE(good_rays >= rays - 2);                   // the handoff band may spoil a ray or two
    // (f) constant per revolution: outside the short step arcs and the merge tapers
    // the stagger holds one exact level, so the vast majority of consecutive points
    // share their z verbatim (a sine would drift at nearly every point).
    size_t same = 0;
    for (size_t i = 1; i < pts.size(); ++i)
        if (std::abs(pts[i].z() - pts[i - 1].z()) <= 2)
            ++same;
    REQUIRE(double(same) >= 0.70 * double(pts.size() - 1));
}

TEST_CASE("Junction ramp holds columns at constant radius between short migrations", "[ContinuousShell]") {
    // The ring <-> spiral continuum. A pure spiral (ramp 100%) migrates continuously, so
    // the interior is one helical thread that can unwind - the printed peel failure. A
    // ramped schedule holds each column at a CONSTANT RADIUS for (1 - ramp) of a turn and
    // concentrates the migration into one short ramp, so each column wraps well past a
    // half turn and is topologically blocked from radial extraction, while the join
    // ridges stay dispersed. Assert on the emitted arm: (a) the migrating fraction of the
    // path tracks the knob, (b) columns hold a constant radius over >180 degrees,
    // (c) coverage stays gap-free (monotone ramp azimuths - non-monotone would double one
    // slot and leave the next empty), and (d) volume is conserved.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const ExPolygon ann = make_annulus(14.0, 10.0 * sp + (0.5 - sp));
    Surfaces surfaces{ Surface(stInternal, ann) };
    const double model_volume = ann.area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;

    // Longest run of consecutive arm points whose radius stays within a quarter spacing,
    // expressed as a fraction of a full turn, plus the fraction of arm length that migrates.
    struct ArmStats { double hold_turn; double migrating; };
    auto arm_stats = [&](int ramp_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        cfg.continuous_shell_interior_spiral.value  = true;
        cfg.continuous_shell_junction_ramp.value    = ramp_pct;
        ExtrusionEntityCollection e;
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);
        std::vector<Point> pts;
        for (const ExtrusionPath *p : paths)
            for (const Point3 &q : p->polyline.points)
                pts.emplace_back(q.x(), q.y());
        REQUIRE(pts.size() > 500);
        Vec2d c(0., 0.);
        for (const Point &p : pts) c += p.cast<double>();
        c /= double(pts.size());
        // Radius and swept angle per step; the two closed surface rings are excluded by
        // dropping steps whose radius sits at the extremes.
        std::vector<double> r(pts.size()), th(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
            const Vec2d d = pts[i].cast<double>() - c;
            r[i]  = d.norm() * SCALING_FACTOR;
            th[i] = std::atan2(d.y(), d.x());
        }
        const double r_lo = *std::min_element(r.begin(), r.end());
        const double r_hi = *std::max_element(r.begin(), r.end());
        double best_turn = 0., run_turn = 0., tot_len = 0.;
        double run_r = r[0];
        std::vector<std::pair<double, double>> steps;   // (radial rate |dr|/dl, dl)
        double tot_dr = 0.;
        for (size_t i = 1; i < pts.size(); ++i) {
            if (r[i] > r_hi - 0.35 * sp || r[i] < r_lo + 0.35 * sp) { run_turn = 0.; run_r = r[i]; continue; }
            double da = th[i] - th[i - 1];
            da = std::atan2(std::sin(da), std::cos(da));
            const double dl = (pts[i] - pts[i - 1]).cast<double>().norm() * SCALING_FACTOR;
            // A column must sweep CONTIGUOUSLY: reject hops and any angular skip, so two
            // separate arcs that merely share a radius (the path is split into sub-paths
            // at every width change) can never be counted as one run.
            if (dl > 0.5 || dl <= 0. || std::abs(da) > 0.15) { run_turn = 0.; run_r = r[i]; continue; }
            const double dr = std::abs(r[i] - r[i - 1]);
            tot_len += dl;
            tot_dr  += dr;
            steps.emplace_back(dr / dl, dl);
            if (std::abs(r[i] - run_r) < 0.25 * sp) {
                run_turn += std::abs(da);
                best_turn = std::max(best_turn, run_turn);
            } else {
                run_turn = 0.;
                run_r    = r[i];
            }
        }
        // Migration concentration, scale-free: the smallest fraction of arm LENGTH that
        // carries 80% of the total radial migration. A pure spiral migrates uniformly, so
        // it needs ~80% of its length; a ramp of f needs only ~f of it.
        std::sort(steps.begin(), steps.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
        double acc_dr = 0., acc_len = 0.;
        for (const auto &[rate, dl] : steps) {
            if (acc_dr >= 0.8 * tot_dr) break;
            acc_dr  += rate * dl;
            acc_len += dl;
        }
        const double mig_len = acc_len;
        double volume = 0.;
        for (const ExtrusionPath *p : paths)
            volume += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
        REQUIRE_THAT(volume, Catch::Matchers::WithinRel(model_volume, 0.10));   // (d)
        return ArmStats{ best_turn / (2.0 * M_PI), tot_len > 0. ? mig_len / tot_len : 1. };
    };

    const ArmStats spiral = arm_stats(100);
    const ArmStats ramped = arm_stats(20);

    // (a) the knob concentrates the migration
    CAPTURE(spiral.migrating, ramped.migrating, spiral.hold_turn, ramped.hold_turn);
    REQUIRE(ramped.migrating < 0.9 * spiral.migrating);
    // (b) THE topological property: a ramped column wraps past a half turn - it cannot be
    // pulled out radially without breaking - where a spiral's "column" is a migrating
    // thread that never holds one radius for anything like that long.
    REQUIRE(ramped.hold_turn > 0.75);
    REQUIRE(ramped.hold_turn > 1.5 * spiral.hold_turn);

    // (c) coverage: every ray still crosses every column, no gap opened by the ramps
    PrintRegionConfig cfg;
    cfg.continuous_shell_adaptive_beading.value = true;
    cfg.continuous_shell_interior_spiral.value  = true;
    cfg.continuous_shell_junction_ramp.value    = 20;
    ExtrusionEntityCollection e;
    const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);
    // Every column is crossed except inside the narrow seam band where the columns hand
    // off; there the transiting bead sits between slots, so a ray may see one fewer
    // crossing. Allow that band and no more, and hold every ray to a gap bound.
    int short_rays = 0;
    for (int i = 0; i < 24; ++i) {
        const double az = -M_PI + (i + 0.5) * (2.0 * M_PI / 24.0);
        const std::vector<double> radii = ray_ring_radii(paths, az);
        CAPTURE(az);
        REQUIRE(radii.size() >= 5);
        if (radii.size() < 6)
            ++short_rays;
        for (size_t j = 1; j < radii.size(); ++j)
            REQUIRE(radii[j] - radii[j - 1] < 2.4 * sp);
    }
    REQUIRE(short_rays <= 3);
}

TEST_CASE("Adaptive beading is safe on generic star tube shapes", "[ContinuousShell]") {
    // The adaptive plan defaults ON, so every star-shaped single-bore section now goes
    // through it (or declines to C-lite). These are INVARIANT tests - they must hold
    // whether the plan engages or declines: the layer builds (no spurious hard-fail),
    // volume is conserved, no ray sees a coverage gap, widths honour the floors, and
    // the width/flow contract holds. Shapes: an elliptical tube, a stadium (slot) tube,
    // and a lumpy three-lobed tube - all star-shaped, none circular.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();

    auto check = [&](const ExPolygon &shape, bool spiral_on, double gap_tol_sp) {
        const double model_volume = shape.area() * SCALING_FACTOR * SCALING_FACTOR * 0.25;
        Surfaces surfaces{ Surface(stInternal, shape) };
        PrintRegionConfig cfg;
        cfg.continuous_shell_interior_spiral.value = spiral_on;   // adaptive on by default
        ExtrusionEntityCollection e;
        const std::vector<ExtrusionPath *> paths = cs_paths(surfaces, flow, cfg, e);
        double volume = 0., total_len = 0.;
        std::map<float, int> width_hist;
        for (const ExtrusionPath *p : paths) {
            REQUIRE(p->width > 0.29f);
            REQUIRE(p->width < 0.80f);
            volume += p->mm3_per_mm * p->polyline.length() * SCALING_FACTOR;
            total_len += p->polyline.length() * SCALING_FACTOR;
            width_hist[std::round(p->width / 0.05f) * 0.05f]++;
            REQUIRE_THAT(double(p->mm3_per_mm), Catch::Matchers::WithinRel(double(flow.with_width(p->width).mm3_per_mm()), 0.001));
        }
        std::string hist;
        for (const auto &kv : width_hist)
            hist += std::to_string(kv.first).substr(0, 4) + "x" + std::to_string(kv.second) + " ";
        CAPTURE(paths.size(), total_len, hist);
        REQUIRE_THAT(volume, Catch::Matchers::WithinRel(model_volume, 0.10));
        // Centroid-based rays: every azimuth must see a gap-free wall.
        double cx = 0., cy = 0.;
        size_t np = 0;
        for (const ExtrusionPath *p : paths)
            for (const Point3 &q : p->polyline.points) { cx += q.x(); cy += q.y(); ++np; }
        cx /= double(np); cy /= double(np);
        for (int i = 0; i < 48; ++i) {
            const double az = -M_PI + (i + 0.5) * (2.0 * M_PI / 48.0);
            std::vector<double> radii;
            for (const ExtrusionPath *p : paths)
                for (const Point3 &q : p->polyline.points) {
                    const double dx = double(q.x()) - cx, dy = double(q.y()) - cy;
                    double d = std::atan2(dy, dx) - az;
                    d = std::atan2(std::sin(d), std::cos(d));
                    // Narrow wedge: on non-circular shapes each ring's radius varies
                    // with theta (dr/dtheta up to a few mm/rad), so a wide wedge smears
                    // discrete rings into one chained cluster.
                    if (std::abs(d) < 0.02)
                        radii.push_back(std::hypot(dx, dy) * SCALING_FACTOR);
                }
            std::sort(radii.begin(), radii.end());
            std::vector<double> clusters;
            for (const double r : radii)
                if (clusters.empty() || r - clusters.back() > 0.15) clusters.push_back(r);
                else clusters.back() = r;
            std::string cl;
            for (const double r : clusters)
                cl += std::to_string(r).substr(0, 5) + " ";
            const double r_lo = radii.empty() ? 0. : radii.front();
            const double r_hi = radii.empty() ? 0. : radii.back();
            CAPTURE(az, clusters.size(), radii.size(), cl, r_lo, r_hi);
            REQUIRE(clusters.size() >= 3);
            for (size_t j = 1; j < clusters.size(); ++j)
                REQUIRE(clusters[j] - clusters[j - 1] < gap_tol_sp * sp);
        }
    };

    SECTION("elliptical tube") {
        Polygon outer, inner;
        for (int i = 0; i < 360; ++i) {
            const double a = 2.0 * M_PI * i / 360.0;
            outer.points.emplace_back(coord_t(std::llround(scale_(14.0) * std::cos(a))),
                                      coord_t(std::llround(scale_(10.0) * std::sin(a))));
            inner.points.emplace_back(coord_t(std::llround(scale_(11.8) * std::cos(a))),
                                      coord_t(std::llround(scale_(7.8) * std::sin(a))));
        }
        ExPolygons ann = diff_ex(Polygons{outer}, Polygons{inner});
        REQUIRE(ann.size() == 1);
        check(ann.front(), false, 2.6);
        check(ann.front(), true, 2.6);                       // and with the spiral requested
    }
    SECTION("stadium (slot) tube") {
        auto stadium = [&](double half_len, double radius) {
            Polygon p;
            for (int i = 0; i <= 90; ++i) {                  // right cap
                const double a = -M_PI / 2 + M_PI * i / 90.0;
                p.points.emplace_back(coord_t(std::llround(scale_(half_len + radius * std::cos(a)))),
                                      coord_t(std::llround(scale_(radius * std::sin(a)))));
            }
            for (int i = 0; i <= 90; ++i) {                  // left cap
                const double a = M_PI / 2 + M_PI * i / 90.0;
                p.points.emplace_back(coord_t(std::llround(scale_(-half_len + radius * std::cos(a)))),
                                      coord_t(std::llround(scale_(radius * std::sin(a)))));
            }
            return p;
        };
        ExPolygons ann = diff_ex(Polygons{stadium(8.0, 8.0)}, Polygons{stadium(8.0, 5.8)});
        REQUIRE(ann.size() == 1);
        check(ann.front(), false, 2.6);
    }
    SECTION("lumpy three-lobed tube") {
        Polygon outer, inner;
        for (int i = 0; i < 480; ++i) {
            const double a = 2.0 * M_PI * i / 480.0;
            const double ro = 11.0 + 1.2 * std::sin(3.0 * a);
            const double ri = ro - 2.3;
            outer.points.emplace_back(coord_t(std::llround(scale_(ro) * std::cos(a))),
                                      coord_t(std::llround(scale_(ro) * std::sin(a))));
            inner.points.emplace_back(coord_t(std::llround(scale_(ri) * std::cos(a))),
                                      coord_t(std::llround(scale_(ri) * std::sin(a))));
        }
        ExPolygons ann = diff_ex(Polygons{outer}, Polygons{inner});
        REQUIRE(ann.size() == 1);
        check(ann.front(), false, 2.8);
    }
}

TEST_CASE("Adaptive beading composes with the brick interlock", "[ContinuousShell]") {
    // The interlock's Z stagger must ride on plan rings exactly as on fixed rings:
    // XY untouched, buried rings offset within +-amp, the join/surface rings pinned flat.
    Flow flow(0.5f, 0.25f, 0.5f);
    const double sp = flow.spacing();
    const ExPolygon ann = make_annulus(10.0, 6.0 * sp + (0.5 - sp));
    Surfaces surfaces{ Surface(stInternal, ann) };

    auto build = [&](double interlock_pct) {
        PrintRegionConfig cfg;
        cfg.continuous_shell_adaptive_beading.value = true;
        cfg.continuous_shell_interlock.value = interlock_pct;
        return cs_signature(surfaces, flow, cfg);
    };
    const CsSig flat = build(0);
    const CsSig locked = build(100);
    REQUIRE(flat.pts.size() == locked.pts.size());
    const coord_t amp = coord_t(std::llround(scale_(0.25 / 4.0)));
    bool any_offset = false;
    for (size_t i = 0; i < flat.pts.size(); ++i) {
        REQUIRE(std::get<0>(flat.pts[i]) == std::get<0>(locked.pts[i]));   // XY identical
        REQUIRE(std::get<1>(flat.pts[i]) == std::get<1>(locked.pts[i]));
        REQUIRE(std::get<2>(flat.pts[i]) == 0);
        const coord_t z = std::get<2>(locked.pts[i]);
        REQUIRE(std::abs(z) <= amp + 2);
        any_offset |= z != 0;
    }
    REQUIRE(any_offset);
}

TEST_CASE("Adaptive beading degrades to C-lite byte-identically out of scope", "[ContinuousShell]") {
    Flow flow(0.5f, 0.25f, 0.5f);

    auto expect_identical = [&](const Surfaces &surfaces, PrintRegionConfig &cfg_on, PrintRegionConfig &cfg_off,
                                size_t layer_index = 3, bool force_solid = false) {
        // The option defaults ON since print validation - pin the states explicitly so
        // this test keeps comparing the adaptive gate against genuine C-lite.
        cfg_on.continuous_shell_adaptive_beading.value  = true;
        cfg_off.continuous_shell_adaptive_beading.value = false;
        const CsSig on  = cs_signature(surfaces, flow, cfg_on,  layer_index, force_solid);
        const CsSig off = cs_signature(surfaces, flow, cfg_off, layer_index, force_solid);
        REQUIRE(on.pts.size() == off.pts.size());
        for (size_t i = 0; i < on.pts.size(); ++i)
            REQUIRE(on.pts[i] == off.pts[i]);
    };

    SECTION("non-star annulus (offset bore) declines to C-lite") {
        // The bore does not enclose the outer centroid, so the angular parametrization is
        // meaningless; the plan declines and the fixed-width rings are byte-identical.
        Polygon outer, inner;
        for (int i = 0; i < 240; ++i) {
            const double a = 2.0 * M_PI * i / 240.0;
            outer.points.emplace_back(coord_t(10000000 * std::cos(a)), coord_t(10000000 * std::sin(a)));
            inner.points.emplace_back(coord_t(5000000 + 3000000 * std::cos(a)), coord_t(3000000 * std::sin(a)));
        }
        ExPolygons ann = diff_ex(Polygons{outer}, Polygons{inner});
        REQUIRE(ann.size() == 1);
        REQUIRE(ann.front().holes.size() == 1);
        Surfaces surfaces{ Surface(stInternal, ann.front()) };
        PrintRegionConfig on, off;
        on.continuous_shell_adaptive_beading.value = true;
        expect_identical(surfaces, on, off);
    }
    SECTION("woven sparse annulus keeps its lattice") {
        const ExPolygon ann = make_annulus(10.0, 4.0);
        Surfaces surfaces{ Surface(stInternal, ann) };
        PrintRegionConfig on, off;
        on.continuous_shell_infill_density.value  = 30;
        off.continuous_shell_infill_density.value = 30;
        on.continuous_shell_adaptive_beading.value = true;
        expect_identical(surfaces, on, off);
    }
    SECTION("solid disc keeps its spiral") {
        Polygon circ;
        for (int i = 0; i < 96; ++i) {
            const double a = 2.0 * M_PI * i / 96.0;
            circ.points.emplace_back(coord_t(8000000 * std::cos(a)), coord_t(8000000 * std::sin(a)));
        }
        Surfaces surfaces{ Surface(stInternal, ExPolygon(circ)) };
        PrintRegionConfig on, off;
        on.continuous_shell_adaptive_beading.value = true;
        expect_identical(surfaces, on, off);
    }
}
