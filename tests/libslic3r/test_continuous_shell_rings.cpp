#include <catch2/catch_all.hpp>

#include "libslic3r/ContinuousShellRings.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/AABBTreeLines.hpp"

#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

namespace {

// The fixtures are the objects with print history: a concentric tube, the eccentric whose
// medial pocket took three rounds to chase, and the crescent whose tips void.
Polygon circle_at(double cx, double cy, double r, int n = 192)
{
    Polygon p;
    for (int i = 0; i < n; ++i) {
        const double a = 2.0 * M_PI * double(i) / double(n);
        p.points.emplace_back(scaled(cx + r * std::cos(a)), scaled(cy + r * std::sin(a)));
    }
    return p;
}

ExPolygon annulus(double r_out, double r_bore, double bore_dx = 0.)
{
    ExPolygon e;
    e.contour = circle_at(0., 0., r_out);
    Polygon h = circle_at(bore_dx, 0., r_bore);
    h.reverse();                                  // a hole winds the other way
    e.holes.push_back(h);
    return e;
}

// Distance from every point of `a` to the SEGMENTS of `b` - the same measurement the shipped
// `ring_gap` makes, so the numbers here are comparable with the campaign's.
std::vector<double> gaps_to(const Polygon &a, const Polygons &bs)
{
    Lines all;
    for (const Polygon &b : bs) {
        const Lines l = b.lines();
        all.insert(all.end(), l.begin(), l.end());
    }
    AABBTreeLines::LinesDistancer<Line> d(all);
    std::vector<double> out;
    out.reserve(a.points.size());
    for (const Point &p : a.points)
        out.push_back(d.distance_from_lines<false>(p) * SCALING_FACTOR);
    return out;
}

std::vector<double> gaps_to(const Polygon &a, const Polygon &b) { return gaps_to(a, Polygons{ b }); }

// A family comes back level by level; past a saddle one level is several loops. Grouping by
// `t` is how a consumer finds "the rings one step inward".
std::vector<Polygons> by_level(const CsRingFamily &fam, std::vector<double> *levels = nullptr)
{
    std::vector<Polygons> out;
    for (const CsRing &r : fam.rings) {
        if (out.empty() || std::abs(r.t - (levels ? levels->back() : 0.)) > 1e-9) {
            out.push_back(Polygons{});
            if (levels) levels->push_back(r.t);
        }
        out.back().push_back(r.poly);
    }
    return out;
}

double percentile(std::vector<double> v, double pct)
{
    if (v.empty()) return 0.;
    const size_t k = std::min(v.size() - 1, size_t(std::llround(pct * double(v.size() - 1))));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

bool self_intersects(const Polygon &p)
{
    // A closed loop that crosses itself does not survive a union with itself as one contour.
    return union_ex(Polygons{ p }).size() != 1;
}

} // namespace

TEST_CASE("Ring field level sets of a concentric annulus are concentric circles", "[ContinuousShellRings]")
{
    const ExPolygon reg = annulus(10.0, 4.0);
    CsRingField f;
    std::string fail;
    INFO(fail);
    REQUIRE(f.build(reg, 0.1, &fail));
    REQUIRE(f.n_bores() == 1);

    // The extent of a concentric annulus is its wall thickness everywhere.
    CHECK_THAT(f.extent_percentile(0.5), WithinAbs(6.0, 0.15));

    CsRingFamily fam;
    REQUIRE(f.family(5, fam, &fail));
    REQUIRE(fam.rings.size() == 5);

    for (size_t j = 0; j < fam.rings.size(); ++j) {
        const double t = double(j + 1) / 6.0;
        const double want = 10.0 - t * 6.0;         // r_out - t * thickness
        double lo = 1e9, hi = -1e9;
        for (const Point &p : fam.rings[j].poly.points) {
            const double r = p.cast<double>().norm() * SCALING_FACTOR;
            lo = std::min(lo, r);
            hi = std::max(hi, r);
        }
        INFO("ring " << j << " t=" << t << " radius " << lo << ".." << hi << " want " << want);
        // 1% of a bead spacing (0.4 mm nozzle) is the tolerance the plan names.
        CHECK(std::abs(0.5 * (lo + hi) - want) < 0.05);
        CHECK(hi - lo < 0.05);
        CHECK_FALSE(self_intersects(fam.rings[j].poly));
    }
}

TEST_CASE("Ring field spaces rings evenly where offsetting does not", "[ContinuousShellRings]")
{
    // THE assertion of this change. Repeated offset_ex on these two sections measures a
    // consecutive-gap p99/median of 4.01 (crescent) and 11.22 (eccentric) in millimetres, and
    // the tall values are LOCAL anomalies - a pocket beside gaps of one slot - which is what
    // made a slot stop being a distance.
    //
    // Interpolated rings are even in the FRACTION, so the gap is W(p)/(n+1) with W the local
    // extent. On a crescent, whose wall runs from 9.6 mm at the back to 0.4 mm at the tips, a
    // millimetre gap that varied 24-fold would be correct; what must not vary is the gap
    // measured against the extent it sits in. So that is what is asserted here, together with
    // the direct no-pocket bound: no gap anywhere may exceed what its own extent allows.
    struct Case { const char *name; ExPolygon reg; };
    std::vector<Case> cases;

    cases.push_back({ "eccentric annulus bore@+2", annulus(10.0, 4.0, 2.0) });

    {   // A crescent: a disc with a large bore pushed nearly to the wall, so the wall runs
        // from 9.6 mm at the back to 0.4 mm at the tips.
        ExPolygon c;
        c.contour = circle_at(0., 0., 10.0);
        Polygon h = circle_at(4.6, 0., 5.0);
        h.reverse();
        c.holes.push_back(h);
        cases.push_back({ "crescent", c });
    }

    const size_t n = 4;
    for (const Case &c : cases) {
        CsRingField f;
        std::string fail;
        REQUIRE(f.build(c.reg, 0.08, &fail));

        CsRingFamily fam;
        REQUIRE(f.family(n, fam, &fail));

        // n LEVELS, not n loops: a level is one loop per bore it wraps, and the crescent's tip
        // is thin enough that the raster cannot always keep its two innermost levels apart -
        // which is the coalescence the degradation ladder exists to catch, not a defect here.
        std::vector<double> levels;
        const std::vector<Polygons> lv = by_level(fam, &levels);
        REQUIRE(lv.size() == n);

        std::vector<double> norm;
        double worst_excess = 0.;
        size_t at = 0;
        for (size_t j = 0; j < lv.size(); ++j)
            for (size_t k = 0; k < lv[j].size(); ++k, ++at) {
                if (j + 1 >= lv.size())
                    continue;
                const CsRing &a = fam.rings[at];
                const std::vector<double> g = gaps_to(a.poly, lv[j + 1]);
                for (size_t i = 0; i < g.size(); ++i) {
                    const double allowed = double(a.extent[i]) / double(n + 1);
                    if (allowed < 1e-6)
                        continue;
                    norm.push_back(g[i] / allowed);
                    worst_excess = std::max(worst_excess, g[i] / allowed);
                }
            }
        REQUIRE_FALSE(norm.empty());
        const double p01 = percentile(norm, 0.01);
        const double p50 = percentile(norm, 0.50);
        const double p99 = percentile(norm, 0.99);
        INFO(c.name << ": gap / (local extent / (n+1)) - p01 " << p01 << ", median " << p50
                    << ", p99 " << p99 << ", max " << worst_excess);

        // Even in the fraction, everywhere.
        CHECK(p50 > 0.85);
        CHECK(p50 < 1.15);
        CHECK(p01 > 0.60);
        CHECK(p99 < 1.40);
        // No pocket: a gap may never exceed the extent it sits in. This is the bound repeated
        // offsetting breaks - it put 11.58 slots between two rings inside a 6 mm wall.
        CHECK(worst_excess < 1.60);

        for (const CsRing &r : fam.rings)
            CHECK_FALSE(self_intersects(r.poly));
    }
}

TEST_CASE("Ring field closes a bore-less region onto its centroid", "[ContinuousShellRings]")
{
    ExPolygon disc;
    disc.contour = circle_at(0., 0., 8.0);

    CsRingField f;
    std::string fail;
    REQUIRE(f.build(disc, 0.1, &fail));
    REQUIRE(f.n_bores() == 0);

    CsRingFamily fam;
    REQUIRE(f.family(4, fam, &fail));
    REQUIRE(fam.rings.size() == 4);

    // Rings shrink monotonically toward the centre and every point is labelled bore-less.
    double prev = 1e9;
    for (const CsRing &r : fam.rings) {
        double hi = 0.;
        for (const Point &p : r.poly.points)
            hi = std::max(hi, p.cast<double>().norm() * SCALING_FACTOR);
        INFO("outer radius " << hi);
        CHECK(hi < prev);
        prev = hi;
        CHECK(std::all_of(r.bore.begin(), r.bore.end(), [](int b) { return b == -1; }));
    }
}

TEST_CASE("Ring field subsections a two-bore region and labels each ring point", "[ContinuousShellRings]")
{
    // Two bores in one slab. Near the wall the level set is one loop around both; past the
    // saddle it is one loop per bore. Both are correct topology for the field, and every
    // point carries the bore it belongs to - which IS the Voronoi cell.
    ExPolygon slab;
    slab.contour = Polygon{ { scaled(-16.), scaled(-9.) }, { scaled(16.), scaled(-9.) },
                            { scaled(16.), scaled(9.) },   { scaled(-16.), scaled(9.) } };
    Polygon a = circle_at(-7.5, 0., 3.0), b = circle_at(7.5, 0., 3.0);
    a.reverse(); b.reverse();
    slab.holes.push_back(a);
    slab.holes.push_back(b);

    CsRingField f;
    std::string fail;
    REQUIRE(f.build(slab, 0.1, &fail));
    REQUIRE(f.n_bores() == 2);

    std::vector<CsRing> near_wall, near_bore;
    const bool got_wall = f.level_set(0.15, near_wall, &fail);
    INFO(fail);
    REQUIRE(got_wall);
    REQUIRE(f.level_set(0.85, near_bore, &fail));

    INFO("loops near the wall: " << near_wall.size() << ", near the bores: " << near_bore.size());
    CHECK(near_wall.size() == 1);
    CHECK(near_bore.size() == 2);

    // The near-wall loop is one curve that visits both cells, so it carries both labels; each
    // near-bore loop sits wholly in one cell.
    const std::vector<int> &lab = near_wall.front().bore;
    CHECK(std::count(lab.begin(), lab.end(), 0) > 0);
    CHECK(std::count(lab.begin(), lab.end(), 1) > 0);
    for (const CsRing &r : near_bore) {
        const int first = r.bore.front();
        CHECK(first >= 0);
        CHECK(std::all_of(r.bore.begin(), r.bore.end(), [first](int x) { return x == first; }));
    }
}

TEST_CASE("Ring field rings join across a Voronoi cut", "[ContinuousShellRings]")
{
    // The design claim that only plastic can finally settle, asserted here in geometry: the
    // cut is the equidistant set of the two bore boundaries, so d_in agrees from either side
    // and u is continuous across it. The single near-wall loop must therefore cross the cut
    // without a gap - and it must do so with UNEQUAL bores too, which is exactly what a
    // centroid-bisector partition would get wrong.
    ExPolygon slab;
    slab.contour = Polygon{ { scaled(-16.), scaled(-9.) }, { scaled(16.), scaled(-9.) },
                            { scaled(16.), scaled(9.) },   { scaled(-16.), scaled(9.) } };
    Polygon a = circle_at(-7.0, 0., 4.0), b = circle_at(8.0, 0., 2.0);
    a.reverse(); b.reverse();
    slab.holes.push_back(a);
    slab.holes.push_back(b);

    CsRingField f;
    std::string fail;
    REQUIRE(f.build(slab, 0.08, &fail));

    std::vector<CsRing> level;
    REQUIRE(f.level_set(0.3, level, &fail));
    REQUIRE(level.size() == 1);

    const CsRing &r = level.front();
    // Walk the loop; wherever the label changes we have crossed the cut. The step across must
    // be no longer than the steps either side of it - i.e. the curve is continuous there.
    size_t crossings = 0;
    double worst = 0.;
    const Points &pts = r.poly.points;
    for (size_t i = 0; i < pts.size(); ++i) {
        const size_t j = (i + 1) % pts.size();
        if (r.bore[i] == r.bore[j])
            continue;
        ++crossings;
        worst = std::max(worst, (pts[j] - pts[i]).cast<double>().norm() * SCALING_FACTOR);
    }
    double maxedge = 0.;
    for (size_t i = 0; i < pts.size(); ++i)
        maxedge = std::max(maxedge, (pts[(i + 1) % pts.size()] - pts[i]).cast<double>().norm()
                                    * SCALING_FACTOR);
    INFO("ring has " << pts.size() << " points, longest edge anywhere " << maxedge
                     << " mm; cut crossings " << crossings
                     << ", longest step across a cut " << worst << " mm");
    CHECK(crossings == 2);                        // one cut, entered and left
    CHECK(worst < 0.4);                           // a bead width: no gap at the cut
}

TEST_CASE("Ring field ring count follows the requested pitch", "[ContinuousShellRings]")
{
    const ExPolygon reg = annulus(10.0, 4.0);     // 6 mm of wall
    CsRingField f;
    std::string fail;
    REQUIRE(f.build(reg, 0.1, &fail));

    // n rings leave n+1 gaps across the extent, so a 1 mm pitch across 6 mm wants 5.
    CHECK(f.ring_count_for_pitch(1.0) == 5);
    CHECK(f.ring_count_for_pitch(2.0) == 2);
    CHECK(f.ring_count_for_pitch(6.0) == 0);      // no room for a buried ring is not a failure

    // And the rings it then builds really do sit at that pitch.
    CsRingFamily fam;
    REQUIRE(f.family(f.ring_count_for_pitch(1.0), fam, &fail));
    const std::vector<double> g = gaps_to(fam.rings.front().poly, fam.rings[1].poly);
    CHECK_THAT(percentile(g, 0.5), WithinAbs(1.0, 0.05));
}

TEST_CASE("Ring field refuses what it cannot fill instead of guessing", "[ContinuousShellRings]")
{
    CsRingField f;
    std::string fail;

    ExPolygon empty;
    CHECK_FALSE(f.build(empty, 0.1, &fail));
    CHECK(fail.find("contour") != std::string::npos);

    const ExPolygon reg = annulus(10.0, 4.0);
    fail.clear();
    CHECK_FALSE(f.build(reg, 0.0, &fail));
    CHECK(fail.find("pitch") != std::string::npos);

    REQUIRE(f.build(reg, 0.1, &fail));
    std::vector<CsRing> out;
    CHECK_FALSE(f.level_set(0.0, out, &fail));
    CHECK_FALSE(f.level_set(1.0, out, &fail));
}

// ---------------------------------------------------------------------------------------
// The ring-field signature table. Deliberately HIDDEN ([.]) and asserted against nothing:
// like the density-0 fence it exists to be diffed across a change, not to encode constants.
//
//     libslic3r_tests.exe "[ContinuousShellRingsSig]" > before.txt
//     libslic3r_tests.exe "[ContinuousShellRingsSig]" > after.txt
//
// It also carries the cost budget, because the grid is the one part of this design that can
// quietly get expensive on a large plate.
TEST_CASE("Continuous shell ring field signature table", "[.][ContinuousShellRingsSig]")
{
    struct Row { const char *name; ExPolygon reg; size_t n; };
    std::vector<Row> rows;
    rows.push_back({ "concentric annulus r10/bore4", annulus(10.0, 4.0), 4 });
    rows.push_back({ "eccentric annulus bore@+2",    annulus(10.0, 4.0, 2.0), 4 });
    {
        ExPolygon c;
        c.contour = circle_at(0., 0., 10.0);
        Polygon h = circle_at(4.6, 0., 5.0);
        h.reverse();
        c.holes.push_back(h);
        rows.push_back({ "crescent bore@+4.6 r5", c, 4 });
    }
    {
        ExPolygon d;
        d.contour = circle_at(0., 0., 8.0);
        rows.push_back({ "solid disc r8 (centroid)", d, 4 });
    }
    {
        ExPolygon slab;
        slab.contour = Polygon{ { scaled(-16.), scaled(-9.) }, { scaled(16.), scaled(-9.) },
                                { scaled(16.), scaled(9.) },   { scaled(-16.), scaled(9.) } };
        Polygon a = circle_at(-7.0, 0., 4.0), b = circle_at(8.0, 0., 2.0);
        a.reverse(); b.reverse();
        slab.holes.push_back(a);
        slab.holes.push_back(b);
        rows.push_back({ "two unequal bores in a slab", slab, 4 });
    }
    {   // A plate-sized region: the cost case.
        ExPolygon big;
        big.contour = circle_at(0., 0., 50.0);
        Polygon h = circle_at(0., 0., 20.0);
        h.reverse();
        big.holes.push_back(h);
        rows.push_back({ "annulus r50/bore20 (cost)", big, 6 });
    }

    std::printf("\n%-30s %5s %7s %7s %7s %6s %6s %7s %8s %7s %7s %8s %8s\n", "fixture", "bores",
                "W p50", "W p75", "W max", "rings", "loops", "pts", "max edge", "gap p50",
                "gap p99", "field ms", "rings ms");
    for (const Row &r : rows) {
        CsRingField f;
        std::string fail;
        const auto t0 = std::chrono::steady_clock::now();
        const bool ok = f.build(r.reg, 0.2, &fail);   // half a bead spacing, 0.4 nozzle
        const auto t1 = std::chrono::steady_clock::now();
        if (! ok) {
            std::printf("%-30s  DECLINED: %s\n", r.name, fail.c_str());
            continue;
        }
        CsRingFamily fam;
        const bool built = f.family(r.n, fam, &fail);
        const auto t2 = std::chrono::steady_clock::now();
        if (! built) {
            std::printf("%-30s  FAILED: %s\n", r.name, fail.c_str());
            continue;
        }

        std::vector<double> levels;
        const std::vector<Polygons> lv = by_level(fam, &levels);
        std::vector<double> norm;
        size_t at = 0;
        for (size_t j = 0; j < lv.size(); ++j) {
            for (size_t k = 0; k < lv[j].size(); ++k, ++at) {
                if (j + 1 >= lv.size())
                    continue;
                const CsRing &a = fam.rings[at];
                const std::vector<double> g = gaps_to(a.poly, lv[j + 1]);
                for (size_t i = 0; i < g.size(); ++i) {
                    const double allowed = double(a.extent[i]) / double(r.n + 1);
                    if (allowed > 1e-6)
                        norm.push_back(g[i] / allowed);
                }
            }
        }
        const double ms_field = std::chrono::duration<double, std::milli>(t1 - t0).count();
        const double ms_rings = std::chrono::duration<double, std::milli>(t2 - t1).count();
        // A long edge inside a ring is a dropped-vertex scar: it is what a cut line or a
        // tangency leaves behind when the corrector cannot place a point.
        size_t pts = 0;
        double maxedge = 0.;
        for (const CsRing &cr : fam.rings) {
            pts += cr.poly.points.size();
            for (size_t i = 0; i < cr.poly.points.size(); ++i) {
                const Point &a = cr.poly.points[i];
                const Point &b = cr.poly.points[(i + 1) % cr.poly.points.size()];
                maxedge = std::max(maxedge, (b - a).cast<double>().norm() * SCALING_FACTOR);
            }
        }
        std::printf("%-30s %5zu %7.3f %7.3f %7.3f %6zu %6zu %7zu %8.3f %7.3f %7.3f %8.1f %8.1f\n",
                    r.name, f.n_bores(), f.extent_percentile(0.5), f.extent_percentile(0.75),
                    f.extent_percentile(1.0), r.n, fam.rings.size(), pts, maxedge,
                    percentile(norm, 0.5), percentile(norm, 0.99), ms_field, ms_rings);
    }
    std::printf("\n");
    SUCCEED("signature table printed");
}
