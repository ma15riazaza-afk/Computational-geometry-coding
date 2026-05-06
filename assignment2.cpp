#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/ch_akl_toussaint.h>

#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>
#include <CGAL/enum.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>

// ============================================================
// Common types
// ============================================================

using Clock = std::chrono::steady_clock;

using EpickKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using EpeckKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using SimpleKernel = CGAL::Simple_cartesian<double>;

enum Orientation {
    LEFT_TURN,
    RIGHT_TURN,
    COLLINEAR
};

const char* orientation_to_string(Orientation o)
{
    switch (o) {
        case LEFT_TURN:  return "LEFT_TURN";
        case RIGHT_TURN: return "RIGHT_TURN";
        case COLLINEAR:  return "COLLINEAR";
        default:         return "UNKNOWN";
    }
}

const char* cgal_orientation_to_string(CGAL::Orientation o)
{
    switch (o) {
        case CGAL::LEFT_TURN:  return "LEFT_TURN";
        case CGAL::RIGHT_TURN: return "RIGHT_TURN";
        case CGAL::COLLINEAR:  return "COLLINEAR";
        default:               return "UNKNOWN";
    }
}

// ============================================================
// Exercises 1 and 3: orientation predicates
// ============================================================

Orientation user_orientation(double px, double py,
                             double qx, double qy,
                             double rx, double ry)
{
    double det = (qx - px) * (ry - py) - (qy - py) * (rx - px);

    if (det > 0.0) return LEFT_TURN;
    if (det < 0.0) return RIGHT_TURN;
    return COLLINEAR;
}

Orientation user_orientation_robust(double px, double py,
                                    double qx, double qy,
                                    double rx, double ry)
{
    using BigFloat = boost::multiprecision::cpp_dec_float_100;

    BigFloat bpx(px);
    BigFloat bpy(py);
    BigFloat bqx(qx);
    BigFloat bqy(qy);
    BigFloat brx(rx);
    BigFloat bry(ry);

    BigFloat det = (bqx - bpx) * (bry - bpy)
                 - (bqy - bpy) * (brx - bpx);

    if (det > 0) return LEFT_TURN;
    if (det < 0) return RIGHT_TURN;
    return COLLINEAR;
}

// ============================================================
// Exercises 1--4
// ============================================================

void run_exercise1()
{
    std::cout << "\n========== Exercise 1 ==========\n";

    double px = 0.0, py = 0.0;
    double qx = 1.0, qy = 0.0;
    double rx = 0.0, ry = 1.0;

    Orientation result = user_orientation(px, py, qx, qy, rx, ry);

    std::cout << "Test points: p=(0,0), q=(1,0), r=(0,1)\n";
    std::cout << "user_orientation: " << orientation_to_string(result) << "\n";
}

void run_exercise2()
{
    std::cout << "\n========== Exercise 2 ==========\n";

    const double px = 0.0, py = 0.0;
    const double qx = std::ldexp(1.0, -537);
    const double qy = std::ldexp(1.0, -537);
    const double rx = 0.0;
    const double ry = std::ldexp(1.0, -538);

    Orientation uo = user_orientation(px, py, qx, qy, rx, ry);
    std::cout << "user_orientation: "
              << orientation_to_string(uo) << '\n';

    {
        using K = CGAL::Simple_cartesian<double>;
        K::Point_2 p(px, py), q(qx, qy), r(rx, ry);

        CGAL::Orientation o = CGAL::orientation(p, q, r);
        std::cout << "Simple_cartesian<double>: "
                  << cgal_orientation_to_string(o) << '\n';
    }

    {
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        K::Point_2 p(px, py), q(qx, qy), r(rx, ry);

        CGAL::Orientation o = CGAL::orientation(p, q, r);
        std::cout << "EPICK: "
                  << cgal_orientation_to_string(o) << '\n';
    }

    {
        using K = CGAL::Exact_predicates_exact_constructions_kernel;
        K::Point_2 p(px, py), q(qx, qy), r(rx, ry);

        CGAL::Orientation o = CGAL::orientation(p, q, r);
        std::cout << "EPECK: "
                  << cgal_orientation_to_string(o) << '\n';
    }
}

void run_exercise3()
{
    std::cout << "\n========== Exercise 3 ==========\n";

    const double px = 0.0, py = 0.0;
    const double qx = std::ldexp(1.0, -537);
    const double qy = std::ldexp(1.0, -537);
    const double rx = 0.0;
    const double ry = std::ldexp(1.0, -538);

    Orientation o1 = user_orientation(px, py, qx, qy, rx, ry);
    Orientation o2 = user_orientation_robust(px, py, qx, qy, rx, ry);

    std::cout << "user_orientation: "
              << orientation_to_string(o1) << '\n';

    std::cout << "user_orientation_robust: "
              << orientation_to_string(o2) << '\n';
}

void run_exercise4()
{
    std::cout << "\n========== Exercise 4 ==========\n";

    const double px = 0.0, py = 0.0;
    const double qx = std::ldexp(1.0, -537);
    const double qy = std::ldexp(1.0, -537);
    const double rx = 0.0;
    const double ry = std::ldexp(1.0, -538);

    std::cout << std::left
              << std::setw(45) << "Implementation"
              << std::setw(20) << "Result"
              << "\n";

    std::cout << std::string(65, '-') << "\n";

    std::cout << std::left
              << std::setw(45) << "user_orientation"
              << std::setw(20)
              << orientation_to_string(user_orientation(px, py, qx, qy, rx, ry))
              << "\n";

    std::cout << std::left
              << std::setw(45) << "user_orientation_robust"
              << std::setw(20)
              << orientation_to_string(user_orientation_robust(px, py, qx, qy, rx, ry))
              << "\n";

    {
        using K = CGAL::Simple_cartesian<double>;
        K::Point_2 p(px, py), q(qx, qy), r(rx, ry);
        std::cout << std::left
                  << std::setw(45) << "Simple_cartesian<double>"
                  << std::setw(20)
                  << cgal_orientation_to_string(CGAL::orientation(p, q, r))
                  << "\n";
    }

    {
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        K::Point_2 p(px, py), q(qx, qy), r(rx, ry);
        std::cout << std::left
                  << std::setw(45) << "EPICK"
                  << std::setw(20)
                  << cgal_orientation_to_string(CGAL::orientation(p, q, r))
                  << "\n";
    }

    {
        using K = CGAL::Exact_predicates_exact_constructions_kernel;
        K::Point_2 p(px, py), q(qx, qy), r(rx, ry);
        std::cout << std::left
                  << std::setw(45) << "EPECK"
                  << std::setw(20)
                  << cgal_orientation_to_string(CGAL::orientation(p, q, r))
                  << "\n";
    }
}

void run_exercise5()
{
    std::cout << "\n========== Exercise 5 (bonus) ==========\n";
    std::cout << std::setprecision(17);

    /*
        Intended mathematical points:

        p = (10^20,     10^20)
        q = (10^20 + 1, 10^20)
        r = (10^20,     10^20 + 1)

        Mathematically, these points form a LEFT_TURN.

        However, in double arithmetic, 10^20 + 1 is rounded to 10^20.
        Therefore, the information is already lost before the robust
        orientation function starts.
    */

    const double base = 1e20;

    const double px = base;
    const double py = base;

    const double qx = base + 1.0;
    const double qy = base;

    const double rx = base;
    const double ry = base + 1.0;

    std::cout << "Intended mathematical result: LEFT_TURN\n\n";

    std::cout << "Stored double coordinates:\n";
    std::cout << "p = (" << px << ", " << py << ")\n";
    std::cout << "q = (" << qx << ", " << qy << ")\n";
    std::cout << "r = (" << rx << ", " << ry << ")\n\n";

    Orientation robust_result =
        user_orientation_robust(px, py, qx, qy, rx, ry);

    std::cout << "user_orientation_robust: "
              << orientation_to_string(robust_result) << '\n';

    {
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;

        K::Point_2 p(px, py);
        K::Point_2 q(qx, qy);
        K::Point_2 r(rx, ry);

        CGAL::Orientation epick_result = CGAL::orientation(p, q, r);

        std::cout << "EPICK: "
                  << cgal_orientation_to_string(epick_result) << '\n';
    }

    std::cout << "\nExplanation: both functions return COLLINEAR because "
              << "the values 10^20 + 1 were rounded to 10^20 when stored "
              << "as double. The lost information cannot be recovered later.\n";
}

// ============================================================
// Exercise 6: segment intersection
// ============================================================

void run_exercise6()
{
    std::cout << "\n========== Exercise 6 ==========\n";
    std::cout << std::setprecision(17);

    {
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;

        K::Point_2 A(0, 0), B(3, 1);
        K::Point_2 C(1, 0), D(1, 3);

        K::Segment_2 s1(A, B);
        K::Segment_2 s2(C, D);

        auto result = CGAL::intersection(s1, s2);

        if (result) {
            if (const K::Point_2* p = std::get_if<K::Point_2>(&*result)) {
                std::cout << "EPICK intersection: " << *p << '\n';
                std::cout << "EPICK x = " << p->x() << '\n';
                std::cout << "EPICK y = " << p->y() << '\n';
            }
        }
    }

    std::cout << '\n';

    {
        using K = CGAL::Exact_predicates_exact_constructions_kernel;

        K::Point_2 A(0, 0), B(3, 1);
        K::Point_2 C(1, 0), D(1, 3);

        K::Segment_2 s1(A, B);
        K::Segment_2 s2(C, D);

        auto result = CGAL::intersection(s1, s2);

        if (result) {
            if (const K::Point_2* p = std::get_if<K::Point_2>(&*result)) {
                auto exact_x = CGAL::exact(p->x());
                auto exact_y = CGAL::exact(p->y());

                std::cout << "EPECK intersection: " << *p << '\n';
                std::cout << "EPECK exact x = " << exact_x << '\n';
                std::cout << "EPECK exact y = " << exact_y << '\n';
                std::cout << "EPECK approx point = ("
                          << CGAL::to_double(p->x()) << ", "
                          << CGAL::to_double(p->y()) << ")\n";
            }
        }
    }
}

// ============================================================
// Exercise 7: point generators
// ============================================================

struct RawPoint {
    double x;
    double y;
};

enum class Distribution {
    Circle,
    Disk,
    Square
};

std::string distribution_name(Distribution d)
{
    switch (d) {
        case Distribution::Circle: return "Circle";
        case Distribution::Disk:   return "Disk";
        case Distribution::Square: return "Square";
    }

    return "Unknown";
}

double pi_value()
{
    return std::acos(-1.0);
}

std::vector<RawPoint> generate_raw_points_on_circle(std::size_t n,
                                                    double R,
                                                    unsigned int seed)
{
    std::vector<RawPoint> points;
    points.reserve(n);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * pi_value());

    for (std::size_t i = 0; i < n; ++i) {
        const double theta = angle_dist(gen);

        points.push_back({
            R * std::cos(theta),
            R * std::sin(theta)
        });
    }

    return points;
}

std::vector<RawPoint> generate_raw_points_in_disk(std::size_t n,
                                                  double R,
                                                  unsigned int seed)
{
    std::vector<RawPoint> points;
    points.reserve(n);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t i = 0; i < n; ++i) {
        const double theta = 2.0 * pi_value() * dist(gen);
        const double radius = R * std::sqrt(dist(gen));

        points.push_back({
            radius * std::cos(theta),
            radius * std::sin(theta)
        });
    }

    return points;
}

std::vector<RawPoint> generate_raw_points_in_square(std::size_t n,
                                                    double R,
                                                    unsigned int seed)
{
    std::vector<RawPoint> points;
    points.reserve(n);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-R / 2.0, R / 2.0);

    for (std::size_t i = 0; i < n; ++i) {
        points.push_back({
            dist(gen),
            dist(gen)
        });
    }

    return points;
}

std::vector<RawPoint> generate_raw_points(Distribution distribution,
                                          std::size_t n,
                                          double R,
                                          unsigned int seed)
{
    switch (distribution) {
        case Distribution::Circle:
            return generate_raw_points_on_circle(n, R, seed);

        case Distribution::Disk:
            return generate_raw_points_in_disk(n, R, seed);

        case Distribution::Square:
            return generate_raw_points_in_square(n, R, seed);
    }

    return {};
}

std::vector<RawPoint> generate_raw_triangle_pointset(std::size_t interior_points_count,
                                                     unsigned int seed)
{
    std::vector<RawPoint> points;
    points.reserve(interior_points_count + 3);

    points.push_back({0.0, 0.0});
    points.push_back({1.0, 0.0});
    points.push_back({0.0, 1.0});

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t i = 0; i < interior_points_count; ++i) {
        double u = dist(gen);
        double v = dist(gen);

        if (u + v > 1.0) {
            u = 1.0 - u;
            v = 1.0 - v;
        }

        points.push_back({u, v});
    }

    return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_2>
convert_to_kernel_points(const std::vector<RawPoint>& raw_points)
{
    using Point = typename Kernel::Point_2;

    std::vector<Point> points;
    points.reserve(raw_points.size());

    for (const RawPoint& p : raw_points) {
        points.emplace_back(p.x, p.y);
    }

    return points;
}

// ============================================================
// Exercise 7: benchmark functions
// ============================================================

template <typename Point>
double benchmark_graham_andrew(const std::vector<Point>& points,
                               std::size_t& hull_size)
{
    std::vector<Point> hull;
    hull.reserve(points.size());

    const auto start = Clock::now();

    CGAL::ch_graham_andrew(
        points.begin(),
        points.end(),
        std::back_inserter(hull)
    );

    const auto end = Clock::now();

    hull_size = hull.size();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

template <typename Point>
double benchmark_jarvis(const std::vector<Point>& points,
                        std::size_t& hull_size)
{
    std::vector<Point> hull;
    hull.reserve(points.size());

    const auto start = Clock::now();

    CGAL::ch_jarvis(
        points.begin(),
        points.end(),
        std::back_inserter(hull)
    );

    const auto end = Clock::now();

    hull_size = hull.size();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

template <typename Point>
double benchmark_akl_toussaint(const std::vector<Point>& points,
                               std::size_t& hull_size)
{
    std::vector<Point> hull;
    hull.reserve(points.size());

    const auto start = Clock::now();

    CGAL::ch_akl_toussaint(
        points.begin(),
        points.end(),
        std::back_inserter(hull)
    );

    const auto end = Clock::now();

    hull_size = hull.size();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

template <typename Kernel>
void run_exercise7a_for_kernel(const std::string& kernel_name,
                               Distribution distribution,
                               std::size_t n,
                               const std::vector<RawPoint>& raw_points)
{
    using Point = typename Kernel::Point_2;

    const std::vector<Point> points =
        convert_to_kernel_points<Kernel>(raw_points);

    std::size_t hull_size_graham = 0;
    std::size_t hull_size_jarvis = 0;
    std::size_t hull_size_akl_toussaint = 0;

    const double graham_time_ms =
        benchmark_graham_andrew(points, hull_size_graham);

    const double jarvis_time_ms =
        benchmark_jarvis(points, hull_size_jarvis);

    const double akl_time_ms =
        benchmark_akl_toussaint(points, hull_size_akl_toussaint);

    if (hull_size_graham != hull_size_jarvis ||
        hull_size_graham != hull_size_akl_toussaint) {
        std::cerr << "Warning: different hull sizes for "
                  << kernel_name << ", "
                  << distribution_name(distribution)
                  << ", n = " << n
                  << ": Graham-Andrew = " << hull_size_graham
                  << ", Jarvis = " << hull_size_jarvis
                  << ", Akl-Toussaint = " << hull_size_akl_toussaint
                  << '\n';
    }

    std::cout
        << std::left
        << std::setw(20) << kernel_name
        << std::setw(15) << distribution_name(distribution)
        << std::setw(15) << points.size()
        << std::setw(15) << hull_size_graham
        << std::setw(20) << graham_time_ms
        << std::setw(20) << jarvis_time_ms
        << std::setw(20) << akl_time_ms
        << "\n";
}

template <typename Function>
double measure_average_time_ms(Function function, int repetitions)
{
    double total_ms = 0.0;

    for (int i = 0; i < repetitions; ++i) {
        const auto start = Clock::now();

        function();

        const auto end = Clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;
        total_ms += duration.count();
    }

    return total_ms / repetitions;
}

template <typename Point>
void run_exercise7b_for_kernel(const std::string& kernel_name,
                               const std::vector<Point>& points,
                               int repetitions)
{
    std::vector<Point> hull_graham;
    std::vector<Point> hull_jarvis;
    std::vector<Point> hull_akl;

    hull_graham.reserve(points.size());
    hull_jarvis.reserve(points.size());
    hull_akl.reserve(points.size());

    const double graham_time_ms = measure_average_time_ms([&]() {
        hull_graham.clear();

        CGAL::ch_graham_andrew(
            points.begin(),
            points.end(),
            std::back_inserter(hull_graham)
        );
    }, repetitions);

    const double jarvis_time_ms = measure_average_time_ms([&]() {
        hull_jarvis.clear();

        CGAL::ch_jarvis(
            points.begin(),
            points.end(),
            std::back_inserter(hull_jarvis)
        );
    }, repetitions);

    const double akl_time_ms = measure_average_time_ms([&]() {
        hull_akl.clear();

        CGAL::ch_akl_toussaint(
            points.begin(),
            points.end(),
            std::back_inserter(hull_akl)
        );
    }, repetitions);

    if (hull_graham.size() != hull_jarvis.size() ||
        hull_graham.size() != hull_akl.size()) {
        std::cerr << "Warning: different hull sizes for "
                  << kernel_name
                  << ": Graham-Andrew = " << hull_graham.size()
                  << ", Jarvis = " << hull_jarvis.size()
                  << ", Akl-Toussaint = " << hull_akl.size()
                  << '\n';
    }

    std::cout
        << std::left
        << std::setw(20) << kernel_name
        << std::setw(15) << points.size()
        << std::setw(15) << hull_graham.size()
        << std::setw(20) << graham_time_ms
        << std::setw(20) << jarvis_time_ms
        << std::setw(20) << akl_time_ms
        << "\n";
}

// ============================================================
// Exercise 7 runners
// ============================================================

void run_exercise7a()
{
    std::cout << "\n========== Exercise 7(a), 7(c), 7(d) ==========\n";
    std::cout << std::fixed << std::setprecision(4);

    const double R = 1.0;

    const std::vector<std::size_t> sizes = {
        1000,
        5000,
        10000,
        50000
    };

    const std::vector<Distribution> distributions = {
        Distribution::Circle,
        Distribution::Disk,
        Distribution::Square
    };

    std::cout
        << std::left
        << std::setw(20) << "Kernel"
        << std::setw(15) << "Distribution"
        << std::setw(15) << "Points"
        << std::setw(15) << "Hull size"
        << std::setw(20) << "Graham-Andrew"
        << std::setw(20) << "Jarvis"
        << std::setw(20) << "Akl-Toussaint"
        << "\n";

    std::cout << std::string(125, '-') << "\n";

    for (Distribution distribution : distributions) {
        for (std::size_t n : sizes) {
            const unsigned int seed =
                1000u
                + static_cast<unsigned int>(distribution) * 100u
                + static_cast<unsigned int>(n);

            const std::vector<RawPoint> raw_points =
                generate_raw_points(distribution, n, R, seed);

            run_exercise7a_for_kernel<EpickKernel>(
                "EPICK",
                distribution,
                n,
                raw_points
            );

            run_exercise7a_for_kernel<SimpleKernel>(
                "Simple Cartesian",
                distribution,
                n,
                raw_points
            );
        }

        std::cout << std::string(125, '-') << "\n";
    }
}

void run_exercise7b()
{
    std::cout << "\n========== Exercise 7(b), 7(c), 7(d) ==========\n";
    std::cout << std::fixed << std::setprecision(4);

    const std::size_t interior_points_count = 50000;
    const int repetitions = 100;
    const unsigned int seed = 42;

    const std::vector<RawPoint> raw_points =
        generate_raw_triangle_pointset(interior_points_count, seed);

    const std::vector<EpickKernel::Point_2> epick_points =
        convert_to_kernel_points<EpickKernel>(raw_points);

    const std::vector<SimpleKernel::Point_2> simple_points =
        convert_to_kernel_points<SimpleKernel>(raw_points);

    std::cout
        << std::left
        << std::setw(20) << "Kernel"
        << std::setw(15) << "Points"
        << std::setw(15) << "Hull size"
        << std::setw(20) << "Graham-Andrew"
        << std::setw(20) << "Jarvis"
        << std::setw(20) << "Akl-Toussaint"
        << "\n";

    std::cout << std::string(110, '-') << "\n";

    run_exercise7b_for_kernel(
        "EPICK",
        epick_points,
        repetitions
    );

    run_exercise7b_for_kernel(
        "Simple Cartesian",
        simple_points,
        repetitions
    );
}

// ============================================================
// main
// ============================================================

int main()
{
    run_exercise1();
    run_exercise2();
    run_exercise3();
    run_exercise4();
    run_exercise5();
    run_exercise6();
    run_exercise7a();
    run_exercise7b();

    return 0;
}