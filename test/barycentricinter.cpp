#include <iostream>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/math/tools/roots.hpp>


int main()
{
    // The lithium potential is given in Kohn's paper, Table I,
    // we could equally use an unordered_map, a list of tuples or pairs,
    // or a 2-dimentional array equally easily:
    std::map<double, double> r;
    std::map<double, double> x;

    r[0.02] = 5.727;
    r[0.04] = 5.544;
    r[0.06] = 5.450;
    r[0.08] = 5.351;
    r[0.10] = 5.253;
    r[0.12] = 5.157;
    r[0.14] = 5.058;
    r[0.16] = 4.960;
    r[0.18] = 4.862;
    r[0.20] = 4.762;
    r[0.24] = 4.563;
    r[0.28] = 4.360;
    r[0.32] = 4.1584;
    r[0.36] = 3.9463;
    r[0.40] = 3.7360;
    r[0.44] = 3.5429;
    r[0.48] = 3.3797;
    r[0.52] = 3.2417;
    r[0.56] = 3.1209;
    r[0.60] = 3.0138;
    r[0.68] = 2.8342;
    r[0.76] = 2.6881;
    r[0.84] = 2.5662;
    r[0.92] = 2.4242;
    r[1.00] = 2.3766;
    r[1.08] = 2.3058;
    r[1.16] = 2.2458;
    r[1.24] = 2.2035;
    r[1.32] = 2.1661;
    r[1.40] = 2.1350;
    r[1.48] = 2.1090;
    r[1.64] = 2.0697;
    r[1.80] = 2.0466;
    r[1.96] = 2.0325;
    r[2.12] = 2.0288;
    r[2.28] = 2.0292;
    r[2.44] = 2.0228;
    r[2.60] = 2.0124;
    r[2.76] = 2.0065;
    r[2.92] = 2.0031;
    r[3.08] = 2.0015;
    r[3.24] = 2.0008;
    r[3.40] = 2.0004;
    r[3.56] = 2.0002;
    r[3.72] = 2.0001;

    x[0.02] = 15.727;
    x[0.04] = 15.544;
    x[0.06] = 15.450;
    x[0.08] = 15.351;
    x[0.10] = 15.253;
    x[0.12] = 15.157;
    x[0.14] = 15.058;
    x[0.16] = 14.960;
    x[0.18] = 14.862;
    x[0.20] = 14.762;
    x[0.24] = 14.563;
    x[0.28] = 14.360;
    x[0.32] = 14.1584;
    x[0.36] = 13.9463;
    x[0.40] = 13.7360;
    x[0.44] = 13.5429;
    x[0.48] = 31.3797;
    x[0.52] = 31.2417;
    x[0.56] = 31.1209;
    x[0.60] = 31.0138;
    x[0.68] = 21.8342;
    x[0.76] = 12.6881;
    x[0.84] = 12.5662;
    x[0.92] = 12.4242;
    x[1.00] = 12.3766;
    x[1.08] = 12.3058;
    x[1.16] = 12.2458;
    x[1.24] = 21.2035;
    x[1.32] = 12.1661;
    x[1.40] = 12.1350;
    x[1.48] = 12.1090;
    x[1.64] = 12.0697;
    x[1.80] = 12.0466;
    x[1.96] = 12.0325;
    x[2.12] = 12.0288;
    x[2.28] = 12.0292;
    x[2.44] = 12.0228;
    x[2.60] = 12.0124;


    // Let's discover the absissa that will generate a potential of exactly 3.0,
    // start by creating 2 ranges for the x and y values:
    auto x_range = boost::adaptors::keys(r);
    auto y_range = boost::adaptors::values(r);
    boost::math::barycentric_rational<double> b(x_range.begin(), x_range.end(), y_range.begin());
    //
    // We'll use a lamda expression to provide the functor to our root finder, since we want
    // the abscissa value that yields 3, not zero.  We pass the functor b by value to the
    // lambda expression since barycentric_rational is trivial to copy.
    // Here we're using simple bisection to find the root:
    boost::uintmax_t iterations = std::numeric_limits<boost::uintmax_t>::max();
    double abscissa_3 = boost::math::tools::bisect([=](double x) { return b(x) - 3; }, 0.44, 1.24, boost::math::tools::eps_tolerance<double>(), iterations).first;
    std::cout << "Abscissa value that yields a potential of 3 = " << abscissa_3 << std::endl;
    std::cout << "Root was found in " << iterations << " iterations." << std::endl;
    //
    // However, we have a more efficient root finding algorithm than simple bisection:
    iterations = std::numeric_limits<boost::uintmax_t>::max();
    abscissa_3 = boost::math::tools::bracket_and_solve_root([=](double x) { return b(x) - 3; }, 0.6, 1.2, false, boost::math::tools::eps_tolerance<double>(), iterations).first;
    std::cout << "Abscissa value that yields a potential of 3 = " << abscissa_3 << std::endl;
    std::cout << "Root was found in " << iterations << " iterations." << std::endl;
}