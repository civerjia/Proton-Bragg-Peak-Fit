 #include "parabolic_cylinder_function.h"
 // translate matlab code to c++
 // author shuang zhou civerjia@gmail.com
 // Reference : E. Cojocaru. Parabolic Cylinder Functions 
 // (https://www.mathworks.com/matlabcentral/fileexchange/22620-parabolic-cylinder-functions)

 double PCF::pu(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute parabolic cylinder function U(a, x)
     //% Input : a-- - Parameter(| a | < 5)
     //    % x-- - Argument(| x | < 5)
     //    % Output : u------U(a, x)
     //    % ================================================================== =
     //    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     //    %This program was inspired by the Matlab program 'specfun' (author
     //        % B.Barrowes) which is a direct conversion by 'f2matlab' (author
     //            % B.Barrowes) of the corresponding Fortran program in
     //    % S.Zhang and J.Jin, 'Computation of Special functions' (Wiley, 1966).
     //    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     //    %E.Cojocaru, January 2009
     //    % Observations, suggestionsand recommendations are welcome at e - mail:
     //    % ecojocaru@theory.nipne.ro
     //% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     std::array<double, 101> c;
     std::array<double, 101> d;
     double eps{ 1e-15 };
     double sqpi{ std::sqrt(M_PI) };
     double c0{ 1.0 }, c1{ a };
     c[0] = a;
     double ck{};
     int m{};
     for (int k1 = 4; k1 <= 200; k1 = k1 + 2)
     {
         m = k1 / 2;
         ck = a * c1 + (double(k1) - 2.0) * (double(k1) - 3.0) * c0 / 4.0;
         c[m - 1] = ck;
         c0 = c1;
         c1 = ck;
     }
     double y1{ 1.0 }, r{ 1.0 }, r1{};
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) - 1.0));
         r1 = c[k - 1] * r;
         y1 = y1 + r1;
         if (std::abs(r1 / y1) <= eps && k > 30)
         {
             break;
         }
     }
     double d1{ 1.0 }, d2{ a };
     d[0] = 1.0;
     d[1] = a;

     for (int k2 = 5; k2 <= 160; k2 = k2 + 2)
     {
         m = ((k2 + 1) / 2);//std::trunc
         double dk{ a * d2 + 0.25 * (double(k2) - 2.0) * (double(k2) - 3.0) * d1 };
         d[m - 1] = dk;
         d1 = d2;
         d2 = dk;
     }
     double y2{ 1.0 };
     r = 1.0;
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) + 1.0));
         r1 = d[k] * r;
         y2 = y2 + r1;
         if (std::abs(r1 / y2) <= eps && k > 30)
         {
             break;
         }
     }
     y2 = x * y2;
     double ar{}, f1{}, f2{}, p0{}, g1{}, g3{}, u{};
     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
     {
         ar = M_PI * (0.25 + a / 2.0);
         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
         u = std::cos(ar) * f1 * y1 - std::sin(ar) * f2 * y2;
     }
     else
     {
         p0 = sqpi / pow(2.0, a / 2.0 + 0.25);
         g1 = std::tgamma(0.25 + a / 2.0);
         g3 = std::tgamma(0.75 + a / 2.0);
         u = p0 * (y1 / g3 - M_SQRT2 * y2 / g1);
     }
     return u;
 }

 double PCF::D(double a, double x)
 {
     // %https://mathworld.wolfram.com/ParabolicCylinderFunction.html
     return PCF::pu(-a - 0.5, x);
 }

 double PCF::dpu(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute derivative of parabolic cylinder functions U(a, x)
     //% Input : a-- - Parameter(| a | < 5)
     //% x-- - Argument(| x | < 5)
     //% Output : du------U'(a,x)                
     //% ================================================================== =
     std::array<double, 101> c;
     std::array<double, 101> d;
     double eps{ 1e-15 };
     double sqpi{ std::sqrt(M_PI) };
     double c0{ 1.0 }, c1{ a };
     c[0] = a;
     double ck{};
     int m{};
     for (int k1 = 4; k1 <= 202; k1 = k1 + 2)
     {
         m = k1 / 2;
         ck = a * c1 + (double(k1) - 2.0) * (double(k1) - 3.0) * c0 / 4.0;
         c[m - 1] = ck;
         c0 = c1;
         c1 = ck;
     }
     double y1{ a }, r{ 1.0 }, r1{};
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) + 1.0));
         r1 = c[k] * r;
         y1 = y1 + r1;
         if (std::abs(r1 / y1) <= eps && k > 30)
         {
             break;
         }
     }
     y1 = x * y1;
     double d1{ 1.0 }, d2{ a };
     d[0] = 1.0;
     d[1] = a;

     for (int k2 = 5; k2 <= 160; k2 = k2 + 2)
     {
         m = ((k2 + 1) / 2);//std::trunc
         double dk{ a * d2 + 0.25 * (double(k2) - 2.0) * (double(k2) - 3.0) * d1 };
         d[m - 1] = dk;
         d1 = d2;
         d2 = dk;
     }
     double y2{ 1.0 };
     r = 1.0;
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) - 1.0));
         r1 = d[k] * r;
         y2 = y2 + r1;
         if (std::abs(r1 / y2) <= eps && k > 30)
         {
             break;
         }
     }
     double ar{}, f1{}, f2{}, p0{}, g1{}, g3{}, du{};
     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
     {
         ar = M_PI * (0.25 + a / 2.0);
         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
         du = std::cos(ar) * f1 * y1 - std::sin(ar) * f2 * y2;
     }
     else
     {
         p0 = sqpi / pow(2.0, a / 2.0 + 0.25);
         g1 = std::tgamma(0.25 + a / 2.0);
         g3 = std::tgamma(0.75 + a / 2.0);
         du = p0 * (y1 / g3 - M_SQRT2 * y2 / g1);
     }
     return du;
 }

 double PCF::dD_dx(double a, double x)
 {
     // %https://mathworld.wolfram.com/ParabolicCylinderFunction.html
     return PCF::dpu(-a - 0.5, x);
 }

 double PCF::pulx(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute parabolic cylinder function U(a, x)
     //% for large argument x(| x| >> |a | and | a | moderate)
     //% Input : x-- - Argument
     //% a-- - Parameter
     //% Output : u-- - U(a, x)
     //% Routine called :
     //% 'pvlx' for computing V(a, x) for large | x |
     //% ================================================================== =
     double eps{ 1e-15 };
     double qe{ std::exp(-x * x / 4.0) };
     double a0{ qe * std::pow(std::abs(x),(-a - 0.5)) };
     double r{ 1.0 }, u{ 1.0 };

     for (int k = 1; k <= 20; ++k)
     {
         r = -r * (2.0 * double(k) + a - 0.5) * (2.0 * double(k) + a - 1.5) / (2.0 * double(k) * x * x);
         u = u + r;
         if (std::abs(r / u) <= eps)
         {
             break;
         }
     }
     u = a0 * u;
     if (x < 0.0)
     {
         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
         {
             u = -u * std::sin(M_PI * a);
         }
         else
         {
             double v{ PCF::pvlx(a, -x) };
             double g0{ std::tgamma(a + 0.5) };
             u = M_PI * v / g0 - std::sin(M_PI * a) * u;
         }
     }

     return u;
 }

 double PCF::dpulx(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute derivative of parabolic cylinder function U(a, x)
     //% for large argument x(| x| >> |a | and | a | moderate)
     //% Input : x-- - Argument
     //% a-- - Parameter
     //% Output : du-- - U'(a,x)
     //% Routine called :
     //% 'dpvlx' for computing V'(a,x) for large |x|
     //% ================================================================== =
     double eps{ 1e-15 };
     double x1{ std::abs(x) };
     double qe{ std::exp(-x1 * x1 / 4.0) };
     double q0{ std::pow(x1, (-a - 0.5)) };
     double dqe{ -x1 * qe / 2.0 };
     double dq0{ (-a - 0.5) * std::pow(x1,(-a - 1.5)) };

     double r1{ 1.0 }, su1{ 1.0 };
     double r2{-1.0 }, su2{ 0.0 };

     for (int k = 1; k <= 20; ++k)
     {
         double ak{ (2.0 * double(k) + a - 0.5) * (2.0 * double(k) + a - 1.5) / (2.0 * double(k)) };
         r1 = -r1 * ak / (x1 * x1);
         r2 = -r2 * ak;
         double s2{ 2.0 * double(k) / (std::pow(x1, 2.0 * double(k) + 1.0)) };
         su1 += r1;
         su2 += r2 * s2;
     }
     double du{};
     du = (qe * dq0 + dqe * q0) * su1 + qe * q0 * su2;
     if (x < 0.0)
     {
         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
         {
             du = du * std::sin(M_PI * a);
         }
         else
         {
             double dv{ PCF::dpvlx(a, -x) };
             double g0{ std::tgamma(a + 0.5) };
             du = -(M_PI * dv / g0 - std::sin(M_PI * a) * du);
         }
     }

     return du;
 }

 double PCF::pv(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute parabolic cylinder function V(a, x)
     //% Input : a-- - Parameter(| a | < 5)
     //% x-- - Argument(| x | < 5)
     //% Output : v------V(a, x)
     //% ================================================================== =
     std::array<double, 101> c;
     std::array<double, 101> d;
     double eps{ 1e-15 };
     double sqpi{ std::sqrt(M_PI) };
     double c0{ 1.0 }, c1{ a };
     c[0] = a;
     double ck{};
     int m{};
     for (int k1 = 4; k1 <= 200; k1 = k1 + 2)
     {
         m = k1 / 2;
         ck = a * c1 + (double(k1) - 2.0) * (double(k1) - 3.0) * c0 / 4.0;
         c[m - 1] = ck;
         c0 = c1;
         c1 = ck;
     }
     double y1{ 1.0 }, r{ 1.0 };
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) - 1.0));
         double r1{ c[k - 1] * r };
         y1 = y1 + r1;
         if (std::abs(r1 / y1) <= eps && k > 30)
         {
             break;
         }
     }
     double d1{ 1.0 }, d2{ a };
     d[0] = 1.0;
     d[1] = a;

     for (int k2 = 5; k2 <= 160; k2 = k2 + 2)
     {
         m = ((k2 + 1) / 2);//std::trunc
         double dk{ a * d2 + 0.25 * (double(k2) - 2.0) * (double(k2) - 3.0) * d1 };
         d[m - 1] = dk;
         d1 = d2;
         d2 = dk;
     }
     double y2{ 1.0 };
     r = 1.0;
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) + 1.0));
         double r1{ d[k] * r };
         y2 = y2 + r1;
         if (std::abs(r1 / y2) <= eps && k > 30)
         {
             break;
         }
     }
     y2 = x * y2;
     double ar{}, f1{}, f2{}, p0{}, g0{}, g1{}, g3{}, v{};
     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
     {
         ar = M_PI * (0.25 + a / 2.0);
         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
         v = (std::sin(ar) * f1 * y1 + std::cos(ar) * f2 * y2) / std::tgamma(0.5 - a);
     }
     else
     {
         double sa{ std::sin(M_PI * a) };
         g0 = std::tgamma(0.5 + a);
         p0 = g0 / (sqpi * pow(2.0, a / 2.0 + 0.25));
         g1 = std::tgamma(0.25 + a / 2.0);
         g3 = std::tgamma(0.75 + a / 2.0);
         v = p0 * (y1 * (1 + sa) / g3 + M_SQRT2 * y2 * (1 - sa) / g1);
     }
     return v;
 }

 double PCF::dpv(double a, double x)
 {
     /*% ================================================================== =
     % Purpose: Compute derivative of parabolic cylinder functions V(a, x)
     % Input : a-- - Parameter(| a | < 5)
     % x-- - Argument(| x | < 7.5)
     % Output : dv------V'(a,x)
     % ================================================================== =*/
     std::array<double, 101> c;
     std::array<double, 101> d;
     double eps{ 1e-15 };
     double sqpi{ std::sqrt(M_PI) };
     double c0{ 1.0 }, c1{ a };
     c[0] = a;
     double ck{};
     int m{};
     for (int k1 = 4; k1 <= 202; k1 = k1 + 2)
     {
         m = k1 / 2;
         ck = a * c1 + (double(k1) - 2.0) * (double(k1) - 3.0) * c0 / 4.0;
         c[m - 1] = ck;
         c0 = c1;
         c1 = ck;
     }
     double y1{ a }, r{ 1.0 }, r1{};
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) + 1.0));
         r1 = c[k] * r;
         y1 = y1 + r1;
         if (std::abs(r1 / y1) <= eps && k > 30)
         {
             break;
         }
     }
     y1 = x * y1;
     double d1{ 1.0 }, d2{ a };
     d[0] = 1.0;
     d[1] = a;

     for (int k2 = 5; k2 <= 160; k2 = k2 + 2)
     {
         m = ((k2 + 1) / 2);//std::trunc
         double dk{ a * d2 + 0.25 * (double(k2) - 2.0) * (double(k2) - 3.0) * d1 };
         d[m - 1] = dk;
         d1 = d2;
         d2 = dk;
     }
     double y2{ 1.0 };
     r = 1.0;
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) - 1.0));
         r1 = d[k] * r;
         y2 = y2 + r1;
         if (std::abs(r1 / y2) <= eps && k > 30)
         {
             break;
         }
     }
     double ar{}, f1{}, f2{}, p0{}, g0{}, g1{}, g3{}, dv{};
     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
     {
         ar = M_PI * (0.25 + a / 2.0);
         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
         dv = (std::sin(ar) * f1 * y1 + std::cos(ar) * f2 * y2) / std::tgamma(0.5 - a);
     }
     else
     {
         double sa{ std::sin(M_PI * a) };
         g0 = std::tgamma(0.5 + a);
         p0 = g0 / (sqpi * pow(2.0, a / 2.0 + 0.25));
         g1 = std::tgamma(0.25 + a / 2.0);
         g3 = std::tgamma(0.75 + a / 2.0);
         dv = p0 * (y1 * (1 + sa) / g3 + M_SQRT2 * y2 * (1 - sa) / g1);
     }
     return dv;
 }

 double PCF::pvlx(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute parabolic cylinder function V(a, x)
     //% for large argument x(| x| >> |a | and | a | moderate)
     //% Input : x-- - Argument
     //% a-- - Parameter
     //% Output : v-- - V(a, x)
     //% Routine called :
     //% 'pulx' for computing U(a, x) for large | x |
     //% ================================================================== =
     double eps{ 1e-15 };
     double x1{ std::abs(x) };
     double qe{ std::exp(x1 * x1 / 4.0) };
     double a0{ std::sqrt(M_2_PI) * qe * std::pow(x1,(a - 0.5)) };
     double r{ 1.0 }, v{ 1.0 };

     for (int k = 1; k <= 20; ++k)
     {
         r = r * (2.0 * double(k) - a - 1.5) * (2.0 * double(k) - a - 0.5) / (2.0 * double(k) * x1 * x1);
         v = v + r;
         if (std::abs(r / v) <= eps)
         {
             break;
         }
     }
     v = a0 * v;
     if (x < 0.0)
     {
         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
         {
             v = v * std::sin(M_PI * a);
         }
         else
         {
             double u{ PCF::pulx(a, -x) };
             double g0{ std::tgamma(a + 0.5) };
             double cp{ std::cos(M_PI * a) * std::cos(M_PI * a) };
             v = cp * g0 * u / M_PI + std::sin(M_PI * a) * v;
         }
     }

     return v;
 }

 double PCF::dpvlx(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute derivative of parabolic cylinder function V(a, x)
     //% for large argument x(| x| >> |a | and | a | moderate)
     //% Input : x-- - Argument
     //% a-- - Parameter
     //% Output : dv-- - V'(a,x)
     //% Routine called :
     //% 'dpulx' for computing U'(a,x) for large |x|
     //% ================================================================== =
     double eps{ 1e-15 };
     double s0{ std::sqrt(M_2_PI) };
     double x1{ std::abs(x) };
     double qe{ std::exp(x1 * x1 / 4.0) };
     double q0{ std::pow(x1, (a - 0.5)) };
     double dqe{ x1 * qe / 2.0 };
     double dq0{ (a - 0.5) * std::pow(x1,(a - 1.5)) };

     double r1{ 1.0 }, sv1{ 1.0 };
     double r2{ 1.0 }, sv2{ 0.0 };

     for (int k = 1; k <= 20; ++k)
     {
         double ak{ (2.0 * double(k) - a - 0.5) * (2.0 * double(k) - a - 1.5) / (2.0 * double(k)) };
         r1 = r1 * ak / (x1 * x1);
         r2 = r2 * ak;
         double s2{ -2.0 * k / (std::pow(x1, 2.0 * double(k) + 1.0)) };
         sv1 += r1;
         sv2 += r2 * s2;
     }
     double dv{};
     dv = s0 * ((qe * dq0 + dqe * q0) * sv1 + qe * q0 * sv2);
     if (x < 0.0)
     {
         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
         {
             dv = -dv * std::sin(M_PI * a);
         }
         else
         {
             double du{ PCF::dpulx(a, -x) };
             double g0{ std::tgamma(a + 0.5) };
             double cp{ std::cos(M_PI * a) * std::cos(M_PI * a) };
             dv = -(cp * g0 * du / M_PI + std::sin(M_PI * a) * dv);
         }
     }

     return dv;
 }

 std::complex<double> PCF::cgamma(double x, double y, int kf)
 {
 //% ================================================================== =
 //% Purpose: Compute complex gamma function Gamma(z) or Ln[Gamma(z)]
 //% Input : x-- - Real part of z
 //% y-- - Imaginary part of z
 //% kf-- - Function code
 //% kf = 0 for Ln[Gamma(z)]
 //% kf = 1 for Gamma(z)
 //% Output : gr-- - Real part of Ln[Gamma(z)] or Gamma(z)
 //% gi-- - Imaginary part of Ln[Gamma(z)] or Gamma(z)
 //% ================================================================== =
     double x1 = 0.0;
     // Bernoulli's numbers B2n divided by [2*n(2*n-1)], n = 1,2,... 
     std::array<double, 10> B = { 1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0,
         -1.0 / 1680.0, 1.0 / 1188.0, -691.0 / 360360.0,
         7.0 / 1092.0, -3617.0 / 122400.0, 43867.0 / 244188.0,
         -174611.0 / 125400.0 };
     // Compute firstly Gamma(| x | , sign(x) * y)
     if (x < 0.0) 
     {
         x1 = x;
         x = -x;
         y = -y;
     }
     double x0 = x;
     int na{};// default value 0
     if (x <= 7.0)
     {
         na = (int)(std::trunc(7.0 - x));
         x0 = x + double(na);
     }
     //% Compute log[Gamma(| x | +na, sign(x) * y)] by using relation for | z|>>1:
     //% log[Gamma(z)] = (z - 1 / 2) * log(z) - z + (1 / 2) * log(2 * pi) + ...
     //%     sum{ B2n / [2 * n * (2 * n - 1) * z ^ (2 * n - 1)] }, where summation is from n = 1 to Inf
     double z1 = std::sqrt(x0 * x0 + y * y);
     double th = std::atan(y / x0);
     double gr = (x0 - 0.5) * std::log(z1) - th * y - x0 + 0.5 * std::log(2.0 * M_PI);
     double gi = th * (x0 - 0.5) + y * std::log(z1) - y;
     for (int k = 1; k <= 10; ++k)
     {
         double t = std::pow(z1 , (1.0 - 2.0 * double(k)));
         gr = gr + B[k-1] * t * std::cos((2.0 * double(k) - 1.0) * th);// with B(n) = B2n / [2 * n * (2 * n - 1)]
         gi = gi - B[k - 1] * t * std::sin((2.0 * double(k) - 1.0) * th);
     }
     //% Compute log[Gamma(| x | , sign(x) * y)] from log[Gamma(| x | +na, sign(x) * y)]
     //% by using recurrence relation : Gamma(z + n) = z * (z + 1)*...(z + n - 1)* Gamma(z)
     double gr1{};
     double gi1{};
     if (x <= 7.0)
     {
         gr1 = 0.0;
         gi1 = 0.0;
         for (int j = 0; j < na; ++j)
         {
             gr1 = gr1 + 0.5 * std::log((x + double(j))* (x + double(j)) + y * y);
             gi1 = gi1 + std::atan(y / (x + double(j)));
         }
         gr = gr - gr1;
         gi = gi - gi1;
     }
     //% If x < 0, compute log[Gamma(z)] by using relation:
     //% Gamma(z)* Gamma(-z) = -pi / [z * sin(pi * z)]
     double th1{}, sr{}, si{}, z2{}, th2{};
     if (x1 < 0.0)
     {
         z1 = std::sqrt(x * x + y * y);
         th1 = std::atan(y / x);
         sr = -std::sin(M_PI * x) * std::cosh(M_PI * y);
         si = -std::cos(M_PI * x) * std::sinh(M_PI * y);
         z2 = std::sqrt(sr * sr + si * si);
         th2 = std::atan(si / sr);
         if (sr < 0.0)
         {
             th2 = M_PI + th2;
         }
         gr = std::log(M_PI / (z1 * z2)) - gr;
         gi = -th1 - th2 - gi;
     }
     //% Compute Gamma(z) from log[Gamma(z)] by using relations:
     //% | Gamma(z) |= exp{ real{log[Gamma(z)]} }; arg[Gamma(z)] = imag{ log[Gamma(z)] }
     if (kf == 1)
     {
         double g0 = std::exp(gr);
         gr = g0 * std::cos(gi);
         gi = g0 * std::sin(gi);
     }
     std::complex<double> z (gr,gi);
     return z;
 }

 double PCF::pw(double a, double x)
 {
 //% ================================================================== =
 //% Purpose: Compute parabolic cylinder function W(a, x)
 //% Input : a-- - Parameter(| a | < 5)
 //% x-- - Argument(| x | < 5)
 //% Output : w------W(a, x)
 //% Routine called :
 //% 'cgamma' for computing complex Gamma function
 //% ================================================================== =
 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     std::array<double, 100> c{};
     std::array<double, 100> d{};
     double eps = 1.0e-15;
     double p0 = std::pow(2.0,-0.75);
     double x1 = 0.25;
     double y1 = a / 2.0;
     std::complex<double> z = PCF::cgamma(x1, y1, 1);
     double g1 = std::abs(z);
     double x2 = 0.75;
     z = PCF::cgamma(x2, y1, 1);
     double g3 = std::abs(z);
     double f1 = std::sqrt(g1 / g3);
     double f2 = sqrt(2.0 * g3 / g1);
     double c0 = 1.0;
     double c1 = a;
     c[0] = a;
     for (int k1 = 4; k1 <= 200; k1 = k1 + 2)
     {
         int m = k1 / 2;
         double ck = a * c1 - (double(k1) - 2.0) * (double(k1) - 3.0) * c0 / 4.0;
         c[m - 1] = ck;
         c0 = c1;
         c1 = ck;
     }
     y1 = 1.0;
     double r = 1.0;
     double r1{};
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) - 1.0));
         r1 = c[k-1] * r;
         y1 = y1 + r1;
         if (std::abs(r1 / y1) <= eps && k > 30)
         {
             break;
         }
     }
     double d1 = 1.0;
     double d2 = a;
     d[0] = 1.0;
     d[1] = a;
     double dk{};
     for (int k2 = 5; k2 <= 160; k2 = k2 + 2)
     {
         int m = (int)(std::trunc((k2 + 1) / 2));
         dk = a * d2 - 0.25 * (double(k2) - 2.0) * (double(k2) - 3.0) * d1;
         d[m-1] = dk;
         d1 = d2;
         d2 = dk;
     }
     double y2 = 1.0;
     r = 1.0;
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) + 1.0));
         r1 = d[k] * r;
         y2 = y2 + r1;
         if (std::abs(r1 / y2) <= eps && k > 30)
         {
             break;
         }
     }
     y2 = x * y2;
     double w = p0 * (f1 * y1 - f2 * y2);
     return w;
 }

 double PCF::dpw(double a, double x)
 {
 //% ================================================================== =
 //% Purpose: Compute derivative of parabolic cylinder functions W(a, x)
 //% Input : a-- - Parameter(| a | < 5)
 //% x-- - Argument(| x | < 5)
 //% Output : dw------W'(a,x)                
 //% Routine called :
 //% 'cgamma' for computing complex Gamma function
 //% ================================================================== =
     std::array<double, 101> c{};
     std::array<double, 101> d{};
     double eps = 1.0e-15;
     double p0 = std::pow(2.0, -0.75);
     double x1 = 0.25;
     double y1 = a / 2.0;
     std::complex<double> z = PCF::cgamma(x1, y1, 1);
     double g1 = std::abs(z);
     double x2 = 0.75;
     z = PCF::cgamma(x2, y1, 1);
     double g3 = std::abs(z);
     double f1 = std::sqrt(g1 / g3);
     double f2 = sqrt(2.0 * g3 / g1);
     double c0 = 1.0;
     double c1 = a;
     c[0] = a;
     double ck{};
     for (int k1 = 4; k1 <= 202; k1 = k1 + 2)
     {
         int m = k1 / 2;
         ck = a * c1 - (double(k1) - 2.0) * (double(k1) - 3.0) * c0 / 4.0;
         c[m-1] = ck;
         c0 = c1;
         c1 = ck;
     }
     y1 = a;
     double r = 1.0;
     double r1{};
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) + 1.0));
         r1 = c[k] * r;
         y1 = y1 + r1;
         if (abs(r1 / y1) <= eps && k > 30)
         {
             break;
         }
     }
     y1 = x * y1;
     double d1 = 1.0;
     double d2 = a;
     d[0] = 1.0;
     d[1] = a;
     double dk{};
     for (int k2 = 5; k2 <= 160; k2 = k2 + 2)
     {
         int m = (int)(std::trunc((k2 + 1) / 2));
         dk = a * d2 - 0.25 * (k2 - 2.0) * (k2 - 3.0) * d1;
         d[m-1] = dk;
         d1 = d2;
         d2 = dk;
     }
     double y2 = 1.0;
     r = 1.0;
     for (int k = 1; k <= 100; ++k)
     {
         r = 0.5 * r * x * x / (double(k) * (2.0 * double(k) - 1.0));
         r1 = d[k] * r;
         y2 = y2 + r1;
         if (abs(r1 / y2) <= eps && k > 30)
         {
             break;
         }
     }
     double dw = p0 * (f1 * y1 - f2 * y2);
     return dw;
 }

 double PCF::pwlx(double a, double x)
 {
 //% ================================================================== =
 //% Purpose: Compute parabolic cylinder function W(a, x)
 //% for large argument x(| x| >> |a | and | a | moderate)
 //% Input : a-- - Parameter
 //% x-- - Argument
 //% Output : w------W(a, x)
 //% Routine called :
 //% 'cgamma' for computing complex Gamma function
 //% ================================================================== =
     std::array<double, 100> u{};
     std::array<double, 100> v{};
     std::complex<double> g0 = PCF::cgamma(0.5, a, 0);
     double phi2 = g0.imag();
     g0 = PCF::cgamma(0.5, a, 1);
     double den = std::norm(g0);
     std::complex<double> gk;
     double uk{}, vk{};
     for (int k = 2; k <= 40; k = k + 2)
     {
         int m = (int)(std::trunc(k / 2));
         gk = PCF::cgamma(double(k) + 0.5, a, 1);
         uk = (gk.real() * g0.real() + gk.imag() * g0.imag()) / den;
         vk = (g0.real() * gk.imag() - gk.real() * g0.imag()) / den;
         u[m-1] = uk;
         v[m-1] = vk;
     }
     double sv1 = v[0] / (2.0 * x * x);
     double su2 = u[0] / (2.0 * x * x);
     int fac = 1;
     double fd{};
     for (int k = 3; k <= 20; k = k + 2)
     {
         fac = -(k) * (k - 1) * fac;
         fd = double(fac) * pow(2.0,double(k)) * pow(x , 2.0*double(k));
         sv1 = sv1 + v[k-1] / fd;
         su2 = su2 + u[k-1] / fd;
     }
     double sv2 = 0.0;
     double su1 = 0.0;
     fac = 1;
     for (int k = 2; k <= 20; k = k + 2)
     {
         fac = -(k) * (k - 1) * fac;
         fd = double(fac) * pow(2.0, double(k)) * pow(x, 2.0 * double(k));
         sv2 = sv2 + v[k-1] / fd;
         su1 = su1 + u[k-1] / fd;
     }
     double s1 = 1 + sv1 - su1;
     double s2 = -sv2 - su2;
     double ea = exp(M_PI * a);
     double sea = sqrt(1.0 + ea * ea);
     double fk = sea - ea;
     double ifk = sea + ea;
     double x1 = abs(x);
     double fa = x1 * x1 / 4.0 - a * log(x1) + M_PI / 4.0 + phi2 / 2.0;
     double w{};
     if (x > 0.0)
     {
         w = sqrt(2.0 * fk / x1) * (s1 * cos(fa) - s2 * sin(fa));
     }
     else if (x < 0.0)
     {
         w = sqrt(-2.0 * ifk / x) * (s1 * sin(fa) + s2 * cos(fa));
     }
     return w;
 }

 double PCF::dpwlx(double a, double x)
 {
     //% ================================================================== =
     //% Purpose: Compute parabolic cylinder function W(a, x)
     //% for large argument x(| x| >> |a | and | a | moderate)
     //% Input : a-- - Parameter
     //% x-- - Argument
     //% Output : w------W(a, x)
     //% Routine called :
     //% 'cgamma' for computing complex Gamma function
     //% ================================================================== =
     std::array<double, 100> u{};
     std::array<double, 100> v{};
     std::complex<double> g0 = PCF::cgamma(0.5, a, 0);
     double phi2 = g0.imag();
     g0 = PCF::cgamma(0.5, a, 1);
     double den = std::norm(g0);
     std::complex<double> gk;
     double uk{}, vk{};
     for (int k = 2; k <= 40; k = k + 2)
     {
         int m = (int)(std::trunc(k / 2));
         gk = PCF::cgamma(double(k) + 0.5, a, 1);
         uk = (gk.real() * g0.real() + gk.imag() * g0.imag()) / den;
         vk = (g0.real() * gk.imag() - gk.real() * g0.imag()) / den;
         u[m - 1] = uk;
         v[m - 1] = vk;
     }
     double x1 = abs(x);
     double sv1 = v[0] / (2.0 * x * x);
     double su2 = u[0] / (2.0 * x * x);
     int fac = 1;
     double fd{}, fdd{};
     double dsv1 = -2.0 * sv1 / x1;
     double dsu2 = -2.0 * su2 / x1;
     for (int k = 3; k <= 20; k = k + 2)
     {
         fac = -(k) * (k - 1) * fac;
         fd = double(fac) * pow(2.0, double(k)) * pow(x, 2.0 * double(k));
         fdd = -2.0 * double(k) / x1;
         sv1 = sv1 + v[k - 1] / fd;
         su2 = su2 + u[k - 1] / fd;
         dsv1 = dsv1 + fdd * v[k-1] / fd;
         dsu2 = dsu2 + fdd * u[k-1] / fd;
     }
     double sv2 = 0.0;
     double su1 = 0.0;
     double dsv2 = 0.0;
     double dsu1 = 0.0;
     fac = 1;
     for (int k = 2; k <= 20; k = k + 2)
     {
         fac = -(k) * (k - 1) * fac;
         fd = double(fac) * pow(2.0, double(k)) * pow(x, 2.0 * double(k));
         fdd = -2.0 * double(k) / x1;
         sv2 = sv2 + v[k - 1] / fd;
         su1 = su1 + u[k - 1] / fd;
         dsv2 = dsv2 + fdd * v[k-1] / fd;
         dsu1 = dsu1 + fdd * u[k-1] / fd;
     }
     double s1 = 1 + sv1 - su1;
     double s2 = -sv2 - su2;
     double ds1 = dsv1 - dsu1;
     double ds2 = -dsv2 - dsu2;
     double ea = exp(M_PI * a);
     double sea = sqrt(1.0 + ea * ea);
     double fk = sea - ea;
     double ifk = sea + ea;
     double fa = x1 * x1 / 4.0 - a * log(x1) + M_PI / 4.0 + phi2 / 2.0;
     double dfa = x1 / 2.0 - a / x1;
     double w{}, dw{};
     if (x > 0.0)
     {
         w = sqrt(2.0 * fk / x1) * (s1 * cos(fa) - s2 * sin(fa));
         dw = -w / (2.0 * x1) + sqrt(2.0 * fk / x1) * ((ds1 - s2 * dfa) * cos(fa) - (s1 * dfa + ds2) * sin(fa));
     }
     else if (x < 0.0)
     {
         w = sqrt(-2.0 * ifk / x) * (s1 * sin(fa) + s2 * cos(fa));
         dw = -(-w / (-2.0 * x) + sqrt(-2.0 * ifk / x) * ((ds1 - s2 * dfa) * sin(fa) + (s1 * dfa + ds2) * cos(fa)));
     }
     return dw;
 }

 double PCF::call_pcf_name(double a, double x, int idx)
 {
     double val{};
     switch (idx)
     {
     case 0:
         val = PCF::pu(a, x);
         break;
     case 1:
         val = PCF::dpu(a, x);
         break;
     case 2:
         val = PCF::pulx(a, x);
         break;
     case 3:
         val = PCF::dpulx(a, x);
         break;
     case 4:
         val = PCF::D(a, x);
         break;
     case 5:
         val = PCF::dD_dx(a, x);
         break;
     case 6:
         val = PCF::pv(a, x);
         break;
     case 7:
         val = PCF::dpv(a, x);
         break;
     case 8:
         val = PCF::pvlx(a, x);
         break;
     case 9:
         val = PCF::dpvlx(a, x);
         break;
     case 10:
         val = PCF::pw(a, x);
         break;
     case 11:
         val = PCF::dpw(a, x);
         break;
     case 12:
         val = PCF::pwlx(a, x);
         break;
     case 13:
         val = PCF::dpwlx(a, x);
         break;
     }
     return val;
 }

 double PCF::call_D(double x, double a)
 {
     return PCF::pu(-a - 0.5, x);
 }

 double PCF::call_dD(double x, double a)
 {
     return PCF::dpu(-a - 0.5, x);
 }


 void PCF::call_pcf_name_array(double a, double* x, int idx, int x_size, double* val)
 {
     double (*fcnPtr)(double, double) { &PCF::pu };
     switch (idx)
     {
     case 0:
         fcnPtr = &PCF::pu;
         break;
     case 1:
         fcnPtr = &PCF::dpu;
         break;
     case 2:
         fcnPtr = &PCF::pulx;
         break;
     case 3:
         fcnPtr = &PCF::dpulx;
         break;
     case 4:
         fcnPtr = &PCF::D;
         break;
     case 5:
         fcnPtr = &PCF::dD_dx;
         break;
     case 6:
         fcnPtr = &PCF::pv;
         break;
     case 7:
         fcnPtr = &PCF::dpv;
         break;
     case 8:
         fcnPtr = &PCF::pvlx;
         break;
     case 9:
         fcnPtr = &PCF::dpvlx;
         break;
     case 10:
         fcnPtr = &PCF::pw;
         break;
     case 11:
         fcnPtr = &PCF::dpw;
         break;
     case 12:
         fcnPtr = &PCF::pwlx;
         break;
     case 13:
         fcnPtr = &PCF::dpwlx;
         break;
     }
 #pragma omp parallel for
     for (int i = 0; i < x_size; ++i)
     {
         *(val + i) = fcnPtr(a, *(x + i));
     }
 }

