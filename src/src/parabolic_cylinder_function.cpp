// #include "parabolic_cylinder_function.h"
// // translate matlab code to c++
// // author shuang zhou civerjia@gmail.com
// // Reference : E. Cojocaru. Parabolic Cylinder Functions 
// // (https://www.mathworks.com/matlabcentral/fileexchange/22620-parabolic-cylinder-functions)
// template <class T>
// T PCF::pu(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute parabolic cylinder function U(a, x)
//     //% Input : a-- - Parameter(| a | < 5)
//     //    % x-- - Argument(| x | < 5)
//     //    % Output : u------U(a, x)
//     //    % ================================================================== =
//     //    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     //    %This program was inspired by the Matlab program 'specfun' (author
//     //        % B.Barrowes) which is a direct conversion by 'f2matlab' (author
//     //            % B.Barrowes) of the corresponding Fortran program in
//     //    % S.Zhang and J.Jin, 'Computation of Special functions' (Wiley, 1966).
//     //    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     //    %E.Cojocaru, January 2009
//     //    % Observations, suggestionsand recommendations are welcome at e - mail:
//     //    % ecojocaru@theory.nipne.ro
//     //% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     std::array<T, 101> c;
//     std::array<T, 101> d;
//     T eps{ 1e-15 };
//     T sqpi{ std::sqrt(M_PI) };
//     T c0{ 1.0 }, c1{ a };
//     c[0] = a;
//     T ck{};
//     long long int m{};
//     for (long long int k1 = 4; k1 <= 200; k1 = k1 + 2)
//     {
//         m = k1 / 2;
//         ck = a * c1 + (T(k1) - 2.0) * (T(k1) - 3.0) * c0 / 4;
//         c[m - 1] = ck;
//         c0 = c1;
//         c1 = ck;
//     }
//     T y1{ 1.0 }, r{ 1.0 }, r1{};
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) - 1.0));
//         r1 = c[k - 1] * r;
//         y1 = y1 + r1;
//         if (std::abs(r1 / y1) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     T d1{ 1.0 }, d2{ a };
//     d[0] = 1.0;
//     d[1] = a;

//     for (long long int k2 = 5; k2 <= 160; k2 = k2 + 2)
//     {
//         m = ((k2 + 1) / 2);//std::trunc
//         T dk{ a * d2 + 0.25 * (T(k2) - 2.0) * (T(k2) - 3.0) * d1 };
//         d[m - 1] = dk;
//         d1 = d2;
//         d2 = dk;
//     }
//     T y2{ 1.0 };
//     r = 1.0;
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) + 1.0));
//         r1 = d[k] * r;
//         y2 = y2 + r1;
//         if (std::abs(r1 / y2) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     y2 = x * y2;
//     T ar{}, f1{}, f2{}, p0{}, g1{}, g3{}, u{};
//     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//     {
//         ar = M_PI * (0.25 + a / 2.0);
//         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
//         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
//         u = std::cos(ar) * f1 * y1 - std::sin(ar) * f2 * y2;
//     }
//     else
//     {
//         p0 = sqpi / pow(2.0, a / 2.0 + 0.25);
//         g1 = std::tgamma(0.25 + a / 2.0);
//         g3 = std::tgamma(0.75 + a / 2.0);
//         u = p0 * (y1 / g3 - M_SQRT2 * y2 / g1);
//     }
//     return u;
// }
// template <class T>
// T PCF::D(T a, T x)
// {
//     // %https://mathworld.wolfram.com/ParabolicCylinderFunction.html
//     return PCF::pu(-a - 0.5, x);
// }
// template <class T>
// T PCF::dpu(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute derivative of parabolic cylinder functions U(a, x)
//     //% Input : a-- - Parameter(| a | < 5)
//     //% x-- - Argument(| x | < 5)
//     //% Output : du------U'(a,x)                
//     //% ================================================================== =
//     std::array<T, 101> c;
//     std::array<T, 101> d;
//     T eps{ 1e-15 };
//     T sqpi{ std::sqrt(M_PI) };
//     T c0{ 1.0 }, c1{ a };
//     c[0] = a;
//     T ck{};
//     long long int m{};
//     for (long long int k1 = 4; k1 <= 202; k1 = k1 + 2)
//     {
//         m = k1 / 2;
//         ck = a * c1 + (T(k1) - 2.0) * (T(k1) - 3.0) * c0 / 4;
//         c[m - 1] = ck;
//         c0 = c1;
//         c1 = ck;
//     }
//     T y1{ a }, r{ 1.0 }, r1{};
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) + 1.0));
//         r1 = c[k] * r;
//         y1 = y1 + r1;
//         if (std::abs(r1 / y1) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     y1 = x * y1;
//     T d1{ 1.0 }, d2{ a };
//     d[0] = 1.0;
//     d[1] = a;

//     for (long long int k2 = 5; k2 <= 160; k2 = k2 + 2)
//     {
//         m = ((k2 + 1) / 2);//std::trunc
//         T dk{ a * d2 + 0.25 * (T(k2) - 2.0) * (T(k2) - 3.0) * d1 };
//         d[m - 1] = dk;
//         d1 = d2;
//         d2 = dk;
//     }
//     T y2{ 1.0 };
//     r = 1.0;
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) - 1.0));
//         r1 = d[k] * r;
//         y2 = y2 + r1;
//         if (std::abs(r1 / y2) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     T ar{}, f1{}, f2{}, p0{}, g1{}, g3{}, du{};
//     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//     {
//         ar = M_PI * (0.25 + a / 2.0);
//         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
//         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
//         du = std::cos(ar) * f1 * y1 - std::sin(ar) * f2 * y2;
//     }
//     else
//     {
//         p0 = sqpi / pow(2.0, a / 2.0 + 0.25);
//         g1 = std::tgamma(0.25 + a / 2.0);
//         g3 = std::tgamma(0.75 + a / 2.0);
//         du = p0 * (y1 / g3 - M_SQRT2 * y2 / g1);
//     }
//     return du;
// }
// template <class T>
// T PCF::dD_dx(T a, T x)
// {
//     // %https://mathworld.wolfram.com/ParabolicCylinderFunction.html
//     return PCF::dpu(-a - 0.5, x);
// }
// template <class T>
// T PCF::pulx(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute parabolic cylinder function U(a, x)
//     //% for large argument x(| x| >> |a | and | a | moderate)
//     //% Input : x-- - Argument
//     //% a-- - Parameter
//     //% Output : u-- - U(a, x)
//     //% Routine called :
//     //% 'pvlx' for computing V(a, x) for large | x |
//     //% ================================================================== =
//     T eps{ 1e-15 };
//     T qe{ std::exp(-x * x / 4.0) };
//     T a0{ qe * std::pow(std::abs(x),(-a - 0.5)) };
//     T r{ 1.0 }, u{ 1.0 };

//     for (long long int k = 1; k <= 20; ++k)
//     {
//         r = -r * (2.0 * T(k) + a - 0.5) * (2.0 * T(k) + a - 1.5) / (2.0 * T(k) * x * x);
//         u = u + r;
//         if (std::abs(r / u) <= eps)
//         {
//             break;
//         }
//     }
//     u = a0 * u;
//     if (x < 0.0)
//     {
//         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//         {
//             u = -u * std::sin(M_PI * a);
//         }
//         else
//         {
//             T v{ PCF::pvlx(a, -x) };
//             T g0{ std::tgamma(a + 0.5) };
//             u = M_PI * v / g0 - std::sin(M_PI * a) * u;
//         }
//     }

//     return u;
// }
// template <class T>
// T PCF::dpulx(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute derivative of parabolic cylinder function U(a, x)
//     //% for large argument x(| x| >> |a | and | a | moderate)
//     //% Input : x-- - Argument
//     //% a-- - Parameter
//     //% Output : du-- - U'(a,x)
//     //% Routine called :
//     //% 'dpvlx' for computing V'(a,x) for large |x|
//     //% ================================================================== =
//     T eps{ 1e-15 };
//     T x1{ std::abs(x) };
//     T qe{ std::exp(-x1 * x1 / 4.0) };
//     T q0{ std::pow(x1, (-a - 0.5)) };
//     T dqe{ -x1 * qe / 2.0 };
//     T dq0{ (-a - 0.5) * std::pow(x1,(-a - 1.5)) };

//     T r1{ 1.0 }, su1{ 1.0 };
//     T r2{-1.0 }, su2{ 0.0 };

//     for (long long int k = 1; k <= 20; ++k)
//     {
//         T ak{ (2.0 * T(k) + a - 0.5) * (2.0 * T(k) + a - 1.5) / (2.0 * T(k)) };
//         r1 = -r1 * ak / (x1 * x1);
//         r2 = -r2 * ak;
//         T s2{ 2.0 * T(k) / (std::pow(x1, 2.0 * T(k) + 1.0)) };
//         su1 += r1;
//         su2 += r2 * s2;
//     }
//     T du{};
//     du = (qe * dq0 + dqe * q0) * su1 + qe * q0 * su2;
//     if (x < 0.0)
//     {
//         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//         {
//             du = du * std::sin(M_PI * a);
//         }
//         else
//         {
//             T dv{ PCF::dpvlx(a, -x) };
//             T g0{ std::tgamma(a + 0.5) };
//             du = -(M_PI * dv / g0 - std::sin(M_PI * a) * du);
//         }
//     }

//     return du;
// }
// template <class T>
// T PCF::pv(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute parabolic cylinder function V(a, x)
//     //% Input : a-- - Parameter(| a | < 5)
//     //% x-- - Argument(| x | < 5)
//     //% Output : v------V(a, x)
//     //% ================================================================== =
//     std::array<T, 101> c;
//     std::array<T, 101> d;
//     T eps{ 1e-15 };
//     T sqpi{ std::sqrt(M_PI) };
//     T c0{ 1.0 }, c1{ a };
//     c[0] = a;
//     T ck{};
//     long long int m{};
//     for (long long int k1 = 4; k1 <= 200; k1 = k1 + 2)
//     {
//         m = k1 / 2;
//         ck = a * c1 + (T(k1) - 2.0) * (T(k1) - 3.0) * c0 / 4.0;
//         c[m - 1] = ck;
//         c0 = c1;
//         c1 = ck;
//     }
//     T y1{ 1.0 }, r{ 1.0 };
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) - 1.0));
//         T r1{ c[k - 1] * r };
//         y1 = y1 + r1;
//         if (std::abs(r1 / y1) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     T d1{ 1.0 }, d2{ a };
//     d[0] = 1.0;
//     d[1] = a;

//     for (long long int k2 = 5; k2 <= 160; k2 = k2 + 2)
//     {
//         m = ((k2 + 1) / 2);//std::trunc
//         T dk{ a * d2 + 0.25 * (T(k2) - 2.0) * (T(k2) - 3.0) * d1 };
//         d[m - 1] = dk;
//         d1 = d2;
//         d2 = dk;
//     }
//     T y2{ 1.0 };
//     r = 1.0;
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) + 1.0));
//         T r1{ d[k] * r };
//         y2 = y2 + r1;
//         if (std::abs(r1 / y2) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     y2 = x * y2;
//     T ar{}, f1{}, f2{}, p0{}, g0{}, g1{}, g3{}, v{};
//     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//     {
//         ar = M_PI * (0.25 + a / 2.0);
//         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
//         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
//         v = (std::sin(ar) * f1 * y1 + std::cos(ar) * f2 * y2) / std::tgamma(0.5 - a);
//     }
//     else
//     {
//         T sa{ std::sin(M_PI * a) };
//         g0 = std::tgamma(0.5 + a);
//         p0 = g0 / (sqpi * pow(2.0, a / 2.0 + 0.25));
//         g1 = std::tgamma(0.25 + a / 2.0);
//         g3 = std::tgamma(0.75 + a / 2.0);
//         v = p0 * (y1 * (1 + sa) / g3 + M_SQRT2 * y2 * (1 - sa) / g1);
//     }
//     return v;
// }
// template <class T>
// T PCF::dpv(T a, T x)
// {
//     /*% ================================================================== =
//     % Purpose: Compute derivative of parabolic cylinder functions V(a, x)
//     % Input : a-- - Parameter(| a | < 5)
//     % x-- - Argument(| x | < 7.5)
//     % Output : dv------V'(a,x)
//     % ================================================================== =*/
//     std::array<T, 101> c;
//     std::array<T, 101> d;
//     T eps{ 1e-15 };
//     T sqpi{ std::sqrt(M_PI) };
//     T c0{ 1.0 }, c1{ a };
//     c[0] = a;
//     T ck{};
//     long long int m{};
//     for (long long int k1 = 4; k1 <= 202; k1 = k1 + 2)
//     {
//         m = k1 / 2;
//         ck = a * c1 + (T(k1) - 2.0) * (T(k1) - 3.0) * c0 / 4.0;
//         c[m - 1] = ck;
//         c0 = c1;
//         c1 = ck;
//     }
//     T y1{ a }, r{ 1.0 }, r1{};
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) + 1.0));
//         r1 = c[k] * r;
//         y1 = y1 + r1;
//         if (std::abs(r1 / y1) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     y1 = x * y1;
//     T d1{ 1.0 }, d2{ a };
//     d[0] = 1.0;
//     d[1] = a;

//     for (long long int k2 = 5; k2 <= 160; k2 = k2 + 2)
//     {
//         m = ((k2 + 1) / 2);//std::trunc
//         T dk{ a * d2 + 0.25 * (T(k2) - 2.0) * (T(k2) - 3.0) * d1 };
//         d[m - 1] = dk;
//         d1 = d2;
//         d2 = dk;
//     }
//     T y2{ 1.0 };
//     r = 1.0;
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) - 1.0));
//         r1 = d[k] * r;
//         y2 = y2 + r1;
//         if (std::abs(r1 / y2) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     T ar{}, f1{}, f2{}, p0{}, g0{}, g1{}, g3{}, dv{};
//     if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//     {
//         ar = M_PI * (0.25 + a / 2.0);
//         f1 = std::tgamma(0.25 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 + 0.25));
//         f2 = std::tgamma(0.75 - a / 2.0) / (sqpi * pow(2.0, a / 2.0 - 0.25));
//         dv = (std::sin(ar) * f1 * y1 + std::cos(ar) * f2 * y2) / std::tgamma(0.5 - a);
//     }
//     else
//     {
//         T sa{ std::sin(M_PI * a) };
//         g0 = std::tgamma(0.5 + a);
//         p0 = g0 / (sqpi * pow(2.0, a / 2.0 + 0.25));
//         g1 = std::tgamma(0.25 + a / 2.0);
//         g3 = std::tgamma(0.75 + a / 2.0);
//         dv = p0 * (y1 * (1 + sa) / g3 + M_SQRT2 * y2 * (1 - sa) / g1);
//     }
//     return dv;
// }
// template <class T>
// T PCF::pvlx(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute parabolic cylinder function V(a, x)
//     //% for large argument x(| x| >> |a | and | a | moderate)
//     //% Input : x-- - Argument
//     //% a-- - Parameter
//     //% Output : v-- - V(a, x)
//     //% Routine called :
//     //% 'pulx' for computing U(a, x) for large | x |
//     //% ================================================================== =
//     T eps{ 1e-15 };
//     T x1{ std::abs(x) };
//     T qe{ std::exp(x1 * x1 / 4.0) };
//     T a0{ std::sqrt(M_2_PI) * qe * std::pow(x1,(a - 0.5)) };
//     T r{ 1.0 }, v{ 1.0 };

//     for (long long int k = 1; k <= 20; ++k)
//     {
//         r = r * (2.0 * T(k) - a - 1.5) * (2.0 * T(k) - a - 0.5) / (2.0 * T(k) * x1 * x1);
//         v = v + r;
//         if (std::abs(r / v) <= eps)
//         {
//             break;
//         }
//     }
//     v = a0 * v;
//     if (x < 0.0)
//     {
//         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//         {
//             v = v * std::sin(M_PI * a);
//         }
//         else
//         {
//             T u{ PCF::pulx(a, -x) };
//             T g0{ std::tgamma(a + 0.5) };
//             T cp{ std::cos(M_PI * a) * std::cos(M_PI * a) };
//             v = cp * g0 * u / M_PI + std::sin(M_PI * a) * v;
//         }
//     }

//     return v;
// }
// template <class T>
// T PCF::dpvlx(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute derivative of parabolic cylinder function V(a, x)
//     //% for large argument x(| x| >> |a | and | a | moderate)
//     //% Input : x-- - Argument
//     //% a-- - Parameter
//     //% Output : dv-- - V'(a,x)
//     //% Routine called :
//     //% 'dpulx' for computing U'(a,x) for large |x|
//     //% ================================================================== =
//     T eps{ 1e-15 };
//     T s0{ std::sqrt(M_2_PI) };
//     T x1{ std::abs(x) };
//     T qe{ std::exp(x1 * x1 / 4.0) };
//     T q0{ std::pow(x1, (a - 0.5)) };
//     T dqe{ x1 * qe / 2.0 };
//     T dq0{ (a - 0.5) * std::pow(x1,(a - 1.5)) };

//     T r1{ 1.0 }, sv1{ 1.0 };
//     T r2{ 1.0 }, sv2{ 0.0 };

//     for (long long int k = 1; k <= 20; ++k)
//     {
//         T ak{ (2.0 * T(k) - a - 0.5) * (2.0 * T(k) - a - 1.5) / (2.0 * T(k)) };
//         r1 = r1 * ak / (x1 * x1);
//         r2 = r2 * ak;
//         T s2{ -2.0 * k / (std::pow(x1, 2.0 * T(k) + 1.0)) };
//         sv1 += r1;
//         sv2 += r2 * s2;
//     }
//     T dv{};
//     dv = s0 * ((qe * dq0 + dqe * q0) * sv1 + qe * q0 * sv2);
//     if (x < 0.0)
//     {
//         if (a < 0 && std::abs(std::trunc(a + 0.5) - (a + 0.5)) <= eps)
//         {
//             dv = -dv * std::sin(M_PI * a);
//         }
//         else
//         {
//             T du{ PCF::dpulx(a, -x) };
//             T g0{ std::tgamma(a + 0.5) };
//             T cp{ std::cos(M_PI * a) * std::cos(M_PI * a) };
//             dv = -(cp * g0 * du / M_PI + std::sin(M_PI * a) * dv);
//         }
//     }

//     return dv;
// }
// template <class T>
// std::complex<T> PCF::cgamma(T x, T y, int kf)
// {
// //% ================================================================== =
// //% Purpose: Compute complex gamma function Gamma(z) or Ln[Gamma(z)]
// //% Input : x-- - Real part of z
// //% y-- - Imaginary part of z
// //% kf-- - Function code
// //% kf = 0 for Ln[Gamma(z)]
// //% kf = 1 for Gamma(z)
// //% Output : gr-- - Real part of Ln[Gamma(z)] or Gamma(z)
// //% gi-- - Imaginary part of Ln[Gamma(z)] or Gamma(z)
// //% ================================================================== =
//     T x1 = 0.0;
//     // Bernoulli's numbers B2n divided by [2*n(2*n-1)], n = 1,2,... 
//     std::array<T, 10> B = { 1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0,
//         -1.0 / 1680.0, 1.0 / 1188.0, -691.0 / 360360.0,
//         7.0 / 1092.0, -3617.0 / 122400.0, 43867.0 / 244188.0,
//         -174611.0 / 125400.0 };
//     // Compute firstly Gamma(| x | , sign(x) * y)
//     if (x < 0.0) 
//     {
//         x1 = x;
//         x = -x;
//         y = -y;
//     }
//     T x0 = x;
//     long long int na{};// default value 0
//     if (x <= 7.0)
//     {
//         na = (long long int)(std::trunc(7.0 - x));
//         x0 = x + T(na);
//     }
//     //% Compute log[Gamma(| x | +na, sign(x) * y)] by using relation for | z|>>1:
//     //% log[Gamma(z)] = (z - 1 / 2) * log(z) - z + (1 / 2) * log(2 * pi) + ...
//     //%     sum{ B2n / [2 * n * (2 * n - 1) * z ^ (2 * n - 1)] }, where summation is from n = 1 to Inf
//     T z1 = std::sqrt(x0 * x0 + y * y);
//     T th = std::atan(y / x0);
//     T gr = (x0 - 0.5) * std::log(z1) - th * y - x0 + 0.5 * std::log(2.0 * M_PI);
//     T gi = th * (x0 - 0.5) + y * std::log(z1) - y;
//     for (long long int k = 1; k <= 10; ++k)
//     {
//         T t = std::pow(z1 , (1.0 - 2.0 * T(k)));
//         gr = gr + B[k-1] * t * std::cos((2.0 * T(k) - 1.0) * th);// with B(n) = B2n / [2 * n * (2 * n - 1)]
//         gi = gi - B[k - 1] * t * std::sin((2.0 * T(k) - 1.0) * th);
//     }
//     //% Compute log[Gamma(| x | , sign(x) * y)] from log[Gamma(| x | +na, sign(x) * y)]
//     //% by using recurrence relation : Gamma(z + n) = z * (z + 1)*...(z + n - 1)* Gamma(z)
//     T gr1{};
//     T gi1{};
//     if (x <= 7.0)
//     {
//         gr1 = 0.0;
//         gi1 = 0.0;
//         for (long long int j = 0; j < na; ++j)
//         {
//             gr1 = gr1 + 0.5 * std::log((x + T(j))* (x + T(j)) + y * y);
//             gi1 = gi1 + std::atan(y / (x + T(j)));
//         }
//         gr = gr - gr1;
//         gi = gi - gi1;
//     }
//     //% If x < 0, compute log[Gamma(z)] by using relation:
//     //% Gamma(z)* Gamma(-z) = -pi / [z * sin(pi * z)]
//     T th1{}, sr{}, si{}, z2{}, th2{};
//     if (x1 < 0.0)
//     {
//         z1 = std::sqrt(x * x + y * y);
//         th1 = std::atan(y / x);
//         sr = -std::sin(M_PI * x) * std::cosh(M_PI * y);
//         si = -std::cos(M_PI * x) * std::sinh(M_PI * y);
//         z2 = std::sqrt(sr * sr + si * si);
//         th2 = std::atan(si / sr);
//         if (sr < 0.0)
//         {
//             th2 = M_PI + th2;
//         }
//         gr = std::log(M_PI / (z1 * z2)) - gr;
//         gi = -th1 - th2 - gi;
//     }
//     //% Compute Gamma(z) from log[Gamma(z)] by using relations:
//     //% | Gamma(z) |= exp{ real{log[Gamma(z)]} }; arg[Gamma(z)] = imag{ log[Gamma(z)] }
//     if (kf == 1)
//     {
//         T g0 = std::exp(gr);
//         gr = g0 * std::cos(gi);
//         gi = g0 * std::sin(gi);
//     }
//     std::complex<T> z (gr,gi);
//     return z;
// }
// template <class T>
// T PCF::pw(T a, T x)
// {
// //% ================================================================== =
// //% Purpose: Compute parabolic cylinder function W(a, x)
// //% Input : a-- - Parameter(| a | < 5)
// //% x-- - Argument(| x | < 5)
// //% Output : w------W(a, x)
// //% Routine called :
// //% 'cgamma' for computing complex Gamma function
// //% ================================================================== =
// //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     std::array<T, 100> c{};
//     std::array<T, 100> d{};
//     T eps = 1.0e-15;
//     T p0 = std::pow(2.0,-0.75);
//     T x1 = 0.25;
//     T y1 = a / 2.0;
//     std::complex<T> z = PCF::cgamma(x1, y1, 1);
//     T g1 = std::abs(z);
//     T x2 = 0.75;
//     z = PCF::cgamma(x2, y1, 1);
//     T g3 = std::abs(z);
//     T f1 = std::sqrt(g1 / g3);
//     T f2 = sqrt(2.0 * g3 / g1);
//     T c0 = 1.0;
//     T c1 = a;
//     c[0] = a;
//     for (long long int k1 = 4; k1 <= 200; k1 = k1 + 2)
//     {
//         long long int m = k1 / 2;
//         T ck = a * c1 - (T(k1) - 2.0) * (T(k1) - 3.0) * c0 / 4.0;
//         c[m - 1] = ck;
//         c0 = c1;
//         c1 = ck;
//     }
//     y1 = 1.0;
//     T r = 1.0;
//     T r1{};
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) - 1.0));
//         r1 = c[k-1] * r;
//         y1 = y1 + r1;
//         if (std::abs(r1 / y1) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     T d1 = 1.0;
//     T d2 = a;
//     d[0] = 1.0;
//     d[1] = a;
//     T dk{};
//     for (long long int k2 = 5; k2 <= 160; k2 = k2 + 2)
//     {
//         long long int m = (long long int)(std::trunc((k2 + 1) / 2));
//         dk = a * d2 - 0.25 * (T(k2) - 2.0) * (T(k2) - 3.0) * d1;
//         d[m-1] = dk;
//         d1 = d2;
//         d2 = dk;
//     }
//     T y2 = 1.0;
//     r = 1.0;
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) + 1.0));
//         r1 = d[k] * r;
//         y2 = y2 + r1;
//         if (std::abs(r1 / y2) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     y2 = x * y2;
//     T w = p0 * (f1 * y1 - f2 * y2);
//     return w;
// }
// template <class T>
// T PCF::dpw(T a, T x)
// {
// //% ================================================================== =
// //% Purpose: Compute derivative of parabolic cylinder functions W(a, x)
// //% Input : a-- - Parameter(| a | < 5)
// //% x-- - Argument(| x | < 5)
// //% Output : dw------W'(a,x)                
// //% Routine called :
// //% 'cgamma' for computing complex Gamma function
// //% ================================================================== =
//     std::array<T, 101> c{};
//     std::array<T, 101> d{};
//     T eps = 1.0e-15;
//     T p0 = std::pow(2.0, -0.75);
//     T x1 = 0.25;
//     T y1 = a / 2.0;
//     std::complex<T> z = PCF::cgamma(x1, y1, 1);
//     T g1 = std::abs(z);
//     T x2 = 0.75;
//     z = PCF::cgamma(x2, y1, 1);
//     T g3 = std::abs(z);
//     T f1 = std::sqrt(g1 / g3);
//     T f2 = sqrt(2.0 * g3 / g1);
//     T c0 = 1.0;
//     T c1 = a;
//     c[0] = a;
//     T ck{};
//     for (long long int k1 = 4; k1 <= 202; k1 = k1 + 2)
//     {
//         long long int m = k1 / 2;
//         ck = a * c1 - (T(k1) - 2.0) * (T(k1) - 3.0) * c0 / 4.0;
//         c[m-1] = ck;
//         c0 = c1;
//         c1 = ck;
//     }
//     y1 = a;
//     T r = 1.0;
//     T r1{};
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) + 1.0));
//         r1 = c[k] * r;
//         y1 = y1 + r1;
//         if (abs(r1 / y1) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     y1 = x * y1;
//     T d1 = 1.0;
//     T d2 = a;
//     d[0] = 1.0;
//     d[1] = a;
//     T dk{};
//     for (long long int k2 = 5; k2 <= 160; k2 = k2 + 2)
//     {
//         long long int m = (long long int)(std::trunc((k2 + 1) / 2));
//         dk = a * d2 - 0.25 * (k2 - 2.0) * (k2 - 3.0) * d1;
//         d[m-1] = dk;
//         d1 = d2;
//         d2 = dk;
//     }
//     T y2 = 1.0;
//     r = 1.0;
//     for (long long int k = 1; k <= 100; ++k)
//     {
//         r = 0.5 * r * x * x / (T(k) * (2.0 * T(k) - 1.0));
//         r1 = d[k] * r;
//         y2 = y2 + r1;
//         if (abs(r1 / y2) <= eps && k > 30)
//         {
//             break;
//         }
//     }
//     T dw = p0 * (f1 * y1 - f2 * y2);
//     return dw;
// }
// template <class T>
// T PCF::pwlx(T a, T x)
// {
// //% ================================================================== =
// //% Purpose: Compute parabolic cylinder function W(a, x)
// //% for large argument x(| x| >> |a | and | a | moderate)
// //% Input : a-- - Parameter
// //% x-- - Argument
// //% Output : w------W(a, x)
// //% Routine called :
// //% 'cgamma' for computing complex Gamma function
// //% ================================================================== =
//     std::array<T, 100> u{};
//     std::array<T, 100> v{};
//     std::complex<T> g0 = PCF::cgamma(0.5, a, 0);
//     T phi2 = g0.imag();
//     g0 = PCF::cgamma(0.5, a, 1);
//     T den = std::norm(g0);
//     std::complex<T> gk;
//     T uk{}, vk{};
//     for (long long int k = 2; k <= 40; k = k + 2)
//     {
//         long long int m = (long long int)(std::trunc(k / 2));
//         gk = PCF::cgamma(T(k) + 0.5, a, 1);
//         uk = (gk.real() * g0.real() + gk.imag() * g0.imag()) / den;
//         vk = (g0.real() * gk.imag() - gk.real() * g0.imag()) / den;
//         u[m-1] = uk;
//         v[m-1] = vk;
//     }
//     T sv1 = v[0] / (2.0 * x * x);
//     T su2 = u[0] / (2.0 * x * x);
//     long long int fac = 1;
//     T fd{};
//     for (long long int k = 3; k <= 20; k = k + 2)
//     {
//         fac = -(k) * (k - 1) * fac;
//         fd = T(fac) * pow(2.0,T(k)) * pow(x , 2.0*T(k));
//         sv1 = sv1 + v[k-1] / fd;
//         su2 = su2 + u[k-1] / fd;
//     }
//     T sv2 = 0.0;
//     T su1 = 0.0;
//     fac = 1;
//     for (long long int k = 2; k <= 20; k = k + 2)
//     {
//         fac = -(k) * (k - 1) * fac;
//         fd = T(fac) * pow(2.0, T(k)) * pow(x, 2.0 * T(k));
//         sv2 = sv2 + v[k-1] / fd;
//         su1 = su1 + u[k-1] / fd;
//     }
//     T s1 = 1 + sv1 - su1;
//     T s2 = -sv2 - su2;
//     T ea = exp(M_PI * a);
//     T sea = sqrt(1.0 + ea * ea);
//     T fk = sea - ea;
//     T ifk = sea + ea;
//     T x1 = abs(x);
//     T fa = x1 * x1 / 4.0 - a * log(x1) + M_PI / 4.0 + phi2 / 2.0;
//     T w{};
//     if (x > 0.0)
//     {
//         w = sqrt(2.0 * fk / x1) * (s1 * cos(fa) - s2 * sin(fa));
//     }
//     else if (x < 0.0)
//     {
//         w = sqrt(-2.0 * ifk / x) * (s1 * sin(fa) + s2 * cos(fa));
//     }
//     return w;
// }
// template <class T>
// T PCF::dpwlx(T a, T x)
// {
//     //% ================================================================== =
//     //% Purpose: Compute parabolic cylinder function W(a, x)
//     //% for large argument x(| x| >> |a | and | a | moderate)
//     //% Input : a-- - Parameter
//     //% x-- - Argument
//     //% Output : w------W(a, x)
//     //% Routine called :
//     //% 'cgamma' for computing complex Gamma function
//     //% ================================================================== =
//     std::array<T, 100> u{};
//     std::array<T, 100> v{};
//     std::complex<T> g0 = PCF::cgamma(0.5, a, 0);
//     T phi2 = g0.imag();
//     g0 = PCF::cgamma(0.5, a, 1);
//     T den = std::norm(g0);
//     std::complex<T> gk;
//     T uk{}, vk{};
//     for (long long int k = 2; k <= 40; k = k + 2)
//     {
//         long long int m = (long long int)(std::trunc(k / 2));
//         gk = PCF::cgamma(T(k) + 0.5, a, 1);
//         uk = (gk.real() * g0.real() + gk.imag() * g0.imag()) / den;
//         vk = (g0.real() * gk.imag() - gk.real() * g0.imag()) / den;
//         u[m - 1] = uk;
//         v[m - 1] = vk;
//     }
//     T x1 = abs(x);
//     T sv1 = v[0] / (2.0 * x * x);
//     T su2 = u[0] / (2.0 * x * x);
//     long long int fac = 1;
//     T fd{}, fdd{};
//     T dsv1 = -2.0 * sv1 / x1;
//     T dsu2 = -2.0 * su2 / x1;
//     for (long long int k = 3; k <= 20; k = k + 2)
//     {
//         fac = -(k) * (k - 1) * fac;
//         fd = T(fac) * pow(2.0, T(k)) * pow(x, 2.0 * T(k));
//         fdd = -2.0 * T(k) / x1;
//         sv1 = sv1 + v[k - 1] / fd;
//         su2 = su2 + u[k - 1] / fd;
//         dsv1 = dsv1 + fdd * v[k-1] / fd;
//         dsu2 = dsu2 + fdd * u[k-1] / fd;
//     }
//     T sv2 = 0.0;
//     T su1 = 0.0;
//     T dsv2 = 0.0;
//     T dsu1 = 0.0;
//     fac = 1;
//     for (long long int k = 2; k <= 20; k = k + 2)
//     {
//         fac = -(k) * (k - 1) * fac;
//         fd = T(fac) * pow(2.0, T(k)) * pow(x, 2.0 * T(k));
//         fdd = -2.0 * T(k) / x1;
//         sv2 = sv2 + v[k - 1] / fd;
//         su1 = su1 + u[k - 1] / fd;
//         dsv2 = dsv2 + fdd * v[k-1] / fd;
//         dsu1 = dsu1 + fdd * u[k-1] / fd;
//     }
//     T s1 = 1 + sv1 - su1;
//     T s2 = -sv2 - su2;
//     T ds1 = dsv1 - dsu1;
//     T ds2 = -dsv2 - dsu2;
//     T ea = exp(M_PI * a);
//     T sea = sqrt(1.0 + ea * ea);
//     T fk = sea - ea;
//     T ifk = sea + ea;
//     T fa = x1 * x1 / 4.0 - a * log(x1) + M_PI / 4.0 + phi2 / 2.0;
//     T dfa = x1 / 2.0 - a / x1;
//     T w{}, dw{};
//     if (x > 0.0)
//     {
//         w = sqrt(2.0 * fk / x1) * (s1 * cos(fa) - s2 * sin(fa));
//         dw = -w / (2.0 * x1) + sqrt(2.0 * fk / x1) * ((ds1 - s2 * dfa) * cos(fa) - (s1 * dfa + ds2) * sin(fa));
//     }
//     else if (x < 0.0)
//     {
//         w = sqrt(-2.0 * ifk / x) * (s1 * sin(fa) + s2 * cos(fa));
//         dw = -(-w / (-2.0 * x) + sqrt(-2.0 * ifk / x) * ((ds1 - s2 * dfa) * sin(fa) + (s1 * dfa + ds2) * cos(fa)));
//     }
//     return dw;
// }
// template <class T>
// T PCF::call_pcf_name(T a, T x, int idx)
// {
//     T val{};
//     switch (idx)
//     {
//     case 0:
//         val = PCF::pu(a, x);
//         break;
//     case 1:
//         val = PCF::dpu(a, x);
//         break;
//     case 2:
//         val = PCF::pulx(a, x);
//         break;
//     case 3:
//         val = PCF::dpulx(a, x);
//         break;
//     case 4:
//         val = PCF::D(a, x);
//         break;
//     case 5:
//         val = PCF::dD_dx(a, x);
//         break;
//     case 6:
//         val = PCF::pv(a, x);
//         break;
//     case 7:
//         val = PCF::dpv(a, x);
//         break;
//     case 8:
//         val = PCF::pvlx(a, x);
//         break;
//     case 9:
//         val = PCF::dpvlx(a, x);
//         break;
//     case 10:
//         val = PCF::pw(a, x);
//         break;
//     case 11:
//         val = PCF::dpw(a, x);
//         break;
//     case 12:
//         val = PCF::pwlx(a, x);
//         break;
//     case 13:
//         val = PCF::dpwlx(a, x);
//         break;
//     }
//     return val;
// }
// template <class T>
// T PCF::call_D(T x, T a)
// {
//     return PCF::pu(-a - 0.5, x);
// }
// template <class T>
// T PCF::call_dD(T x, T a)
// {
//     return PCF::dpu(-a - 0.5, x);
// }

// template <class T>
// void PCF::call_pcf_name_array(T a, T* x, int idx, int x_size, T* val)
// {
//     T (*fcnPtr)(T, T) { &PCF::pu };
//     switch (idx)
//     {
//     case 0:
//         fcnPtr = &PCF::pu;
//         break;
//     case 1:
//         fcnPtr = &PCF::dpu;
//         break;
//     case 2:
//         fcnPtr = &PCF::pulx;
//         break;
//     case 3:
//         fcnPtr = &PCF::dpulx;
//         break;
//     case 4:
//         fcnPtr = &PCF::D;
//         break;
//     case 5:
//         fcnPtr = &PCF::dD_dx;
//         break;
//     case 6:
//         fcnPtr = &PCF::pv;
//         break;
//     case 7:
//         fcnPtr = &PCF::dpv;
//         break;
//     case 8:
//         fcnPtr = &PCF::pvlx;
//         break;
//     case 9:
//         fcnPtr = &PCF::dpvlx;
//         break;
//     case 10:
//         fcnPtr = &PCF::pw;
//         break;
//     case 11:
//         fcnPtr = &PCF::dpw;
//         break;
//     case 12:
//         fcnPtr = &PCF::pwlx;
//         break;
//     case 13:
//         fcnPtr = &PCF::dpwlx;
//         break;
//     }
// #pragma omp parallel for
//     for (int i = 0; i < x_size; ++i)
//     {
//         *(val + i) = fcnPtr(a, *(x + i));
//     }
// }
