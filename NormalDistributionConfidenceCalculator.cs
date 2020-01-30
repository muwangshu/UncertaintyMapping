using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SegmentGeneration
{
    public static class NormalDistributionConfidenceCalculator
    {
        /// <summary>
        /// 
        /// </summary>
        public static double InverseNormalDistribution(double probability, double min, double max)
        {
            double x = 0;
            double a = 0;
            double b = 1;

            double precision = Math.Pow(10, -3);

            while ((b - a) > precision)
            {
                x = (a + b) / 2;
                if (NormInv(x) > probability)
                {
                    b = x;
                }
                else
                {
                    a = x;
                }
            }

            if ((max > 0) && (min > 0))
            {
                x = x * (max - min) + min;
            }
            return x;
        }

        /// <summary>
        /// Returns the cumulative density function evaluated at A given value.
        /// </summary>
        /// <param name="x">A position on the x-axis.</param>
        /// <param name="mean"></param>
        /// <param name="sigma"></param>
        /// <returns>The cumulative density function evaluated at <C>x</C>.</returns>
        /// <remarks>The value of the cumulative density function at A point <C>x</C> is
        /// probability that the value of A random variable having this normal density is
        /// less than or equal to <C>x</C>.
        /// </remarks>
        public static double NormalDistribution(double x, double mean, double sigma)
        {
            // This algorithm is ported from dcdflib:
            // Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
            // Package of Special Function Routines and Test Drivers"
            // acm Transactions on Mathematical Software. 19, 22-32.
            int i;
            double del, xden, xnum, xsq;
            double result, ccum;
            double arg = (x - mean) / sigma;
            const double sixten = 1.60e0;
            const double sqrpi = 3.9894228040143267794e-1;
            const double thrsh = 0.66291e0;
            const double root32 = 5.656854248e0;
            const double zero = 0.0e0;
            const double min = Double.Epsilon;
            double z = arg;
            double y = Math.Abs(z);
            const double half = 0.5e0;
            const double one = 1.0e0;

            double[] a =
                {
                2.2352520354606839287e00, 1.6102823106855587881e02, 1.0676894854603709582e03,
                1.8154981253343561249e04, 6.5682337918207449113e-2
            };

            double[] b =
                {
                4.7202581904688241870e01, 9.7609855173777669322e02, 1.0260932208618978205e04,
                4.5507789335026729956e04
            };

            double[] c =
                {
                3.9894151208813466764e-1, 8.8831497943883759412e00, 9.3506656132177855979e01,
                5.9727027639480026226e02, 2.4945375852903726711e03, 6.8481904505362823326e03,
                1.1602651437647350124e04, 9.8427148383839780218e03, 1.0765576773720192317e-8
            };

            double[] d =
                {
                2.2266688044328115691e01, 2.3538790178262499861e02, 1.5193775994075548050e03,
                6.4855582982667607550e03, 1.8615571640885098091e04, 3.4900952721145977266e04,
                3.8912003286093271411e04, 1.9685429676859990727e04
            };
            double[] p =
                {
                2.1589853405795699e-1, 1.274011611602473639e-1, 2.2235277870649807e-2,
                1.421619193227893466e-3, 2.9112874951168792e-5, 2.307344176494017303e-2
            };


            double[] q =
                {
                1.28426009614491121e00, 4.68238212480865118e-1, 6.59881378689285515e-2,
                3.78239633202758244e-3, 7.29751555083966205e-5
            };
            if (y <= thrsh)
            {
                //
                // Evaluate  anorm  for  |X| <= 0.66291
                //
                xsq = zero;
                if (y > double.Epsilon) xsq = z * z;
                xnum = a[4] * xsq;
                xden = xsq;
                for (i = 0; i < 3; i++)
                {
                    xnum = (xnum + a[i]) * xsq;
                    xden = (xden + b[i]) * xsq;
                }
                result = z * (xnum + a[3]) / (xden + b[3]);
                double temp = result;
                result = half + temp;
            }

            //
            // Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
            //
            else if (y <= root32)
            {
                xnum = c[8] * y;
                xden = y;
                for (i = 0; i < 7; i++)
                {
                    xnum = (xnum + c[i]) * y;
                    xden = (xden + d[i]) * y;
                }
                result = (xnum + c[7]) / (xden + d[7]);
                xsq = Math.Floor(y * sixten) / sixten;
                del = (y - xsq) * (y + xsq);
                result = Math.Exp(-(xsq * xsq * half)) * Math.Exp(-(del * half)) * result;
                ccum = one - result;
                if (z > zero)
                {
                    result = ccum;
                }
            }

            //
            // Evaluate  anorm  for |X| > sqrt(32)
            //
            else
            {
                xsq = one / (z * z);
                xnum = p[5] * xsq;
                xden = xsq;
                for (i = 0; i < 4; i++)
                {
                    xnum = (xnum + p[i]) * xsq;
                    xden = (xden + q[i]) * xsq;
                }
                result = xsq * (xnum + p[4]) / (xden + q[4]);
                result = (sqrpi - result) / y;
                xsq = Math.Floor(z * sixten) / sixten;
                del = (z - xsq) * (z + xsq);
                result = Math.Exp(-(xsq * xsq * half)) * Math.Exp(-(del * half)) * result;
                ccum = one - result;
                if (z > zero)
                {
                    result = ccum;
                }
            }

            if (result < min)
                result = 0.0e0;
            return result;
        }

        /// <summary>
        /// Given a probability, a mean, and a standard deviation, an x value can be calculated.
        /// </summary>
        /// <returns></returns>
        public static double NormInv(double probability)
        {
            const double a1 = -39.6968302866538;
            const double a2 = 220.946098424521;
            const double a3 = -275.928510446969;
            const double a4 = 138.357751867269;
            const double a5 = -30.6647980661472;
            const double a6 = 2.50662827745924;

            const double b1 = -54.4760987982241;
            const double b2 = 161.585836858041;
            const double b3 = -155.698979859887;
            const double b4 = 66.8013118877197;
            const double b5 = -13.2806815528857;

            const double c1 = -7.78489400243029E-03;
            const double c2 = -0.322396458041136;
            const double c3 = -2.40075827716184;
            const double c4 = -2.54973253934373;
            const double c5 = 4.37466414146497;
            const double c6 = 2.93816398269878;

            const double d1 = 7.78469570904146E-03;
            const double d2 = 0.32246712907004;
            const double d3 = 2.445134137143;
            const double d4 = 3.75440866190742;

            //Define break-points
            // using Epsilon is wrong; see link above for reference to 0.02425 value
            //const double pLow = double.Epsilon;
            const double pLow = 0.02425;

            const double pHigh = 1 - pLow;

            //Define work variables
            double q;
            double result = 0;

            // if argument out of bounds.
            // set it to a value within desired precision.
            if (probability <= 0)
                probability = pLow;

            if (probability >= 1)
                probability = pHigh;

            if (probability < pLow)
            {
                //Rational approximation for lower region
                q = Math.Sqrt(-2 * Math.Log(probability));
                result = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
            }
            else if (probability <= pHigh)
            {
                //Rational approximation for lower region
                q = probability - 0.5;
                double r = q * q;
                result = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
                         (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
            }
            else if (probability < 1)
            {
                //Rational approximation for upper region
                q = Math.Sqrt(-2 * Math.Log(1 - probability));
                result = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
            }

            return result;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="probability"></param>
        /// <param name="mean"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double NormInv(double probability, double mean, double sigma)
        {
            double x = NormInv(probability);
            return sigma * x + mean;
        }


        /// <summary>
        /// Returns the value of the gaussian error function at <paramref name="x"/>.
        /// </summary>
        public static double Erf(double x)
        {
            /*
            Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
            *
            * Developed at SunPro, a Sun Microsystems, Inc. business.
            * Permission to use, copy, modify, and distribute this
            * software is freely granted, provided that this notice
            * is preserved.
            */

            #region Constants

            const double tiny = 1e-300;
            const double erx = 8.45062911510467529297e-01;

            // Coefficients for approximation to erf on [0, 0.84375]
            const double efx = 1.28379167095512586316e-01; /* 0x3FC06EBA; 0x8214DB69 */
            const double efx8 = 1.02703333676410069053e+00; /* 0x3FF06EBA; 0x8214DB69 */
            const double pp0 = 1.28379167095512558561e-01; /* 0x3FC06EBA; 0x8214DB68 */
            const double pp1 = -3.25042107247001499370e-01; /* 0xBFD4CD7D; 0x691CB913 */
            const double pp2 = -2.84817495755985104766e-02; /* 0xBF9D2A51; 0xDBD7194F */
            const double pp3 = -5.77027029648944159157e-03; /* 0xBF77A291; 0x236668E4 */
            const double pp4 = -2.37630166566501626084e-05; /* 0xBEF8EAD6; 0x120016AC */
            const double qq1 = 3.97917223959155352819e-01; /* 0x3FD97779; 0xCDDADC09 */
            const double qq2 = 6.50222499887672944485e-02; /* 0x3FB0A54C; 0x5536CEBA */
            const double qq3 = 5.08130628187576562776e-03; /* 0x3F74D022; 0xC4D36B0F */
            const double qq4 = 1.32494738004321644526e-04; /* 0x3F215DC9; 0x221C1A10 */
            const double qq5 = -3.96022827877536812320e-06; /* 0xBED09C43; 0x42A26120 */

            // Coefficients for approximation to erf in [0.84375, 1.25]
            const double pa0 = -2.36211856075265944077e-03; /* 0xBF6359B8; 0xBEF77538 */
            const double pa1 = 4.14856118683748331666e-01; /* 0x3FDA8D00; 0xAD92B34D */
            const double pa2 = -3.72207876035701323847e-01; /* 0xBFD7D240; 0xFBB8C3F1 */
            const double pa3 = 3.18346619901161753674e-01; /* 0x3FD45FCA; 0x805120E4 */
            const double pa4 = -1.10894694282396677476e-01; /* 0xBFBC6398; 0x3D3E28EC */
            const double pa5 = 3.54783043256182359371e-02; /* 0x3FA22A36; 0x599795EB */
            const double pa6 = -2.16637559486879084300e-03; /* 0xBF61BF38; 0x0A96073F */
            const double qa1 = 1.06420880400844228286e-01; /* 0x3FBB3E66; 0x18EEE323 */
            const double qa2 = 5.40397917702171048937e-01; /* 0x3FE14AF0; 0x92EB6F33 */
            const double qa3 = 7.18286544141962662868e-02; /* 0x3FB2635C; 0xD99FE9A7 */
            const double qa4 = 1.26171219808761642112e-01; /* 0x3FC02660; 0xE763351F */
            const double qa5 = 1.36370839120290507362e-02; /* 0x3F8BEDC2; 0x6B51DD1C */
            const double qa6 = 1.19844998467991074170e-02; /* 0x3F888B54; 0x5735151D */

            // Coefficients for approximation to erfc in [1.25, 1/0.35]
            const double ra0 = -9.86494403484714822705e-03; /* 0xBF843412; 0x600D6435 */
            const double ra1 = -6.93858572707181764372e-01; /* 0xBFE63416; 0xE4BA7360 */
            const double ra2 = -1.05586262253232909814e+01; /* 0xC0251E04; 0x41B0E726 */
            const double ra3 = -6.23753324503260060396e+01; /* 0xC04F300A; 0xE4CBA38D */
            const double ra4 = -1.62396669462573470355e+02; /* 0xC0644CB1; 0x84282266 */
            const double ra5 = -1.84605092906711035994e+02; /* 0xC067135C; 0xEBCCABB2 */
            const double ra6 = -8.12874355063065934246e+01; /* 0xC0545265; 0x57E4D2F2 */
            const double ra7 = -9.81432934416914548592e+00; /* 0xC023A0EF; 0xC69AC25C */
            const double sa1 = 1.96512716674392571292e+01; /* 0x4033A6B9; 0xBD707687 */
            const double sa2 = 1.37657754143519042600e+02; /* 0x4061350C; 0x526AE721 */
            const double sa3 = 4.34565877475229228821e+02; /* 0x407B290D; 0xD58A1A71 */
            const double sa4 = 6.45387271733267880336e+02; /* 0x40842B19; 0x21EC2868 */
            const double sa5 = 4.29008140027567833386e+02; /* 0x407AD021; 0x57700314 */
            const double sa6 = 1.08635005541779435134e+02; /* 0x405B28A3; 0xEE48AE2C */
            const double sa7 = 6.57024977031928170135e+00; /* 0x401A47EF; 0x8E484A93 */
            const double sa8 = -6.04244152148580987438e-02; /* 0xBFAEEFF2; 0xEE749A62 */

            // Coefficients for approximation to erfc in [1/0.35, 28]
            const double rb0 = -9.86494292470009928597e-03; /* 0xBF843412; 0x39E86F4A */
            const double rb1 = -7.99283237680523006574e-01; /* 0xBFE993BA; 0x70C285DE */
            const double rb2 = -1.77579549177547519889e+01; /* 0xC031C209; 0x555F995A */
            const double rb3 = -1.60636384855821916062e+02; /* 0xC064145D; 0x43C5ED98 */
            const double rb4 = -6.37566443368389627722e+02; /* 0xC083EC88; 0x1375F228 */
            const double rb5 = -1.02509513161107724954e+03; /* 0xC0900461; 0x6A2E5992 */
            const double rb6 = -4.83519191608651397019e+02; /* 0xC07E384E; 0x9BDC383F */
            const double sb1 = 3.03380607434824582924e+01; /* 0x403E568B; 0x261D5190 */
            const double sb2 = 3.25792512996573918826e+02; /* 0x40745CAE; 0x221B9F0A */
            const double sb3 = 1.53672958608443695994e+03; /* 0x409802EB; 0x189D5118 */
            const double sb4 = 3.19985821950859553908e+03; /* 0x40A8FFB7; 0x688C246A */
            const double sb5 = 2.55305040643316442583e+03; /* 0x40A3F219; 0xCEDF3BE6 */
            const double sb6 = 4.74528541206955367215e+02; /* 0x407DA874; 0xE79FE763 */
            const double sb7 = -2.24409524465858183362e+01; /* 0xC03670E2; 0x42712D62 */

            #endregion

            if (double.IsNaN(x))
                return double.NaN;

            if (double.IsNegativeInfinity(x))
                return -1.0;

            if (double.IsPositiveInfinity(x))
                return 1.0;

            int n0, hx, ix, i;
            double R, S, P, Q, s, y, z, r;
            unsafe
            {
                double one = 1.0;
                n0 = ((*(int*)&one) >> 29) ^ 1;
                hx = *(n0 + (int*)&x);
            }
            ix = hx & 0x7FFFFFFF;

            if (ix < 0x3FEB0000) // |x| < 0.84375
            {
                if (ix < 0x3E300000) // |x| < 2**-28
                {
                    if (ix < 0x00800000)
                        return 0.125 * (8.0 * x + efx8 * x); // avoid underflow
                    return x + efx * x;
                }
                z = x * x;
                r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
                s = 1.0 + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
                y = r / s;
                return x + x * y;
            }
            if (ix < 0x3FF40000) // 0.84375 <= |x| < 1.25
            {
                s = Math.Abs(x) - 1.0;
                P = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
                Q = 1.0 + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
                if (hx >= 0)
                    return erx + P / Q;
                else
                    return -erx - P / Q;
            }
            if (ix >= 0x40180000) // inf > |x| >= 6
            {
                if (hx >= 0)
                    return 1.0 - tiny;
                else
                    return tiny - 1.0;
            }
            x = Math.Abs(x);
            s = 1.0 / (x * x);
            if (ix < 0x4006DB6E) // |x| < 1/0.35
            {
                R = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
                S = 1.0 + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 + s * sa8)))))));
            }
            else // |x| >= 1/0.35
            {
                R = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
                S = 1.0 + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
            }
            z = x;
            unsafe { *(1 - n0 + (int*)&z) = 0; }
            r = Math.Exp(-z * z - 0.5625) * Math.Exp((z - x) * (z + x) + R / S);
            if (hx >= 0)
                return 1.0 - r / x;
            else
                return r / x - 1.0;
        }

        /// <summary>
        /// Returns the value of the complementary error function at <paramref name="x"/>.
        /// </summary>
        public static double Erfc(double x)
        {
            /*
            Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
            *
            * Developed at SunPro, a Sun Microsystems, Inc. business.
            * Permission to use, copy, modify, and distribute this
            * software is freely granted, provided that this notice
            * is preserved.
            */

            #region Constants

            const double tiny = 1e-300;
            const double erx = 8.45062911510467529297e-01;

            // Coefficients for approximation to erf on [0, 0.84375]
            const double efx = 1.28379167095512586316e-01; /* 0x3FC06EBA; 0x8214DB69 */
            const double efx8 = 1.02703333676410069053e+00; /* 0x3FF06EBA; 0x8214DB69 */
            const double pp0 = 1.28379167095512558561e-01; /* 0x3FC06EBA; 0x8214DB68 */
            const double pp1 = -3.25042107247001499370e-01; /* 0xBFD4CD7D; 0x691CB913 */
            const double pp2 = -2.84817495755985104766e-02; /* 0xBF9D2A51; 0xDBD7194F */
            const double pp3 = -5.77027029648944159157e-03; /* 0xBF77A291; 0x236668E4 */
            const double pp4 = -2.37630166566501626084e-05; /* 0xBEF8EAD6; 0x120016AC */
            const double qq1 = 3.97917223959155352819e-01; /* 0x3FD97779; 0xCDDADC09 */
            const double qq2 = 6.50222499887672944485e-02; /* 0x3FB0A54C; 0x5536CEBA */
            const double qq3 = 5.08130628187576562776e-03; /* 0x3F74D022; 0xC4D36B0F */
            const double qq4 = 1.32494738004321644526e-04; /* 0x3F215DC9; 0x221C1A10 */
            const double qq5 = -3.96022827877536812320e-06; /* 0xBED09C43; 0x42A26120 */

            // Coefficients for approximation to erf in [0.84375, 1.25]
            const double pa0 = -2.36211856075265944077e-03; /* 0xBF6359B8; 0xBEF77538 */
            const double pa1 = 4.14856118683748331666e-01; /* 0x3FDA8D00; 0xAD92B34D */
            const double pa2 = -3.72207876035701323847e-01; /* 0xBFD7D240; 0xFBB8C3F1 */
            const double pa3 = 3.18346619901161753674e-01; /* 0x3FD45FCA; 0x805120E4 */
            const double pa4 = -1.10894694282396677476e-01; /* 0xBFBC6398; 0x3D3E28EC */
            const double pa5 = 3.54783043256182359371e-02; /* 0x3FA22A36; 0x599795EB */
            const double pa6 = -2.16637559486879084300e-03; /* 0xBF61BF38; 0x0A96073F */
            const double qa1 = 1.06420880400844228286e-01; /* 0x3FBB3E66; 0x18EEE323 */
            const double qa2 = 5.40397917702171048937e-01; /* 0x3FE14AF0; 0x92EB6F33 */
            const double qa3 = 7.18286544141962662868e-02; /* 0x3FB2635C; 0xD99FE9A7 */
            const double qa4 = 1.26171219808761642112e-01; /* 0x3FC02660; 0xE763351F */
            const double qa5 = 1.36370839120290507362e-02; /* 0x3F8BEDC2; 0x6B51DD1C */
            const double qa6 = 1.19844998467991074170e-02; /* 0x3F888B54; 0x5735151D */

            // Coefficients for approximation to erfc in [1.25, 1/0.35]
            const double ra0 = -9.86494403484714822705e-03; /* 0xBF843412; 0x600D6435 */
            const double ra1 = -6.93858572707181764372e-01; /* 0xBFE63416; 0xE4BA7360 */
            const double ra2 = -1.05586262253232909814e+01; /* 0xC0251E04; 0x41B0E726 */
            const double ra3 = -6.23753324503260060396e+01; /* 0xC04F300A; 0xE4CBA38D */
            const double ra4 = -1.62396669462573470355e+02; /* 0xC0644CB1; 0x84282266 */
            const double ra5 = -1.84605092906711035994e+02; /* 0xC067135C; 0xEBCCABB2 */
            const double ra6 = -8.12874355063065934246e+01; /* 0xC0545265; 0x57E4D2F2 */
            const double ra7 = -9.81432934416914548592e+00; /* 0xC023A0EF; 0xC69AC25C */
            const double sa1 = 1.96512716674392571292e+01; /* 0x4033A6B9; 0xBD707687 */
            const double sa2 = 1.37657754143519042600e+02; /* 0x4061350C; 0x526AE721 */
            const double sa3 = 4.34565877475229228821e+02; /* 0x407B290D; 0xD58A1A71 */
            const double sa4 = 6.45387271733267880336e+02; /* 0x40842B19; 0x21EC2868 */
            const double sa5 = 4.29008140027567833386e+02; /* 0x407AD021; 0x57700314 */
            const double sa6 = 1.08635005541779435134e+02; /* 0x405B28A3; 0xEE48AE2C */
            const double sa7 = 6.57024977031928170135e+00; /* 0x401A47EF; 0x8E484A93 */
            const double sa8 = -6.04244152148580987438e-02; /* 0xBFAEEFF2; 0xEE749A62 */

            // Coefficients for approximation to erfc in [1/0.35, 28]
            const double rb0 = -9.86494292470009928597e-03; /* 0xBF843412; 0x39E86F4A */
            const double rb1 = -7.99283237680523006574e-01; /* 0xBFE993BA; 0x70C285DE */
            const double rb2 = -1.77579549177547519889e+01; /* 0xC031C209; 0x555F995A */
            const double rb3 = -1.60636384855821916062e+02; /* 0xC064145D; 0x43C5ED98 */
            const double rb4 = -6.37566443368389627722e+02; /* 0xC083EC88; 0x1375F228 */
            const double rb5 = -1.02509513161107724954e+03; /* 0xC0900461; 0x6A2E5992 */
            const double rb6 = -4.83519191608651397019e+02; /* 0xC07E384E; 0x9BDC383F */
            const double sb1 = 3.03380607434824582924e+01; /* 0x403E568B; 0x261D5190 */
            const double sb2 = 3.25792512996573918826e+02; /* 0x40745CAE; 0x221B9F0A */
            const double sb3 = 1.53672958608443695994e+03; /* 0x409802EB; 0x189D5118 */
            const double sb4 = 3.19985821950859553908e+03; /* 0x40A8FFB7; 0x688C246A */
            const double sb5 = 2.55305040643316442583e+03; /* 0x40A3F219; 0xCEDF3BE6 */
            const double sb6 = 4.74528541206955367215e+02; /* 0x407DA874; 0xE79FE763 */
            const double sb7 = -2.24409524465858183362e+01; /* 0xC03670E2; 0x42712D62 */

            #endregion

            if (double.IsNaN(x))
                return double.NaN;

            if (double.IsNegativeInfinity(x))
                return 2.0;

            if (double.IsPositiveInfinity(x))
                return 0.0;

            int n0, hx, ix;
            double R, S, P, Q, s, y, z, r;
            unsafe
            {
                double one = 1.0;
                n0 = ((*(int*)&one) >> 29) ^ 1;
                hx = *(n0 + (int*)&x);
            }
            ix = hx & 0x7FFFFFFF;

            if (ix < 0x3FEB0000) // |x| < 0.84375
            {
                if (ix < 0x3C700000) // |x| < 2**-56
                    return 1.0 - x;
                z = x * x;
                r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
                s = 1.0 + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
                y = r / s;
                if (hx < 0x3FD00000) // x < 1/4
                    return 1.0 - (x + x * y);
                else
                {
                    r = x * y;
                    r += (x - 0.5);
                    return 0.5 - r;
                }
            }
            if (ix < 0x3FF40000) // 0.84375 <= |x| < 1.25
            {
                s = Math.Abs(x) - 1.0;
                P = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
                Q = 1.0 + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
                if (hx >= 0)
                {
                    z = 1.0 - erx;
                    return z - P / Q;
                }
                else
                {
                    z = erx + P / Q;
                    return 1.0 + z;
                }
            }
            if (ix < 0x403C0000) // |x| < 28
            {
                x = Math.Abs(x);
                s = 1.0 / (x * x);
                if (ix < 0x4006DB6D) // |x| < 1/.35 ~ 2.857143
                {
                    R = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
                    S = 1.0 + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 + s * sa8)))))));
                }
                else // |x| >= 1/.35 ~ 2.857143
                {
                    if (hx < 0 && ix >= 0x40180000)
                        return 2.0 - tiny; // x < -6
                    R = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
                    S = 1.0 + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
                }
                z = x;
                unsafe { *(1 - n0 + (int*)&z) = 0; }
                r = Math.Exp(-z * z - 0.5625) *
                Math.Exp((z - x) * (z + x) + R / S);
                if (hx > 0)
                    return r / x;
                else
                    return 2.0 - r / x;
            }
            else
            {
                if (hx > 0)
                    return tiny * tiny;
                else
                    return 2.0 - tiny;
            }
        }

        public static double CumDensity(double z)
        {
            return 0.5 * (1.0 + Erf(z/Math.Sqrt(2)));
        }

        public static double CumDensity2(double z)
        {
            double p = 0.3275911;
            double a1 = 0.254829592;
            double a2 = -0.284496736;
            double a3 = 1.421413741;
            double a4 = -1.453152027;
            double a5 = 1.061405429;

            int sign;
            if (z < 0.0)
                sign = -1;
            else
                sign = 1;

            double x = Math.Abs(z) / Math.Sqrt(2.0);
            double t = 1.0 / (1.0 + p * x);
            double erf = 1.0 - (((((a5 * t + a4) * t) + a3)
              * t + a2) * t + a1) * t * Math.Exp(-x * x);
            return 0.5 * (1.0 + sign * erf);
        }
    }
}
