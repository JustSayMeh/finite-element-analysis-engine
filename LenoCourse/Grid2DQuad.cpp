#include "Grid2DQuad.h"
void fill_localmatrixM(double ms[81], double hr, double hz, double rk, double lambda)
{
	ms[0] = (0.00222222222222073 * hr * hr * hz + 0.0177777777777751 * hr * hz * rk) * lambda;
	ms[9] = ms[1] = (-1.58391818179856e-15 * hr * hr * hz + 0.00888888888889199 * hr * hz * rk) * lambda;
	ms[10] = (0.0355555555555584 * hr * hr * hz + 0.071111111111107 * hr * hz * rk) * lambda;
	ms[18] = ms[2] = (-0.00222222222222093 * hr * hr * hz - 0.004444444444446 * hr * hz * rk) * lambda;
	ms[19] = ms[11] = (0.00888888888889061 * hr * hr * hz + 0.00888888888889042 * hr * hz * rk) * lambda;
	ms[20] = (0.0155555555555558 * hr * hr * hz + 0.0177777777777767 * hr * hz * rk) * lambda;
	ms[27] = ms[3] = (0.00111111111111275 * hr * hr * hz + 0.00888888888888784 * hr * hz * rk) * lambda;
	ms[28] = ms[12] = (1.09542005096349e-15 * hr * hr * hz + 0.00444444444444582 * hr * hz * rk) * lambda;
	ms[29] = ms[21] = (-0.00111111111111057 * hr * hr * hz - 0.00222222222222291 * hr * hz * rk) * lambda;
	ms[30] = (0.00888888888888868 * hr * hr * hz + 0.0711111111111109 * hr * hz * rk) * lambda;
	ms[36] = ms[4] = (1.09542005096349e-15 * hr * hr * hz + 0.00444444444444582 * hr * hz * rk) * lambda;
	ms[37] = ms[13] = (0.0177777777777792 * hr * hr * hz + 0.0355555555555539 * hr * hz * rk) * lambda;
	ms[38] = ms[22] = (0.00444444444444292 * hr * hr * hz + 0.00444444444444513 * hr * hz * rk) * lambda;
	ms[39] = ms[31] = (-1.1842378929335e-16 * hr * hr * hz + 0.0355555555555556 * hr * hz * rk) * lambda;
	ms[40] = (32 * hr * hr * hz / 225 + 64 * hr * hz * rk / 225) * lambda;
	ms[45] = ms[5] = (-0.00111111111111057 * hr * hr * hz - 0.00222222222222291 * hr * hz * rk) * lambda;
	ms[46] = ms[14] = (0.00444444444444292 * hr * hr * hz + 0.00444444444444513 * hr * hz * rk) * lambda;
	ms[47] = ms[23] = (0.00777777777777788 * hr * hr * hz + 0.00888888888888892 * hr * hz * rk) * lambda;
	ms[48] = ms[32] = (-0.00888888888888891 * hr * hr * hz - 0.0177777777777778 * hr * hz * rk) * lambda;
	ms[49] = ms[41] = (0.0355555555555558 * hr * hr * hz + 0.0355555555555555 * hr * hz * rk) * lambda;
	ms[50] = (0.0622222222222222 * hr * hr * hz + 0.0711111111111111 * hr * hz * rk) * lambda;
	ms[54] = ms[6] = (-0.000555555555554987 * hr * hr * hz - 0.00444444444444481 * hr * hz * rk) * lambda;
	ms[55] = ms[15] = (3.77475828372553e-16 * hr * hr * hz - 0.00222222222222179 * hr * hz * rk) * lambda;
	ms[56] = ms[24] = (0.000555555555555741 * hr * hr * hz + 0.00111111111111089 * hr * hz * rk) * lambda;
	ms[57] = ms[33] = (0.00111111111111208 * hr * hr * hz + 0.00888888888889 * hr * hz * rk) * lambda;
	ms[58] = ms[42] = (6.51330841113425e-16 * hr * hr * hz + 0.00444444444444313 * hr * hz * rk) * lambda;
	ms[59] = ms[51] = (-0.00111111111111167 * hr * hr * hz - 0.00222222222222156 * hr * hz * rk) * lambda;
	ms[60] = (0.00222222222222217 * hr * hr * hz + 0.0177777777777777 * hr * hz * rk) * lambda;
	ms[63] = ms[7] = (3.77475828372553e-16 * hr * hr * hz - 0.00222222222222179 * hr * hz * rk) * lambda;
	ms[64] = ms[16] = (-0.00888888888888845 * hr * hr * hz - 0.0177777777777784 * hr * hz * rk) * lambda;
	ms[65] = ms[25] = (-0.00222222222222276 * hr * hr * hz - 0.00222222222222201 * hr * hz * rk) * lambda;
	ms[66] = ms[34] = (6.51330841113425e-16 * hr * hr * hz + 0.00444444444444313 * hr * hz * rk) * lambda;
	ms[67] = ms[43] = (0.0177777777777764 * hr * hr * hz + 0.0355555555555573 * hr * hz * rk) * lambda;
	ms[68] = ms[52] = (0.00444444444444336 * hr * hr * hz + 0.00444444444444378 * hr * hz * rk) * lambda;
	ms[69] = ms[61] = (-2.96059473233375e-17 * hr * hr * hz + 0.00888888888888892 * hr * hz * rk) * lambda;
	ms[70] = (0.0355555555555555 * hr * hr * hz + 0.0711111111111111 * hr * hz * rk) * lambda;
	ms[72] = ms[8] = (0.000555555555555741 * hr * hr * hz + 0.00111111111111089 * hr * hz * rk) * lambda;
	ms[73] = ms[17] = (-0.00222222222222276 * hr * hr * hz - 0.00222222222222201 * hr * hz * rk) * lambda;
	ms[74] = ms[26] = (-0.00388888888888886 * hr * hr * hz - 0.00444444444444446 * hr * hz * rk) * lambda;
	ms[75] = ms[35] = (-0.00111111111111167 * hr * hr * hz - 0.00222222222222156 * hr * hz * rk) * lambda;
	ms[76] = ms[44] = (0.00444444444444336 * hr * hr * hz + 0.00444444444444378 * hr * hz * rk) * lambda;
	ms[77] = ms[53] = (0.00777777777777766 * hr * hr * hz + 0.00888888888888933 * hr * hz * rk) * lambda;
	ms[78] = ms[62] = (-0.00222222222222223 * hr * hr * hz - 0.00444444444444446 * hr * hz * rk) * lambda;
	ms[79] = ms[71] = (0.00888888888888895 * hr * hr * hz + 0.00888888888888888 * hr * hz * rk) * lambda;
	ms[80] = (0.0155555555555555 * hr * hr * hz + 0.0177777777777778 * hr * hz * rk) * lambda;
}
void fill_localmatrixG(double ms[81], double hr, double hz, double rk, double lambda)
{
	ms[0] = (0.5 * (-0.399999999999999 * hr * hr - 3.2 * hr * rk) / hz + 1.0 * (0.149999999999999 * hr * hr + 1.2 * hr * rk) / hz + 0.333333333333333 * (0.266666666666659 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.25 * (-6.0 * hr * hz - 28.0 * hz * rk) / hr + 0.5 * (-3.0 * hr * hz - 14.0 * hz * rk) / hr + 1.0 * (0.5 * hr * hz + 2.33333333333333 * hz * rk) / hr + 0.2 * (2.0 * hr * hz + 9.33333333333333 * hz * rk) / hr + 0.333333333333333 * (6.5 * hr * hz + 30.3333333333333 * hz * rk) / hr) * lambda;
	ms[9] = ms[1] = (0.155555555555563 * hr * rk / hz + 0.333333333333333 * (-8.66666666666669 * hr * hz - 34.6666666666667 * hz * rk) / hr + 0.2 * (-2.66666666666667 * hr * hz - 10.6666666666667 * hz * rk) / hr + 1.0 * (-0.666666666666668 * hr * hz - 2.66666666666667 * hz * rk) / hr + 0.5 * (4.0 * hr * hz + 16.0 * hz * rk) / hr + 0.25 * (8.0 * hr * hz + 32.0 * hz * rk) / hr) * lambda;
	ms[10] = (0.5 * (-6.39999999999998 * hr * hr - 12.8 * hr * rk) / hz + 1.0 * (2.4 * hr * hr + 4.8 * hr * rk) / hz + 0.333333333333333 * (4.26666666666666 * hr * hr + 8.53333333333333 * hr * rk) / hz + 0.25 * (-32.0 * hr * hz - 64.0 * hz * rk) / hr + 0.5 * (-16.0 * hr * hz - 32.0 * hz * rk) / hr + 1.0 * (2.66666666666667 * hr * hz + 5.33333333333333 * hz * rk) / hr + 0.2 * (10.6666666666667 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (34.6666666666667 * hr * hz + 69.3333333333333 * hz * rk) / hr) * lambda;
	ms[18] = ms[2] = (0.333333333333333 * (-0.266666666666668 * hr * hr - 0.533333333333335 * hr * rk) / hz + 1.0 * (-0.15 * hr * hr - 0.300000000000001 * hr * rk) / hz + 0.5 * (0.400000000000006 * hr * hr + 0.799999999999997 * hr * rk) / hz + 0.25 * (-2.0 * hr * hz - 4.0 * hz * rk) / hr + 0.5 * (-1.0 * hr * hz - 2.0 * hz * rk) / hr + 1.0 * (0.166666666666667 * hr * hz + 0.333333333333333 * hz * rk) / hr + 0.2 * (0.666666666666668 * hr * hz + 1.33333333333333 * hz * rk) / hr + 0.333333333333333 * (2.16666666666667 * hr * hz + 4.33333333333333 * hz * rk) / hr) * lambda;
	ms[19] = ms[11] = (0.5 * (-1.6 * hr * hr - 1.59999999999999 * hr * rk) / hz + 1.0 * (0.600000000000001 * hr * hr + 0.6 * hr * rk) / hz + 0.333333333333333 * (1.06666666666667 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.333333333333333 * (-26.0 * hr * hz - 34.6666666666667 * hz * rk) / hr + 0.2 * (-8.0 * hr * hz - 10.6666666666667 * hz * rk) / hr + 1.0 * (-2.0 * hr * hz - 2.66666666666667 * hz * rk) / hr + 0.5 * (12.0 * hr * hz + 16.0 * hz * rk) / hr + 0.25 * (24.0 * hr * hz + 32.0 * hz * rk) / hr) * lambda;
	ms[20] = (0.5 * (-2.8 * hr * hr - 3.2 * hr * rk) / hz + 1.0 * (1.05 * hr * hr + 1.2 * hr * rk) / hz + 0.333333333333333 * (1.86666666666667 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.25 * (-22.0 * hr * hz - 28.0 * hz * rk) / hr + 0.5 * (-11.0 * hr * hz - 14.0 * hz * rk) / hr + 1.0 * (1.83333333333333 * hr * hz + 2.33333333333333 * hz * rk) / hr + 0.2 * (7.33333333333333 * hr * hz + 9.33333333333333 * hz * rk) / hr + 0.333333333333333 * (23.8333333333333 * hr * hz + 30.3333333333333 * hz * rk) / hr) * lambda;
	ms[27] = ms[3] = (0.333333333333333 * (-0.533333333333317 * hr * hr - 4.26666666666665 * hr * rk) / hz + 1.0 * (-0.199999999999999 * hr * hr - 1.6 * hr * rk) / hz + 0.5 * (0.666666666666657 * hr * hr + 5.33333333333331 * hr * rk) / hz + 0.333333333333333 * (-8.0 * hr * hz - 37.3333333333333 * hz * rk) / hr + 0.2 * (-4.0 * hr * hz - 18.6666666666667 * hz * rk) / hr + 0.5 * (2.0 * hr * hz + 9.33333333333333 * hz * rk) / hr + 0.25 * (10.0 * hr * hz + 46.6666666666667 * hz * rk) / hr) * lambda;
	ms[28] = ms[12] = (-0.177777777777772 * hr * rk / hz + 0.25 * (-13.3333333333334 * hr * hz - 53.3333333333333 * hz * rk) / hr + 0.5 * (-2.66666666666667 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.2 * (5.33333333333334 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (10.6666666666667 * hr * hz + 42.6666666666667 * hz * rk) / hr) * lambda;
	ms[29] = ms[21] = (0.5 * (-0.666666666666668 * hr * hr - 1.33333333333334 * hr * rk) / hz + 1.0 * (0.200000000000003 * hr * hr + 0.399999999999999 * hr * rk) / hz + 0.333333333333333 * (0.533333333333337 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.333333333333333 * (-2.66666666666667 * hr * hz - 5.33333333333333 * hz * rk) / hr + 0.2 * (-1.33333333333334 * hr * hz - 2.66666666666666 * hz * rk) / hr + 0.5 * (0.666666666666668 * hr * hz + 1.33333333333333 * hz * rk) / hr + 0.25 * (3.33333333333334 * hr * hz + 6.66666666666666 * hz * rk) / hr) * lambda;
	ms[30] = (0.5 * (-1.06666666666663 * hr * hr - 8.5333333333333 * hr * rk) / hz + 1.0 * (0.266666666666659 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.333333333333333 * (1.06666666666663 * hr * hr + 8.5333333333333 * hr * rk) / hz + 0.25 * (-16.0 * hr * hz - 74.6666666666667 * hz * rk) / hr + 0.533333333333333 * (8.0 * hr * hz + 37.3333333333333 * hz * rk) / hr) * lambda;
	ms[36] = ms[4] = (-0.177777777777772 * hr * rk / hz + 0.25 * (-13.3333333333334 * hr * hz - 53.3333333333333 * hz * rk) / hr + 0.5 * (-2.66666666666667 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.2 * (5.33333333333334 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (10.6666666666667 * hr * hz + 42.6666666666667 * hz * rk) / hr) * lambda;
	ms[37] = ms[13] = (0.333333333333333 * (-8.53333333333332 * hr * hr - 17.0666666666667 * hr * rk) / hz + 1.0 * (-3.19999999999999 * hr * hr - 6.40000000000001 * hr * rk) / hz + 0.5 * (10.6666666666667 * hr * hr + 21.3333333333333 * hr * rk) / hz + 0.333333333333333 * (-42.6666666666667 * hr * hz - 85.3333333333333 * hz * rk) / hr + 0.2 * (-21.3333333333333 * hr * hz - 42.6666666666667 * hz * rk) / hr + 0.5 * (10.6666666666667 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.25 * (53.3333333333334 * hr * hz + 106.666666666667 * hz * rk) / hr) * lambda;
	ms[38] = ms[22] = (0.333333333333333 * (-2.13333333333335 * hr * hr - 2.13333333333333 * hr * rk) / hz + 1.0 * (-0.800000000000001 * hr * hr - 0.799999999999997 * hr * rk) / hz + 0.5 * (2.66666666666667 * hr * hr + 2.66666666666667 * hr * rk) / hz + 0.25 * (-40.0 * hr * hz - 53.3333333333333 * hz * rk) / hr + 0.5 * (-8.0 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.2 * (16.0 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (32.0 * hr * hz + 42.6666666666667 * hz * rk) / hr) * lambda;
	ms[39] = ms[31] = (0.355555555555559 * hr * rk / hz + 0.533333333333333 * (-10.6666666666667 * hr * hz - 42.6666666666667 * hz * rk) / hr + 0.25 * (21.3333333333334 * hr * hz + 85.3333333333333 * hz * rk) / hr) * lambda;
	ms[40] = ((-128 * hr * hr - 256 * hr * rk) / (15 * hz) + (64 * hr * hr + 128 * hr * rk) / (15 * hz) + (256 * hr * hr + 512 * hr * rk) / (45 * hz) + (-64 * hr * hz - 128 * hz * rk) / (3 * hr) + 8 * (128 * hr * hz + 256 * hz * rk) / (45 * hr)) * lambda;
	ms[45] = ms[5] = (0.5 * (-0.666666666666668 * hr * hr - 1.33333333333334 * hr * rk) / hz + 1.0 * (0.200000000000003 * hr * hr + 0.399999999999999 * hr * rk) / hz + 0.333333333333333 * (0.533333333333337 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.333333333333333 * (-2.66666666666667 * hr * hz - 5.33333333333333 * hz * rk) / hr + 0.2 * (-1.33333333333334 * hr * hz - 2.66666666666666 * hz * rk) / hr + 0.5 * (0.666666666666668 * hr * hz + 1.33333333333333 * hz * rk) / hr + 0.25 * (3.33333333333334 * hr * hz + 6.66666666666666 * hz * rk) / hr) * lambda;
	ms[46] = ms[14] = (0.333333333333333 * (-2.13333333333335 * hr * hr - 2.13333333333333 * hr * rk) / hz + 1.0 * (-0.800000000000001 * hr * hr - 0.799999999999997 * hr * rk) / hz + 0.5 * (2.66666666666667 * hr * hr + 2.66666666666667 * hr * rk) / hz + 0.25 * (-40.0 * hr * hz - 53.3333333333333 * hz * rk) / hr + 0.5 * (-8.0 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.2 * (16.0 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (32.0 * hr * hz + 42.6666666666667 * hz * rk) / hr) * lambda;
	ms[47] = ms[23] = (0.333333333333333 * (-3.73333333333333 * hr * hr - 4.26666666666667 * hr * rk) / hz + 1.0 * (-1.4 * hr * hr - 1.6 * hr * rk) / hz + 0.5 * (4.66666666666666 * hr * hr + 5.33333333333333 * hr * rk) / hz + 0.333333333333333 * (-29.3333333333333 * hr * hz - 37.3333333333333 * hz * rk) / hr + 0.2 * (-14.6666666666667 * hr * hz - 18.6666666666667 * hz * rk) / hr + 0.5 * (7.33333333333333 * hr * hz + 9.33333333333333 * hz * rk) / hr + 0.25 * (36.6666666666667 * hr * hz + 46.6666666666667 * hz * rk) / hr) * lambda;
	ms[48] = ms[32] = (0.333333333333333 * (-1.06666666666667 * hr * hr - 2.13333333333334 * hr * rk) / hz + 1.0 * (-0.266666666666668 * hr * hr - 0.533333333333335 * hr * rk) / hz + 0.5 * (1.06666666666667 * hr * hr + 2.13333333333334 * hr * rk) / hz + 0.25 * (-5.33333333333334 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.533333333333333 * (2.66666666666667 * hr * hz + 5.33333333333333 * hz * rk) / hr) * lambda;
	ms[49] = ms[41] = (0.5 * (-4.26666666666669 * hr * hr - 4.26666666666667 * hr * rk) / hz + 1.0 * (1.06666666666667 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.333333333333333 * (4.26666666666669 * hr * hr + 4.26666666666667 * hr * rk) / hz + 0.533333333333333 * (-32.0 * hr * hz - 42.6666666666667 * hz * rk) / hr + 0.25 * (64.0 * hr * hz + 85.3333333333333 * hz * rk) / hr) * lambda;
	ms[50] = (0.5 * (-7.46666666666666 * hr * hr - 8.53333333333333 * hr * rk) / hz + 1.0 * (1.86666666666667 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.333333333333333 * (7.46666666666666 * hr * hr + 8.53333333333333 * hr * rk) / hz + 0.25 * (-58.6666666666667 * hr * hz - 74.6666666666667 * hz * rk) / hr + 0.533333333333333 * (29.3333333333333 * hr * hz + 37.3333333333333 * hz * rk) / hr) * lambda;
	ms[54] = ms[6] = (0.5 * (-0.266666666666659 * hr * hr - 2.13333333333333 * hr * rk) / hz + 1.0 * (0.0499999999999998 * hr * hr + 0.4 * hr * rk) / hz + 0.333333333333333 * (0.266666666666659 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.25 * (-4.0 * hr * hz - 18.6666666666667 * hz * rk) / hr + 0.5 * (-0.5 * hr * hz - 2.33333333333333 * hz * rk) / hr + 0.2 * (2.0 * hr * hz + 9.33333333333333 * hz * rk) / hr + 0.333333333333333 * (2.5 * hr * hz + 11.6666666666667 * hz * rk) / hr) * lambda;
	ms[55] = ms[15] = (0.0222222222222197 * hr * rk / hz + 0.333333333333333 * (-3.33333333333334 * hr * hz - 13.3333333333333 * hz * rk) / hr + 0.2 * (-2.66666666666667 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.5 * (0.666666666666668 * hr * hz + 2.66666666666667 * hz * rk) / hr + 0.25 * (5.33333333333334 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[56] = ms[24] = (0.333333333333333 * (-0.266666666666668 * hr * hr - 0.533333333333335 * hr * rk) / hz + 1.0 * (-0.0500000000000007 * hr * hr - 0.0999999999999996 * hr * rk) / hz + 0.5 * (0.266666666666668 * hr * hr + 0.533333333333335 * hr * rk) / hz + 0.25 * (-1.33333333333334 * hr * hz - 2.66666666666666 * hz * rk) / hr + 0.5 * (-0.166666666666667 * hr * hz - 0.333333333333333 * hz * rk) / hr + 0.2 * (0.666666666666668 * hr * hz + 1.33333333333333 * hz * rk) / hr + 0.333333333333333 * (0.833333333333336 * hr * hz + 1.66666666666666 * hz * rk) / hr) * lambda;
	ms[57] = ms[33] = (0.333333333333333 * (-0.533333333333317 * hr * hr - 4.26666666666665 * hr * rk) / hz + 1.0 * (-0.0666666666666647 * hr * hr - 0.533333333333331 * hr * rk) / hz + 0.5 * (0.399999999999999 * hr * hr + 3.2 * hr * rk) / hz + 0.2 * (-4.0 * hr * hz - 18.6666666666667 * hz * rk) / hr + 0.333333333333333 * (-2.0 * hr * hz - 9.33333333333333 * hz * rk) / hr + 0.25 * (6.0 * hr * hz + 28.0 * hz * rk) / hr) * lambda;
	ms[58] = ms[42] = (-0.17777777777779 * hr * rk / hz + 0.25 * (-8.0 * hr * hz - 32.0 * hz * rk) / hr + 0.333333333333333 * (2.66666666666667 * hr * hz + 10.6666666666667 * hz * rk) / hr + 0.2 * (5.33333333333334 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[59] = ms[51] = (0.5 * (-0.400000000000006 * hr * hr - 0.799999999999997 * hr * rk) / hz + 1.0 * (0.0666666666666671 * hr * hr + 0.133333333333334 * hr * rk) / hz + 0.333333333333333 * (0.533333333333337 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.2 * (-1.33333333333334 * hr * hz - 2.66666666666666 * hz * rk) / hr + 0.333333333333333 * (-0.666666666666668 * hr * hz - 1.33333333333333 * hz * rk) / hr + 0.25 * (2.0 * hr * hz + 4.0 * hz * rk) / hr) * lambda;
	ms[60] = (0.5 * (-0.133333333333329 * hr * hr - 1.06666666666666 * hr * rk) / hz + 1.0 * (0.0166666666666662 * hr * hr + 0.133333333333333 * hr * rk) / hz + 0.333333333333333 * (0.266666666666659 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.25 * (-2.0 * hr * hz - 9.33333333333333 * hz * rk) / hr + 0.333333333333333 * (0.5 * hr * hz + 2.33333333333333 * hz * rk) / hr + 0.2 * (2.0 * hr * hz + 9.33333333333333 * hz * rk) / hr) * lambda;
	ms[63] = ms[7] = (0.0222222222222197 * hr * rk / hz + 0.333333333333333 * (-3.33333333333334 * hr * hz - 13.3333333333333 * hz * rk) / hr + 0.2 * (-2.66666666666667 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.5 * (0.666666666666668 * hr * hz + 2.66666666666667 * hz * rk) / hr + 0.25 * (5.33333333333334 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[64] = ms[16] = (0.5 * (-4.26666666666666 * hr * hr - 8.53333333333333 * hr * rk) / hz + 1.0 * (0.799999999999997 * hr * hr + 1.6 * hr * rk) / hz + 0.333333333333333 * (4.26666666666666 * hr * hr + 8.53333333333333 * hr * rk) / hz + 0.25 * (-21.3333333333333 * hr * hz - 42.6666666666667 * hz * rk) / hr + 0.5 * (-2.66666666666667 * hr * hz - 5.33333333333333 * hz * rk) / hr + 0.2 * (10.6666666666667 * hr * hz + 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (13.3333333333333 * hr * hz + 26.6666666666667 * hz * rk) / hr) * lambda;
	ms[65] = ms[25] = (0.5 * (-1.06666666666667 * hr * hr - 1.06666666666667 * hr * rk) / hz + 1.0 * (0.2 * hr * hr + 0.199999999999999 * hr * rk) / hz + 0.333333333333333 * (1.06666666666667 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.333333333333333 * (-10.0 * hr * hz - 13.3333333333333 * hz * rk) / hr + 0.2 * (-8.0 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.5 * (2.0 * hr * hz + 2.66666666666667 * hz * rk) / hr + 0.25 * (16.0 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[66] = ms[34] = (-0.17777777777779 * hr * rk / hz + 0.25 * (-8.0 * hr * hz - 32.0 * hz * rk) / hr + 0.333333333333333 * (2.66666666666667 * hr * hz + 10.6666666666667 * hz * rk) / hr + 0.2 * (5.33333333333334 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[67] = ms[43] = (0.333333333333333 * (-8.53333333333332 * hr * hr - 17.0666666666667 * hr * rk) / hz + 1.0 * (-1.06666666666666 * hr * hr - 2.13333333333333 * hr * rk) / hz + 0.5 * (6.39999999999998 * hr * hr + 12.8 * hr * rk) / hz + 0.2 * (-21.3333333333333 * hr * hz - 42.6666666666667 * hz * rk) / hr + 0.333333333333333 * (-10.6666666666667 * hr * hz - 21.3333333333333 * hz * rk) / hr + 0.25 * (32.0 * hr * hz + 64.0 * hz * rk) / hr) * lambda;
	ms[68] = ms[52] = (0.333333333333333 * (-2.13333333333335 * hr * hr - 2.13333333333333 * hr * rk) / hz + 1.0 * (-0.266666666666668 * hr * hr - 0.266666666666667 * hr * rk) / hz + 0.5 * (1.6 * hr * hr + 1.59999999999999 * hr * rk) / hz + 0.25 * (-24.0 * hr * hz - 32.0 * hz * rk) / hr + 0.333333333333333 * (8.0 * hr * hz + 10.6666666666667 * hz * rk) / hr + 0.2 * (16.0 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[69] = ms[61] = (0.155555555555557 * hr * rk / hz + 0.2 * (-2.66666666666667 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.333333333333333 * (-0.666666666666668 * hr * hz - 2.66666666666667 * hz * rk) / hr + 0.25 * (2.66666666666667 * hr * hz + 10.6666666666667 * hz * rk) / hr) * lambda;
	ms[70] = (0.5 * (-2.13333333333333 * hr * hr - 4.26666666666667 * hr * rk) / hz + 1.0 * (0.266666666666666 * hr * hr + 0.533333333333333 * hr * rk) / hz + 0.333333333333333 * (4.26666666666666 * hr * hr + 8.53333333333333 * hr * rk) / hz + 0.25 * (-10.6666666666667 * hr * hz - 21.3333333333333 * hz * rk) / hr + 0.333333333333333 * (2.66666666666667 * hr * hz + 5.33333333333333 * hz * rk) / hr + 0.2 * (10.6666666666667 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[72] = ms[8] = (0.333333333333333 * (-0.266666666666668 * hr * hr - 0.533333333333335 * hr * rk) / hz + 1.0 * (-0.0500000000000007 * hr * hr - 0.0999999999999996 * hr * rk) / hz + 0.5 * (0.266666666666668 * hr * hr + 0.533333333333335 * hr * rk) / hz + 0.25 * (-1.33333333333334 * hr * hz - 2.66666666666666 * hz * rk) / hr + 0.5 * (-0.166666666666667 * hr * hz - 0.333333333333333 * hz * rk) / hr + 0.2 * (0.666666666666668 * hr * hz + 1.33333333333333 * hz * rk) / hr + 0.333333333333333 * (0.833333333333336 * hr * hz + 1.66666666666666 * hz * rk) / hr) * lambda;
	ms[73] = ms[17] = (0.5 * (-1.06666666666667 * hr * hr - 1.06666666666667 * hr * rk) / hz + 1.0 * (0.2 * hr * hr + 0.199999999999999 * hr * rk) / hz + 0.333333333333333 * (1.06666666666667 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.333333333333333 * (-10.0 * hr * hz - 13.3333333333333 * hz * rk) / hr + 0.2 * (-8.0 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.5 * (2.0 * hr * hz + 2.66666666666667 * hz * rk) / hr + 0.25 * (16.0 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[74] = ms[26] = (0.5 * (-1.86666666666667 * hr * hr - 2.13333333333333 * hr * rk) / hz + 1.0 * (0.35 * hr * hr + 0.4 * hr * rk) / hz + 0.333333333333333 * (1.86666666666667 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.25 * (-14.6666666666667 * hr * hz - 18.6666666666667 * hz * rk) / hr + 0.5 * (-1.83333333333333 * hr * hz - 2.33333333333333 * hz * rk) / hr + 0.2 * (7.33333333333333 * hr * hz + 9.33333333333333 * hz * rk) / hr + 0.333333333333333 * (9.16666666666667 * hr * hz + 11.6666666666667 * hz * rk) / hr) * lambda;
	ms[75] = ms[35] = (0.5 * (-0.400000000000006 * hr * hr - 0.799999999999997 * hr * rk) / hz + 1.0 * (0.0666666666666671 * hr * hr + 0.133333333333334 * hr * rk) / hz + 0.333333333333333 * (0.533333333333337 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.2 * (-1.33333333333334 * hr * hz - 2.66666666666666 * hz * rk) / hr + 0.333333333333333 * (-0.666666666666668 * hr * hz - 1.33333333333333 * hz * rk) / hr + 0.25 * (2.0 * hr * hz + 4.0 * hz * rk) / hr) * lambda;
	ms[76] = ms[44] = (0.333333333333333 * (-2.13333333333335 * hr * hr - 2.13333333333333 * hr * rk) / hz + 1.0 * (-0.266666666666668 * hr * hr - 0.266666666666667 * hr * rk) / hz + 0.5 * (1.6 * hr * hr + 1.59999999999999 * hr * rk) / hz + 0.25 * (-24.0 * hr * hz - 32.0 * hz * rk) / hr + 0.333333333333333 * (8.0 * hr * hz + 10.6666666666667 * hz * rk) / hr + 0.2 * (16.0 * hr * hz + 21.3333333333333 * hz * rk) / hr) * lambda;
	ms[77] = ms[53] = (0.333333333333333 * (-3.73333333333333 * hr * hr - 4.26666666666667 * hr * rk) / hz + 1.0 * (-0.466666666666666 * hr * hr - 0.533333333333333 * hr * rk) / hz + 0.5 * (2.8 * hr * hr + 3.2 * hr * rk) / hz + 0.2 * (-14.6666666666667 * hr * hz - 18.6666666666667 * hz * rk) / hr + 0.333333333333333 * (-7.33333333333333 * hr * hz - 9.33333333333333 * hz * rk) / hr + 0.25 * (22.0 * hr * hz + 28.0 * hz * rk) / hr) * lambda;
	ms[78] = ms[62] = (0.333333333333333 * (-0.266666666666668 * hr * hr - 0.533333333333335 * hr * rk) / hz + 1.0 * (-0.0166666666666668 * hr * hr - 0.0333333333333334 * hr * rk) / hz + 0.5 * (0.133333333333334 * hr * hr + 0.266666666666667 * hr * rk) / hz + 0.25 * (-0.666666666666668 * hr * hz - 1.33333333333333 * hz * rk) / hr + 0.333333333333333 * (0.166666666666667 * hr * hz + 0.333333333333333 * hz * rk) / hr + 0.2 * (0.666666666666668 * hr * hz + 1.33333333333333 * hz * rk) / hr) * lambda;
	ms[79] = ms[71] = (0.5 * (-0.533333333333337 * hr * hr - 0.533333333333333 * hr * rk) / hz + 1.0 * (0.0666666666666671 * hr * hr + 0.0666666666666667 * hr * rk) / hz + 0.333333333333333 * (1.06666666666667 * hr * hr + 1.06666666666667 * hr * rk) / hz + 0.2 * (-8.0 * hr * hz - 10.6666666666667 * hz * rk) / hr + 0.333333333333333 * (-2.0 * hr * hz - 2.66666666666667 * hz * rk) / hr + 0.25 * (8.0 * hr * hz + 10.6666666666667 * hz * rk) / hr) * lambda;
	ms[80] = (0.5 * (-0.933333333333333 * hr * hr - 1.06666666666667 * hr * rk) / hz + 1.0 * (0.116666666666667 * hr * hr + 0.133333333333333 * hr * rk) / hz + 0.333333333333333 * (1.86666666666667 * hr * hr + 2.13333333333333 * hr * rk) / hz + 0.25 * (-7.33333333333333 * hr * hz - 9.33333333333333 * hz * rk) / hr + 0.333333333333333 * (1.83333333333333 * hr * hz + 2.33333333333333 * hz * rk) / hr + 0.2 * (7.33333333333333 * hr * hz + 9.33333333333333 * hz * rk) / hr) * lambda;
}

void fill_localmatrixM1D(double ms[9], double hr, double rp)
{
	ms[0] = (0.133333333333333 * hr) * rp;
	ms[3] = ms[1] = (0.0666666666666673 * hr)* rp;
	ms[4] = (8 * hr / 15)* rp;
	ms[6] = ms[2] = (-0.0333333333333332 * hr)* rp;
	ms[7] = ms[5] = (0.0666666666666667 * hr)* rp;
	ms[8] = (0.133333333333333 * hr) * rp;
}



//void fill_localmatrixM(double ms[81], double hr, double hz, double rk, double lambda)
//{
//	ms[0] = (0.0177777777777751 * hr * hz) * lambda;
//	ms[9] = ms[1] = (0.008888888888892 * hr * hz) * lambda;
//	ms[10] = (0.0711111111111071 * hr * hz) * lambda;
//	ms[18] = ms[2] = (-0.004444444444446 * hr * hz) * lambda;
//	ms[19] = ms[11] = (0.00888888888889044 * hr * hz) * lambda;
//	ms[20] = (0.0177777777777768 * hr * hz) * lambda;
//	ms[27] = ms[3] = (0.00888888888888784 * hr * hz) * lambda;
//	ms[28] = ms[12] = (0.00444444444444582 * hr * hz) * lambda;
//	ms[29] = ms[21] = (-0.0022222222222229 * hr * hz) * lambda;
//	ms[30] = (0.0711111111111109 * hr * hz) * lambda;
//	ms[36] = ms[4] = (0.00444444444444582 * hr * hz) * lambda;
//	ms[37] = ms[13] = (0.0355555555555539 * hr * hz) * lambda;
//	ms[38] = ms[22] = (0.00444444444444513 * hr * hz) * lambda;
//	ms[39] = ms[31] = (0.0355555555555557 * hr * hz) * lambda;
//	ms[40] = (64 * hr * hz / 225) * lambda;
//	ms[45] = ms[5] = (-0.0022222222222229 * hr * hz) * lambda;
//	ms[46] = ms[14] = (0.00444444444444513 * hr * hz) * lambda;
//	ms[47] = ms[23] = (0.00888888888888892 * hr * hz) * lambda;
//	ms[48] = ms[32] = (-0.0177777777777778 * hr * hz) * lambda;
//	ms[49] = ms[41] = (0.0355555555555556 * hr * hz) * lambda;
//	ms[50] = (0.0711111111111111 * hr * hz) * lambda;
//	ms[54] = ms[6] = (-0.00444444444444481 * hr * hz) * lambda;
//	ms[55] = ms[15] = (-0.00222222222222179 * hr * hz) * lambda;
//	ms[56] = ms[24] = (0.00111111111111089 * hr * hz) * lambda;
//	ms[57] = ms[33] = (0.00888888888889 * hr * hz) * lambda;
//	ms[58] = ms[42] = (0.00444444444444313 * hr * hz) * lambda;
//	ms[59] = ms[51] = (-0.00222222222222157 * hr * hz) * lambda;
//	ms[60] = (0.0177777777777777 * hr * hz) * lambda;
//	ms[63] = ms[7] = (-0.00222222222222179 * hr * hz) * lambda;
//	ms[64] = ms[16] = (-0.0177777777777783 * hr * hz) * lambda;
//	ms[65] = ms[25] = (-0.00222222222222201 * hr * hz) * lambda;
//	ms[66] = ms[34] = (0.00444444444444313 * hr * hz) * lambda;
//	ms[67] = ms[43] = (0.0355555555555573 * hr * hz) * lambda;
//	ms[68] = ms[52] = (0.00444444444444378 * hr * hz) * lambda;
//	ms[69] = ms[61] = (0.00888888888888892 * hr * hz) * lambda;
//	ms[70] = (0.0711111111111111 * hr * hz) * lambda;
//	ms[72] = ms[8] = (0.00111111111111089 * hr * hz) * lambda;
//	ms[73] = ms[17] = (-0.00222222222222201 * hr * hz) * lambda;
//	ms[74] = ms[26] = (-0.00444444444444447 * hr * hz) * lambda;
//	ms[75] = ms[35] = (-0.00222222222222157 * hr * hz) * lambda;
//	ms[76] = ms[44] = (0.00444444444444378 * hr * hz) * lambda;
//	ms[77] = ms[53] = (0.00888888888888933 * hr * hz) * lambda;
//	ms[78] = ms[62] = (-0.00444444444444446 * hr * hz) * lambda;
//	ms[79] = ms[71] = (0.00888888888888889 * hr * hz) * lambda;
//	ms[80] = (0.0177777777777778 * hr * hz) * lambda;
//}
//void fill_localmatrixG(double ms[81], double hr, double hz, double rk, double lambda)
//{
//	ms[0] = (0.31111111111111 * hr / hz + 0.31111111111111 * hz / hr) * lambda;
//	ms[9] = ms[1] = (0.155555555555557 * hr / hz - 0.355555555555554 * hz / hr) * lambda;
//	ms[10] = (1.24444444444444 * hr / hz + 0.711111111111109 * hz / hr) * lambda;
//	ms[18] = ms[2] = (-0.0777777777777775 * hr / hz + 0.0444444444444441 * hz / hr) * lambda;
//	ms[19] = ms[11] = (0.155555555555555 * hr / hz - 0.355555555555554 * hz / hr) * lambda;
//	ms[20] = (0.311111111111111 * hr / hz + 0.31111111111111 * hz / hr) * lambda;
//	ms[27] = ms[3] = (-0.355555555555554 * hr / hz + 0.155555555555557 * hz / hr) * lambda;
//	ms[28] = ms[12] = (-0.177777777777779 * hr / hz - 0.177777777777779 * hz / hr) * lambda;
//	ms[29] = ms[21] = (0.0888888888888886 * hr / hz + 0.0222222222222224 * hz / hr) * lambda;
//	ms[30] = (0.711111111111109 * hr / hz + 1.24444444444444 * hz / hr) * lambda;
//	ms[36] = ms[4] = (-0.177777777777779 * hr / hz - 0.177777777777779 * hz / hr) * lambda;
//	ms[37] = ms[13] = (-1.42222222222222 * hr / hz + 0.355555555555559 * hz / hr) * lambda;
//	ms[38] = ms[22] = (-0.177777777777778 * hr / hz - 0.177777777777779 * hz / hr) * lambda;
//	ms[39] = ms[31] = (0.355555555555559 * hr / hz - 1.42222222222222 * hz / hr) * lambda;
//	ms[40] = (128 * hr / (45 * hz) + 128 * hz / (45 * hr)) * lambda;
//	ms[45] = ms[5] = (0.0888888888888886 * hr / hz + 0.0222222222222224 * hz / hr) * lambda;
//	ms[46] = ms[14] = (-0.177777777777778 * hr / hz - 0.177777777777779 * hz / hr) * lambda;
//	ms[47] = ms[23] = (-0.355555555555555 * hr / hz + 0.155555555555557 * hz / hr) * lambda;
//	ms[48] = ms[32] = (-0.177777777777777 * hr / hz + 0.177777777777777 * hz / hr) * lambda;
//	ms[49] = ms[41] = (0.355555555555555 * hr / hz - 1.42222222222222 * hz / hr) * lambda;
//	ms[50] = (0.711111111111111 * hr / hz + 1.24444444444444 * hz / hr) * lambda;
//	ms[54] = ms[6] = (0.0444444444444441 * hr / hz - 0.0777777777777775 * hz / hr) * lambda;
//	ms[55] = ms[15] = (0.0222222222222224 * hr / hz + 0.0888888888888886 * hz / hr) * lambda;
//	ms[56] = ms[24] = (-0.011111111111111 * hr / hz - 0.011111111111111 * hz / hr) * lambda;
//	ms[57] = ms[33] = (-0.355555555555554 * hr / hz + 0.155555555555555 * hz / hr) * lambda;
//	ms[58] = ms[42] = (-0.177777777777779 * hr / hz - 0.177777777777778 * hz / hr) * lambda;
//	ms[59] = ms[51] = (0.0888888888888886 * hr / hz + 0.0222222222222221 * hz / hr) * lambda;
//	ms[60] = (0.31111111111111 * hr / hz + 0.311111111111111 * hz / hr) * lambda;
//	ms[63] = ms[7] = (0.0222222222222224 * hr / hz + 0.0888888888888886 * hz / hr) * lambda;
//	ms[64] = ms[16] = (0.177777777777777 * hr / hz - 0.177777777777777 * hz / hr) * lambda;
//	ms[65] = ms[25] = (0.0222222222222221 * hr / hz + 0.0888888888888886 * hz / hr) * lambda;
//	ms[66] = ms[34] = (-0.177777777777779 * hr / hz - 0.177777777777778 * hz / hr) * lambda;
//	ms[67] = ms[43] = (-1.42222222222222 * hr / hz + 0.355555555555555 * hz / hr) * lambda;
//	ms[68] = ms[52] = (-0.177777777777778 * hr / hz - 0.177777777777778 * hz / hr) * lambda;
//	ms[69] = ms[61] = (0.155555555555557 * hr / hz - 0.355555555555555 * hz / hr) * lambda;
//	ms[70] = (1.24444444444444 * hr / hz + 0.711111111111111 * hz / hr) * lambda;
//	ms[72] = ms[8] = (-0.011111111111111 * hr / hz - 0.011111111111111 * hz / hr) * lambda;
//	ms[73] = ms[17] = (0.0222222222222221 * hr / hz + 0.0888888888888886 * hz / hr) * lambda;
//	ms[74] = ms[26] = (0.0444444444444443 * hr / hz - 0.0777777777777775 * hz / hr) * lambda;
//	ms[75] = ms[35] = (0.0888888888888886 * hr / hz + 0.0222222222222221 * hz / hr) * lambda;
//	ms[76] = ms[44] = (-0.177777777777778 * hr / hz - 0.177777777777778 * hz / hr) * lambda;
//	ms[77] = ms[53] = (-0.355555555555555 * hr / hz + 0.155555555555556 * hz / hr) * lambda;
//	ms[78] = ms[62] = (-0.0777777777777775 * hr / hz + 0.0444444444444443 * hz / hr) * lambda;
//	ms[79] = ms[71] = (0.155555555555556 * hr / hz - 0.355555555555555 * hz / hr) * lambda;
//	ms[80] = (0.311111111111111 * hr / hz + 0.311111111111111 * hz / hr) * lambda;
//}
//
//void fill_localmatrixM1D(double ms[9], double hr, double rp)
//{
//	ms[0] = (0.133333333333333 * hr);
//	ms[3] = ms[1] = (0.0666666666666673 * hr) ;
//	ms[4] = (8 * hr / 15) ;
//	ms[6] = ms[2] = (-0.0333333333333332 * hr) ;
//	ms[7] = ms[5] = (0.0666666666666667 * hr) ;
//	ms[8] = (0.133333333333333 * hr) ;
//}

void Grid2DQuad::buildMatrix()
{

	bf.clear();
	for (int j = 0; j < nodes.size(); ++j)
	{
		bf.push_back(0);
	}


	for (int i = 0; i < elems.size(); ++i)
	{

		double A[81], M[81], b[9] = {0}, llh[2] = { 0 }, lambda = 0;//llh-������� ��������� �� ��������(x,y,z)
		Element* th = elems[i];//�������� �������� �������
		int n = th->nodes.size();
		llh[0] = abs(nodes[th->nodes[0]]->coords[0] - nodes[th->nodes[2]]->coords[0]);
		llh[1] = abs(nodes[th->nodes[0]]->coords[1] - nodes[th->nodes[6]]->coords[1]);
		for (int j = 0; j < n; ++j)
		{
			lambda += nodes[th->nodes[j]]->params[0];
		}
		lambda /= n;


		fill_localmatrixG(A, llh[0], llh[1], nodes[th->nodes[0]]->coords[0], lambda);
		fill_localmatrixM(M, llh[0], llh[1], nodes[th->nodes[0]]->coords[0], 1);

		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++)
				b[j] += M[j * n + k] * F[th->nodes[k]];
		}


		for (int j = 0; j < n; j++)
		{

			for (int k = 0; k < n; k++)
			{
				printf("%lf, ", A[j * n + k]);
			}
			printf(";\n");
		}

		//printf("\n");
		/*	for (int k = 0; k < 4; k++)
			{
				printf("%lf, ", F(nodes[th->nodes[k]]->r, nodes[th->nodes[k]]->z));
			}*/
	//	printf("\n");
		for (int j = 0; j < n; ++j)
		{
			bf[th->nodes[j]] += b[j];
			for (int jj = 0; jj < n; ++jj)
			{
				if (elems[i]->nodes[jj] >= th->nodes[j])
					continue;
				int indx = j * n + jj;
				int s = ig[elems[i]->nodes[j]];
				int e = ig[elems[i]->nodes[j] + 1];
				for (; s < e && jg[s] != elems[i]->nodes[jj]; ++s)
					;
				if (s != e)
				{
					al[s] += A[indx];
					au[s] += A[indx];
				}
				//AG[elems[i]->nodes[j] * nodes.size() + elems[i]->nodes[jj]] += A[indx];
			}
			diag[elems[i]->nodes[j]] += A[j * n + j];
		}
	}
	//for (int i = 0; i < nodes.size(); ++i)
	//{
	//	for (int j = 0; j < nodes.size(); ++j)
	//		printf("%.14lf,", AG[i * nodes.size() + j]);
	//	printf(";");
	//}
	return;
}


void Grid2DQuad::secondBoundary()
{
	double M[9];
	double h, rp;
	for (int i = 0; i < secondB.size(); i++)
	{
		Element* el = secondB[i];

		h = -nodes[el->nodes[0]]->coords[0] + nodes[el->nodes[el->nodes.size() - 1]]->coords[0];
		rp = nodes[el->nodes[0]]->coords[0];
		if (h == 0)
		{
			h = -nodes[el->nodes[0]]->coords[1] + nodes[el->nodes[el->nodes.size() - 1]]->coords[1];
		}
			
		fill_localmatrixM1D(M, h, 0.5);
		for (int j = 0; j < el->nodes.size(); j++)
		{
			int node_index = el->nodes[j];
			bf[node_index] += M[j * 3] * el->parameters[0] + M[j * 3 + 1] * el->parameters[1] + M[j * 3 + 2] * el->parameters[2];
		}
	}
};
