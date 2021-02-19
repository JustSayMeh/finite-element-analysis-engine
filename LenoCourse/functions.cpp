#include <math.h>
double F(double r, double z)
{
   return 0;
}


//


//void fill_local_matrixM1d(double mt[4], double hr, double rp) {
//	mt[0] = hr / 3;
//	mt[1] = mt[2] = hr / 6;
//	mt[3] = hr / 3;
//}
//void fill_localmatrixM(double mt[16], double hr, double hz, double rp) {
//	mt[0] = hr * hz / 9;
//	mt[1] = mt[4] = hr * hz / 18;
//	mt[2] = mt[8] = hr * hz / 18;
//	mt[3] = mt[12] = hr * hz / 36;
//	mt[5] = hr * hz / 9;
//	mt[6] = mt[9] = hr * hz / 36;
//	mt[7] = mt[13] = hr * hz / 18;
//	mt[10] = hr * hz / 9;
//	mt[11] = mt[14] = hr * hz / 18;
//	mt[15] = hr * hz / 9;
//}
//void fill_localmatrixG(double mt[16], double hr, double hz, double rp) {
//	mt[0] = hr / (3 * hz) + hz / (3 * hr);
//	mt[1] = mt[4] = hr / (6 * hz) - hz / (3 * hr);
//	mt[2] = mt[8] = -hr / (3 * hz) + hz / (6 * hr);
//	mt[3] = mt[12] = -hr / (6 * hz) - hz / (6 * hr);
//	mt[5] = hr / (3 * hz) + hz / (3 * hr);
//	mt[6] = mt[9] = -hr / (6 * hz) - hz / (6 * hr);
//	mt[7] = mt[13] = -hr / (3 * hz) + hz / (6 * hr);
//	mt[10] = hr / (3 * hz) + hz / (3 * hr);
//	mt[11] = mt[14] = hr / (6 * hz) - hz / (3 * hr);
//	mt[15] = hr / (3 * hz) + hz / (3 * hr);
//}
