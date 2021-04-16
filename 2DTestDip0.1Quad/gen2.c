#include <stdio.h>
#include <math.h>
double U(double r, double z){
    return -log(r) + log(3);
}
int main(){
    int n = 0, j = 0;
    double rms[1000], zms[1000];
    double r = 1, z = 1;
    double dr = 0.125, dz = 1;
    for (; z <= 3; z += dz)
    {
        j += 1;
        for (r = 1; r <= 3; r += dr)
        {
            rms[n] = r;
            zms[n] = z;
            n += 1;
          
            printf("%lf %lf\n", r, z);
        }
    }
    int rcount = n / j;
    printf("%d %lf", n, n * 1.0 / j);
    printf("\n-------------------------------\n");
    for (int i  = 0 ; i < n; i++)
    {
        if (i % rcount == 0 || i % rcount == rcount - 1 || i < rcount || i > n - rcount)
            printf("%d %.15lf\n", i, U(rms[i], zms[i]));
    }
    printf("\n-------------------------------\n");
    for (int i = 0;  i < n - rcount; i++)
    {
        if (i % rcount == rcount - 1)
            continue;
        printf("%d %d %d %d\n", i, i + 1, i + rcount, i + rcount + 1);
    }
    printf("\n-------------------------------\n");
    for (int i = 0 ; i < n; i++)
    {
        printf("%.15lf;\n", U(rms[i], zms[i]));
    }
}