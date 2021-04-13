#include <stdio.h>
double U(double r, double z){
    return r*r*r;
}
int main(){
    int n;
    double r, z;
    FILE * fl = fopen("nodes.txt", "r");
    fscanf(fl, "%d", &n);
    for (int i = 0; i < n; ++i)
    {
        fscanf(fl, "%lf%lf", &r, &z);
        printf("%d %lf;\n", i, U(r, z));
    }

}