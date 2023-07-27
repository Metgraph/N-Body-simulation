#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DIMENSION 3

void get_energy(double pos[], double vel[], double mass[], double G, int sz, double *KE, double *PE) {

    // calculate KE
    double local_KE=0;
    for (int i = 0; i < sz; i++) {
        // row = i-th body vel
        double *row = &vel[i * DIMENSION];
        double local_sum = 0;
        for (int axis = 0; axis < DIMENSION; axis++) {
            local_sum += row[axis] * row[axis];
        }
        //every cell is multiplied by corrispetive mass and summed in local_KE
        local_KE += mass[i] * local_sum;
    }
    *KE=local_KE*0.5;

    //calculate PE
    double local_PE = 0;

    //calculate only the upper triangle in matrix instead of calculate whole matrix and take only the
    //upper triangle values like using np.triu
    for (int i = 0; i < sz; i++) {
        for (int j = i; j < sz; j++) {
            double sum=0;
            for (int axis = 0; axis < DIMENSION; axis++) {
                // d_value is the value of d{axis} matrix in position [i,j]
                double d_value =
                    pos[j * DIMENSION + axis] - pos[i * DIMENSION + axis];
                sum += d_value * d_value;
            }
            //calculate a cell of -(mass*mass.T)*inv_rl
            double val = sqrt(sum);
            val = val > 0 ? 1.0 / val : val;
            val *= -mass[i] * mass[j];
            //every cell is summed in local_PE
            local_PE += val;
        }
    }

    *PE=local_PE*G;
}

#ifdef TEST_GET_ENERGY
int main() {
    double vel[] = {5,65,105, -18,9,-89, -64, 23, 32, 45,45,45, -29, -71, -61};
    double pos[]={90,44,-32,33,100,9,-84,-4,0,103,182,71,9,2,-103};

    double mass[] = {103, 50,88, 34, 59};
    double KE, PE;
    get_energy(pos, vel, mass, 0.7, sizeof(mass) / sizeof(double), &KE, &PE);
    printf("KE: %.15lf, PE: %.15lf\n", KE, PE);
}
#endif