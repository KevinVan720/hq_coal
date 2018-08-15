#include "sample.h"

int main(int argc, char **argv)
{
    Coal_Sampler sampler;
    sampler.SetHydro(0., 0., 0., 165.);

    for (int i = 0; i < 100; i++)
    {
        double p_x = 200 * i;
        double p_y = 0.;
        double p_z = 0.;
        FourVector p_h({p_x, p_y, p_z, sqrt(1270 * 1270 + p_x * p_x + p_y * p_y + p_z * p_z)});
        cout << "--------------------------" << endl;
        cout <<"heavy quark momentum: " << p_x << " " << p_y << " " << p_z << endl;
        double rst = sampler.recomb_prob(4, p_h);
        cout << "total probability: " << rst << endl;
        //cout << p_x <<" "<< sampler.mc_integral(1, 5, p_h) << endl;

    }

    return 0;
}