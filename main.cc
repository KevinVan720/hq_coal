#include "sample.h"

int main(int argc, char **argv)
{
    Coal_Sampler sampler(0, 0, 0, 165);

    for (int i = 0; i < 100; i++)
    {
        double p_x = 200 * i;
        FourVector p_h({p_x, 0, 0, sqrt(1270*1270+p_x*p_x)});
        cout << p_x <<" "<< sampler.recomb_prob(4,p_h) << endl;
        //cout << p_x <<" "<< sampler.mc_integral(1, 3,  5, p_h) << endl;

    }

    return 0;
}