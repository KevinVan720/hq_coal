#include "math.h"
#include "FourVector.h"
#include <map>
#include <random>
#include <chrono>
#include <utility>

using namespace Jetscape;

#define min(x,y) x > y ? y : x
#define max(x,y) x > y ? x : y

FourVector boost(const FourVector& p, double vx, double vy, double vz)
{
  FourVector prest;
  double px=p.x();
  double py=p.y();
  double pz=p.z();
  double p0=p.t();
  double prx,pry,prz,pr0;
  double beta = sqrt( vx*vx + vy*vy + vz*vz );
  double gamma;
  double cosPhi;
          // Set momentum in fluid cell's frame
          // 1: for brick
  if (beta < 1e-10)
    {
      gamma = 1.;
      cosPhi = 1.;
      prest = p;
    }
          // 2: for evolving medium
  else
    {
      
      gamma  = 1./sqrt( 1. - beta*beta );
      cosPhi = ( px*vx + py*vy + pz*vz )/( p0*beta );

              // boost particle to the local rest frame of fluid cell
       pr0 = p0*gamma*( 1. - beta*cosPhi );

       prx = -vx*gamma*p0
             + (1.+(gamma-1.)*vx*vx/(beta*beta))*px
             + (gamma-1.)*vx*vy/(beta*beta)*py
             + (gamma-1.)*vx*vz/(beta*beta)*pz;
       pry = -vy*gamma*p0
             + (1.+(gamma-1.)*vy*vy/(beta*beta))*py
             + (gamma-1.)*vx*vy/(beta*beta)*px
             + (gamma-1.)*vy*vz/(beta*beta)*pz;
       prz = -vz*gamma*p0
             + (1.+(gamma-1.)*vz*vz/(beta*beta))*pz
             + (gamma-1.)*vx*vz/(beta*beta)*px
             + (gamma-1.)*vy*vz/(beta*beta)*py;

       prest = FourVector ( prx, pry, prz, pr0 );

    }
  return prest;
}

class Coal_Sampler
{
  private:
    int mc_max_iter;
    int mc_burnin;
    int jump_time;
    double pl_max;
    std::map<int, double> mass_dict;
    std::map<int, double> meson_omega_dict;
    std::map<int, double> baryon_omega_dict;
    std::map<int, double> gM_dict;
    std::map<int, double> gB_dict;
    std::map<int, int> meson_dict;
    std::map<int, std::pair<int,int>> baryon_dict;
    std::vector<int> charm_meson_list;
    std::vector<int> bottom_meson_list;
    std::vector<int> charm_baryon_list;
    std::vector<int> bottom_baryon_list;
    unsigned seed;
    int gq;
    int gg;
    double v_x;
    double v_y;
    double v_z;
    double T;
    //std::map<int, double> meson_prob;
    //std::map<int, double> baryon_prob;

  public:
    Coal_Sampler(double fv_x, double fv_y, double fv_z, double Temp): v_x(fv_x), v_y(fv_y), v_z(fv_z), T(Temp)
    {
        //MC parameters
        seed = std::chrono::system_clock::now().time_since_epoch().count();
        mc_max_iter = 10000000;
        mc_burnin = 100000;
        jump_time = 100;
        pl_max = 10. * T;
        //model parameters
        gq = 6;
        gg = 16;
        meson_omega_dict[4] = 330;
        meson_omega_dict[5] = 330;
        baryon_omega_dict[4] = 430;
        baryon_omega_dict[5] = 410;
        charm_meson_list = {411, 421, 413, 423, 431, 433};
        bottom_meson_list = {511, 521, 513, 523, 531, 533};
        charm_baryon_list = {4122, 4222, 4112, 4212, 4224, 4114, 4214, 4232, 4132, 4324, 4314, 4332, 4334};
        bottom_baryon_list = {5122, 5232, 5132, 5332};
        mass_dict[1] = 300;
        mass_dict[2] = 300;
        mass_dict[3] = 475;
        mass_dict[4] = 1270;
        mass_dict[5] = 4200;
        mass_dict[21] = 600;
        mass_dict[411] = 1870;
        mass_dict[421] = 1865;
        mass_dict[413] = 2010;
        mass_dict[423] = 2007;
        mass_dict[431] = 1968;
        mass_dict[433] = 2112;
        mass_dict[511] = 5280;
        mass_dict[521] = 5279;
        mass_dict[513] = 5325;
        mass_dict[523] = 5325;
        mass_dict[531] = 5367;
        mass_dict[533] = 5415;
        mass_dict[4122] = 2286;
        mass_dict[4222] = 2454;
        mass_dict[4112] = 2454;
        mass_dict[4212] = 2453;
        mass_dict[4224] = 2518;
        mass_dict[4114] = 2518;
        mass_dict[4214] = 2517;
        mass_dict[4232] = 2468;
        mass_dict[4132] = 2471;
        mass_dict[4324] = 2646;
        mass_dict[4314] = 2646;
        mass_dict[4332] = 2695;
        mass_dict[4334] = 2766;
        mass_dict[5122] = 5619;
        mass_dict[5232] = 5792;
        mass_dict[5132] = 5794;
        mass_dict[5332] = 6046;
        gM_dict[411] = 1 / 36.;
        gM_dict[421] = 1 / 36.;
        gM_dict[413] = 1 / 12.;
        gM_dict[423] = 1 / 12.;
        gM_dict[431] = 1 / 36.;
        gM_dict[433] = 1 / 12.;
        gM_dict[511] = 1 / 36.;
        gM_dict[521] = 1 / 36.;
        gM_dict[513] = 1 / 12.;
        gM_dict[523] = 1 / 12.;
        gM_dict[531] = 1 / 36.;
        gM_dict[533] = 1 / 12.;
        gB_dict[4122] = 1 / 108.;
        gB_dict[4222] = 1 / 36.;
        gB_dict[4112] = 1 / 36.;
        gB_dict[4212] = 1 / 36.;
        gB_dict[4224] = 1 / 18.;
        gB_dict[4114] = 1 / 18.;
        gB_dict[4214] = 1 / 18.;
        gB_dict[4232] = 1 / 54.;
        gB_dict[4132] = 1 / 54.;
        gB_dict[4324] = 1 / 18.;
        gB_dict[4314] = 1 / 18.;
        gB_dict[4332] = 1 / 108.;
        gB_dict[4334] = 1 / 36.;
        gB_dict[5122] = 1 / 108.;
        gB_dict[5232] = 1 / 54.;
        gB_dict[5132] = 1 / 54.;
        gB_dict[5332] = 1 / 108.;

        meson_dict[411] = 2;
        meson_dict[421] = 1;
        meson_dict[413] = 2;
        meson_dict[423] = 1;
        meson_dict[431] = 3;
        meson_dict[433] = 3;
        meson_dict[511] = 2;
        meson_dict[521] = 1;
        meson_dict[513] = 2;
        meson_dict[523] = 1;
        meson_dict[531] = 3;
        meson_dict[533] = 3;
        baryon_dict[4122] = std::make_pair(1,2);
        baryon_dict[4222] = std::make_pair(1,1);
        baryon_dict[4112] = std::make_pair(2,2);
        baryon_dict[4212] = std::make_pair(1,2);
        baryon_dict[4224] = std::make_pair(1,1);
        baryon_dict[4114] = std::make_pair(2,2);
        baryon_dict[4214] = std::make_pair(1,2);
        baryon_dict[4232] = std::make_pair(1,3);
        baryon_dict[4132] = std::make_pair(2,3);
        baryon_dict[4324] = std::make_pair(1,3);
        baryon_dict[4314] = std::make_pair(2,3);
        baryon_dict[4332] = std::make_pair(3,3);
        baryon_dict[4334] = std::make_pair(3,3);
        baryon_dict[5122] = std::make_pair(1,2);
        baryon_dict[5232] = std::make_pair(1,3);
        baryon_dict[5132] = std::make_pair(2,3);
        baryon_dict[5332] = std::make_pair(3,3);

    }

    double q2(FourVector p1, FourVector p2);
    double q12(FourVector p1, FourVector p2, FourVector p3);
    double q22(FourVector p1, FourVector p2, FourVector p3);
    double Wigner(FourVector p1, FourVector p2, double sigma);
    double Wigner(FourVector p1, FourVector p2, FourVector p3, double sigma1, double sigma2);
    double FD_dist(double E, double T);
    double BE_dist(double E, double T);
    double prob2meson(int light_id, FourVector &p_l, int heavy_id, FourVector &p_h);
    double prob2baryon(int light_id1, FourVector &p_l1, int light_id2, FourVector &p_l2, int heavy_id, FourVector &p_h);
    double alpha_meson(int light_id, FourVector &p_l, FourVector &p_lnew, int heavy_id, FourVector &p_h);
    double alpha_baryon(int light_id1, FourVector &p_l1, FourVector &p_l1new, int light_id2, FourVector &p_l2, FourVector &p_l2new, int heavy_id, FourVector &p_h);
    double mc_integral(int light_id, int heavy_id, FourVector &p_h);
    double mc_integral(int light_id1, int light_id2, int heavy_id, FourVector &p_h);
    FourVector mc_sample(int light_id, int heavy_id, FourVector &p_h);
    FourVector mc_sample(int light_id1, int light_id2, int heavy_id, FourVector &p_h);
    double recomb_prob(int heavy_id, FourVector &p_h);
};
