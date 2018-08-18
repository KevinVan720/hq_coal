#ifndef SAMPLE_H
#define SAMPLE_H

#include <cmath>
#include "FourVector.h"
#include <map>
#include <random>
#include <chrono>
#include <utility>
#include "H5Cpp.h"
#include <algorithm>
#include <fstream>
//#include <boost/filesystem.hpp>

using namespace Jetscape;

#define min(x,y) x > y ? y : x
#define max(x,y) x > y ? x : y

class Coal_Sampler
{
  private:
    int mc_max_iter;
    int mc_burnin;
    int jump_time;
    double T_start;
    size_t T_N;
    size_t p_N;
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
    
    std::string fname;

  public:
    Coal_Sampler()
    {
      cout << "Initialize..." << endl;
      //MC parameters
      fname = "recomb_table.h5";
      seed = std::chrono::system_clock::now().time_since_epoch().count();
      mc_max_iter = 10000000;
      mc_burnin = 50000;
      jump_time = 100;
      T_start = 150;
      T_N =10;
      p_N = 50;
      //model parameters
      gq = 6;
      gg = 16;
      //omega's
      meson_omega_dict[4] = 230;
      meson_omega_dict[5] = 110;
      baryon_omega_dict[4] = 230;
      baryon_omega_dict[5] = 110;
      //meson and baryons considered in this model
      charm_meson_list = {411, 421, 413, 423, 431, 433};
      bottom_meson_list = {511, 521, 513, 523, 531, 533};
      charm_baryon_list = {4122, 4222, 4112, 4212, 4224, 4114, 4214, 4232, 4132, 4324, 4314, 4332, 4334};
      bottom_baryon_list = {5122, 5222, 5112, 5212, 5224, 5114, 5214, 5232, 5132, 5324, 5314, 5332, 5334};
      //mass of various particles
      mass_dict[1] = 300;
      mass_dict[2] = 300;
      mass_dict[3] = 475;
      mass_dict[4] = 1270;
      mass_dict[5] = 4190;
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
      mass_dict[5222] = 5811;
      mass_dict[5112] = 5815;
      mass_dict[5212] = 5811;
      mass_dict[5224] = 5832;
      mass_dict[5114] = 5835;
      mass_dict[5214] = 5832;
      mass_dict[5232] = 5792;
      mass_dict[5132] = 5794;
      mass_dict[5324] = 5949;
      mass_dict[5314] = 5949;
      mass_dict[5332] = 6046;
      mass_dict[5334] = 6046;
      //g factor of mesons and baryons
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
      gB_dict[4222] = 1 / 36. * 1 / 2.;
      gB_dict[4112] = 1 / 36. * 1 / 2.;
      gB_dict[4212] = 1 / 36.;
      gB_dict[4224] = 1 / 18. * 1 / 2.;
      gB_dict[4114] = 1 / 18. * 1 / 2.;
      gB_dict[4214] = 1 / 18.;
      gB_dict[4232] = 1 / 54.;
      gB_dict[4132] = 1 / 54.;
      gB_dict[4324] = 1 / 18.;
      gB_dict[4314] = 1 / 18.;
      gB_dict[4332] = 1 / 108. * 1 / 2.;
      gB_dict[4334] = 1 / 36. * 1 / 2.;

      gB_dict[5122] = 1 / 108.;
      gB_dict[5222] = 1 / 36. * 1 / 2.;
      gB_dict[5112] = 1 / 36. * 1 / 2.;
      gB_dict[5212] = 1 / 36.;
      gB_dict[5224] = 1 / 18. * 1 / 2.;
      gB_dict[5114] = 1 / 18. * 1 / 2.;
      gB_dict[5214] = 1 / 18.;
      gB_dict[5232] = 1 / 54.;
      gB_dict[5132] = 1 / 54.;
      gB_dict[5324] = 1 / 18.;
      gB_dict[5314] = 1 / 18.;
      gB_dict[5332] = 1 / 108. * 1 / 2.;
      gB_dict[5334] = 1 / 36. * 1 / 2.;
      //constituent light quark id of meson and baryons
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

      baryon_dict[4122] = std::make_pair(1, 2);
      baryon_dict[4222] = std::make_pair(1, 1);
      baryon_dict[4112] = std::make_pair(2, 2);
      baryon_dict[4212] = std::make_pair(1, 2);
      baryon_dict[4224] = std::make_pair(1, 1);
      baryon_dict[4114] = std::make_pair(2, 2);
      baryon_dict[4214] = std::make_pair(1, 2);
      baryon_dict[4232] = std::make_pair(1, 3);
      baryon_dict[4132] = std::make_pair(2, 3);
      baryon_dict[4324] = std::make_pair(1, 3);
      baryon_dict[4314] = std::make_pair(2, 3);
      baryon_dict[4332] = std::make_pair(3, 3);
      baryon_dict[4334] = std::make_pair(3, 3);

      baryon_dict[5122] = std::make_pair(1, 2);
      baryon_dict[5222] = std::make_pair(1, 1);
      baryon_dict[5112] = std::make_pair(2, 2);
      baryon_dict[5212] = std::make_pair(1, 2);
      baryon_dict[5224] = std::make_pair(1, 1);
      baryon_dict[5114] = std::make_pair(2, 2);
      baryon_dict[5214] = std::make_pair(1, 2);
      baryon_dict[5232] = std::make_pair(1, 3);
      baryon_dict[5132] = std::make_pair(2, 3);
      baryon_dict[5324] = std::make_pair(1, 3);
      baryon_dict[5314] = std::make_pair(2, 3);
      baryon_dict[5332] = std::make_pair(3, 3);
      baryon_dict[5334] = std::make_pair(3, 3);

    }
    
    void SetHydro(double vx, double vy, double vz, double Thydro)
    {
      v_x = vx;
      v_y = vy;
      v_z = vz;
      T = Thydro;
      pl_max = 10. * T;
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
    void SaveTable(int pid);
    void ReadTable(int pid);
    void GenerateTable();
};

#endif
