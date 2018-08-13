#include "sample.h"

double Coal_Sampler::q2(FourVector p1, FourVector p2)
{
    double v_cm_x = (p1.x() + p2.x()) / (p1.t() + p2.t());
    double v_cm_y = (p1.y() + p2.y()) / (p1.t() + p2.t());
    double v_cm_z = (p1.z() + p2.z()) / (p1.t() + p2.t());
    p1 = boost(p1, v_cm_x, v_cm_y, v_cm_z);
    p2 = boost(p2, v_cm_x, v_cm_y, v_cm_z);
    double q_x = (p2.t() * p1.x() - p1.t() * p2.x()) / (p1.t() + p2.t());
    double q_y = (p2.t() * p1.y() - p1.t() * p2.y()) / (p1.t() + p2.t());
    double q_z = (p2.t() * p1.z() - p1.t() * p2.z()) / (p1.t() + p2.t());
    return (q_x * q_x + q_y * q_y + q_z * q_z);
}

double Coal_Sampler::q12(FourVector p1, FourVector p2, FourVector p3)
{
    double v_cm_x = (p1.x() + p2.x() + p3.x()) / (p1.t() + p2.t() + p3.t());
    double v_cm_y = (p1.y() + p2.y() + p3.y()) / (p1.t() + p2.t() + p3.t());
    double v_cm_z = (p1.z() + p2.z() + p3.z()) / (p1.t() + p2.t() + p3.t());
    p1 = boost(p1, v_cm_x, v_cm_y, v_cm_z);
    p2 = boost(p2, v_cm_x, v_cm_y, v_cm_z);
    p3 = boost(p3, v_cm_x, v_cm_y, v_cm_z);
    double q1_x = (p2.t() * p1.x() - p1.t() * p2.x()) / (p1.t() + p2.t());
    double q1_y = (p2.t() * p1.y() - p1.t() * p2.y()) / (p1.t() + p2.t());
    double q1_z = (p2.t() * p1.z() - p1.t() * p2.z()) / (p1.t() + p2.t());
    return (q1_x * q1_x + q1_y * q1_y + q1_z * q1_z);
}

double Coal_Sampler::q22(FourVector p1, FourVector p2, FourVector p3)
{
    double v_cm_x = (p1.x() + p2.x() + p3.x()) / (p1.t() + p2.t() + p3.t());
    double v_cm_y = (p1.y() + p2.y() + p3.y()) / (p1.t() + p2.t() + p3.t());
    double v_cm_z = (p1.z() + p2.z() + p3.z()) / (p1.t() + p2.t() + p3.t());
    p1 = boost(p1, v_cm_x, v_cm_y, v_cm_z);
    p2 = boost(p2, v_cm_x, v_cm_y, v_cm_z);
    p3 = boost(p3, v_cm_x, v_cm_y, v_cm_z);
    double q2_x = (p3.t() * (p1.x() + p2.x()) - (p1.t() + p2.t()) * p3.x()) / (p1.t() + p2.t() + p3.t());
    double q2_y = (p3.t() * (p1.y() + p2.y()) - (p1.t() + p2.t()) * p3.y()) / (p1.t() + p2.t() + p3.t());
    double q2_z = (p3.t() * (p1.z() + p2.z()) - (p1.t() + p2.t()) * p3.z()) / (p1.t() + p2.t() + p3.t());
    return (q2_x * q2_x + q2_y * q2_y + q2_z * q2_z);
}

double Coal_Sampler::Wigner(FourVector p1, FourVector p2, double sigma)
{
    return pow(2. * sqrt(pi) * sigma, 3) * exp(-q2(p1, p2) * pow(sigma, 2));
}

double Coal_Sampler::Wigner(FourVector p1, FourVector p2, FourVector p3, double sigma1, double sigma2)
{
    return pow(2. * sqrt(pi) * sigma1, 3) * pow(2. * sqrt(pi) * sigma2, 3) * exp(-q12(p1, p2, p3) * pow(sigma1, 2) - q22(p1, p2, p3) * pow(sigma2, 2));
}

double Coal_Sampler::FD_dist(double E, double T)
{
    return 1. / (exp(E / T) + 1.);
}

double Coal_Sampler::BE_dist(double E, double T)
{
    return 1. / (exp(E / T) - 1.);
}

double Coal_Sampler::prob2meson(int light_id, FourVector &p_l, int heavy_id, FourVector &p_h)
{
    double m_l = mass_dict[light_id];
    double m_h = mass_dict[heavy_id];
    double mu = m_l * m_h / (m_l + m_h);
    double sigma = 1. / sqrt(mu * meson_omega_dict[heavy_id]);
    double p = sqrt(pow(p_l.x(), 2) + pow(p_l.y(), 2) + pow(p_l.z(), 2));
    return (
               gq * FD_dist(p_l.t(), T) +
               8. / 3. * gg * BE_dist(sqrt(pow(2. * p, 2) + pow(mass_dict[21], 2)), T)) *
           Wigner(p_h, p_l, sigma);
    ;
}

double Coal_Sampler::alpha_meson(int light_id, FourVector &p_l, FourVector &p_lnew, int heavy_id, FourVector &p_h)
{
    double m_l = mass_dict[light_id];
    double m_h = mass_dict[heavy_id];
    double mu = m_l * m_h / (m_l + m_h);
    double sigma = 1. / sqrt(mu * meson_omega_dict[heavy_id]);
    double p = sqrt(pow(p_l.x(), 2) + pow(p_l.y(), 2) + pow(p_l.z(), 2));
    double p_new = sqrt(pow(p_lnew.x(), 2) + pow(p_lnew.y(), 2) + pow(p_lnew.z(), 2));
    double ratio = (gq * FD_dist(p_lnew.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p_new, 2) + pow(mass_dict[21], 2)), T)) /
                   (gq * FD_dist(p_l.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p, 2) + pow(mass_dict[21], 2)), T));

    ratio *= exp(-(q2(p_lnew, p_h) - q2(p_l, p_h)) * pow(sigma, 2));
    return ratio;
}

double Coal_Sampler::alpha_baryon(int light_id1, FourVector &p_l1, FourVector &p_l1new, int light_id2, FourVector &p_l2, FourVector &p_l2new, int heavy_id, FourVector &p_h)
{
    double m_l1 = mass_dict[light_id1];
    double m_l2 = mass_dict[light_id2];
    double m_h = mass_dict[heavy_id];
    double mu1 = m_l1 * m_l2 / (m_l1 + m_l2);
    double mu2 = (m_l1 + m_l2) * m_h / (m_l1 + m_l2 + m_h);
    double sigma1 = 1. / sqrt(mu1 * baryon_omega_dict[heavy_id]);
    double sigma2 = 1. / sqrt(mu2 * baryon_omega_dict[heavy_id]);
    double p1 = sqrt(pow(p_l1.x(), 2) + pow(p_l1.y(), 2) + pow(p_l1.z(), 2));
    double p2 = sqrt(pow(p_l2.x(), 2) + pow(p_l2.y(), 2) + pow(p_l2.z(), 2));
    double p1new = sqrt(pow(p_l1new.x(), 2) + pow(p_l1new.y(), 2) + pow(p_l1new.z(), 2));
    double p2new = sqrt(pow(p_l2new.x(), 2) + pow(p_l2new.y(), 2) + pow(p_l2new.z(), 2));
    double ratio = (gq * FD_dist(p_l1new.t(), T) * gq * FD_dist(p_l2new.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p1new, 2) + pow(mass_dict[21], 2)), T) * gq * FD_dist(p_l2new.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p2new, 2) + pow(mass_dict[21], 2)), T) * gq * FD_dist(p_l1new.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p1new, 2) + pow(mass_dict[21], 2)), T) * 8. / 3. * gg * BE_dist(sqrt(pow(2. * p2new, 2) + pow(mass_dict[21], 2)), T)) /
                   (gq * FD_dist(p_l1.t(), T) * gq * FD_dist(p_l2.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p1, 2) + pow(mass_dict[21], 2)), T) * gq * FD_dist(p_l2.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p2, 2) + pow(mass_dict[21], 2)), T) * gq * FD_dist(p_l1.t(), T) +
                    8. / 3. * gg * BE_dist(sqrt(pow(2. * p1, 2) + pow(mass_dict[21], 2)), T) * 8. / 3. * gg * BE_dist(sqrt(pow(2. * p2, 2) + pow(mass_dict[21], 2)), T));
    double exponent = (-q12(p_l1new, p_l2new, p_h) * pow(sigma1, 2) - q22(p_l1new, p_l2new, p_h) * pow(sigma2, 2)) -
                      (-q12(p_l1, p_l2, p_h) * pow(sigma1, 2) - q22(p_l1, p_l2, p_h) * pow(sigma2, 2));
    ratio *= exp(exponent);
    return ratio;
}

double Coal_Sampler::prob2baryon(int light_id1, FourVector &p_l1, int light_id2, FourVector &p_l2, int heavy_id, FourVector &p_h)
{
    double m_l1 = mass_dict[light_id1];
    double m_l2 = mass_dict[light_id2];
    double m_h = mass_dict[heavy_id];
    double mu1 = m_l1 * m_l2 / (m_l1 + m_l2);
    double mu2 = (m_l1 + m_l2) * m_h / (m_l1 + m_l2 + m_h);
    double sigma1 = 1. / sqrt(mu1 * baryon_omega_dict[heavy_id]);
    double sigma2 = 1. / sqrt(mu2 * baryon_omega_dict[heavy_id]);
    double p1 = sqrt(pow(p_l1.x(), 2) + pow(p_l1.y(), 2) + pow(p_l1.z(), 2));
    double p2 = sqrt(pow(p_l2.x(), 2) + pow(p_l2.y(), 2) + pow(p_l2.z(), 2));    
    return (
               gq * FD_dist(p_l1.t(), T) * gq * FD_dist(p_l2.t(), T) +
               8. / 3. * gg * BE_dist(sqrt(pow(2. * p1, 2) + pow(mass_dict[21], 2)), T) * gq * FD_dist(p_l2.t(), T) +
               8. / 3. * gg * BE_dist(sqrt(pow(2. * p2, 2) + pow(mass_dict[21], 2)), T) * gq * FD_dist(p_l1.t(), T) +
               8. / 3. * gg * BE_dist(sqrt(pow(2. * p1, 2) + pow(mass_dict[21], 2)), T) * 8. / 3. * gg * BE_dist(sqrt(pow(2. * p2, 2) + pow(mass_dict[21], 2)), T)) *
           Wigner(p_l1, p_l2, p_h, sigma1, sigma2);
           
}


double Coal_Sampler::mc_integral(int light_id, int heavy_id, FourVector &p_h)
{
    double sum = 0.;
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    double p_xstart=0.,p_ystart=0.,p_zstart=0.;
    for (int i = 0; i < 20;i++)
    {
        FourVector p_l_start = mc_sample(light_id, heavy_id, p_h);
        p_xstart += p_l_start.x();
        p_ystart += p_l_start.y();
        p_zstart += p_l_start.z();
    }
    p_xstart /= 20.;
    p_ystart /= 20.;
    p_zstart /= 20.;
    std::uniform_real_distribution<double> uniform_distx(p_xstart - 0.5 * pl_max, p_xstart + 0.5 * pl_max);
    std::uniform_real_distribution<double> uniform_disty(p_ystart - 0.5 * pl_max, p_ystart + 0.5 * pl_max);
    std::uniform_real_distribution<double> uniform_distz(p_zstart - 0.5 * pl_max, p_zstart + 0.5 * pl_max);
    double m_l = mass_dict[light_id];
    FourVector p_h_fluid = boost(p_h, v_x, v_y, v_z);
    for (int i = 0; i < mc_max_iter; i++)
    {
        double p_x = uniform_distx(generator);
        double p_y = uniform_disty(generator);
        double p_z = uniform_distz(generator);
        double p = sqrt(p_x * p_x + p_y * p_y + p_z * p_z);
        double E_l = sqrt(p * p + m_l * m_l);
        FourVector p_l(p_x, p_y, p_z, E_l);
        sum += prob2meson(light_id, p_l, heavy_id, p_h_fluid);
    }
    return sum / mc_max_iter * pow(pl_max, 3) / pow(2. * pi, 3);
}


FourVector Coal_Sampler::mc_sample(int light_id, int heavy_id, FourVector &p_h)
{
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> uniform_dist(-pl_max, pl_max);
    std::uniform_real_distribution<double> ZeroOne_dist(0., 1.);
    int iter_n = mc_burnin + mc_burnin * ZeroOne_dist(generator);
    FourVector p_h_fluid = boost(p_h, v_x, v_y, v_z);
    double m_l = mass_dict[light_id];
    double p_x = uniform_dist(generator);
    double p_y = uniform_dist(generator);
    double p_z = uniform_dist(generator);
    double p = sqrt(p_x * p_x + p_y * p_y + p_z * p_z);
    double E_l = sqrt(p * p + m_l * m_l);
    double jump_width1 = max(p_h_fluid.x(), pl_max);
    double jump_width2 = max(p_h_fluid.y(), pl_max);
    double jump_width3 = max(p_h_fluid.z(), pl_max);
    int jump_count = 0;
    FourVector p_l(p_x, p_y, p_z, E_l);
    double alpha;
    for (int i = 0; i < 2 * mc_burnin; i++)
    {
        double tempp_x = p_l.x() + jump_width1 * (ZeroOne_dist(generator) - 0.5);
        double tempp_y = p_l.y() + jump_width2 * (ZeroOne_dist(generator) - 0.5);
        double tempp_z = p_l.z() + jump_width3 * (ZeroOne_dist(generator) - 0.5);
        p = sqrt(tempp_x * tempp_x + tempp_y * tempp_y + tempp_z * tempp_z);
        E_l = sqrt(p * p + m_l * m_l);
        FourVector tempp_l(tempp_x, tempp_y, tempp_z, E_l);
        alpha = alpha_meson(light_id, p_l, tempp_l, heavy_id, p_h_fluid);

        double u = ZeroOne_dist(generator);
        if (alpha >= 1. || u < alpha)
        {
            p_l = tempp_l;
            jump_count++;
        }
        if (i >= iter_n && jump_count > jump_time)
        {
            break;
        }
    }
    return p_l;
}

double Coal_Sampler::mc_integral(int light_id1, int light_id2, int heavy_id, FourVector &p_h)
{
    double sum = 0.;
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    double p_xstart=0.,p_ystart=0.,p_zstart=0.;
    for (int i = 0; i < 20;i++)
    {
        FourVector p_l_start = mc_sample(light_id1, light_id2, heavy_id, p_h);
        p_xstart += p_l_start.x();
        p_ystart += p_l_start.y();
        p_zstart += p_l_start.z();
    }
    p_xstart /= 20.;
    p_ystart /= 20.;
    p_zstart /= 20.;
    std::uniform_real_distribution<double> uniform_distx(p_xstart - 0.5 * pl_max, p_xstart + 0.5 * pl_max);
    std::uniform_real_distribution<double> uniform_disty(p_ystart - 0.5 * pl_max, p_ystart + 0.5 * pl_max);
    std::uniform_real_distribution<double> uniform_distz(p_zstart - 0.5 * pl_max, p_zstart + 0.5 * pl_max);
    FourVector p_h_fluid = boost(p_h, v_x, v_y, v_z);
    double m_l1 = mass_dict[light_id1];
    double m_l2 = mass_dict[light_id2];
    for (int i = 0; i < mc_max_iter; i++)
    {
        double p1_x = uniform_distx(generator);
        double p1_y = uniform_disty(generator);
        double p1_z = uniform_distz(generator);
        double p1 = sqrt(p1_x * p1_x + p1_y * p1_y + p1_z * p1_z);
        double p2_x = uniform_distx(generator);
        double p2_y = uniform_disty(generator);
        double p2_z = uniform_distz(generator);
        double p2 = sqrt(p2_x * p2_x + p2_y * p2_y + p2_z * p2_z);
        double E_l1 = sqrt(p1 * p1 + m_l1 * m_l1);
        double E_l2 = sqrt(p2 * p2 + m_l2 * m_l2);
        FourVector p_l1(p1_x, p1_y, p1_z, E_l1);
        FourVector p_l2(p2_x, p2_y, p2_z, E_l2);

        sum += prob2baryon(light_id1, p_l1, light_id2, p_l2, heavy_id, p_h_fluid);
    }
    return sum / mc_max_iter * pow(pl_max, 6) / pow(2. * pi, 6);
}


FourVector Coal_Sampler::mc_sample(int light_id1, int light_id2, int heavy_id, FourVector &p_h)
{
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> uniform_dist(-pl_max, pl_max);
    std::uniform_real_distribution<double> ZeroOne_dist(0., 1.);
    int iter_n = mc_burnin + mc_burnin * ZeroOne_dist(generator);
    FourVector p_h_fluid = boost(p_h, v_x, v_y, v_z);
    double m_l1 = mass_dict[light_id1];
    double p_x1 = uniform_dist(generator);
    double p_y1 = uniform_dist(generator);
    double p_z1 = uniform_dist(generator);
    double p1 = sqrt(p_x1 * p_x1 + p_y1 * p_y1 + p_z1 * p_z1);
    double E_l1 = sqrt(p1 * p1 + m_l1 * m_l1);
    double m_l2 = mass_dict[light_id2];
    double p_x2 = uniform_dist(generator);
    double p_y2 = uniform_dist(generator);
    double p_z2 = uniform_dist(generator);
    double p2 = sqrt(p_x2 * p_x2 + p_y2 * p_y2 + p_z2 * p_z2);
    double E_l2 = sqrt(p2 * p2 + m_l2 * m_l2);
    FourVector p_l1(p_x1, p_y1, p_z1, E_l1);
    FourVector p_l2(p_x2, p_y2, p_z2, E_l2);
    double jump_width1 = max(p_h_fluid.x(), pl_max);
    double jump_width2 = max(p_h_fluid.y(), pl_max);
    double jump_width3 = max(p_h_fluid.z(), pl_max);
    int jump_count = 0;
    double alpha;
    for (int i = 0; i < 2 * mc_burnin; i++)
    {
        double tempp_x1 = p_l1.x() + jump_width1 * (ZeroOne_dist(generator) - 0.5);
        double tempp_y1 = p_l1.y() + jump_width2 * (ZeroOne_dist(generator) - 0.5);
        double tempp_z1 = p_l1.z() + jump_width3 * (ZeroOne_dist(generator) - 0.5);
        double tempp_x2 = p_l2.x() + jump_width1 * (ZeroOne_dist(generator) - 0.5);
        double tempp_y2 = p_l2.y() + jump_width2 * (ZeroOne_dist(generator) - 0.5);
        double tempp_z2 = p_l2.z() + jump_width3 * (ZeroOne_dist(generator) - 0.5);
        p1 = sqrt(tempp_x1 * tempp_x1 + tempp_y1 * tempp_y1 + tempp_z1 * tempp_z1);
        p2 = sqrt(tempp_x2 * tempp_x2 + tempp_y2 * tempp_y2 + tempp_z2 * tempp_z2);
        E_l1 = sqrt(p1 * p1 + m_l1 * m_l1);
        E_l2 = sqrt(p2 * p2 + m_l2 * m_l2);
        FourVector tempp_l1(tempp_x1, tempp_y1, tempp_z1, E_l1);
        FourVector tempp_l2(tempp_x2, tempp_y2, tempp_z2, E_l2);
        alpha = alpha_baryon(light_id1, p_l1, tempp_l1, light_id2, p_l2, tempp_l2, heavy_id, p_h_fluid);

        double u = ZeroOne_dist(generator);
        if (alpha >= 1. || u < alpha)
        {
            p_l1 = tempp_l1;
            p_l2 = tempp_l2;
            jump_count++;
        }
        if (i >= iter_n && jump_count > jump_time)
        {
            break;
        }
    }
    return p_l1;
}

double Coal_Sampler::recomb_prob(int heavy_id, FourVector &p_h)
{
    std::vector<int> meson_list;
    std::vector<int> baryon_list;
    if (heavy_id == 4)
    {
        meson_list = charm_meson_list;
        baryon_list = charm_baryon_list;
    }
    else if(heavy_id==5)
    {
        meson_list = bottom_meson_list;
        baryon_list = bottom_baryon_list;
    }
    else
    {
        std::cerr << "not a heavy quark!" << endl;
    }
    std::vector<double> meson_prob;
    std::vector<double> baryon_prob;
    double tot_prob = 0.;
    for (auto pid : meson_list)
    {
        double temp = gM_dict[pid] * mc_integral(meson_dict[pid], heavy_id, p_h);
        meson_prob.push_back(temp);
        cout << pid <<" " << temp <<endl;
        tot_prob += temp;
    }
    for(auto pid : baryon_list)
    {
        double temp = gB_dict[pid] * mc_integral(baryon_dict[pid].first, baryon_dict[pid].second, heavy_id, p_h);
        baryon_prob.push_back(temp);
        cout << pid <<" " << temp <<endl;
        tot_prob += temp;
    }
    return tot_prob;
}

int main(int argc, char **argv)
{
    Coal_Sampler sampler(0, 0, 0, 165);

    
    for (int i = 0; i < 100; i++)
    {
        double p_x = 200 * i;
        FourVector p_h({p_x, 0, 0, sqrt(4200*4200+p_x*p_x)});
        cout << p_x <<" "<< sampler.recomb_prob(5,p_h) << endl;
        //cout << p_x <<" "<< sampler.mc_integral(1, 3,  5, p_h) << endl;

    }
    
    /*
    for (int i = 0; i < 100; i++)
    {
        double p_xh = strtod(argv[1],NULL);
        double p_yh = 0.;
        double p_zh = -30000.;
        double E_h = sqrt(pow(p_xh, 2) + pow(p_yh, 2) + pow(p_zh, 2) + pow(4200, 2));
        FourVector p_h(p_xh, p_yh, p_zh, E_h);
        FourVector p_l = sampler.mc_sample(1, 1, 5, p_h);
        cout << p_l.x() << " " << p_l.y() << " " << p_l.z() << " " << p_l.t() << endl;
    }
    */


    return 0;
}
