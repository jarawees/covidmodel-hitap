// reporter.h

#ifndef REPORTER_H
#define REPORTER_H

struct Parameters;

#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <limits>
#include <initializer_list>
#include <string>

using namespace std;

enum ReportIndex
{
    riS = 0, //1
    riSv_l,
    riSv_m,
    riSv_h,
    riE, //5
    riEv_l, 
    riEv_m,
    riEv_h,
    riIp, 
    riIp_l,//10
    riIp_m,
    riIp_h,
    riIs, 
    riIs_l,
    riIs_m, //15
    riIs_h,
    riIa, 
    riIa_l,
    riIa_m,
    riIa_h,//20 
    riR,
    riRv_l,
    riRv_m,
    riRv_h,
    ricases, //25 
    ricases_reported, 
    risubclinical,
    rilambda,
    rilambdav_l,
    rilambdav_m, //30
    rilambdav_h

};

const vector<string> ref_col_names = {
    "S",  //1
    "Sv_l",
    "Sv_m",
    "Sv_h",
    "E", //5
    "Ev_l", 
    "Ev_m",
    "Ev_h",
    "Ip", 
    "Ip_l",//10
    "Ip_m",
    "Ip_h",
    "Is", 
    "Is_l",
    "Is_m", //15
    "Is_h",
    "Ia", 
    "Ia_l",
    "Ia_m",
    "Ia_h",//20 
    "R",
    "Rv_l",
    "Rv_m",
    "Rv_h",
    "cases", //25 
    "cases_reported",
    "subclinical",
    "foi",
    "foiv_l",
    "foiv_m",//30
    "foiv_h"
};

// For reporting results
class Reporter
{
public:
    Reporter(Parameters& P);

    // Access / modify data
    double& operator()(double t, unsigned int p, unsigned int a, unsigned int c)
    {
        unsigned int row = (unsigned int)(t - t0) * n_populations * n_age_groups + p * n_age_groups + a;
        return data[c][row];
    }

    // Access data, summed over populations and groups
    double operator()(string compartment, double t, initializer_list<unsigned int> p, initializer_list<unsigned int> a);

    // Access / modify observer data
    double& Obs(double t, unsigned int p, unsigned int a, unsigned int c)
    {
        if (c >= obs.size())
            obs.resize(c + 1, vector<double>(n_times * n_populations * n_age_groups, 0.));
        unsigned int row = (unsigned int)(t - t0) * n_populations * n_age_groups + p * n_age_groups + a;
        return obs[c][row];
    }

    // Save data to file
    void Save(string basename, unsigned long int seed);

//private:
    double t0;
    unsigned int n_times;
    unsigned int n_populations;
    unsigned int n_age_groups;
    vector<string> col_names;
    unsigned int user_defined_offset;

    vector<vector<double>> data;
    vector<vector<double>> obs;
    string csv;
};

#endif 
