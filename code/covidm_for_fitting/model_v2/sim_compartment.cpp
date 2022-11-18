// sim_compartment.cpp

#include "sim_compartment.h"
#include "parameters.h"
#include "reporter.h"
#include "randomizer.h"
#include "user_defined.h"
//
// MODEL DYNAMICS
//

Population::Population(Parameters& P, unsigned int pindex)
 : seed_row(0), p(pindex)
{
    // Set up built-in compartments
    N  = P.pop[p].size;
    S  = N;
    Sv_l = vector<double>(S.size(), 0.);
    Sv_m = vector<double>(S.size(), 0.);
    Sv_h = vector<double>(S.size(), 0.);
    
    E  = vector<Compartment>(S.size());
    Ev_l = vector<Compartment>(S.size());
    Ev_m = vector<Compartment>(S.size());
    Ev_h = vector<Compartment>(S.size());

    Ip = vector<Compartment>(S.size());
    Ia = vector<Compartment>(S.size());
    Is = vector<Compartment>(S.size());
    
    Ip_l = vector<Compartment>(S.size());
    Ia_l = vector<Compartment>(S.size());
    Is_l = vector<Compartment>(S.size());
    
    Ip_m = vector<Compartment>(S.size());
    Ia_m = vector<Compartment>(S.size());
    Is_m = vector<Compartment>(S.size());
    
    Ip_h = vector<Compartment>(S.size());
    Ia_h = vector<Compartment>(S.size());
    Is_h = vector<Compartment>(S.size());
    
    C  = vector<Compartment>(S.size());

    R  = vector<double>(S.size(), 0.);
    Rv_l = vector<double>(S.size(), 0.);
    Rv_m = vector<double>(S.size(), 0.);
    Rv_h = vector<double>(S.size(), 0.);

    // Initial immunity
    for (unsigned int a = 0; a < S.size(); ++a) {
        double imm = 0;
        if (P.deterministic) imm = S[a] * P.pop[p].imm0[a];
        else                 imm = (unsigned int)(S[a] * P.pop[p].imm0[a] + 0.5);
        S[a] -= imm;
        R[a] += imm;
    }

    pc = vector<vector<Compartment>>(
        P.processes.state_count, vector<Compartment>(S.size())
    );
    pci = vector<double>(pc.size(), 0.);
    pco = vector<double>(pc.size(), 0.);

}

// Do seeding and calculate contagiousness
void Population::Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag)
{
    auto add = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        E[age].Add(P, Rand, n, P.pop[p].dE);
    };

    // Do seeding
    while (seed_row < P.pop[p].seed_times.size() && t >= P.pop[p].seed_times[seed_row])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add(a, 1);
                    break;
                }
            }
        }
        ++seed_row;
    }

    // Calculate contagiousness from this population
    for (unsigned int a = 0; a < contag.size(); ++a)
        contag[a] = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip[a].Size() + P.pop[p].fIa[a] * Ia[a].Size() + P.pop[p].fIs[a] * Is[a].Size() + P.pop[p].fIp[a] * Ip_l[a].Size() + P.pop[p].fIa[a] * Ia_l[a].Size() + P.pop[p].fIs[a] * Is_l[a].Size() + P.pop[p].fIp[a] * Ip_m[a].Size() + P.pop[p].fIa[a] * Ia_m[a].Size() + P.pop[p].fIs[a] * Is_m[a].Size() + P.pop[p].fIp[a] * Ip_h[a].Size() + P.pop[p].fIa[a] * Ia_h[a].Size() + P.pop[p].fIs[a] * Is_h[a].Size()) / N[a];
}

// Execute one time step's events
void Population::Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, Reporter& rep)
{
    // Calculate force of infection in this compartment
    lambda.assign(infec.size(), 0.0);

    // NB, this now does not include susceptibility term
    for (unsigned int a = 0; a < lambda.size(); ++a)
        for (unsigned int b = 0; b < lambda.size(); ++b)
            lambda[a] += P.pop[p].cm(a,b) * infec[b];

    // Account for seasonality
    if (P.pop[p].season_A != 0)
    {
        double f = 1.0 + P.pop[p].season_A * cos(2. * M_PI * (t - P.pop[p].season_phi) / P.pop[p].season_T);
        for (unsigned int a = 0; a < lambda.size(); ++a)
            lambda[a] = lambda[a] * f;
    }

    // TODO: applying u later changes the meaning of omega term
    // Account for importation
    for (unsigned int a = 0; a < lambda.size(); ++a)
        lambda[a] += P.pop[p].omega[a];

    // TODO: construction + internal testing in these helpers every Tick invocation is code smell
    // maybe compiler is smart enough to get it? seems unlikely
    // Helpers
    auto multinomial = [&](double n, vector<double>& p, vector<double>& nd_out, vector<unsigned int>& ni_out) {
        nd_out.resize(p.size(), 0.);
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = n * p[i];
        }
        else
        {
            ni_out.resize(p.size(), 0);
            Rand.Multinomial(n, p, ni_out);
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = ni_out[i];
        }
    };

    auto poisson = [&](double l) {
        if (P.deterministic)
            return l;
        else
            return (double)Rand.Poisson(l);
    };

    auto binomial = [&](double n, double p) {
        if (P.deterministic)
            return n * p;
        else
            return (double)Rand.Binomial(n, p);
    };

    auto num = [&](double n) {
        if (P.deterministic)
            return n;
        else
            return round(n);
    };

    // Do state transitions and reporting for each age group
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        // 0. Report prevalences
        if (t == (int)t)
        {
            // TODO magic number code smell
            // Built-in states
            rep(t, p, a, riS)    = S[a];
            rep(t, p, a, riSv_l)  = Sv_l[a];
            rep(t, p, a, riSv_m)  = Sv_m[a];
            rep(t, p, a, riSv_h)  = Sv_h[a];
            
            rep(t, p, a, riE)    = E[a].Size();
            rep(t, p, a, riEv_l) = Ev_l[a].Size();
            rep(t, p, a, riEv_m) = Ev_m[a].Size();
            rep(t, p, a, riEv_h) = Ev_h[a].Size();
            
            rep(t, p, a, riIp) = Ip[a].Size();
            rep(t, p, a, riIs) = Is[a].Size();
            rep(t, p, a, riIa) = Ia[a].Size();
            
            rep(t, p, a, riIp_l) = Ip_l[a].Size();
            rep(t, p, a, riIs_l) = Is_l[a].Size();
            rep(t, p, a, riIa_l) = Ia_l[a].Size();
            
            rep(t, p, a, riIp_m) = Ip_m[a].Size();
            rep(t, p, a, riIs_m) = Is_m[a].Size();
            rep(t, p, a, riIa_m) = Ia_m[a].Size();
            
            rep(t, p, a, riIp_h) = Ip_h[a].Size();
            rep(t, p, a, riIs_h) = Is_h[a].Size();
            rep(t, p, a, riIa_h) = Ia_h[a].Size();
            
            rep(t, p, a, riR)    = R[a];
            rep(t, p, a, riRv_l) = Rv_l[a];
            rep(t, p, a, riRv_m) = Rv_m[a];
            rep(t, p, a, riRv_h) = Rv_h[a];
            
            rep(t, p, a, rilambda)    = P.pop[p].u[a]*lambda[a];
            rep(t, p, a, rilambdav_l) = P.pop[p].uv_l[a]*lambda[a];
            rep(t, p, a, rilambdav_m) = P.pop[p].uv_m[a]*lambda[a];
            rep(t, p, a, rilambdav_h) = P.pop[p].uv_h[a]*lambda[a];    
            
            // User-specified processes
            for (size_t i=0;
                 i < P.processes.prevalence_states.size();
                 i++
            ) rep(t, p, a, rep.user_defined_offset + i) = pc[P.processes.prevalence_states[i]][a].Size();

        }

        // 1. Built-in states
        double primary_dose_eligible   = S[a] + E[a].Size() + Ip[a].Size() + Ia[a].Size() + R[a];
        double booster_dose_eligible_l = Sv_l[a] + Ev_l[a].Size() + Ip_l[a].Size() + Ia_l[a].Size() + Rv_l[a];
        double booster_dose_eligible_m = Sv_m[a] + Ev_m[a].Size() + Ip_m[a].Size() + Ia_m[a].Size() + Rv_m[a];
        double booster_dose_eligible_h = Sv_h[a] + Ev_h[a].Size() + Ip_h[a].Size() + Ia_h[a].Size() + Rv_h[a];
        double booster_dose_eligible = booster_dose_eligible_l + booster_dose_eligible_m + booster_dose_eligible_h;
        
        // initial vaccination campaign, primarys doses
        // (2-4) S -> Sv_m; S -> Sv_l; S -> Sv_h
        // (46-48) R -> Rv_m; R -> Rv_l; R -> Rv_h
        // min S[a] potentially problematic, used twice when doses > humans
        double nS_Sv_l = min(num(P.pop[p].ev_p[a] * S[a] * P.time_step * P.pop[p].v_p_2l[a]),                            num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * (S[a] / primary_dose_eligible) * P.time_step * P.pop[p].v_p_2l[a]));
        double nS_Sv_m = min(num(P.pop[p].ev_p[a] * S[a] * P.time_step * P.pop[p].v_p_2m[a]),                            num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * (S[a] / primary_dose_eligible) * P.time_step * P.pop[p].v_p_2m[a]));
        double nS_Sv_h = min(num(P.pop[p].ev_p[a] * S[a] * P.time_step * (1 - P.pop[p].v_p_2l[a] - P.pop[p].v_p_2m[a])), num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * (S[a] / primary_dose_eligible) * P.time_step * (1 - P.pop[p].v_p_2l[a] - P.pop[p].v_p_2m[a])));
        
        double nR_Rv_l = min(num(P.pop[p].ev_p[a] * R[a] * P.time_step * P.pop[p].v_p_2l[a]),                            num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * (R[a] / primary_dose_eligible) * P.time_step * P.pop[p].v_p_2l[a]));
        double nR_Rv_m = min(num(P.pop[p].ev_p[a] * R[a] * P.time_step * P.pop[p].v_p_2m[a]),                            num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * (R[a] / primary_dose_eligible) * P.time_step * P.pop[p].v_p_2m[a]));
        double nR_Rv_h = min(num(P.pop[p].ev_p[a] * R[a] * P.time_step * (1 - P.pop[p].v_p_2l[a] - P.pop[p].v_p_2m[a])), num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * (R[a] / primary_dose_eligible) * P.time_step * (1 - P.pop[p].v_p_2l[a] - P.pop[p].v_p_2m[a])));
        
        S[a]    -=  nS_Sv_l + nS_Sv_m + nS_Sv_h;
        Sv_l[a] += nS_Sv_l;
        Sv_m[a] += nS_Sv_m;
        Sv_h[a] += nS_Sv_h;
        
        R[a]    -= nR_Rv_l + nR_Rv_m + nR_Rv_h;
        Rv_l[a] += nR_Rv_l;
        Rv_m[a] += nR_Rv_m;
        Rv_h[a] += nR_Rv_h;
        
        // booster vaccination campaign, booster doses
        // (7-9) Sv_l -> Sv_m; Sv_m -> Sv_h; Sv_l -> Sv_h 
        // (41-43) Rv_l -> Rv_m; Rv_m -> Rv_h; Rv_l -> Rv_h
        double nSv_l_Sv_m = min(Sv_l[a] * P.pop[p].v_b_l2m[a],       num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * (Sv_l[a] / booster_dose_eligible) * P.time_step * P.pop[p].v_b_l2m[a]));
        double nSv_l_Sv_h = min(Sv_l[a] * (1 - P.pop[p].v_b_l2m[a]), num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * (Sv_l[a] / booster_dose_eligible) * P.time_step * (1-P.pop[p].v_b_l2m[a])));
        double nSv_m_Sv_h = min(Sv_m[a],                             num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * (Sv_m[a] / booster_dose_eligible) * P.time_step));
        
        
        double nRv_l_Rv_m = min(Rv_l[a] * P.pop[p].v_b_l2m[a],       num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * (Rv_l[a] / booster_dose_eligible) * P.time_step * P.pop[p].v_b_l2m[a]));
        double nRv_l_Rv_h = min(Rv_l[a] * (1 - P.pop[p].v_b_l2m[a]), num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * (Rv_l[a] / booster_dose_eligible) * P.time_step * (1-P.pop[p].v_b_l2m[a])));
        double nRv_m_Rv_h = min(Rv_m[a],                             num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * (Rv_m[a] / booster_dose_eligible) * P.time_step));
        
        Sv_l[a] -= nSv_l_Sv_m + nSv_l_Sv_h;
        Sv_m[a] += nSv_l_Sv_m - nSv_m_Sv_h;
        Sv_h[a] += nSv_l_Sv_h + nSv_m_Sv_h;
        
        Rv_l[a] -= nRv_l_Rv_m + nRv_l_Rv_h;
        Rv_m[a] += nRv_l_Rv_m - nRv_m_Rv_h;
        Rv_h[a] += nRv_l_Rv_h + nRv_m_Rv_h;
        
        // waning of vaccine induced protection
        // (5-6) Sv_m -> Sv_l; Sv_h -> Sv_m
        // (44-45) Rv_m -> Rv_l; Rv_h -> Rv_m`
        double nSv_m_Sv_l = binomial(Sv_m[a], 1.0 - exp(-P.pop[p].wv_m2l[a] * P.time_step));
        double nSv_h_Sv_m = binomial(Sv_h[a], 1.0 - exp(-P.pop[p].wv_h2m[a] * P.time_step));
        double nRv_m_Rv_l = binomial(Rv_m[a], 1.0 - exp(-P.pop[p].wv_m2l[a] * P.time_step));
        double nRv_h_Rv_m = binomial(Rv_h[a], 1.0 - exp(-P.pop[p].wv_h2m[a] * P.time_step));
        
        Sv_l[a] += nSv_m_Sv_l;
        Sv_m[a] += nSv_h_Sv_m - nSv_m_Sv_l;
        Sv_h[a] -= nSv_h_Sv_m;
        
        Rv_l[a] += nRv_m_Rv_l;
        Rv_m[a] += nRv_h_Rv_m - nRv_m_Rv_l;
        Rv_h[a] -= nRv_h_Rv_m;
        
        // waning of infection induced protection
        // (29) R -> S
        // (32) Rv_l -> Sv_l
        // (35) Rv_m -> Sv_m
        // (38) Rv_h -> Sv_h
        
        double nR_S = binomial(R[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        double nRv_l_Sv_l = binomial(Rv_l[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        double nRv_m_Sv_m = binomial(Rv_m[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        double nRv_h_Sv_h = binomial(Rv_h[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        
        R[a] -= nR_S;
        S[a] += nR_S;

        Rv_l[a] -= nRv_l_Sv_l;
        Sv_l[a] += nRv_l_Sv_l;
        
        Rv_m[a] -= nRv_m_Sv_m;
        Sv_m[a] += nRv_m_Sv_m;
        
        Rv_h[a] -= nRv_h_Sv_h;
        Sv_h[a] += nRv_h_Sv_h;
        
        // Infecting
        // (1) S -> E and (31) R -> E
        double nS_E = binomial(S[a], 1.0 - exp(-P.pop[p].u[a]*lambda[a] * P.time_step));
        double nR_E = binomial(R[a], 1.0 - exp(-P.pop[p].ur[a]*lambda[a] * P.time_step));
        S[a] -= nS_E;
        R[a] -= nR_E;
        E[a].Add(P, Rand, nS_E + nR_E, P.pop[p].dE);
        
        // (10) Sv_l -> Ev_l and (34) Rv_l -> Ev_l
        double nSv_l_Ev_l = binomial(Sv_l[a], 1.0 - exp(-P.pop[p].uv_l[a]*lambda[a] * P.time_step));
        double nRv_l_Ev_l = binomial(Rv_l[a], 1.0 - exp(-P.pop[p].uvr_l[a]*lambda[a] * P.time_step));
        Sv_l[a] -= nSv_l_Ev_l;
        Rv_l[a] -= nRv_l_Ev_l;
        Ev_l[a].Add(P, Rand, nSv_l_Ev_l + nRv_l_Ev_l, P.pop[p].dEv_l);

        // (11) Sv_m -> Ev_m and (37) Rv_m -> Ev_m
        double nSv_m_Ev_m = binomial(Sv_m[a], 1.0 - exp(-P.pop[p].uv_m[a]*lambda[a] * P.time_step));
        double nRv_m_Ev_m = binomial(Rv_m[a], 1.0 - exp(-P.pop[p].uvr_m[a]*lambda[a] * P.time_step));
        Sv_m[a] -= nSv_m_Ev_m;
        Rv_m[a] -= nRv_m_Ev_m;
        Ev_m[a].Add(P, Rand, nSv_m_Ev_m + nRv_m_Ev_m, P.pop[p].dEv_m);
        
        // (12) Sv_h -> Ev_h and (40) Rv_h -> Ev_h
        double nSv_h_Ev_h = binomial(Sv_h[a], 1.0 - exp(-P.pop[p].uv_h[a]*lambda[a] * P.time_step));
        double nRv_h_Ev_h = binomial(Rv_h[a], 1.0 - exp(-P.pop[p].uvr_h[a]*lambda[a] * P.time_step));
        Sv_h[a] -= nSv_h_Ev_h;
        Rv_h[a] -= nRv_h_Ev_h;
        Ev_h[a].Add(P, Rand, nSv_h_Ev_h + nRv_h_Ev_h, P.pop[p].dEv_h);
    
        // progressing from infection
        // (13)-(14) E -> Ip and E -> Ia
        double nE_Ipa = E[a].Mature();
        double nE_Ip = binomial(nE_Ipa, P.pop[p].y[a]);
        double nE_Ia = nE_Ipa - nE_Ip;
        Ip[a].Add(P, Rand, nE_Ip, P.pop[p].dIp);
        Ia[a].Add(P, Rand, nE_Ia, P.pop[p].dIa);
        
        // (15)-(16) Ev_l -> Ip_l and Ev_l -> Is_l 
        double nEv_l_Ipa = Ev_l[a].Mature();
        double nEv_l_Ip = binomial(nEv_l_Ipa, P.pop[p].yv_l[a]);
        double nEv_l_Ia = nEv_l_Ipa - nEv_l_Ip;
        Ip_l[a].Add(P, Rand, nEv_l_Ip, P.pop[p].dIp_l);
        Ia_l[a].Add(P, Rand, nEv_l_Ia, P.pop[p].dIa_l);
        
        // (17)-(18) Ev_m -> Ip_m and Ev_m -> Is_m
        double nEv_m_Ipa = Ev_m[a].Mature();
        double nEv_m_Ip = binomial(nEv_m_Ipa, P.pop[p].yv_m[a]);
        double nEv_m_Ia = nEv_m_Ipa - nEv_m_Ip;
        Ip_m[a].Add(P, Rand, nEv_m_Ip, P.pop[p].dIp_m);
        Ia_m[a].Add(P, Rand, nEv_m_Ia, P.pop[p].dIa_m);
        
        // (19)-(20) Ev_h -> Ip_h and Ev_h -> Is_h
        double nEv_h_Ipa = Ev_h[a].Mature();
        double nEv_h_Ip = binomial(nEv_h_Ipa, P.pop[p].yv_h[a]);
        double nEv_h_Ia = nEv_h_Ipa - nEv_h_Ip;
        Ip_h[a].Add(P, Rand, nEv_h_Ip, P.pop[p].dIp_h);
        Ia_h[a].Add(P, Rand, nEv_h_Ia, P.pop[p].dIa_h);
        
        // progressing from presym to sym
        // (21) Ip -> Is
        double nIp_Is = Ip[a].Mature();
        Is[a].Add(P, Rand, nIp_Is, P.pop[p].dIs);
        
        // (22) Ip_l -> Is_l
        double nIp_l_Is_l = Ip_l[a].Mature();
        Is_l[a].Add(P, Rand, nIp_l_Is_l, P.pop[p].dIs_l);
        
        // (23) Ip_m -> Is_m
        double nIp_m_Is_m = Ip_m[a].Mature();
        Is_m[a].Add(P, Rand, nIp_m_Is_m, P.pop[p].dIs_m);
        
        // (24) Ip_h -> Is_h
        double nIp_h_Is_h = Ip_h[a].Mature();
        Is_h[a].Add(P, Rand, nIp_h_Is_h, P.pop[p].dIs_h);
        
        // Reported cases
        double n_to_report = binomial(nIp_Is + nIp_l_Is_l + nIp_m_Is_m + nIp_h_Is_h, P.pop[p].rho[a]);
        C[a].Add(P, Rand, n_to_report, P.pop[p].dC);
        double n_reported = C[a].Mature();
        
        // removal
        // (30) Is -> R
        double nIs_R = Is[a].Mature();
        R[a] += nIs_R;

        // (25) Ia -> R
        double nIa_R = Ia[a].Mature();
        R[a] += nIa_R;
        
        // (33) Is_l -> R_l
        double nIs_l_Rv_l = Is_l[a].Mature();
        Rv_l[a] += nIs_l_Rv_l;
        
        // (26) Ia_l -> R_l
        double nIa_l_Rv_l = Ia_l[a].Mature();
        Rv_l[a] += nIa_l_Rv_l;
        
        // (36) Is_m -> R_m
        double nIs_m_Rv_m = Is_m[a].Mature();
        Rv_m[a] += nIs_m_Rv_m;
        
        // (27) Ia_m -> R_m
        double nIa_m_Rv_m = Ia_m[a].Mature();
        Rv_m[a] += nIa_m_Rv_m;
        
        // (39) Is_h -> R_h
        double nIs_h_Rv_h = Is_h[a].Mature();
        Rv_h[a] += nIs_h_Rv_h;
        
        // (28) Ia_h -> R_h
        double nIa_h_Rv_h = Ia_h[a].Mature();
        Rv_h[a] += nIa_h_Rv_h;
        
        // 2. User-specified processes
        // assert: processes are ordered such that when iterating
        // all sources have received all their inputs from prior processes
        // so:
        // if a negative input exists for a needed source, it can be matured from pc at that time
        // after the loop, all the pcos that are still negative can be matured from the pcs
        fill(pco.begin(), pco.end(), -1.);
        fill(pci.begin(), pci.end(), 0.);

        for (auto& process : P.processes.flows)
        {
            // Determine number of individuals entering the process
            // only if you have to called it as an src
            // the labels of these variables are in processes_spec.h
            double n_entering = 0.;
            switch (process.source_id)
            {
                // 1. all infections
            case src_newE_all:
                n_entering = nS_E + nR_E + nSv_l_Ev_l + nRv_l_Ev_l + nSv_m_Ev_m + nRv_m_Ev_m + nSv_h_Ev_h + nRv_h_Ev_h; break;
                // 2. infections with no vaccine effect
            case src_newE:
                n_entering = nS_E + nR_E; break;
                // 3. infections with low vaccine effect
            case src_newEv_l:
                n_entering = nSv_l_Ev_l + nRv_l_Ev_l; break;
                // 4. infections with medium vaccine effect
            case src_newEv_m:
                n_entering = nSv_m_Ev_m + nRv_m_Ev_m; break;
                // 5. infections with high vaccine effect
            case src_newEv_h:
                n_entering = nSv_h_Ev_h + nRv_h_Ev_h; break;
                // 6. cases reported
            case srcCasesReported:
                n_entering = n_to_report; break;
                
                default:
                    if (pco[process.source_id] < 0) {
                        pco[process.source_id] = pc[process.source_id][a].Mature();
                    }
                    n_entering = pco[process.source_id];
                    // TODO: replace below with some other kind of validation re processes ordering?
                    // if (n_entering < 0)
                    //     throw logic_error("Process sourced from unset user process:" + process.source_name + " - have user processes been specified in the right order?");
                    // break;
            }

            multinomial(n_entering, process.prob[a], nd_out, ni_out);

            // add to the relevant process compartments

            // Seed and mature this process's compartments
            unsigned int c = 0;
            for (unsigned int compartment_id : process.sink_ids)
            {
                // TODO should be irrelevant?
                if (compartment_id != Null)
                {
                    pc[compartment_id][a].Add(P, Rand, nd_out[c], process.delays[c]);
                    pci[compartment_id] += nd_out[c];
                }
                ++c;
            }
        }

        // mature all compartments not yet updated
        for (size_t i = 0; i < pco.size(); i++) if (pco[i] < 0) pco[i] = pc[i][a].Mature();

        // 3. Report incidence / outcidence
        // Built-in states
        rep(t, p, a, ricases) += nIp_Is + nIp_l_Is_l + nIp_m_Is_m + nIp_h_Is_h;
        rep(t, p, a, ricases_reported) += n_reported;
        rep(t, p, a, risubclinical) += nE_Ia + nEv_l_Ia + nEv_m_Ia + nEv_h_Ia;

        // User-specified incidence + outcidence flows
        for (size_t i=0;
             i < P.processes.incidence_states.size();
             i++
        ) rep(
            t, p, a,
            rep.user_defined_offset + P.processes.inc_offset + i
        ) += pci[P.processes.incidence_states[i]];

        for (size_t i=0;
             i < P.processes.outcidence_states.size();
             i++
        ) rep(
            t, p, a,
            rep.user_defined_offset + P.processes.out_offset + i
        ) += pco[P.processes.outcidence_states[i]];

    }

    // Births, deaths, aging
    double Ntot = accumulate(N.begin(), N.end(), 0.0);
    for (unsigned int a = N.size() - 1; ; --a)
    {
        // Births
        double B = poisson(Ntot * (exp(P.pop[p].B[a] * P.time_step) - 1.));

        // Deaths
        double death_prob = 1.0 - exp(-P.pop[p].D[a] * P.time_step);
        double DS   = binomial(S[a],        death_prob);
        double DSv_l  = binomial(Sv_l[a],        death_prob);
        double DSv_m = binomial(Sv_m[a],        death_prob);
        double DSv_h = binomial(Sv_h[a],        death_prob);

        double DR    = binomial(R[a],        death_prob);
        double DRv_l = binomial(Rv_l[a],        death_prob);
        double DRv_m = binomial(Rv_m[a],        death_prob);
        double DRv_h = binomial(Rv_h[a],        death_prob);

        // Changes
        N[a] += B;
        S[a] += B;

        S[a]  -= DS;
        Sv_l[a] -= DSv_l;
        Sv_m[a] -= DSv_m;
        Sv_h[a] -= DSv_h;
        
        R[a]  -= DR;
        Rv_l[a] -= DRv_l;
        Rv_m[a] -= DRv_m;
        Rv_h[a] -= DRv_h;

        // NB - RemoveProb also takes care of removal, so no -= needed
        double DE  = E[a] .RemoveProb(P, Rand, death_prob);
        double DEv_l = Ev_l[a].RemoveProb(P, Rand, death_prob);
        double DEv_m = Ev_m[a].RemoveProb(P, Rand, death_prob);
        double DEv_h = Ev_h[a].RemoveProb(P, Rand, death_prob);
        
        double DIp   = Ip[a].RemoveProb(P, Rand, death_prob);
        double DIp_l = Ip_l[a].RemoveProb(P, Rand, death_prob);
        double DIp_m = Ip_m[a].RemoveProb(P, Rand, death_prob);
        double DIp_h = Ip_h[a].RemoveProb(P, Rand, death_prob);
        
        double DIa = Ia[a].RemoveProb(P, Rand, death_prob);
        double DIa_l = Ia_l[a].RemoveProb(P, Rand, death_prob);
        double DIa_m = Ia_m[a].RemoveProb(P, Rand, death_prob);
        double DIa_h = Ia_h[a].RemoveProb(P, Rand, death_prob);
        
        double DIs = Is[a].RemoveProb(P, Rand, death_prob);
        double DIs_l = Is_l[a].RemoveProb(P, Rand, death_prob);
        double DIs_m = Is_m[a].RemoveProb(P, Rand, death_prob);
        double DIs_h = Is_h[a].RemoveProb(P, Rand, death_prob);

        N[a]  -=
            DS + DSv_l + DSv_m + DSv_h +
            DE + DEv_l + DEv_m + DEv_h +
            DIp + DIp_l + DIp_m + DIp_h +
            DIa + DIa_l + DIa_m + DIa_h +
            DIs + DIs_l + DIs_m + DIs_h +
            DR  + DRv_l + DRv_m + DRv_h;

        // Agings
        if (a != lambda.size() - 1)
        {
            double age_prob = 1.0 - exp(-P.pop[p].A[a] * P.time_step);
            double AS    = binomial(S[a],        age_prob);
            double ASv_l = binomial(Sv_l[a],        age_prob);
            double ASv_m = binomial(Sv_m[a],        age_prob);
            double ASv_h = binomial(Sv_h[a],        age_prob);
            
            
            double AR    = binomial(R[a],        age_prob);
            double ARv_l = binomial(Rv_l[a],        age_prob);
            double ARv_m = binomial(Rv_m[a],        age_prob);
            double ARv_h = binomial(Rv_h[a],        age_prob);

            S[a]        -= AS;
            S[a + 1]    += AS;
            Sv_l[a]     -= ASv_l;
            Sv_l[a + 1] += ASv_l;
            Sv_m[a]     -= ASv_m;
            Sv_m[a + 1] += ASv_m;
            Sv_h[a]     -= ASv_h;
            Sv_h[a + 1] += ASv_h;

            R[a]        -= AR;
            R[a + 1]    += AR;
            Rv_l[a]     -= ARv_l;
            Rv_l[a + 1] += ARv_l;
            Rv_m[a]     -= ARv_m;
            Rv_m[a + 1] += ARv_m;
            Rv_h[a]     -= ARv_h;
            Rv_h[a + 1] += ARv_h;
            
            double AE    = E[a] .MoveProb(E [a + 1], P, Rand, age_prob);
            double AEv_l = Ev_l[a].MoveProb(Ev_l[a + 1], P, Rand, age_prob);
            double AEv_m = Ev_m[a].MoveProb(Ev_m[a + 1], P, Rand, age_prob);
            double AEv_h = Ev_h[a].MoveProb(Ev_h[a + 1], P, Rand, age_prob);
            
            double AIp   = Ip[a].MoveProb(Ip[a + 1], P, Rand, age_prob);
            double AIp_l = Ip_l[a].MoveProb(Ip_l[a + 1], P, Rand, age_prob);
            double AIp_m = Ip_m[a].MoveProb(Ip_m[a + 1], P, Rand, age_prob);
            double AIp_h = Ip_m[a].MoveProb(Ip_h[a + 1], P, Rand, age_prob);
            
            double AIa   = Ia[a].MoveProb(Ia[a + 1], P, Rand, age_prob);
            double AIa_l = Ia_l[a].MoveProb(Ia_l[a + 1], P, Rand, age_prob);
            double AIa_m = Ia_m[a].MoveProb(Ia_m[a + 1], P, Rand, age_prob);
            double AIa_h = Ia_m[a].MoveProb(Ia_h[a + 1], P, Rand, age_prob);
            
            double AIs   = Is[a].MoveProb(Is[a + 1], P, Rand, age_prob);
            double AIs_l = Is_l[a].MoveProb(Is_l[a + 1], P, Rand, age_prob);
            double AIs_m = Is_m[a].MoveProb(Is_m[a + 1], P, Rand, age_prob);
            double AIs_h = Is_h[a].MoveProb(Is_h[a + 1], P, Rand, age_prob);

            N[a]      -= AS + ASv_l + ASv_m + ASv_h + AE + AEv_l + AEv_m + AEv_h + AIp + AIp_l + AIp_m + AIp_h + AIa + AIa_l + AIa_m + AIa_h + AIs + AIs_l + AIs_m + AIs_h + AR + ARv_l + ARv_m + ARv_h;
            N[a + 1]  += AS + ASv_l + ASv_m + ASv_h + AE + AEv_l + AEv_m + AEv_h + AIp + AIp_l + AIp_m + AIp_h + AIa + AIa_l + AIa_m + AIa_h + AIs + AIs_l + AIs_m + AIs_h + AR + ARv_l + ARv_m + ARv_h;
        }

        if (a == 0)
            break;
    }
}

// Print full population details
void Population::DebugPrint() const
{
    auto vecprint = [&](const vector<double>& vec, string name) {
        cout << name;
        for (auto& v : vec)
            cout << " " << v;
        cout << "\n";
    };

    auto comprint = [&](const vector<Compartment>& comp, string name) {
        cout << name;
        for (unsigned int c = 0; c < comp.size(); ++c) {
            cout << "element " << c << "\n";
            comp[c].DebugPrint();
        }
    };

    vecprint(lambda, "lambda");
    vecprint(N, "N");

    vecprint(S, "S");
    vecprint(Sv_l, "Sv_l");
    vecprint(Sv_m, "Sv_m");
    
    vecprint(R, "R");
    vecprint(Rv_l, "Rv_l");
    vecprint(Rv_m, "Rv_m");
    
    comprint(E, "E");
    comprint(Ev_l, "Ev_l");
    comprint(Ev_m, "Ev_m");
    
    comprint(Ip, "Ip");
    comprint(Ip_l, "Ip_l");
    comprint(Ip_m, "Ip_m");
    
    comprint(Ia, "Ia");
    comprint(Ia_l, "Ia_l");
    comprint(Ia_m, "Ia_m");
    
    comprint(Is, "Is");
    comprint(Is_l, "Is_l");
    comprint(Is_m, "Is_m");
    
    comprint(C, "C");
    cout << "seed_row " << seed_row << " p " << p << "\n";
    for (auto& c : pc)
        comprint(c, "User");

    cout << "\n\n";
}


Metapopulation::Metapopulation(Parameters& P)
{
    P.changes.Capture(P);

    for (unsigned int i = 0; i < P.pop.size(); ++i)
        pops.push_back(Population(P, i));
}

// Execute one time step's events
bool Metapopulation::Tick(Parameters& P, Randomizer& Rand, double t, unsigned int ts, Reporter& rep)
{
    // Apply any changes to parameters
    P.changes.Apply(P, t);

    unsigned int n_ages = P.pop[0].size.size();

    // Calculate contagiousness from each population
    // NOTE -- 'contag' subscripted first by j, then by a.
    // It's the effective number of infectious individuals FROM subpop j of age a.
    contag.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int j = 0; j < pops.size(); ++j)
        pops[j].Contagiousness(P, Rand, t, contag[j]);

    // note -- 'infec' subscripted first by i, then by a
    // It's the effective number of infectious individuals who are CURRENTLY IN subpop i of age a.
    infec.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int i = 0; i < pops.size(); ++i)
        for (unsigned int j = 0; j < pops.size(); ++j)
            for (unsigned int a = 0; a < n_ages; ++a)
                infec[i][a] += P.travel(j, i) * contag[j][a] * (j != i ? P.pop[j].tau[a] : 1.0);

    // Update populations
    //#pragma omp parallel for schedule(dynamic) reduction(&&:keep_going)
    for (unsigned int i = 0; i < pops.size(); ++i)
        pops[i].Tick(P, Rand, t, infec[i], rep);

    // Run observer at the last time step of each day.
    if (t + P.time_step == int(t + P.time_step))
        return CppObserver(P, Rand, rep, (int)t, x);

    return true;
}

void Metapopulation::Run(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double> x_fit)
{
    x = x_fit;

    #ifdef _OPENMP
    omp_set_num_threads(6);
    #endif

    // Run simulation
    unsigned int time_steps = (1 + P.time1 - P.time0) / P.time_step;
    for (unsigned int ts = 0; ts < time_steps; ++ts)
    {
        if (!Tick(P, Rand, P.time0 + ts * P.time_step, ts, rep))
            break;
    }
}
