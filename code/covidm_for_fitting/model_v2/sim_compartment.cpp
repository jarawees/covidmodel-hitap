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
    
    E  = vector<Compartment>(S.size());
    // breakthrough infections
    Ev_l = vector<Compartment>(S.size());
    Ev_m = vector<Compartment>(S.size());

    Ip = vector<Compartment>(S.size());
    Ia = vector<Compartment>(S.size());
    Is = vector<Compartment>(S.size());
    Ip_l = vector<Compartment>(S.size());
    Ia_l = vector<Compartment>(S.size());
    Is_l = vector<Compartment>(S.size());
    Ip_m = vector<Compartment>(S.size());
    Ia_m = vector<Compartment>(S.size());
    Is_m = vector<Compartment>(S.size());
    
    C  = vector<Compartment>(S.size());
    // 1dosed waned
    // Sw = vector<double>(S.size(), 0.);
    // vaccinated & R (1 and 2 dose)
    R  = vector<double>(S.size(), 0.);
    Rl = vector<double>(S.size(), 0.);
    Rm = vector<double>(S.size(), 0.);

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
        contag[a] = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip[a].Size() + P.pop[p].fIa[a] * Ia[a].Size() + P.pop[p].fIs[a] * Is[a].Size() + P.pop[p].fIp[a] * Ip_l[a].Size() + P.pop[p].fIa[a] * Ia_l[a].Size() + P.pop[p].fIs[a] * Is_l[a].Size() + P.pop[p].fIp[a] * Ip_m[a].Size() + P.pop[p].fIa_m[a] * Ia[a].Size() + P.pop[p].fIs_m[a] * Is[a].Size()) / N[a];
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
            // does ordering matter in this section?
            rep(t, p, a, riS)    = S[a];
            rep(t, p, a, riSv_l)  = Sv_l[a];
            rep(t, p, a, riSv_m)  = Sv_m[a];
            // rep(t, p, a, riSw) = Sw[a];
            
            rep(t, p, a, riE)    = E[a].Size();
            rep(t, p, a, riEv_l) = Ev_l[a].Size();
            rep(t, p, a, riEv_m) = Ev_m[a].Size();
            
            rep(t, p, a, riIp) = Ip[a].Size();
            rep(t, p, a, riIs) = Is[a].Size();
            rep(t, p, a, riIa) = Ia[a].Size();
            
            rep(t, p, a, riIp_l) = Ip_l[a].Size();
            rep(t, p, a, riIs_l) = Is_l[a].Size();
            rep(t, p, a, riIa_l) = Ia_l[a].Size();
            
            rep(t, p, a, riIp_m) = Ip_m[a].Size();
            rep(t, p, a, riIs_m) = Is_m[a].Size();
            rep(t, p, a, riIa_m) = Ia_m[a].Size();
            
            rep(t, p, a, riR)    = R[a];
            rep(t, p, a, riRv_l) = Rv_l[a];
            rep(t, p, a, riRv_m) = Rv_m[a];
            
            rep(t, p, a, rilambda)    = P.pop[p].u[a]*lambda[a];
            rep(t, p, a, rilambdav_l) = P.pop[p].uv_l[a]*lambda[a];
            rep(t, p, a, rilambdav_m) = P.pop[p].uv_m[a]*lambda[a];

            // User-specified processes
            for (size_t i=0;
                 i < P.processes.prevalence_states.size();
                 i++
            ) rep(t, p, a, rep.user_defined_offset + i) = pc[P.processes.prevalence_states[i]][a].Size();

        }

        // 1. Built-in states
        // Vaccination and waning of natural immunity and vaccine protection
        // when figuring out distro of second doses, only giving to those
        // having received 1st dose (and not subsequently infectious)
        // double primary_dose_pop = Sv[a]+Sw[a]+Rv[a]+Ev[a].Size();
        // when giving first doses, only to population that hasn't already received
        // first or second dose
        // N[a] - first_dose_pop - Sv2[a] - Rv2[a];
        
        double primary_dose_eligible = S[a] + E[a].Size() + Ip[a].Size() + Is[a].Size() + R[a]
        double booster_dose_eligible_l = Sv_l[a] + Ev_l.Size() + Ip_l[a].Size() + Is_l[a].Size() + Rv_l[a]
        double booster_dose_eligible_m = Sv_m[a] + Ev_m.Size() + Ip_m[a].Size() + Is_m[a].Size() + Rv_m[a]
        double booster_dose_eligible = booster_dose_eligible_l + booster_dose_eligible_m

        // initial vaccination campaign, primary doses
        
        // (2) S -> Sv_l
        double nS_Sv_l = min(S[a], num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * S[a] * 0.5 / primary_dose_eligible * P.time_step));
        S[a] -= nS_Sv_l;
        Sv_l[a] += nS_Sv_l;
        
        // (3) S -> Sv_m
        double nS_Sv_m = min(S[a], num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * S[a] * 0.5 / primary_dose_eligible * P.time_step));
        S[a] -= nS_Sv_m;
        Sv_m[a] += nS_Sv_m;
        
        // (46) R -> Rv_l
        double nR_Rv_l = min(R[a], num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * R[a] * 0.5 / primary_dose_eligible * P.time_step));
        R[a] -= nR_Rv_l;
        Rv_l[a] += nR_Rv_l;
        
        // (47) R -> Rv_m
        double nR_Rv_m = min(R[a], num(P.pop[p].v_p[a] * P.pop[p].ev_p[a] * R[a] * 0.5 / primary_dose_eligible * P.time_step));
        R[a] -= nR_Rv_m;
        Rv_m[a] += nR_Rv_m;
          
        // double nS_Sv = min(S[a], num(P.pop[p].v[a] * P.pop[p].ev[a] * S[a] / first_dose_eligible * P.time_step));
        // double nR_Rv = min(R[a], num(P.pop[p].v[a] * P.pop[p].ev[a] * R[a] / first_dose_eligible * P.time_step));

        // waning of vaccine induced protection
        // (5) Sv_m -> Sv_l
        double nSv_m_Sv_l = binomial(Sv_m[a], 1.0 - exp(-P.pop[p].wv_ml[a] * P.time_step))
        Sv_m[a] -= nSv_m_Sv_l;
        Sv_l[a] += nSv_m_Sv_l;
        
        // (44) Rv_m -> Rv_l
        double nRv_m_Rv_l = binomial(Rv_m[a], 1.0 - exp(-P.pop[p].wv_ml[a] * P.time_step))
        Rv_m[a] -= nRv_m_Rv_l;
        Rv_l[a] += nRv_m_Rv_l;
        
        // double nSv_Sw  = binomial(Sv[a], 1.0 - exp(-P.pop[p].wv[a] * P.time_step));
        
        // booster vaccination campaign, booster doses
        // (8) Sv_l -> Sv_m
        double nSv_l_Sv_m = min(Sv_l[a], num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * Sv_l[a] / booster_dose_eligible * P.time_step));
        Sv_l[a] -= nSv_l_Sv_m;
        Sv_m[a] += nSv_l_Sv_m;
        
        // (41) Rv_l -> Rv_m
        double nRv_l_Rv_m = min(Rv_l[a], num(P.pop[p].v_b[a] * P.pop[p].ev_b[a] * Rv_l[a] / booster_dose_eligible * P.time_step));
        Rv_l[a] -= nRv_l_Rv_m;
        Rv_m[a] += nRv_l_Rv_m;
        
        // double nSw_Sv2 = min(Sw[a], num(P.pop[p].v2[a] * P.pop[p].ev2[a] * Sw[a] / first_dose_pop * P.time_step));
        // double nRv_Rv2 = min(Rv[a], num(P.pop[p].v2[a] * P.pop[p].ev2[a] * Rv[a] / first_dose_pop * P.time_step));
        
        // is += - and -= the same?
        // S[a]   += -nS_Sv;
        // Sv[a]  += nS_Sv - nSv_Sv2;
        // Sv2[a] += nSv_Sv2 + nSw_Sv2;
        // Sw[a]  += -nSw_Sv2;

        // waning of infection induced protection
        // (29) R -> S
        double nR_S = binomial(R[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        R[a] -= nR_S;
        S[a] += nR_S;
        
        // (32) Rv_l -> Sv_l
        double nRv_l_Sv_l = binomial(Rv_l[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        Rv_l[a] -= nRv_l_Sv_l;
        Sv_l[a] += nRv_l_Sv_l;
        
        // (35) Rv_l -> Sv_l
        double nRv_m_Sv_m = binomial(Rv_m[a], 1.0 - exp(-P.pop[p].wn[a] * P.time_step));
        Rv_m[a] -= nRv_m_Sv_m;
        Sv_m[a] += nRv_m_Sv_m;
        
        // double nRv_RSv = binomial(Rv[a], 1.0 - exp(-(P.pop[p].wv[a] + P.pop[p].wn[a]) * P.time_step));
        // double nRv_R = binomial(nRv_RSv, P.pop[p].wn[a] + P.pop[p].wv[a] == 0 ? 0 : P.pop[p].wv[a] / (P.pop[p].wn[a] + P.pop[p].wv[a]));
        // double nRv_Sv = nRv_RSv - nRv_R;
        // Rv[a]  += nR_Rv - nRv_Rv2;
        // Rv2[a] += nRv_Rv2;

        // infection things
        // S -> E; Sv -> Ev; Sv2 -> Ev2; Sw -> E
        double nS_E = binomial(S[a], 1.0 - exp(-P.pop[p].u[a]*lambda[a] * P.time_step));
        S[a] -= nS_E;
        double nSw_E = binomial(Sw[a], 1.0 - exp(-P.pop[p].u[a]*lambda[a] * P.time_step));
        Sw[a] -= nSw_E;
        
        E[a].Add(P, Rand, nS_E + nSw_E, P.pop[p].dE);

        double nSv_Ev = binomial(Sv[a], 1.0 - exp(-P.pop[p].uv[a]*lambda[a] * P.time_step));
        Sv[a] -= nSv_Ev;
        Ev[a].Add(P, Rand, nSv_Ev, P.pop[p].dEv);

        double nSv2_Ev2 = binomial(Sv2[a], 1.0 - exp(-P.pop[p].uv2[a]*lambda[a] * P.time_step));
        Sv2[a] -= nSv2_Ev2;
        Ev2[a].Add(P, Rand, nSv2_Ev2, P.pop[p].dEv2);

        // E -> Ip/Ia
        double nE_Ipa = E[a].Mature();
        double nE_Ip = binomial(nE_Ipa, P.pop[p].y[a]);
        double nE_Ia = nE_Ipa - nE_Ip;

        double nEv_Ipa = Ev[a].Mature();
        //double nEv_Ip = binomial(nEv_Ipa*0.62, P.pop[p].yv[a]);
        //double nEv_Ia = nEv_Ipa*0.62 - nEv_Ip;
        //double nEv_Rv  = nEv_Ipa*0.38;
        
        double nEv_Ip = binomial(nEv_Ipa, P.pop[p].yv[a]);
        double nEv_Ia = nEv_Ipa - nEv_Ip;

        double nEv2_Ipa = Ev2[a].Mature();
        // double nEv2_Ip = binomial(nEv2_Ipa*0.5, P.pop[p].yv2[a]);
        // double nEv2_Ia = nEv2_Ipa*0.5 - nEv2_Ip;
        // double nEv2_Rv2  = nEv2_Ipa*0.5;
        
        double nEv2_Ip = binomial(nEv2_Ipa, P.pop[p].yv2[a]);
        double nEv2_Ia = nEv2_Ipa - nEv2_Ip;

        Ip[a].Add(P, Rand, nE_Ip + nEv_Ip + nEv2_Ip, P.pop[p].dIp);
        Ia[a].Add(P, Rand, nE_Ia + nEv_Ia + nEv2_Ia, P.pop[p].dIa);
        // Rv[a]  += nEv_Rv;
        // Rv2[a]  += nEv2_Rv2;

        // Ip -> Is -- also, true case onsets
        double nIp_Is = Ip[a].Mature();
        Is[a].Add(P, Rand, nIp_Is, P.pop[p].dIs);

        // Reported cases
        double n_to_report = binomial(nIp_Is, P.pop[p].rho[a]);
        C[a].Add(P, Rand, n_to_report, P.pop[p].dC);
        double n_reported = C[a].Mature();

        // N.B. assuming that infected vaccines do not become Rv

        // Is -> R
        double nIs_R = Is[a].Mature();
        R[a] += nIs_R;

        // Ia -> R
        double nIa_R = Ia[a].Mature();
        R[a] += nIa_R;

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
            double n_entering = 0.;
            switch (process.source_id)
            {
                case srcS:
                    n_entering = nS_E; break;
                case srcE:
                    n_entering = nE_Ipa; break;
                case srcEv:
                    n_entering = nEv_Ipa; break;
                case srcEv2:
                    n_entering = nEv2_Ipa; break;
                case srcEp:
                    n_entering = nE_Ip; break;
                case srcEvp:
                    n_entering = nEv_Ip; break;
                case srcEv2p:
                    n_entering = nEv2_Ip; break;
                case srcEa:
                    n_entering = nE_Ia; break;
                case srcEva:
                    n_entering = nEv_Ia; break;
                case srcEv2a:
                    n_entering = nEv2_Ia; break;
                case srcIp:
                    n_entering = nIp_Is; break;
                case srcIs:
                    n_entering = nIs_R; break;
                case srcIa:
                    n_entering = nIa_R; break;
                case srcI:
                    n_entering = nIs_R + nIa_R; break;
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
        rep(t, p, a, ricases) += nIp_Is;
        rep(t, p, a, ricases_reported) += n_reported;
        rep(t, p, a, risubclinical) += nE_Ia;

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
        double DS  = binomial(S [a],        death_prob);
        double DSv = binomial(Sv[a],        death_prob);
        double DSv2 = binomial(Sv2[a],        death_prob);
        double DSw = binomial(Sw[a],        death_prob);

        double DR  = binomial(R [a],        death_prob);
        double DRv = binomial(Rv[a],        death_prob);
        double DRv2 = binomial(Rv2[a],        death_prob);

        // Changes
        N[a] += B;
        S[a] += B;

        S[a]  -= DS;
        Sv[a] -= DSv;
        Sv2[a] -= DSv2;
        Sw[a] -= DSw;
        R[a]  -= DR;
        Rv[a] -= DRv;
        Rv2[a] -= DRv2;

        // NB - RemoveProb also takes care of removal, so no -= needed
        double DE  = E[a] .RemoveProb(P, Rand, death_prob);
        double DEv = Ev[a].RemoveProb(P, Rand, death_prob);
        double DEv2 = Ev2[a].RemoveProb(P, Rand, death_prob);
        double DIp = Ip[a].RemoveProb(P, Rand, death_prob);
        double DIa = Ia[a].RemoveProb(P, Rand, death_prob);
        double DIs = Is[a].RemoveProb(P, Rand, death_prob);


        N[a]  -=
            DS + DSv + DSv2 + DSw +
            DE + DEv + DEv2 +
            DIp + DIa + DIs +
            DR + DRv + DRv2;

        // Agings
        if (a != lambda.size() - 1)
        {
            double age_prob = 1.0 - exp(-P.pop[p].A[a] * P.time_step);
            double AS  = binomial(S [a],        age_prob);
            double ASv = binomial(Sv[a],        age_prob);
            double ASv2 = binomial(Sv2[a],        age_prob);
            double ASw = binomial(Sw[a],        age_prob);
            double AR  = binomial(R [a],        age_prob);
            double ARv = binomial(Rv[a],        age_prob);
            double ARv2 = binomial(Rv2[a],        age_prob);

            S[a]      -= AS;
            S[a + 1]  += AS;
            Sv[a]     -= ASv;
            Sv[a + 1] += ASv;
            Sv2[a]     -= ASv2;
            Sv2[a + 1] += ASv2;
            Sw[a]     -= ASw;
            Sw[a + 1] += ASw;
            R[a]      -= AR;
            R[a + 1]  += AR;
            Rv[a]     -= ARv;
            Rv[a + 1] += ARv;
            Rv2[a]     -= ARv2;
            Rv2[a + 1] += ARv2;

            double AE  = E[a] .MoveProb(E [a + 1], P, Rand, age_prob);
            double AEv = Ev[a].MoveProb(Ev[a + 1], P, Rand, age_prob);
            double AEv2 = Ev2[a].MoveProb(Ev2[a + 1], P, Rand, age_prob);
            double AIp = Ip[a].MoveProb(Ip[a + 1], P, Rand, age_prob);
            double AIa = Ia[a].MoveProb(Ia[a + 1], P, Rand, age_prob);
            double AIs = Is[a].MoveProb(Is[a + 1], P, Rand, age_prob);

            N[a]      -= AS + ASv + ASv2 + ASw + AE + AEv + AEv2 + AIp + AIa + AIs + AR + ARv + ARv2;
            N[a + 1]  += AS + ASv + ASv2 + ASw + AE + AEv + AEv2 + AIp + AIa + AIs + AR + ARv + ARv2;
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
    vecprint(R, "R");
    vecprint(Sv, "Sv");
    vecprint(Sv2, "Sv2");
    vecprint(Sw, "Sw");
    vecprint(Rv, "Rv");
    vecprint(Rv2, "Rv2");
    comprint(E, "E");
    comprint(Ev, "Ev");
    comprint(Ev2, "Ev2");
    comprint(Ip, "Ip");
    comprint(Ia, "Ia");
    comprint(Is, "Is");
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
