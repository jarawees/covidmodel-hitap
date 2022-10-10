// changes.cpp

#include "changes.h"
#include "parameters.h"
#include "helper.h"
#include <stdexcept>

#define PVector(x)   x
#define PDiscrete(x) x.weights

#define ParamCapture(param, get_vec) \
    else if (param_name == #param) { \
        if (pops.empty()) { \
            for (auto& pp : P.pop) { \
                param_ptr.push_back(&get_vec(pp.param)); \
                param_orig.push_back(get_vec(pp.param)); \
            } \
        } else { \
            for (auto& pi : pops) { \
                if (pi >= P.pop.size()) \
                    throw logic_error("Referring to nonexistent population in Change::Capture."); \
                param_ptr.push_back(&get_vec(P.pop[pi].param)); \
                param_orig.push_back(get_vec(P.pop[pi].param)); \
            } \
        } \
    }

// Construct a change impacting parameter pname in populations po of parameters P;
// apply value v with mode m at times t
Change::Change(Parameters& P, vector<unsigned int>& po, vector<unsigned int>& ru, string pname, 
    Mode m, vector<double>& t, vector<vector<double>>& v)
 : mode(m), times(t), values(v), param_name(pname), current(-1), pops(po), runs(ru)
{
    if (times.size() != values.size())
        throw logic_error("Change: times and values must be same length.");
}

// Capture parameters to change
void Change::Capture(Parameters& P)
{
    current = -1;
    if (false) {}
    ParamCapture(dE,      PDiscrete)
    ParamCapture(dEv_l,      PDiscrete)
    ParamCapture(dEv_m,      PDiscrete)
    ParamCapture(dEv_h,      PDiscrete)
    ParamCapture(dIp,     PDiscrete)
    ParamCapture(dIa,     PDiscrete)
    ParamCapture(dIs,     PDiscrete)
    ParamCapture(dIp_l,     PDiscrete)
    ParamCapture(dIa_l,     PDiscrete)
    ParamCapture(dIs_l,     PDiscrete)
    ParamCapture(dIp_m,     PDiscrete)
    ParamCapture(dIa_m,     PDiscrete)
    ParamCapture(dIs_m,     PDiscrete)
    ParamCapture(dIp_h,     PDiscrete)
    ParamCapture(dIa_h,     PDiscrete)
    ParamCapture(dIs_h,     PDiscrete)
    ParamCapture(dC,      PDiscrete)
    ParamCapture(contact, PVector)
    ParamCapture(contact_mult,    PVector)
    ParamCapture(contact_lowerto, PVector)
    ParamCapture(u,       PVector)
    ParamCapture(uv_l,       PVector)
    ParamCapture(uv_m,       PVector)
    ParamCapture(uv_h,       PVector)
    ParamCapture(ur,       PVector)
    ParamCapture(uvr_l,       PVector)
    ParamCapture(uvr_m,       PVector)
    ParamCapture(uvr_h,       PVector)
    ParamCapture(fIp,     PVector)
    ParamCapture(fIa,     PVector)
    ParamCapture(fIs,     PVector)
    ParamCapture(y,       PVector)
    ParamCapture(yv_l,       PVector)
    ParamCapture(yv_m,       PVector)
    ParamCapture(yv_h,       PVector)
    ParamCapture(omega,   PVector)
    ParamCapture(rho,     PVector)
    ParamCapture(tau,     PVector)
    // primary dose, initial vaccination campaign
    ParamCapture(v_p,       PVector)
    ParamCapture(ev_p,      PVector)
    // booster dose, boosting vaccination campaign
    ParamCapture(v_b,       PVector)
    ParamCapture(ev_b,      PVector)
    // transition probability from none to low by primary doses
    ParamCapture(v_p_2l,      PVector)
    // transition probability from none to medium by primary doses
    ParamCapture(v_p_2m,      PVector)
    // transition probability from low to medium by booster doses
    ParamCapture(v_b_l2m,      PVector)
    // waning of infection induced immunity
    ParamCapture(wn,      PVector)
    // waning of vaccine effects from medium to low
    ParamCapture(wv_m2l,      PVector)
    // waning of vaccine effects from high to medium
    ParamCapture(wv_h2m,      PVector)
    ParamCapture(A,       PVector)
    ParamCapture(B,       PVector)
    ParamCapture(D,       PVector)

    else if (param_name == "travel") {
        param_ptr.push_back(&P.travel.x);
        param_orig.push_back(P.travel.x);
    }
    else
        throw logic_error("Unsupported parameter change for " + param_name);
}

// Return true if parameters will change at time t
bool Change::Update(double t)
{
    bool refresh = false;
    while (current + 1 < times.size() && t >= times[current + 1])
    {
        ++current;
        refresh = true;
    }
    return refresh;
}

// Apply parameter changes
void Change::Apply(double t)
{
    auto op = [&](vector<double>* ptr) {
        if (ptr->size() == values[current].size() || values[current].empty())
        {
            if (!values[current].empty())
            {
                for (unsigned int j = 0; j < ptr->size(); ++j)
                {
                    switch (mode)
                    {
                        case Assign:    (*ptr)[j] = values[current][j]; break;
                        case Add:       (*ptr)[j] += values[current][j]; break;
                        case Multiply:  (*ptr)[j] *= values[current][j]; break;
                        case LowerTo:   (*ptr)[j] = min((*ptr)[j], values[current][j]); break;
                        case RaiseTo:   (*ptr)[j] = max((*ptr)[j], values[current][j]); break;
                    }
                }
            }
        }
        else
        {
            cout << current << "\n";
            cout << ptr->size() << "\n";
            cout << values[current].size() << "\n";
            throw logic_error("Change::Apply: change value has different size from parameter.");
        }
    };

    if (current >= 0)
        for (auto ptr : param_ptr)
            op(ptr);
}

// Reset all linked parameters to their original values
void Change::Reset()
{
    for (unsigned int i = 0; i < param_ptr.size(); ++i)
        *param_ptr[i] = param_orig[i];
}

// Capture parameters to change
void ChangeSet::Capture(Parameters& P)
{
    for (auto& c : ch)
        c.Capture(P);
}


// Apply any needed changes
void ChangeSet::Apply(Parameters& P, double t)
{
    bool refresh = false;

    // Update all changes
    for (auto& c : ch)
    {
        if (c.Update(t))
        {
            refresh = true;
        }
    }

    if (refresh)
    {
        for (auto& c : ch)
            c.Reset();
        for (auto& c : ch)
            c.Apply(t);
        for (auto& pp : P.pop) { // TODO -- slightly wasteful; Change could keep track if linked parameter sets need refreshing; then again in most cases probably needed
            pp.needs_recalc = true;
            pp.Recalculate();
        }
    }
}

