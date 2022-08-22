//
// PROCESS
//

#ifndef PROCESS_SPEC_H
#define PROCESS_SPEC_H

#include "helper.h"
#include <algorithm>
#include <vector>

// why do we not need all compartments here?
enum SourceID
{
    srcS = 1000000,
    srcR,
    srcSv_l,
    srcRv_l,
    srcSv_m,
    srcRv_m,
    srcE,
    srcE_Ip,
    srcE_Ia,
    srcEv_l,
    srcEv_l_Ip,
    srcEv_l_Ia,
    srcEv_m,
    srcEv_m_Ip,
    srcEv_m_Ia,
    srcIp_Is,
    srcIp_l_Is_l,
    srcIp_m_Is_m,
    srcIa_R,
    srcIa_l_Rv_l,
    srcIa_m_Rv_m,
    srcIs_R,
    srcIs_l_Rv_l,
    srcIs_m_Rv_m,
    srcCasesReported,
    
};

const map<string,SourceID> processSourceMap = {
    {"S", srcS},
    {"R", srcR},
    {"Sv_l", srcSv_l},
    {"Rv_l", srcRv_l},
    {"Sv_m", srcSv_m},
    {"Rv_m", srcRv_m},
    {"E", srcE},
    {"E_Ip", srcE_Ip},
    {"E_Ia", srcE_Ia},
    {"Ev_l", srcEv_l},
    {"Ev_l_Ip", srcEv_l_Ip},
    {"Ev_l_Ia", srcEv_l_Ia},
    {"Ev_m", srcEv_l},
    {"Ev_m_Ip", srcEv_m_Ip},
    {"Ev_m_Ia", srcEv_m_Ia},
    {"Ip_Is",srcIp_Is},
    {"Ip_l_Is_l",srcIp_l_Is_l},
    {"Ip_m_Is_m",srcIp_m_Is_m},
    {"Ia_R",srcIa_R},
    {"Ia_l_Rv_l",srcIa_l_Rv_l},
    {"Ia_m_Rv_m",srcIa_m_Rv_m},
    {"Is_R",srcIs_R},
    {"Is_l_Rv_l",srcIs_l_Rv_l},
    {"Is_m_Rv_m",srcIs_m_Rv_m},
    {"cases_reported", srcCasesReported}
};

const unsigned int Null = 999999;
class Discrete;

// plain-old object container for Processes
// when initially constructed, source and sink ids
// not specified; should only be added when building a
// ProcessList
struct ProcessSpec
{
    // name of the source; must be either in processSourceMap keys
    // *or* the names of some other process
    // used to identify appropriate source_id
    string source_name;
    
    // exits from compartment source_id enter this process
    // source_id must be in SourceID *or* some other process sink_ids
    unsigned int source_id;
    
    // name of the sinks; cannot be in processSourceMap keys
    // used to identify appropriate sink_ids
    // special case: the last name may be "null"
    vector<string> names;
    
    // the ids for end points
    vector<unsigned int> sink_ids;
    
    string type;            // ignored for now - multinomial or dirichlet multinomial

    // reporting mode of sub-processes: empty, "i", "o", "p", or a combination of these
    vector<string> report;
    
    // probability by group of entering each sub-process from the source above;
    // indexed by group then by subprocess
    // e.g. prob[a] { p outcome 1, p 2, ..., p n} for age a 
    vector<vector<double>> prob;
    vector<Discrete> delays;    // delays for each sub-process

    // 
    // 
    // vector<unsigned int> p_cols;
    // vector<unsigned int> p_ids;
    // vector<unsigned int> i_cols;
    // vector<unsigned int> i_ids;
    // vector<unsigned int> o_cols;
    // vector<unsigned int> o_ids;     // prevalence, incidence, and outcidence data columns and sub-process identifiers for this process
};

class ProcessList {
public:
    vector<ProcessSpec> flows;
    // the size of user process containers to use
    unsigned int state_count = 0;
    vector<string> state_names;
    
    // the sink ids (i.e. indices in the process containers)
    // corresponding to prevalence, incidence, and outcidence for recording
    
    // when doing column numbered reporting
    // prevalence columns = base_col_offset + 0:prevalence_states.size()-1
    // incidence column  = base_col_offset + prevalence_states.size() + ...
    // incidence column  = base_col_offset + prevalence_states.size() + incidince_states.size()
    vector<unsigned int> prevalence_states;
    vector<unsigned int> incidence_states;
    vector<unsigned int> outcidence_states;
    unsigned int recording_count = 0, inc_offset = 0, out_offset = 0;
    
    void Update(vector<ProcessSpec>& ps);
    
};

#endif
