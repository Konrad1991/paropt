//*********************************************************************
//
//  file: state.h
//
//  contents:
//    - Definition of ODE state_t class
//
//**********************************************************************

// Precompile rule: This header file is only included once...KindOfGuard

#include "header.hpp"
// ================================================================= //

class HarvestStateClass {
private:
  std::vector<double> m_harvest_state;
  std::vector<double> m_harvest_time;

public:
  HarvestStateClass (std::vector<double> t_harvest_time, std::vector<double> t_harvest_state);
  HarvestStateClass (double t_harvest_state);
  ~HarvestStateClass ();

  std::vector<double> getHarvestState (void) const;
  double getHarvestState (int i) const;
  std::vector<double> getHarvestTime () const;
  double getHarvestTime (int i) const;
  int getHarvestStateLength () const;
};

/* ================================================================= //

class InitStatesClass {
private:
  std::vector<double> m_init_states;

public:
  InitStateClass (std::vector<double> t_init_states);
  ~InitStateClass ();

  int getNoStates ();
  std::vector<double> getInitStates();
};

// ================================================================= */

class HarvestStateOrderClass {
private:
  std::list<HarvestStateClass> m_HarvestStateList;

public:
  HarvestStateOrderClass(std::list<HarvestStateClass> t_HarvestStateList);
  ~HarvestStateOrderClass ();

  int lengthcount ();
  void cut_idx(std::vector<int> &t_cut_idx_vec);
  void get_harvest_time_combi(std::vector<double> &t_time_combi_vec);
  void get_harvest_state_combi(std::vector<double> &t_harvest_state_combi_vec);

};

// ================================================================= //
