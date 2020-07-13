/********************************************************************
//
//  file: state.cpp
//
//  contents:
//    - Definition of ODE state_t class
//
*********************************************************************/

#include "state.hpp"


// HarvestStateClass
// ========================================================================
// ========================================================================
// ========================================================================
HarvestStateClass::HarvestStateClass (
  std::vector<double> t_harvest_time, std::vector<double> t_harvest_state
){
  m_harvest_time.resize(t_harvest_time.size());
  m_harvest_state.resize(t_harvest_state.size());
  m_harvest_time = t_harvest_time;
  m_harvest_state = t_harvest_state;
}
// ================================================================= //
HarvestStateClass::HarvestStateClass (double t_harvest_state){
  m_harvest_time.resize(1);
  m_harvest_state.resize(1);
  m_harvest_time[0] = 0.;
  m_harvest_state[0] = t_harvest_state;
}
// ================================================================= //
HarvestStateClass::~HarvestStateClass(){
  //std::cout << "HarvestStateClass Destructor ~HarvestStateClass" << '\n';
}
// ================================================================= //
std::vector<double> HarvestStateClass::getHarvestState (void) const {
  return m_harvest_state;
}
// ================================================================= //
double HarvestStateClass::getHarvestState (int i) const {
  return m_harvest_state[i];
}
// ================================================================= //
std::vector<double> HarvestStateClass::getHarvestTime (void) const {
  return m_harvest_time;
}
// ================================================================= //
double HarvestStateClass::getHarvestTime (int i) const {
  return m_harvest_time[i];
}
// ================================================================= //
int HarvestStateClass::getHarvestStateLength (void) const {
  return m_harvest_state.size();
}

// HarvestStateOrderClass
// ========================================================================
// ========================================================================
// ========================================================================
HarvestStateOrderClass::HarvestStateOrderClass(std::list<HarvestStateClass> t_HarvestStateList){
  m_HarvestStateList = t_HarvestStateList;
}
// ================================================================= //
HarvestStateOrderClass::~HarvestStateOrderClass(){
  //std::cout << "HarvestStateOrderClass Destructor ~HarvestStateOrderClass" << '\n';
}
// ================================================================= //
int HarvestStateOrderClass::lengthcount (void) {
  std::list<HarvestStateClass>::iterator it;
  int t_length_no = 0;
  for (it = m_HarvestStateList.begin(); it != m_HarvestStateList.end(); ++it){
    t_length_no += it->getHarvestStateLength();
  }
  return t_length_no;
}
// ================================================================= //
void HarvestStateOrderClass::cut_idx(
  std::vector<int> &t_cut_idx_vec
){
  std::list<HarvestStateClass>::iterator it;
  int tmpCount = 0;
  for (it = m_HarvestStateList.begin(); it != m_HarvestStateList.end(); ++it){
    ++tmpCount;
  }
  t_cut_idx_vec.resize(tmpCount);
  tmpCount=0;
  for (it = m_HarvestStateList.begin(); it != m_HarvestStateList.end(); ++it){
    t_cut_idx_vec[tmpCount] = it->getHarvestStateLength();
    ++tmpCount;
  }
}
// ================================================================= //
void HarvestStateOrderClass::get_harvest_time_combi(
  std::vector<double> &t_time_combi_vec
){
  std::list<HarvestStateClass>::iterator it;
  int lengthcount = HarvestStateOrderClass::lengthcount();
  //
  t_time_combi_vec.resize(lengthcount);
  //
  int Cc=0;
  //
  for (it = m_HarvestStateList.begin(); it != m_HarvestStateList.end(); ++it){
    for (int iC = 0; iC < it->getHarvestStateLength(); ++iC) {
      t_time_combi_vec[Cc] = it->getHarvestTime(iC);
      Cc+=1;
    }
  }
}
// ================================================================= //
void HarvestStateOrderClass::get_harvest_state_combi(
  std::vector<double> &t_harvest_state_combi_vec
){
  std::list<HarvestStateClass>::iterator it;
  //int lengthcount = HarvestStateOrderClass::lengthcount();
  //
  t_harvest_state_combi_vec.resize(HarvestStateOrderClass::lengthcount());
  //
  int Cc=0;
  //
  for (it = m_HarvestStateList.begin(); it != m_HarvestStateList.end(); ++it){
    for (int iC = 0; iC < it->getHarvestStateLength(); ++iC) {
      t_harvest_state_combi_vec[Cc] = it->getHarvestState(iC);
      Cc+=1;
    }
  }
}
