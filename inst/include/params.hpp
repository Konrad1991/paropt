#ifndef PARAMS
#define PARAMS
//*********************************************************************
//
//  file: params.h
//
//  contents:
//    - Global Params (declared as extern)
//
//**********************************************************************

// Precompile rule: This header file is only included once...KindOfGuard
// ================================================================= //
#include "header.hpp"

class ParamClass {
private:
  int m_no_spl_pts;
  // Number of sample points, i.e. number of measurement times
  std::vector<double> m_spl_pts_arr;
  // vector filled with parameterpoints at each sample point respectivly
  std::vector<double> m_time_vec;
  // respective time pts required for the interpolation
  std::vector<double> m_lb_arr;
  // vector filled with lower bounds of the parameterpoints at each sample point respectivly
  std::vector<double> m_ub_arr;
  // vector filled with upper bounds of the parameterpoints at each sample point respectivly
public:
  //
  ParamClass( // double = km; only solve system
    double t_par_val
  );
  ParamClass( // vector = vmax; only solve system
    int t_no_spl_pts,
    std::vector<double> t_time_vec,
    std::vector<double> t_par_vec
  );
  //
  ParamClass ( // upper and lower bound for km; optimize system
      double t_low_bound,
      double t_up_bound
  );
  ParamClass ( // upper and lower bound for vmax; optimize system
      int t_no_spl_pts,
      std::vector<double> t_time_vec,
      std::vector<double> t_low_bound,
      std::vector<double> t_up_bound
  );
  ParamClass ( // upper, lower bound as well as start values for km
    double t_par_val,
    double t_low_bound,
    double t_up_bound
  );
  ParamClass ( // upper, lower bound as well as start values for vmax
    int t_no_spl_pts,
    std::vector<double> t_time_vec,
    std::vector<double> t_par_vec,
    std::vector<double> t_low_bound,
    std::vector<double> t_up_bound
  );
  //
  ~ParamClass () {}
  //
  std::vector<double> get_time_points() const {return m_time_vec;}
  double get_time_points(int i) const {return m_time_vec[i];}
  int get_no_spl_pts() const {return m_no_spl_pts;}
  double get_spl_pts_arr(int m) const {return m_spl_pts_arr[m];}
  double get_low_bound(int m) const {return m_lb_arr[m];}
  double get_up_bound(int m) const {return m_ub_arr[m];}
    //
};


class ParamOrderClass {
  private:
    std::list<ParamClass> m_ParamList;

  public:
    ParamOrderClass (std::list<ParamClass> t_ParamList);
    ~ParamOrderClass () {}
    void c_count(int &t_col_no);
    void cut_idx(std::vector<int> &t_cut_idx_vec);
    void get_time_combi(std::vector<double> &t_time_combi_vec);
    void get_param_combi(std::vector<double> &t_param_combi_vec);
    void get_lb_combi(std::vector<double> &t_lb_combi_vec);
    void get_up_combi(std::vector<double> &t_ub_combi_vec);
};

#endif // PARAMS
