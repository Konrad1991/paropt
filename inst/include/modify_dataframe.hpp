#ifndef MODIFYDATAFRAME
#define MODIFYDATAFRAME

#include "header.hpp"

typedef std::vector<double> VD;
typedef std::vector<int> VI;
typedef std::vector<std::vector<double> > MD;
typedef std::vector<std::vector<int> > MI;
typedef std::vector<std::string> VS;
typedef Rcpp::DataFrame DF;

enum class IMPORT_PARAMETER {
  UNDEFINED,
  Error,
  SUCCESS
};

enum class TIME {
  Error,
  UNDEFINED,
  NOT_TIME_AS_NAME,
  CONTAINS_NA,
  SUCCESS
};

enum class EC1 {
  UNDEFINED,
  C_UgL, // Column_Upper greater lower
  C_LgU,
  R_UgL,
  R_LgU,
  SUCCESS
};

enum class EC2 {
  UNDEFINED,
  SUCCESS,
  Error
};

enum class EC3 {
  SUCCESS,
  Error
};


enum class IMPORT_STATES {
  UNDEFINED,
  Error,
  SUCCESS
};

void DF_to_MD (DF x, MD &r, VS &rs);

void remove_NA2 (MD &Res, MD &Inp, int NROW, int NCOL);

enum TIME CTC (VD T, std::string HEADS);

enum EC1 CHECK1 (MD L, MD U, int &LINE);

enum EC2 CHECK2 (VS L, VS U, int &LINE);

enum EC3 CHECK3 (MI T_L, MI T_U, VI L, VI U);

enum IMPORT_PARAMETER ip (DF lb, DF ub,
  VI &params_cut_idx_vec, VD &params_time_combi_vec, VD &param_combi_lb, VD &param_combi_ub, VS &header_parameter);

enum IMPORT_PARAMETER ip_start (DF Start,
    VI &params_cut_idx_vec, VD &params_time_combi_vec, VD &param_combi, VS &header_parameter);

enum IMPORT_STATES Import_states(DF Start,
VI &hs_cut_idx_vec,
VD &hs_time_combi_vec,
VD &hs_harvest_state_combi_vec,
      VS &headers);

#endif // MODIFYDATAFRAME
