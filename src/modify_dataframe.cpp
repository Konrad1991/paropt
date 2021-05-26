#include "modify_dataframe.hpp"
#include "state.hpp"
#include "params.hpp"


#define NA std::nan("l")

typedef std::vector<double> VD;
typedef std::vector<int> VI;
typedef std::vector<std::vector<double> > MD;
typedef std::vector<std::vector<int> > MI;
typedef std::vector<std::string> VS;
typedef Rcpp::DataFrame DF;

/*
Dataframe to MD
*/
void DF_to_MD (DF x, MD &r, VS &rs) {

  /*
  Extract Names
  */
  int NCOL = x.size();
  Rcpp::CharacterVector NAMES(NCOL);
  NAMES = x.names();
  rs.resize(NCOL);
  for(int i = 0; i < NCOL; i++) {
    rs[i] = NAMES(i);
  }

  /*
  Extract Data
  */
  int NROW = x.nrows();
  Rcpp::NumericMatrix M(NROW, NCOL);
  r.resize(NCOL);
  for(int i = 0; i < NCOL; i++) {
    r[i].resize(NROW);
  }

  for(int i = 0; i < NCOL; i++) {
    Rcpp::NumericVector temp = x[i];
    Rcpp::NumericMatrix::Column col = M(Rcpp::_, i);
    col = temp;
  }

  for(int i = 0; i < NCOL; i++) {
    for(int j = 0; j < NROW; j++) {
      r[i][j] = M(j,i); //M(i , j);
    }
  }

}

/*
Remove NAs
*/
void remove_NA2 (MD &Res, MD &Inp, int NROW, int NCOL) {

  Res.resize(NCOL);
  VI R_min_NA(NCOL);

  for(int i = 0; i < NCOL; i++) {
    for(int j = 0; j < NROW; j++) {
      if (!std::isnan(Inp[i][j])) {
        R_min_NA[i] = R_min_NA[i] + 1;
      }
    }
  }

  for(int i = 0; i < NCOL; i++) {
    Res[i].resize(R_min_NA[i]);
  }

  int p = 0;
  for(int i = 0; i < NCOL; i++) {
    for(int j = 0; j < NROW; j++) {
      if (!std::isnan(Inp[i][j])) {
          Res[i][p] = Inp[i][j];
          p++;
      }
    }
    p = 0;
  }

}

/*
check time column
*/
enum TIME CTC (VD T, std::string HEADS) {

  enum TIME res = TIME::UNDEFINED;

  bool name = false;
  bool values_all_numeric = false;

  std::string time = HEADS;

  unsigned int counter = 0;
  if(time == "time") {
    name = true;
  } else {
		Rcpp::Rcerr << "Error: First column has to be the Name time" << std::endl;
  }

  for(size_t i = 0; i < T.size(); i++) {
    if(std::isnan(T[i])) {
      Rcpp::Rcerr << "Error: Time vector is not allowed to contain NAs" << std::endl;
      res = TIME::CONTAINS_NA;
  		break;
    }
    counter++;
  }

  values_all_numeric = counter == T.size() ? true : false;

  if(name == true && values_all_numeric == true) {
    res = TIME::SUCCESS;
  } else if(name == false && values_all_numeric == true) {
    res = TIME::NOT_TIME_AS_NAME;
  } else if(name == true && values_all_numeric == false) {
    res = TIME::CONTAINS_NA;
  } else if (name == false && values_all_numeric == false) {
    res = TIME::Error;
  }

  return res;
}

/*
Error checker Nr.1
*/
enum EC1 CHECK1 (MD L, MD U, int &LINE) {
  enum EC1 res = EC1::UNDEFINED;

  if(L.size() > U.size()) {
    res = EC1::C_LgU;
  } else if (U.size() > L.size()) {
    res = EC1::C_UgL;
  }

  if(res == EC1::C_UgL || res == EC1::C_LgU) {
    return res;
  }

  int NCOL = L.size();
  for(int i = 0; i < NCOL; i++) {

    if(L[i].size() > U[i].size()) {
      res = EC1::R_LgU;
      LINE = i;
    } else if (U[i].size() > L[i].size()) {
      res = EC1::R_UgL;
      LINE = i;
    }

  }

  return res;
}

/*
Error checker Nr.2
*/
enum EC2 CHECK2 (VS L, VS U, int &LINE) {

  enum EC2 res = EC2::UNDEFINED;

  for(unsigned int i = 0; i < L.size(); i++) {

    if (L[i] != U[i]) {
      LINE = i;
      res = EC2::Error;
      return res;
    }
  }

  return res;
}

/*
Error checker Nr.3
*/
enum EC3 CHECK3 (MI T_L, MI T_U, VI L, VI U) {

  enum EC3 res = EC3::Error;

  if(T_L.size() != T_U.size()) {
    Rcpp::stop("\nError: Different number of time_points between startvalues, lower bounds and upper bounds");
  }

  for(unsigned int i = 0; i < T_L.size(); i++) {

    if(T_L[i].size() != T_U[i].size()) {
      Rcpp::Rcerr << "In column:  " << i << std::endl;
      Rcpp::stop("\nError: Different number of time_points between lower bounds and upper bounds");
    }

  }

  for(unsigned int i = 0; i < L.size(); i++) {

    if(L[i] != U[i]) {
      Rcpp::Rcerr << "In column:  " << i << std::endl;
      Rcpp::stop("\nError: NA values not at the same position");
    }

  }

  res = EC3::SUCCESS;

  return res;
}


/*
Import parameter
*/
enum IMPORT_PARAMETER ip (DF lb, DF ub,
  VI &params_cut_idx_vec, VD &params_time_combi_vec, VD &param_combi_lb, VD &param_combi_ub, VS &header_parameter) {

  enum IMPORT_PARAMETER ret = IMPORT_PARAMETER::UNDEFINED;

  // extract names Ls, Us and data L and U
  MD L; MD U;
  VS Ls; VS Us;
  DF_to_MD(lb, L, Ls);
  DF_to_MD(ub, U, Us);

  header_parameter.resize(L.size());
  for(unsigned int i = 0; i < L.size(); i++) {
    header_parameter[i] = Ls[i];
  }

  // remove NA in data
  MD Lf; MD Uf;
  remove_NA2(Lf, L, lb.nrows(), lb.size());
  remove_NA2(Uf, U, ub.nrows(), ub.size());

  // check time column
  enum TIME TL = CTC(Lf[0], Ls[0]);
  enum TIME TU = CTC(Uf[0], Us[0]);
  if(TL != TIME::SUCCESS) {
    Rcpp::stop("\nError: Time column of lower bounds is not correct");
  } else if ( TU != TIME::SUCCESS) {
    Rcpp::stop("\nError: Time column of upper bounds is not correct");
  }

  // Error checks
  int line_Error = -1;
  enum EC1 e1 = CHECK1(Lf, Uf, line_Error);
  if(e1 == EC1::C_UgL) {
    Rcpp::stop("\nError: More columns in dataframe for upper bounds than in dataframes containing lower bounds");
  } else if (e1 == EC1::C_LgU) {
    Rcpp::stop("\nError: More columns in dataframe for lower bounds than in dataframes containing upper bounds");
  } else if (e1 == EC1::R_UgL) {
    Rcpp::Rcerr << "In column:  " << line_Error << std::endl;
    Rcpp::stop("\nError: More rows in dataframe for upper bounds than in dataframes containing lower bounds");
  } else if (e1 == EC1::R_LgU) {
    Rcpp::Rcerr << "In column:  " << line_Error << std::endl;
    Rcpp::stop("\nError: More rows in dataframe for lower bounds than in dataframes containing upper bounds");
  }

  line_Error = -1;
  enum EC2 e2 = CHECK2(Ls, Us, line_Error);
  if(e2 == EC2::Error) {
    Rcpp::Rcerr << "In column:  " << std::endl;
    Rcpp::stop("\nError: Headers are different!");
  }

  // extract time column
  VD T_L(lb.nrows()); VD T_U(ub.nrows());

  for(int i = 0; i < lb.nrows(); i++) {
    T_L[i] = Lf[0][i];
    T_U[i] = Uf[0][i];
  }

  // extract number of rows without NA and timepoints1
  VI NROW_L_MIN_NA(lb.size());
  VI NROW_U_MIN_NA(ub.size());
  MI TP_L(lb.size()); // TP = time points
  MI TP_U(ub.size());

  int COLS = lb.size();
  int ROWS = lb.nrows();

  for(int i = 0; i < COLS; i++) {
    for(int j = 0; j < ROWS; j++) {

      if(!std::isnan(L[i][j])) {
        NROW_L_MIN_NA[i] = NROW_L_MIN_NA[i] + 1;
        TP_L[i].push_back(j);
      }

      if(!std::isnan(U[i][j])) {
        NROW_U_MIN_NA[i] = NROW_U_MIN_NA[i] + 1;
        TP_U[i].push_back(j);
      }

    }
  }
  CHECK3(TP_L, TP_U, NROW_L_MIN_NA, NROW_U_MIN_NA);



  // fill vectors
  std::list<ParamClass> paramlist;
  for(int i = 1; i < COLS; i++) {
    std::vector<double> temp_lb(NROW_L_MIN_NA[i]);
    std::vector<double> temp_ub(NROW_L_MIN_NA[i]);
    std::vector<double> temp_time(NROW_L_MIN_NA[i]);
    for(int j = 0; j < NROW_L_MIN_NA[i]; j++) {
      temp_lb[j] = Lf[i][j];
      temp_ub[j] = Uf[i][j];
      temp_time[j] = L[0][TP_L[i][j]];
    }
    if(temp_time.size() == 1) {
      Rcpp::Rcerr << Ls[i] << ":" << "\t" << "Parameter is const" << std::endl;
      ParamClass temp(temp_lb[0], temp_ub[0]);
      paramlist.push_back(temp);
    } else {
      Rcpp::Rcerr << Ls[i] << ":" << "\t" << "Parameter is variabel" << std::endl;
      ParamClass temp(temp_time.size(), temp_time, temp_lb, temp_ub);
      paramlist.push_back(temp);
    }
  }

  ParamOrderClass ParamOrder(paramlist);
  ParamOrder.cut_idx(params_cut_idx_vec);
  ParamOrder.get_time_combi(params_time_combi_vec);
  ParamOrder.get_lb_combi(param_combi_lb);
  ParamOrder.get_up_combi(param_combi_ub);

  ret = IMPORT_PARAMETER::SUCCESS;
  return ret;
}


/*
Import start parameter
*/
enum IMPORT_PARAMETER ip_start (DF Start,
  VI &params_cut_idx_vec, VD &params_time_combi_vec, VD &param_combi, VS &header_parameter) {

  enum IMPORT_PARAMETER ret = IMPORT_PARAMETER::UNDEFINED;

  // extract names N and data D
  MD D;
  VS N;
  DF_to_MD(Start, D, N);

  header_parameter.resize(Start.size());
  for(unsigned int i = 0; i < N.size(); i++) {
    header_parameter[i] = N[i];
  }

  // remove NA in data
  MD Df;
  remove_NA2(Df, D, Start.nrows(), Start.size());

  // check time column
  enum TIME Enum_TIME = CTC(Df[0], N[0]);
  if(Enum_TIME != TIME::SUCCESS) {
    Rcpp::stop("\nError: Time column of parameters is not correct");
  }

  // extract time column
  VD T_C(Start.nrows());
  for(int i = 0; i < Start.nrows(); i++) {
    T_C[i] = Df[0][i];
  }

  // extract number of rows without NA and timepoints1
  VI NROW_MIN_NA(Start.size());
  MI TP(Start.size()); // TP = time points

  int COLS = Start.size();
  int ROWS = Start.nrows();

  for(int i = 0; i < COLS; i++) {
    for(int j = 0; j < ROWS; j++) {

      if(!std::isnan(D[i][j])) {
        NROW_MIN_NA[i] = NROW_MIN_NA[i] + 1;
        TP[i].push_back(j);
      }

    }
  }

  std::list<ParamClass> paramlist;
  for(int i = 1; i < COLS; i++) {
    std::vector<double> temp_param(NROW_MIN_NA[i]);
    std::vector<double> temp_time(NROW_MIN_NA[i]);
    for(int j = 0; j < NROW_MIN_NA[i]; j++) {
      temp_param[j] = Df[i][j];;
      temp_time[j] = D[0][TP[i][j]];
    }
    if(temp_time.size() == 1) {
      Rcpp::Rcerr << N[i] << ":" << "\t" << "Parameter is const" << std::endl;
      ParamClass temp(temp_param[0]);
      paramlist.push_back(temp);
    } else {
      Rcpp::Rcerr << N[i] << ":" << "\t" << "Parameter is variabel" << std::endl;
      ParamClass temp(temp_time.size(), temp_time, temp_param);
      paramlist.push_back(temp);
    }
  }

  ParamOrderClass ParamOrder(paramlist);
  ParamOrder.cut_idx(params_cut_idx_vec);
  ParamOrder.get_time_combi(params_time_combi_vec);
  ParamOrder.get_param_combi(param_combi);

  ret = IMPORT_PARAMETER::SUCCESS;
  return ret;
}


enum IMPORT_STATES Import_states(DF Start,
VI &hs_cut_idx_vec,
VD &hs_time_combi_vec,
VD &hs_harvest_state_combi_vec,
  VS &headers) {

  enum IMPORT_STATES ret = IMPORT_STATES::UNDEFINED;

  // extract names N and data D
  MD Df;
  VS N;
  DF_to_MD(Start, Df, N);

  headers.resize(N.size());
  for(unsigned int i = 0; i < N.size(); i++) {
    headers[i] = N[i];
  }

  int ncol = Start.size();
  //int nrwo = Start.nrows();

  for(int i = 1; i < (ncol); i++) {
    if(std::isnan(Df[i][0])) {
      Rcpp::stop("\nError: first line of states contains NA");
    }
  }

  // check time column
  enum TIME Enum_TIME = CTC(Df[0], N[0]);
  if(Enum_TIME != TIME::SUCCESS) {
    Rcpp::stop("\nError: Time column of states is not correct");
  }

  // extract time vector
  std::vector<double> time(Df[0].size());

  for(size_t i = 0; i < time.size(); i++) {
    time[i] = Df[0][i];
  }

  int size_of_states = time.size();
  // define combination of state and time and store it in std::list
  std::list<HarvestStateClass> paramlist;
  for(int i = 1; i < ncol; i++) {
    std::vector<double> temp_state(size_of_states);
    std::vector<double> temp_time(size_of_states);
    for(int j = 0; j < size_of_states; j++) {
      temp_state[j] = Df[i][j];
      temp_time[j] = Df[0][j];
    }
      Rcpp::Rcerr << headers[i] << std::endl;
      HarvestStateClass temp(temp_time, temp_state);
      paramlist.push_back(temp);
  }

  HarvestStateOrderClass HarvestOrder(paramlist);
  HarvestOrder.cut_idx(hs_cut_idx_vec);
  HarvestOrder.get_harvest_time_combi(hs_time_combi_vec);
  HarvestOrder.get_harvest_state_combi(hs_harvest_state_combi_vec);

  ret = IMPORT_STATES::SUCCESS;

  return ret;
}
