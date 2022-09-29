#include <algorithm>
#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;
using namespace std;

// needed for weighted medians
static NumericVector dwrk, dwrk_x, dwrk_w;
static IntegerVector iwrk;

// for debugging purposes
void printMat(const NumericMatrix &mat) {
  for(int i = 0; i < mat.nrow(); i++) {
    for(int j = 0; j < mat(i, _).size(); j++) {
      Rcout << mat(i, j) << ", ";
    }
    Rcout << '\n';
  }
}

static NumericMatrix
  cmeans_setup(int nr_objects, int nr_centers, int dist_metric)
  {
    if(dist_metric == 1) {
      // Needed for weighted medians. (manhattan)
      dwrk_x = NumericVector(nr_objects);
      dwrk_w = NumericVector(nr_objects);
      dwrk = NumericVector(nr_objects);
      iwrk = IntegerVector(nr_objects);
    }
    return NumericMatrix(nr_objects, nr_centers);
  }

/* Update the dissimilarities/distances (between objects and prototypes) for a
 single object (i.e., a single row of the dissimilarity matrix. */
static void
  object_dissimilarities(const NumericMatrix& feature_mat, const NumericMatrix& centers,
                         const int nr_objects, const int nr_features, const int nr_centers,
                         const int dist_metric, const int object_nr, NumericMatrix& dist_mat, 
                         const LogicalMatrix & missing_vals)
  {
    double sum, v;
    for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
      sum = 0;
      v = 0;
      for(int feature_nr = 0; feature_nr < nr_features; feature_nr++) {
        if(!missing_vals(object_nr, feature_nr)) {
          // feature_entry - center_entry
          v = feature_mat(object_nr, feature_nr) - centers(center_nr, feature_nr);
        }
        // euclidean-distance
        if(dist_metric == 0)
          sum += v * v;
        // manhattan-distance
        else if(dist_metric == 1)
          sum += fabs(v);
        
      }
      // save distance from object to center in distance matrix
      dist_mat(object_nr, center_nr) = sum;
    }
  }

static void
  cmeans_dissimilarities(const NumericMatrix& feature_mat, const NumericMatrix& centers,
                         const int nr_objects, const int nr_features, 
                         const int nr_centers, const int dist_metric, 
                         NumericMatrix& dist_mat, const LogicalMatrix & missing_vals)
  {
    
    // loop over all objects
    for(int object_nr = 0; object_nr < nr_objects; object_nr++) {
      object_dissimilarities(feature_mat, centers, nr_objects, nr_features,
                             nr_centers, dist_metric, object_nr, dist_mat, missing_vals);
    }
  }

static void
  object_memberships(const NumericMatrix & dist_mat,
                     const int nr_centers, const NumericVector & fuzz,
                     const int object_nr, NumericMatrix & membership_mat)
  {
    double sum, v;
    
    int n_of_zeroes = 0;
    for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
      if(dist_mat(object_nr, center_nr) == 0)
        n_of_zeroes++;
    }
    if(n_of_zeroes > 0) {
      v = 1 / n_of_zeroes;
      for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
        // update membership_mat only with v and on positions which are 0 in dist_mat
        membership_mat(object_nr, center_nr) = 
          dist_mat(object_nr, center_nr) == 0 ? v: 0;
      }
    }
    else {
      /* Use the assumption that in general, pow() is more
       expensive than subscripting. */
      sum = 0;
      for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
        // check if individual fuzzifier values
        v = pow(dist_mat(object_nr, center_nr), -1.0/(fuzz[object_nr]-1.0));
        sum += v;
        membership_mat(object_nr, center_nr) = v;
      }
      sum = 1/sum;
      for(int center_nr = 0; center_nr < nr_centers; center_nr++)
        membership_mat(object_nr, center_nr) *= sum;
    }
    
    sum = 0;
    double epsilon = 0.00001;
    for(int center_nr = 0; center_nr < nr_centers; center_nr++)
      sum +=  membership_mat(object_nr, center_nr);
  }


static void
  cmeans_memberships(const NumericMatrix & dist_mat,
                     const int nr_objects, const int nr_centers,
                     const NumericVector & exponent, NumericMatrix & membership_mat)
  {
    // loop over all objects
    for(int object_nr = 0; object_nr < nr_objects; object_nr++) {
      object_memberships(dist_mat, nr_centers, exponent,
                         object_nr, membership_mat);
    }
  }

static double
  cmeans_error_fn(const NumericMatrix & membership_mat, const NumericMatrix & dist_mat,
                  NumericVector weight, const int nr_objects, 
                  const int nr_centers, const NumericVector & fuzz)
  {
    double sum;
    
    sum = 0;
    double currFuzzVal = fuzz[0];
    for(int object_nr = 0; object_nr < nr_objects; object_nr++) {
      currFuzzVal = fuzz[object_nr];  
      for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
        sum += weight[object_nr] * pow(membership_mat(object_nr, center_nr), currFuzzVal)
        * dist_mat(object_nr, center_nr);
      }
    }
    return(sum);
  }

static void
  rsort_with_index(NumericVector & feature_values, IntegerVector & iwrk, 
                   const int len)
  {
    // using lambda function to sort index-vector based on values in feature-vector
    std::stable_sort(iwrk.begin(), iwrk.end(), [&](int i, int j) {
      return feature_values[i] < feature_values[j];
    });
    
    // sort feature vector in ascending order
    feature_values.sort();
  }

static double
  cmeans_weighted_median(NumericVector & feature_values, 
                         NumericVector &  weight, int len)
  {
    double val, mval, cumsum_w, cumsum_w_x, marg;
    
    /* Sort feature_values. */
    for(int i = 0; i < len; i++)
      iwrk[i] = i;
    
    rsort_with_index(feature_values, iwrk, len);
    
    /* Permute weight using iwrk, and normalize. */
    double sum = 0;    
    for(int i = 0; i < len; i++) {
      dwrk[i] = weight[iwrk[i]];
      sum += dwrk[i];
    }
    for(int i = 0; i < len; i++) {
      weight[i] = dwrk[i] / sum;
    }
    
    cumsum_w = cumsum_w_x = 0;
    mval = R_PosInf;
    // first value is returned, if val never less than mval
    marg = feature_values[0];
    
    for(int i = 0; i < len; i++) {
      cumsum_w += weight[i];
      cumsum_w_x += weight[i] * feature_values[i];
      val = feature_values[i] * (cumsum_w - .5) - cumsum_w_x;
      if(val < mval) {
        marg = feature_values[i];
        mval = val;
      }
    }
    return(marg);
  }

static void 
  euclidean_prototypes(const NumericMatrix & feature_mat, const NumericMatrix & membership_mat,
                       const NumericVector & weight, const int nr_objects,
                       const int nr_features, const int nr_centers, const NumericVector & fuzz,
                       const int dist_metric, NumericMatrix & centers, const NumericVector & ratio_missing_vals, 
                       const LogicalMatrix & missing_vals, const double weight_missing) 
  {
    double currFuzz = fuzz[0];
    double k = weight_missing;
    
    for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
      for(int feature_nr = 0; feature_nr < nr_features; feature_nr++) {
        // reset centers
        centers(center_nr, feature_nr) = 0;
      }
      double sum = 0;
      for(int object_nr = 0; object_nr < nr_objects; object_nr++) {
        currFuzz = fuzz[object_nr]; 
        double v = weight[object_nr]* pow(membership_mat(object_nr, center_nr), currFuzz) * (1 - (1 - ratio_missing_vals(object_nr)) * k);
        // denominator: sum of all membership values ^ fuzz * weight
        sum += v;
        for(int feature_nr = 0; feature_nr < nr_features; feature_nr++) {
          // only allow calculation if current value is not missing
          if(!missing_vals(object_nr, feature_nr)) {
            // numerator: sum(denominator * feature values)
            centers(center_nr, feature_nr) += v * feature_mat(object_nr, feature_nr);
          }
        }
        
      }
      sum = 1.0/sum;
      for(int feature_nr = 0; feature_nr < nr_features; feature_nr++) {
        // new center: numerator / denominator
        centers(center_nr, feature_nr) *= sum;
      }
    }
  }



static void 
  manhattan_prototypes(const NumericMatrix & feature_mat, const NumericMatrix & membership_mat,
                       const NumericVector & weight, const int nr_objects,
                       const int nr_features, const int nr_centers, const NumericVector & fuzz,
                       const int dist_metric, NumericMatrix & centers)
  {
    double currFuzz = fuzz[0];
    
    for(int center_nr = 0; center_nr < nr_centers; center_nr++) {
      for(int feature_nr = 0; feature_nr < nr_features; feature_nr++) {
        for(int object_nr = 0; object_nr < nr_objects; object_nr++) {
          currFuzz = fuzz[object_nr]; 
          dwrk_x[object_nr] = feature_mat(object_nr, feature_nr);
          dwrk_w[object_nr] = weight[object_nr] * 
            pow(membership_mat(object_nr, center_nr), currFuzz);
        }
        centers(center_nr, feature_nr) =
          cmeans_weighted_median(dwrk_x, dwrk_w, nr_objects);
      }
    }
  }


static void
  cmeans_prototypes(const NumericMatrix & feature_mat, NumericMatrix & membership_mat,
                    const NumericVector & weight, const int nr_objects,
                    const int nr_features, const int nr_centers, const NumericVector & fuzz,
                    const int dist_metric, NumericMatrix & centers, const NumericVector & ratio_missing_vals,  
                    const LogicalMatrix & missing_vals, const double weight_missing)
  {
    if(dist_metric == 0) {
      /* Euclidean: weighted means. */
      euclidean_prototypes(feature_mat, membership_mat, weight, nr_objects,
                           nr_features, nr_centers, fuzz, dist_metric, centers, ratio_missing_vals, missing_vals, weight_missing);
    }
    else {
      /* Manhattan: weighted medians. */
      manhattan_prototypes(feature_mat, membership_mat, weight, nr_objects,
                           nr_features, nr_centers, fuzz, dist_metric, centers);
    }
  }

// [[Rcpp::export]]
void fill_missing_vals_and_ratio(const NumericMatrix & feature_mat, LogicalMatrix & missing_vals, 
                                 NumericVector & ratio_missing_vals, double missing_value = NA_REAL) {
  
  // declare function pointer for checking if value is missing value based on the parameter missing_value
  bool (*missing_val_check)(double val1, double val2);
  // assign comparator function to function pointer based on missing_value
  if(NumericVector::is_na(missing_value)) {
    missing_val_check = [](double val1, double val2){return NumericVector::is_na(val2);};
  } else {
    missing_val_check = [](double val1, double val2){return val1 == val2;};  
  }
  
  double ratio = 0;
  double nr_missing_values = 0;
  
  for(int i = 0; i < feature_mat.nrow(); i++) {
    for(int j = 0; j < feature_mat.ncol(); j++) {
      if((*missing_val_check)(missing_value, feature_mat(i,j))) {
        missing_vals(i, j) = true;
        nr_missing_values++;
      }
    }
    // calculate ratio and save in vector
    ratio = nr_missing_values/feature_mat.ncol();
    ratio_missing_vals(i) = ratio;
    
    // reset ratio and nr_missing_values for new object
    ratio = 0;
    nr_missing_values = 0;
  }
}




// attribute tells the compiler to make this function invisible in R
// returns ermin(fitness), because the conversion from R primitive datatype to C++ double doesn't allow changes directly
// [[Rcpp::export]]
double c_plusplus_means(const NumericMatrix & feature_mat, NumericMatrix & centers,
                        NumericVector & weight, NumericVector & fuzz, int dist_metric, int iter_max, double rel_tol,
                        int verbose, NumericMatrix & membership_mat, double ermin, IntegerVector & iter, double missing_value = NA_REAL, 
                        double weight_missing = 0) {

  // check for user interrupts
  Rcpp::checkUserInterrupt();

  int nr_objects = feature_mat.nrow();
  int nr_centers = centers.nrow();
  // for symmetric matrix
  int nr_features = feature_mat.ncol();
  
  LogicalMatrix missing_vals(nr_objects, nr_features);
  NumericVector ratio_missing_vals(nr_objects);
  fill_missing_vals_and_ratio(feature_mat, missing_vals, ratio_missing_vals, missing_value);

  double old_fitness, new_fitness;
  
  NumericMatrix dist_mat = cmeans_setup(nr_objects, nr_centers, dist_metric);
  
  // calculate distances based on dist_metric(euclidean/manhattan)
  cmeans_dissimilarities(feature_mat, centers, nr_objects, nr_features, 
                         nr_centers, dist_metric, dist_mat, missing_vals);
  
  // calculate memberships based on distances -> save in membership_mat
  cmeans_memberships(dist_mat, nr_objects, nr_centers, fuzz, membership_mat);
  
  // calculate fitness: J(c, m) -> see literature
  old_fitness = new_fitness = cmeans_error_fn(membership_mat, dist_mat, weight,
                                              nr_objects, nr_centers, fuzz);
  
  iter[0] = 0;
  
  // iterate iter_max-times to find cluster centers
  while(iter[0]++ < iter_max) {
    cmeans_prototypes(feature_mat, membership_mat, weight, nr_objects, nr_features,
                      nr_centers, fuzz, dist_metric, centers, ratio_missing_vals, missing_vals, weight_missing);
    cmeans_dissimilarities(feature_mat, centers, nr_objects, nr_features, 
                           nr_centers, dist_metric, dist_mat, missing_vals);
    cmeans_memberships(dist_mat, nr_objects, nr_centers, fuzz, membership_mat);

    new_fitness = cmeans_error_fn(membership_mat, dist_mat, weight, 
                                  nr_objects, nr_centers, fuzz);

    if(fabs(old_fitness - new_fitness) < rel_tol * (old_fitness + rel_tol)) {
      if(verbose)
        Rcout << "Iteration: " << iter << " converged with fitness: " 
              << new_fitness << "\n";
        ermin = new_fitness;
        break;
    }
    else {
      if(verbose) {
        ermin = cmeans_error_fn(membership_mat, dist_mat, weight, nr_objects,
                                nr_centers, fuzz);
        Rcout << "Iteration: " << iter << " with fitness: " << new_fitness << "\n";
      }
      old_fitness = new_fitness;
    }
  }
  
  ermin = new_fitness;
  return ermin;
}
