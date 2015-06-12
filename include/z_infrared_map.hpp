//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

// A set of functions for calculation of electric field / frequency maps
// as outlined in the papers of Skinner et al. See, for example, S.M. Gruenbaum,
// et al. J. Chem. Theory Comput. 2013 9 (7), 3109. These maps are necessary
// for calculation of vibrational spectra using the methodology  outline in the
// aforementioned papers.

// Abbreviations:
//   FTIR = fourier-transform infrared
//   SFG = sum-frequency generation

#include <armadillo>
#include "z_frequency_map.hpp"
#include "z_histogram.hpp"
#include "z_molecule_group.hpp"
#include "z_cx_tcf.hpp"

#ifndef _Z_INFRARED_MAP_HPP_
#define _Z_INFRARED_MAP_HPP_

class InfraredMap:
  public FrequencyMap {
 public:
  inline InfraredMap(const MoleculeGroup& group, const Chromophore chromophore,
                      const double timestep, const int correlation_length = 200)
      : FrequencyMap(group, chromophore, timestep, correlation_length) {
    mu_01_ = arma::zeros<arma::cube>(DIMS, steps_guess_, num_chromophores());
  }

 private:
  arma::cube mu_01_;

  inline void ResizeArrays() {
    omega_01_.resize(num_chromophores(), steps_guess_);
    mu_01_.resize(DIMS, steps_guess_, num_chromophores());
  }

  inline double MuOrAlphaProduct(const int chromophore, const int corr_start,
                                 const int i_corr) {
    return arma::dot(mu_01_.slice(chromophore).col(corr_start),
                     mu_01_.slice(chromophore).col(corr_start+i_corr));
  }

  // Calculates omega and mu/alpha using the frequency maps.
  void UseMapping(const int chromophore,
                  const arma::rowvec& box);
};
#endif
