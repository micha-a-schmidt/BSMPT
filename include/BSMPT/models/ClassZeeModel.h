// Copyright (C) 2018  Philipp Basler and Margarete MA~\274hlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete MA~\274hlleitner and Jonas MA~\274ller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#pragma once

#include <string> // for string
#include <vector> // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_ZeeModel class
 * Template for implementing a new model
 */
class Class_ZeeModel : public Class_Potential_Origin
{
public:
  Class_ZeeModel();
  virtual ~Class_ZeeModel();

  // Add here your parameters for the Lagrangian as well as for the counterterm
  // potential Add here your variables in which you will save the Debye
  // correction factors

double mu1s, mu2s, mu3s, muhs, lh, l1, l2, l3, l4, l5r, l5i, l8, l9, l6r, l6i, l7r, l7i, l10r, l10i, mu; 
double dT1ar, dT1br, dT2ar, dT2br, dTr, dT1ai, dT1bi, dT2ai, dT2bi, dTi; 
double yt, g, gp;
 
  // double ms, lambda;

  // double dms, dlambda, dT, yt, g;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
