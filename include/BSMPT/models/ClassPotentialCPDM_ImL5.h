#ifndef CLASSPOTENTIALCPDM_IML5_H
#define CLASSPOTENTIALCPDM_IML5_H

/*
 * ClassTemplate.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
  * @file
  */

#pragma once

#include <string>                               // for string
#include <vector>                               // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT{
namespace Models{

/**
 * @brief The Class_Template class
 * Template for implementing a new model
 */
class Class_CPDM_ImL5 : public Class_Potential_Origin
{
public:
  Class_CPDM_ImL5 ();
  virtual
  ~Class_CPDM_ImL5 ();


  double tCT_free= 0;
  double tCT_dep = 0;
  double m11sqrt,m22sqrt,mSsqrt,ReA,ImA,L1,L2,L3,L4,L5,L6,L7,L8;

  double CTm11sqrt,CTm22sqrt,CTmSsqrt,CTReA,CTImA,CTL1,CTL2,CTL3,CTL4,CTL5,CTL6,CTL7,CTL8,CTImL5;
  double dT1,dT2,dTCP,dTCB,dTs;




  void ReadAndSet(const std::string& linestr, std::vector<double>& par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double>& par) override;
  void set_CT_Pot_Par(const std::vector<double>& par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double>& v) const override;
  double VCounterSimplified(const std::vector<double>& v) const override;
  void Debugging(const std::vector<double>& input, std::vector<double>& output) const override;
};

}
}



#endif // CLASSPOTENTIALCPDM_IML5_H
