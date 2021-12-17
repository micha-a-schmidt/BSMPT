// Copyright (C) 2018  Philipp Basler and Margarete MA~\274hlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete MA~\274hlleitner and Jonas
// MA~\274ller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include <algorithm>              // for max, copy
#include <iostream>               // for operator<<, endl, basic_o...
#include <memory>                 // for allocator_traits<>::value...
#include <stddef.h>               // for std::size_t

#include <BSMPT/models/ClassZeeModel.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

namespace BSMPT
{
namespace Models
{
/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of
 * Lagrangian parameters AFTER using the tadpole conditions), nParCT (number of
 * counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_ZeeModel::Class_ZeeModel()
{
  Model =
      ModelID::ModelIDs::TEMPLATE; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 10;               // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 20; // number of parameters in the tree-Level Lagrangian
  nParCT = 10; // number of parameters in the counterterm potential

  nVEV = 10; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
VevOrder[0]=0;
VevOrder[1]=1;
VevOrder[2]=2;
VevOrder[3]=3;
VevOrder[4]=4;
VevOrder[5]=5;
VevOrder[6]=6;
VevOrder[7]=7;
VevOrder[8]=8;
VevOrder[9]=9;


  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_ZeeModel::~Class_ZeeModel()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_ZeeModel::addLegendCT() const
{
  std::vector<std::string> labels;
labels.push_back("dT1ar");
labels.push_back("dT1br");
labels.push_back("dT2ar");
labels.push_back("dT2br");
labels.push_back("dTr");
labels.push_back("dT1ai");
labels.push_back("dT1bi");
labels.push_back("dT2ai");
labels.push_back("dT2bi");
labels.push_back("dTi");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_ZeeModel::addLegendTemp() const
{
//UPDATE LATER
  std::vector<std::string> labels;
  labels.push_back("T_c"); // Label for the critical temperature
  labels.push_back("v_c"); // Label for the critical vev
  labels.push_back("v_c/T_c"); // Label for v_c/T_c, you could use xi_c also for example
  // out += "Your VEV order"; // Now you have to put the label for your vevs
  labels.push_back("omega");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_ZeeModel::addLegendTripleCouplings() const
{
// DO NOT IMPLEMENT
   std::vector<std::string> labels;
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_ZeeModel::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
labels.push_back("v1ar");
labels.push_back("v1br");
labels.push_back("v2ar");
labels.push_back("v2br");
labels.push_back("vr");
labels.push_back("v1ai");
labels.push_back("v1bi");
labels.push_back("v2ai");
labels.push_back("v2bi");
labels.push_back("vi");

  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_ZeeModel::ReadAndSet(const std::string &linestr,
                                std::vector<double> &par)
{

  std::stringstream ss(linestr);
  double tmp;


  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k=0; k<20; k++)
  {
    ss >> par[k];
  }
  set_gen(par); // This you have to call so that everything will be set
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_ZeeModel::set_gen(const std::vector<double> &par)
{
// UPDATE CHECK
mu1s = par[0];
mu2s = par[1];
mu3s = par[2];
muhs = par[3];
lh = par[4];
l1 = par[5];
l2 = par[6];
l3 = par[7];
l4 = par[8];
l5r = par[9];
l5i = par[10];
l8 = par[11];
l9 = par[12];
l6r = par[13];
l6i = par[14];
l7r = par[15];
l7i = par[16];
l10r = par[17];
l10i = par[18];
mu = par[19];

//  ms     = par[0]; // Class member is set accordingly to the input parameters
//  lambda = par[1]; // Class member is set accordingly to the input parameters
  g      = C_g;    // SM SU (2) gauge coupling --> SMparam .h
  yt = std::sqrt(2) / C_vev0 * C_MassTop; // Top Yukawa coupling --> SMparam .h
  scale = C_vev0; // Renormalisation scale is set to the SM VEV
  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // Here you have to set the vector vevTreeMin. The vector vevTree will then be
  // set by the function MinimizeOrderVEV
  vevTreeMin[0] = C_vev0;
  vevTree       = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_ZeeModel::set_CT_Pot_Par(const std::vector<double> &par)
{
Curvature_Higgs_CT_L1[0]=par[0];
Curvature_Higgs_CT_L1[1]=par[1];
Curvature_Higgs_CT_L1[2]=par[2];
Curvature_Higgs_CT_L1[3]=par[3];
Curvature_Higgs_CT_L1[4]=par[4];
Curvature_Higgs_CT_L1[5]=par[5];
Curvature_Higgs_CT_L1[6]=par[6];
Curvature_Higgs_CT_L1[7]=par[7];
Curvature_Higgs_CT_L1[8]=par[8];
Curvature_Higgs_CT_L1[9]=par[9];

}

/**
 * console output of all Parameters
 */
void Class_ZeeModel::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
ss << "mu1s = " << mu1s << std::endl << "mu2s = " << mu2s << std::endl << "mu3s = " << mu3s << std::endl << "muhs = " << muhs << std::endl << "lh = " << lh << std::endl << "l1 = " << l1 << std::endl << "l2 = " << l2 << std::endl << "l3 = " << l3 << std::endl << "l4 = " << l4 << std::endl << "l5r = " << l5r << std::endl << "l5i = " << l5i << std::endl << "l8 = " << l8 << std::endl << "l9 = " << l9 << std::endl << "l6r = " << l6r << std::endl << "l6i = " << l6i << std::endl << "l7r = " << l7r << std::endl << "l7i = " << l7i << std::endl << "l10r = " << l10r << std::endl << "l10i = " << l10i << std::endl << "mu = " << mu << std::endl  ;
//  ss << "lambda = " << lambda << std::endl << "m^2 = " << ms << std::endl;

  ss << "The counterterm parameters are : " << std::endl;
ss << "dT1ar = " << dT1ar << std::endl << "dT1br = " << dT1br << std::endl << "dT2ar = " << dT2ar << std::endl << "dT2br = " << dT2br << std::endl << "dTr = " << dTr << std::endl << "dT1ai = " << dT1ai << std::endl << "dT1bi = " << dT1bi << std::endl << "dT2ai = " << dT2ai << std::endl << "dT2bi = " << dT2bi << std::endl << "dTi = " << dTi << std::endl  ;
//  ss << "dT = " << dT << std::endl
//     << "dlambda = " << dlambda << std::endl
//     << "dm^2 = " << dms << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_ZeeModel::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // Here you have to use your formulae for the counterterm scheme
parCT.push_back(-NablaWeinberg(0));
parCT.push_back(-NablaWeinberg(1));
parCT.push_back(-NablaWeinberg(2));
parCT.push_back(-NablaWeinberg(3));
parCT.push_back(-NablaWeinberg(4));
parCT.push_back(-NablaWeinberg(5));
parCT.push_back(-NablaWeinberg(6));
parCT.push_back(-NablaWeinberg(7));
parCT.push_back(-NablaWeinberg(8));
parCT.push_back(-NablaWeinberg(9));

  // double t = 0;
  // parCT.push_back(t); // dT
  // parCT.push_back(3.0 * t / std::pow(C_vev0, 3) +
  //                 3.0 / std::pow(C_vev0, 3) * NablaWeinberg(0) -
  //                 3.0 / std::pow(C_vev0, 2) * HesseWeinberg(0, 0)); // dlambda
  // parCT.push_back(-3.0 / (2 * C_vev0) * NablaWeinberg(0) +
  //                 1.0 / 2.0 * HesseWeinberg(0, 0) -
  //                 3.0 * t / (2 * C_vev0)); // dms

  return parCT;
}

void Class_ZeeModel::TripleHiggsCouplings()
{
// DO NOT IMPLEMENT...NOT NECESSARY FOR US
  }

void Class_ZeeModel::SetCurvatureArrays()
{
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

Curvature_Higgs_L2[0][0]=mu1s;
Curvature_Higgs_L2[0][2]=mu3s;
Curvature_Higgs_L2[1][1]=mu1s;
Curvature_Higgs_L2[1][3]=mu3s;
Curvature_Higgs_L2[2][0]=mu3s;
Curvature_Higgs_L2[2][2]=mu2s;
Curvature_Higgs_L2[3][1]=mu3s;
Curvature_Higgs_L2[3][3]=mu2s;
Curvature_Higgs_L2[4][4]=muhs;
Curvature_Higgs_L2[5][5]=mu1s;
Curvature_Higgs_L2[5][7]=mu3s;
Curvature_Higgs_L2[6][6]=mu1s;
Curvature_Higgs_L2[6][8]=mu3s;
Curvature_Higgs_L2[7][5]=mu3s;
Curvature_Higgs_L2[7][7]=mu2s;
Curvature_Higgs_L2[8][6]=mu3s;
Curvature_Higgs_L2[8][8]=mu2s;
Curvature_Higgs_L2[9][9]=muhs;
Curvature_Higgs_L3[0][3][4]=mu/std::sqrt(2);
Curvature_Higgs_L3[0][4][3]=mu/std::sqrt(2);
Curvature_Higgs_L3[0][8][9]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[0][9][8]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[1][2][4]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[1][4][2]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[1][7][9]=mu/std::sqrt(2);
Curvature_Higgs_L3[1][9][7]=mu/std::sqrt(2);
Curvature_Higgs_L3[2][1][4]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[2][4][1]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[2][6][9]=mu/std::sqrt(2);
Curvature_Higgs_L3[2][9][6]=mu/std::sqrt(2);
Curvature_Higgs_L3[3][0][4]=mu/std::sqrt(2);
Curvature_Higgs_L3[3][4][0]=mu/std::sqrt(2);
Curvature_Higgs_L3[3][5][9]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[3][9][5]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[4][0][3]=mu/std::sqrt(2);
Curvature_Higgs_L3[4][1][2]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[4][2][1]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[4][3][0]=mu/std::sqrt(2);
Curvature_Higgs_L3[4][5][8]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[4][6][7]=mu/std::sqrt(2);
Curvature_Higgs_L3[4][7][6]=mu/std::sqrt(2);
Curvature_Higgs_L3[4][8][5]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[5][3][9]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[5][4][8]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[5][8][4]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[5][9][3]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[6][2][9]=mu/std::sqrt(2);
Curvature_Higgs_L3[6][4][7]=mu/std::sqrt(2);
Curvature_Higgs_L3[6][7][4]=mu/std::sqrt(2);
Curvature_Higgs_L3[6][9][2]=mu/std::sqrt(2);
Curvature_Higgs_L3[7][1][9]=mu/std::sqrt(2);
Curvature_Higgs_L3[7][4][6]=mu/std::sqrt(2);
Curvature_Higgs_L3[7][6][4]=mu/std::sqrt(2);
Curvature_Higgs_L3[7][9][1]=mu/std::sqrt(2);
Curvature_Higgs_L3[8][0][9]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[8][4][5]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[8][5][4]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[8][9][0]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[9][0][8]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[9][1][7]=mu/std::sqrt(2);
Curvature_Higgs_L3[9][2][6]=mu/std::sqrt(2);
Curvature_Higgs_L3[9][3][5]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[9][5][3]=-(mu/std::sqrt(2));
Curvature_Higgs_L3[9][6][2]=mu/std::sqrt(2);
Curvature_Higgs_L3[9][7][1]=mu/std::sqrt(2);
Curvature_Higgs_L3[9][8][0]=-(mu/std::sqrt(2));
Curvature_Higgs_L4[0][0][0][0]=3*l1;
Curvature_Higgs_L4[0][0][0][2]=3*l6r;
Curvature_Higgs_L4[0][0][0][7]=-3*l6i;
Curvature_Higgs_L4[0][0][1][1]=l1;
Curvature_Higgs_L4[0][0][1][3]=l6r;
Curvature_Higgs_L4[0][0][1][8]=-l6i;
Curvature_Higgs_L4[0][0][2][0]=3*l6r;
Curvature_Higgs_L4[0][0][2][2]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[0][0][2][5]=l6i;
Curvature_Higgs_L4[0][0][2][7]=-l5i;
Curvature_Higgs_L4[0][0][3][1]=l6r;
Curvature_Higgs_L4[0][0][3][3]=l3;
Curvature_Higgs_L4[0][0][3][6]=l6i;
Curvature_Higgs_L4[0][0][4][4]=l8;
Curvature_Higgs_L4[0][0][5][2]=l6i;
Curvature_Higgs_L4[0][0][5][5]=l1;
Curvature_Higgs_L4[0][0][5][7]=l6r;
Curvature_Higgs_L4[0][0][6][3]=l6i;
Curvature_Higgs_L4[0][0][6][6]=l1;
Curvature_Higgs_L4[0][0][6][8]=l6r;
Curvature_Higgs_L4[0][0][7][0]=-3*l6i;
Curvature_Higgs_L4[0][0][7][2]=-l5i;
Curvature_Higgs_L4[0][0][7][5]=l6r;
Curvature_Higgs_L4[0][0][7][7]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[0][0][8][1]=-l6i;
Curvature_Higgs_L4[0][0][8][6]=l6r;
Curvature_Higgs_L4[0][0][8][8]=l3;
Curvature_Higgs_L4[0][0][9][9]=l8;
Curvature_Higgs_L4[0][1][0][1]=l1;
Curvature_Higgs_L4[0][1][0][3]=l6r;
Curvature_Higgs_L4[0][1][0][8]=-l6i;
Curvature_Higgs_L4[0][1][1][0]=l1;
Curvature_Higgs_L4[0][1][1][2]=l6r;
Curvature_Higgs_L4[0][1][1][7]=-l6i;
Curvature_Higgs_L4[0][1][2][1]=l6r;
Curvature_Higgs_L4[0][1][2][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][1][2][8]=-0.5*l5i;
Curvature_Higgs_L4[0][1][3][0]=l6r;
Curvature_Higgs_L4[0][1][3][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][1][3][7]=-0.5*l5i;
Curvature_Higgs_L4[0][1][7][1]=-l6i;
Curvature_Higgs_L4[0][1][7][3]=-0.5*l5i;
Curvature_Higgs_L4[0][1][7][8]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[0][1][8][0]=-l6i;
Curvature_Higgs_L4[0][1][8][2]=-0.5*l5i;
Curvature_Higgs_L4[0][1][8][7]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[0][2][0][0]=3*l6r;
Curvature_Higgs_L4[0][2][0][2]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[0][2][0][5]=l6i;
Curvature_Higgs_L4[0][2][0][7]=-l5i;
Curvature_Higgs_L4[0][2][1][1]=l6r;
Curvature_Higgs_L4[0][2][1][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][2][1][8]=-0.5*l5i;
Curvature_Higgs_L4[0][2][2][0]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[0][2][2][2]=3*l7r;
Curvature_Higgs_L4[0][2][2][5]=l5i;
Curvature_Higgs_L4[0][2][2][7]=-l7i;
Curvature_Higgs_L4[0][2][3][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][2][3][3]=l7r;
Curvature_Higgs_L4[0][2][3][6]=l5i/2.;
Curvature_Higgs_L4[0][2][4][4]=l10r;
Curvature_Higgs_L4[0][2][5][0]=l6i;
Curvature_Higgs_L4[0][2][5][2]=l5i;
Curvature_Higgs_L4[0][2][5][5]=l6r;
Curvature_Higgs_L4[0][2][5][7]=l5r;
Curvature_Higgs_L4[0][2][6][3]=l5i/2.;
Curvature_Higgs_L4[0][2][6][6]=l6r;
Curvature_Higgs_L4[0][2][6][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][2][7][0]=-l5i;
Curvature_Higgs_L4[0][2][7][2]=-l7i;
Curvature_Higgs_L4[0][2][7][5]=l5r;
Curvature_Higgs_L4[0][2][7][7]=l7r;
Curvature_Higgs_L4[0][2][8][1]=-0.5*l5i;
Curvature_Higgs_L4[0][2][8][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][2][8][8]=l7r;
Curvature_Higgs_L4[0][2][9][9]=l10r;
Curvature_Higgs_L4[0][3][0][1]=l6r;
Curvature_Higgs_L4[0][3][0][3]=l3;
Curvature_Higgs_L4[0][3][0][6]=l6i;
Curvature_Higgs_L4[0][3][1][0]=l6r;
Curvature_Higgs_L4[0][3][1][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][3][1][7]=-0.5*l5i;
Curvature_Higgs_L4[0][3][2][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][3][2][3]=l7r;
Curvature_Higgs_L4[0][3][2][6]=l5i/2.;
Curvature_Higgs_L4[0][3][3][0]=l3;
Curvature_Higgs_L4[0][3][3][2]=l7r;
Curvature_Higgs_L4[0][3][3][7]=-l7i;
Curvature_Higgs_L4[0][3][6][0]=l6i;
Curvature_Higgs_L4[0][3][6][2]=l5i/2.;
Curvature_Higgs_L4[0][3][6][7]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][3][7][1]=-0.5*l5i;
Curvature_Higgs_L4[0][3][7][3]=-l7i;
Curvature_Higgs_L4[0][3][7][6]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][4][0][4]=l8;
Curvature_Higgs_L4[0][4][2][4]=l10r;
Curvature_Higgs_L4[0][4][4][0]=l8;
Curvature_Higgs_L4[0][4][4][2]=l10r;
Curvature_Higgs_L4[0][4][4][7]=-l10i;
Curvature_Higgs_L4[0][4][7][4]=-l10i;
Curvature_Higgs_L4[0][5][0][2]=l6i;
Curvature_Higgs_L4[0][5][0][5]=l1;
Curvature_Higgs_L4[0][5][0][7]=l6r;
Curvature_Higgs_L4[0][5][2][0]=l6i;
Curvature_Higgs_L4[0][5][2][2]=l5i;
Curvature_Higgs_L4[0][5][2][5]=l6r;
Curvature_Higgs_L4[0][5][2][7]=l5r;
Curvature_Higgs_L4[0][5][5][0]=l1;
Curvature_Higgs_L4[0][5][5][2]=l6r;
Curvature_Higgs_L4[0][5][5][7]=-l6i;
Curvature_Higgs_L4[0][5][7][0]=l6r;
Curvature_Higgs_L4[0][5][7][2]=l5r;
Curvature_Higgs_L4[0][5][7][5]=-l6i;
Curvature_Higgs_L4[0][5][7][7]=-l5i;
Curvature_Higgs_L4[0][6][0][3]=l6i;
Curvature_Higgs_L4[0][6][0][6]=l1;
Curvature_Higgs_L4[0][6][0][8]=l6r;
Curvature_Higgs_L4[0][6][2][3]=l5i/2.;
Curvature_Higgs_L4[0][6][2][6]=l6r;
Curvature_Higgs_L4[0][6][2][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][6][3][0]=l6i;
Curvature_Higgs_L4[0][6][3][2]=l5i/2.;
Curvature_Higgs_L4[0][6][3][7]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][6][6][0]=l1;
Curvature_Higgs_L4[0][6][6][2]=l6r;
Curvature_Higgs_L4[0][6][6][7]=-l6i;
Curvature_Higgs_L4[0][6][7][3]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][6][7][6]=-l6i;
Curvature_Higgs_L4[0][6][7][8]=-0.5*l5i;
Curvature_Higgs_L4[0][6][8][0]=l6r;
Curvature_Higgs_L4[0][6][8][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][6][8][7]=-0.5*l5i;
Curvature_Higgs_L4[0][7][0][0]=-3*l6i;
Curvature_Higgs_L4[0][7][0][2]=-l5i;
Curvature_Higgs_L4[0][7][0][5]=l6r;
Curvature_Higgs_L4[0][7][0][7]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[0][7][1][1]=-l6i;
Curvature_Higgs_L4[0][7][1][3]=-0.5*l5i;
Curvature_Higgs_L4[0][7][1][8]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[0][7][2][0]=-l5i;
Curvature_Higgs_L4[0][7][2][2]=-l7i;
Curvature_Higgs_L4[0][7][2][5]=l5r;
Curvature_Higgs_L4[0][7][2][7]=l7r;
Curvature_Higgs_L4[0][7][3][1]=-0.5*l5i;
Curvature_Higgs_L4[0][7][3][3]=-l7i;
Curvature_Higgs_L4[0][7][3][6]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][7][4][4]=-l10i;
Curvature_Higgs_L4[0][7][5][0]=l6r;
Curvature_Higgs_L4[0][7][5][2]=l5r;
Curvature_Higgs_L4[0][7][5][5]=-l6i;
Curvature_Higgs_L4[0][7][5][7]=-l5i;
Curvature_Higgs_L4[0][7][6][3]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][7][6][6]=-l6i;
Curvature_Higgs_L4[0][7][6][8]=-0.5*l5i;
Curvature_Higgs_L4[0][7][7][0]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[0][7][7][2]=l7r;
Curvature_Higgs_L4[0][7][7][5]=-l5i;
Curvature_Higgs_L4[0][7][7][7]=-3*l7i;
Curvature_Higgs_L4[0][7][8][1]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[0][7][8][6]=-0.5*l5i;
Curvature_Higgs_L4[0][7][8][8]=-l7i;
Curvature_Higgs_L4[0][7][9][9]=-l10i;
Curvature_Higgs_L4[0][8][0][1]=-l6i;
Curvature_Higgs_L4[0][8][0][6]=l6r;
Curvature_Higgs_L4[0][8][0][8]=l3;
Curvature_Higgs_L4[0][8][1][0]=-l6i;
Curvature_Higgs_L4[0][8][1][2]=-0.5*l5i;
Curvature_Higgs_L4[0][8][1][7]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[0][8][2][1]=-0.5*l5i;
Curvature_Higgs_L4[0][8][2][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][8][2][8]=l7r;
Curvature_Higgs_L4[0][8][6][0]=l6r;
Curvature_Higgs_L4[0][8][6][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[0][8][6][7]=-0.5*l5i;
Curvature_Higgs_L4[0][8][7][1]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[0][8][7][6]=-0.5*l5i;
Curvature_Higgs_L4[0][8][7][8]=-l7i;
Curvature_Higgs_L4[0][8][8][0]=l3;
Curvature_Higgs_L4[0][8][8][2]=l7r;
Curvature_Higgs_L4[0][8][8][7]=-l7i;
Curvature_Higgs_L4[0][9][0][9]=l8;
Curvature_Higgs_L4[0][9][2][9]=l10r;
Curvature_Higgs_L4[0][9][7][9]=-l10i;
Curvature_Higgs_L4[0][9][9][0]=l8;
Curvature_Higgs_L4[0][9][9][2]=l10r;
Curvature_Higgs_L4[0][9][9][7]=-l10i;
Curvature_Higgs_L4[1][0][0][1]=l1;
Curvature_Higgs_L4[1][0][0][3]=l6r;
Curvature_Higgs_L4[1][0][0][8]=-l6i;
Curvature_Higgs_L4[1][0][1][0]=l1;
Curvature_Higgs_L4[1][0][1][2]=l6r;
Curvature_Higgs_L4[1][0][1][7]=-l6i;
Curvature_Higgs_L4[1][0][2][1]=l6r;
Curvature_Higgs_L4[1][0][2][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][0][2][8]=-0.5*l5i;
Curvature_Higgs_L4[1][0][3][0]=l6r;
Curvature_Higgs_L4[1][0][3][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][0][3][7]=-0.5*l5i;
Curvature_Higgs_L4[1][0][7][1]=-l6i;
Curvature_Higgs_L4[1][0][7][3]=-0.5*l5i;
Curvature_Higgs_L4[1][0][7][8]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[1][0][8][0]=-l6i;
Curvature_Higgs_L4[1][0][8][2]=-0.5*l5i;
Curvature_Higgs_L4[1][0][8][7]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[1][1][0][0]=l1;
Curvature_Higgs_L4[1][1][0][2]=l6r;
Curvature_Higgs_L4[1][1][0][7]=-l6i;
Curvature_Higgs_L4[1][1][1][1]=3*l1;
Curvature_Higgs_L4[1][1][1][3]=3*l6r;
Curvature_Higgs_L4[1][1][1][8]=-3*l6i;
Curvature_Higgs_L4[1][1][2][0]=l6r;
Curvature_Higgs_L4[1][1][2][2]=l3;
Curvature_Higgs_L4[1][1][2][5]=l6i;
Curvature_Higgs_L4[1][1][3][1]=3*l6r;
Curvature_Higgs_L4[1][1][3][3]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[1][1][3][6]=l6i;
Curvature_Higgs_L4[1][1][3][8]=-l5i;
Curvature_Higgs_L4[1][1][4][4]=l8;
Curvature_Higgs_L4[1][1][5][2]=l6i;
Curvature_Higgs_L4[1][1][5][5]=l1;
Curvature_Higgs_L4[1][1][5][7]=l6r;
Curvature_Higgs_L4[1][1][6][3]=l6i;
Curvature_Higgs_L4[1][1][6][6]=l1;
Curvature_Higgs_L4[1][1][6][8]=l6r;
Curvature_Higgs_L4[1][1][7][0]=-l6i;
Curvature_Higgs_L4[1][1][7][5]=l6r;
Curvature_Higgs_L4[1][1][7][7]=l3;
Curvature_Higgs_L4[1][1][8][1]=-3*l6i;
Curvature_Higgs_L4[1][1][8][3]=-l5i;
Curvature_Higgs_L4[1][1][8][6]=l6r;
Curvature_Higgs_L4[1][1][8][8]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[1][1][9][9]=l8;
Curvature_Higgs_L4[1][2][0][1]=l6r;
Curvature_Higgs_L4[1][2][0][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][2][0][8]=-0.5*l5i;
Curvature_Higgs_L4[1][2][1][0]=l6r;
Curvature_Higgs_L4[1][2][1][2]=l3;
Curvature_Higgs_L4[1][2][1][5]=l6i;
Curvature_Higgs_L4[1][2][2][1]=l3;
Curvature_Higgs_L4[1][2][2][3]=l7r;
Curvature_Higgs_L4[1][2][2][8]=-l7i;
Curvature_Higgs_L4[1][2][3][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][2][3][2]=l7r;
Curvature_Higgs_L4[1][2][3][5]=l5i/2.;
Curvature_Higgs_L4[1][2][5][1]=l6i;
Curvature_Higgs_L4[1][2][5][3]=l5i/2.;
Curvature_Higgs_L4[1][2][5][8]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][2][8][0]=-0.5*l5i;
Curvature_Higgs_L4[1][2][8][2]=-l7i;
Curvature_Higgs_L4[1][2][8][5]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][3][0][0]=l6r;
Curvature_Higgs_L4[1][3][0][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][3][0][7]=-0.5*l5i;
Curvature_Higgs_L4[1][3][1][1]=3*l6r;
Curvature_Higgs_L4[1][3][1][3]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[1][3][1][6]=l6i;
Curvature_Higgs_L4[1][3][1][8]=-l5i;
Curvature_Higgs_L4[1][3][2][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][3][2][2]=l7r;
Curvature_Higgs_L4[1][3][2][5]=l5i/2.;
Curvature_Higgs_L4[1][3][3][1]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[1][3][3][3]=3*l7r;
Curvature_Higgs_L4[1][3][3][6]=l5i;
Curvature_Higgs_L4[1][3][3][8]=-l7i;
Curvature_Higgs_L4[1][3][4][4]=l10r;
Curvature_Higgs_L4[1][3][5][2]=l5i/2.;
Curvature_Higgs_L4[1][3][5][5]=l6r;
Curvature_Higgs_L4[1][3][5][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][3][6][1]=l6i;
Curvature_Higgs_L4[1][3][6][3]=l5i;
Curvature_Higgs_L4[1][3][6][6]=l6r;
Curvature_Higgs_L4[1][3][6][8]=l5r;
Curvature_Higgs_L4[1][3][7][0]=-0.5*l5i;
Curvature_Higgs_L4[1][3][7][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][3][7][7]=l7r;
Curvature_Higgs_L4[1][3][8][1]=-l5i;
Curvature_Higgs_L4[1][3][8][3]=-l7i;
Curvature_Higgs_L4[1][3][8][6]=l5r;
Curvature_Higgs_L4[1][3][8][8]=l7r;
Curvature_Higgs_L4[1][3][9][9]=l10r;
Curvature_Higgs_L4[1][4][1][4]=l8;
Curvature_Higgs_L4[1][4][3][4]=l10r;
Curvature_Higgs_L4[1][4][4][1]=l8;
Curvature_Higgs_L4[1][4][4][3]=l10r;
Curvature_Higgs_L4[1][4][4][8]=-l10i;
Curvature_Higgs_L4[1][4][8][4]=-l10i;
Curvature_Higgs_L4[1][5][1][2]=l6i;
Curvature_Higgs_L4[1][5][1][5]=l1;
Curvature_Higgs_L4[1][5][1][7]=l6r;
Curvature_Higgs_L4[1][5][2][1]=l6i;
Curvature_Higgs_L4[1][5][2][3]=l5i/2.;
Curvature_Higgs_L4[1][5][2][8]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][5][3][2]=l5i/2.;
Curvature_Higgs_L4[1][5][3][5]=l6r;
Curvature_Higgs_L4[1][5][3][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][5][5][1]=l1;
Curvature_Higgs_L4[1][5][5][3]=l6r;
Curvature_Higgs_L4[1][5][5][8]=-l6i;
Curvature_Higgs_L4[1][5][7][1]=l6r;
Curvature_Higgs_L4[1][5][7][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][5][7][8]=-0.5*l5i;
Curvature_Higgs_L4[1][5][8][2]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][5][8][5]=-l6i;
Curvature_Higgs_L4[1][5][8][7]=-0.5*l5i;
Curvature_Higgs_L4[1][6][1][3]=l6i;
Curvature_Higgs_L4[1][6][1][6]=l1;
Curvature_Higgs_L4[1][6][1][8]=l6r;
Curvature_Higgs_L4[1][6][3][1]=l6i;
Curvature_Higgs_L4[1][6][3][3]=l5i;
Curvature_Higgs_L4[1][6][3][6]=l6r;
Curvature_Higgs_L4[1][6][3][8]=l5r;
Curvature_Higgs_L4[1][6][6][1]=l1;
Curvature_Higgs_L4[1][6][6][3]=l6r;
Curvature_Higgs_L4[1][6][6][8]=-l6i;
Curvature_Higgs_L4[1][6][8][1]=l6r;
Curvature_Higgs_L4[1][6][8][3]=l5r;
Curvature_Higgs_L4[1][6][8][6]=-l6i;
Curvature_Higgs_L4[1][6][8][8]=-l5i;
Curvature_Higgs_L4[1][7][0][1]=-l6i;
Curvature_Higgs_L4[1][7][0][3]=-0.5*l5i;
Curvature_Higgs_L4[1][7][0][8]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[1][7][1][0]=-l6i;
Curvature_Higgs_L4[1][7][1][5]=l6r;
Curvature_Higgs_L4[1][7][1][7]=l3;
Curvature_Higgs_L4[1][7][3][0]=-0.5*l5i;
Curvature_Higgs_L4[1][7][3][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][7][3][7]=l7r;
Curvature_Higgs_L4[1][7][5][1]=l6r;
Curvature_Higgs_L4[1][7][5][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][7][5][8]=-0.5*l5i;
Curvature_Higgs_L4[1][7][7][1]=l3;
Curvature_Higgs_L4[1][7][7][3]=l7r;
Curvature_Higgs_L4[1][7][7][8]=-l7i;
Curvature_Higgs_L4[1][7][8][0]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[1][7][8][5]=-0.5*l5i;
Curvature_Higgs_L4[1][7][8][7]=-l7i;
Curvature_Higgs_L4[1][8][0][0]=-l6i;
Curvature_Higgs_L4[1][8][0][2]=-0.5*l5i;
Curvature_Higgs_L4[1][8][0][7]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[1][8][1][1]=-3*l6i;
Curvature_Higgs_L4[1][8][1][3]=-l5i;
Curvature_Higgs_L4[1][8][1][6]=l6r;
Curvature_Higgs_L4[1][8][1][8]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[1][8][2][0]=-0.5*l5i;
Curvature_Higgs_L4[1][8][2][2]=-l7i;
Curvature_Higgs_L4[1][8][2][5]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][8][3][1]=-l5i;
Curvature_Higgs_L4[1][8][3][3]=-l7i;
Curvature_Higgs_L4[1][8][3][6]=l5r;
Curvature_Higgs_L4[1][8][3][8]=l7r;
Curvature_Higgs_L4[1][8][4][4]=-l10i;
Curvature_Higgs_L4[1][8][5][2]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[1][8][5][5]=-l6i;
Curvature_Higgs_L4[1][8][5][7]=-0.5*l5i;
Curvature_Higgs_L4[1][8][6][1]=l6r;
Curvature_Higgs_L4[1][8][6][3]=l5r;
Curvature_Higgs_L4[1][8][6][6]=-l6i;
Curvature_Higgs_L4[1][8][6][8]=-l5i;
Curvature_Higgs_L4[1][8][7][0]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[1][8][7][5]=-0.5*l5i;
Curvature_Higgs_L4[1][8][7][7]=-l7i;
Curvature_Higgs_L4[1][8][8][1]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[1][8][8][3]=l7r;
Curvature_Higgs_L4[1][8][8][6]=-l5i;
Curvature_Higgs_L4[1][8][8][8]=-3*l7i;
Curvature_Higgs_L4[1][8][9][9]=-l10i;
Curvature_Higgs_L4[1][9][1][9]=l8;
Curvature_Higgs_L4[1][9][3][9]=l10r;
Curvature_Higgs_L4[1][9][8][9]=-l10i;
Curvature_Higgs_L4[1][9][9][1]=l8;
Curvature_Higgs_L4[1][9][9][3]=l10r;
Curvature_Higgs_L4[1][9][9][8]=-l10i;
Curvature_Higgs_L4[2][0][0][0]=3*l6r;
Curvature_Higgs_L4[2][0][0][2]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[2][0][0][5]=l6i;
Curvature_Higgs_L4[2][0][0][7]=-l5i;
Curvature_Higgs_L4[2][0][1][1]=l6r;
Curvature_Higgs_L4[2][0][1][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][0][1][8]=-0.5*l5i;
Curvature_Higgs_L4[2][0][2][0]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[2][0][2][2]=3*l7r;
Curvature_Higgs_L4[2][0][2][5]=l5i;
Curvature_Higgs_L4[2][0][2][7]=-l7i;
Curvature_Higgs_L4[2][0][3][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][0][3][3]=l7r;
Curvature_Higgs_L4[2][0][3][6]=l5i/2.;
Curvature_Higgs_L4[2][0][4][4]=l10r;
Curvature_Higgs_L4[2][0][5][0]=l6i;
Curvature_Higgs_L4[2][0][5][2]=l5i;
Curvature_Higgs_L4[2][0][5][5]=l6r;
Curvature_Higgs_L4[2][0][5][7]=l5r;
Curvature_Higgs_L4[2][0][6][3]=l5i/2.;
Curvature_Higgs_L4[2][0][6][6]=l6r;
Curvature_Higgs_L4[2][0][6][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][0][7][0]=-l5i;
Curvature_Higgs_L4[2][0][7][2]=-l7i;
Curvature_Higgs_L4[2][0][7][5]=l5r;
Curvature_Higgs_L4[2][0][7][7]=l7r;
Curvature_Higgs_L4[2][0][8][1]=-0.5*l5i;
Curvature_Higgs_L4[2][0][8][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][0][8][8]=l7r;
Curvature_Higgs_L4[2][0][9][9]=l10r;
Curvature_Higgs_L4[2][1][0][1]=l6r;
Curvature_Higgs_L4[2][1][0][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][1][0][8]=-0.5*l5i;
Curvature_Higgs_L4[2][1][1][0]=l6r;
Curvature_Higgs_L4[2][1][1][2]=l3;
Curvature_Higgs_L4[2][1][1][5]=l6i;
Curvature_Higgs_L4[2][1][2][1]=l3;
Curvature_Higgs_L4[2][1][2][3]=l7r;
Curvature_Higgs_L4[2][1][2][8]=-l7i;
Curvature_Higgs_L4[2][1][3][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][1][3][2]=l7r;
Curvature_Higgs_L4[2][1][3][5]=l5i/2.;
Curvature_Higgs_L4[2][1][5][1]=l6i;
Curvature_Higgs_L4[2][1][5][3]=l5i/2.;
Curvature_Higgs_L4[2][1][5][8]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][1][8][0]=-0.5*l5i;
Curvature_Higgs_L4[2][1][8][2]=-l7i;
Curvature_Higgs_L4[2][1][8][5]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][2][0][0]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[2][2][0][2]=3*l7r;
Curvature_Higgs_L4[2][2][0][5]=l5i;
Curvature_Higgs_L4[2][2][0][7]=-l7i;
Curvature_Higgs_L4[2][2][1][1]=l3;
Curvature_Higgs_L4[2][2][1][3]=l7r;
Curvature_Higgs_L4[2][2][1][8]=-l7i;
Curvature_Higgs_L4[2][2][2][0]=3*l7r;
Curvature_Higgs_L4[2][2][2][2]=3*l2;
Curvature_Higgs_L4[2][2][2][5]=3*l7i;
Curvature_Higgs_L4[2][2][3][1]=l7r;
Curvature_Higgs_L4[2][2][3][3]=l2;
Curvature_Higgs_L4[2][2][3][6]=l7i;
Curvature_Higgs_L4[2][2][4][4]=l9;
Curvature_Higgs_L4[2][2][5][0]=l5i;
Curvature_Higgs_L4[2][2][5][2]=3*l7i;
Curvature_Higgs_L4[2][2][5][5]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[2][2][5][7]=l7r;
Curvature_Higgs_L4[2][2][6][3]=l7i;
Curvature_Higgs_L4[2][2][6][6]=l3;
Curvature_Higgs_L4[2][2][6][8]=l7r;
Curvature_Higgs_L4[2][2][7][0]=-l7i;
Curvature_Higgs_L4[2][2][7][5]=l7r;
Curvature_Higgs_L4[2][2][7][7]=l2;
Curvature_Higgs_L4[2][2][8][1]=-l7i;
Curvature_Higgs_L4[2][2][8][6]=l7r;
Curvature_Higgs_L4[2][2][8][8]=l2;
Curvature_Higgs_L4[2][2][9][9]=l9;
Curvature_Higgs_L4[2][3][0][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][3][0][3]=l7r;
Curvature_Higgs_L4[2][3][0][6]=l5i/2.;
Curvature_Higgs_L4[2][3][1][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][3][1][2]=l7r;
Curvature_Higgs_L4[2][3][1][5]=l5i/2.;
Curvature_Higgs_L4[2][3][2][1]=l7r;
Curvature_Higgs_L4[2][3][2][3]=l2;
Curvature_Higgs_L4[2][3][2][6]=l7i;
Curvature_Higgs_L4[2][3][3][0]=l7r;
Curvature_Higgs_L4[2][3][3][2]=l2;
Curvature_Higgs_L4[2][3][3][5]=l7i;
Curvature_Higgs_L4[2][3][5][1]=l5i/2.;
Curvature_Higgs_L4[2][3][5][3]=l7i;
Curvature_Higgs_L4[2][3][5][6]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[2][3][6][0]=l5i/2.;
Curvature_Higgs_L4[2][3][6][2]=l7i;
Curvature_Higgs_L4[2][3][6][5]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[2][4][0][4]=l10r;
Curvature_Higgs_L4[2][4][2][4]=l9;
Curvature_Higgs_L4[2][4][4][0]=l10r;
Curvature_Higgs_L4[2][4][4][2]=l9;
Curvature_Higgs_L4[2][4][4][5]=l10i;
Curvature_Higgs_L4[2][4][5][4]=l10i;
Curvature_Higgs_L4[2][5][0][0]=l6i;
Curvature_Higgs_L4[2][5][0][2]=l5i;
Curvature_Higgs_L4[2][5][0][5]=l6r;
Curvature_Higgs_L4[2][5][0][7]=l5r;
Curvature_Higgs_L4[2][5][1][1]=l6i;
Curvature_Higgs_L4[2][5][1][3]=l5i/2.;
Curvature_Higgs_L4[2][5][1][8]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][5][2][0]=l5i;
Curvature_Higgs_L4[2][5][2][2]=3*l7i;
Curvature_Higgs_L4[2][5][2][5]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[2][5][2][7]=l7r;
Curvature_Higgs_L4[2][5][3][1]=l5i/2.;
Curvature_Higgs_L4[2][5][3][3]=l7i;
Curvature_Higgs_L4[2][5][3][6]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[2][5][4][4]=l10i;
Curvature_Higgs_L4[2][5][5][0]=l6r;
Curvature_Higgs_L4[2][5][5][2]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[2][5][5][5]=3*l6i;
Curvature_Higgs_L4[2][5][5][7]=l5i;
Curvature_Higgs_L4[2][5][6][3]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[2][5][6][6]=l6i;
Curvature_Higgs_L4[2][5][6][8]=l5i/2.;
Curvature_Higgs_L4[2][5][7][0]=l5r;
Curvature_Higgs_L4[2][5][7][2]=l7r;
Curvature_Higgs_L4[2][5][7][5]=l5i;
Curvature_Higgs_L4[2][5][7][7]=l7i;
Curvature_Higgs_L4[2][5][8][1]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][5][8][6]=l5i/2.;
Curvature_Higgs_L4[2][5][8][8]=l7i;
Curvature_Higgs_L4[2][5][9][9]=l10i;
Curvature_Higgs_L4[2][6][0][3]=l5i/2.;
Curvature_Higgs_L4[2][6][0][6]=l6r;
Curvature_Higgs_L4[2][6][0][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][6][2][3]=l7i;
Curvature_Higgs_L4[2][6][2][6]=l3;
Curvature_Higgs_L4[2][6][2][8]=l7r;
Curvature_Higgs_L4[2][6][3][0]=l5i/2.;
Curvature_Higgs_L4[2][6][3][2]=l7i;
Curvature_Higgs_L4[2][6][3][5]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[2][6][5][3]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[2][6][5][6]=l6i;
Curvature_Higgs_L4[2][6][5][8]=l5i/2.;
Curvature_Higgs_L4[2][6][6][0]=l6r;
Curvature_Higgs_L4[2][6][6][2]=l3;
Curvature_Higgs_L4[2][6][6][5]=l6i;
Curvature_Higgs_L4[2][6][8][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][6][8][2]=l7r;
Curvature_Higgs_L4[2][6][8][5]=l5i/2.;
Curvature_Higgs_L4[2][7][0][0]=-l5i;
Curvature_Higgs_L4[2][7][0][2]=-l7i;
Curvature_Higgs_L4[2][7][0][5]=l5r;
Curvature_Higgs_L4[2][7][0][7]=l7r;
Curvature_Higgs_L4[2][7][2][0]=-l7i;
Curvature_Higgs_L4[2][7][2][5]=l7r;
Curvature_Higgs_L4[2][7][2][7]=l2;
Curvature_Higgs_L4[2][7][5][0]=l5r;
Curvature_Higgs_L4[2][7][5][2]=l7r;
Curvature_Higgs_L4[2][7][5][5]=l5i;
Curvature_Higgs_L4[2][7][5][7]=l7i;
Curvature_Higgs_L4[2][7][7][0]=l7r;
Curvature_Higgs_L4[2][7][7][2]=l2;
Curvature_Higgs_L4[2][7][7][5]=l7i;
Curvature_Higgs_L4[2][8][0][1]=-0.5*l5i;
Curvature_Higgs_L4[2][8][0][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][8][0][8]=l7r;
Curvature_Higgs_L4[2][8][1][0]=-0.5*l5i;
Curvature_Higgs_L4[2][8][1][2]=-l7i;
Curvature_Higgs_L4[2][8][1][5]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][8][2][1]=-l7i;
Curvature_Higgs_L4[2][8][2][6]=l7r;
Curvature_Higgs_L4[2][8][2][8]=l2;
Curvature_Higgs_L4[2][8][5][1]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][8][5][6]=l5i/2.;
Curvature_Higgs_L4[2][8][5][8]=l7i;
Curvature_Higgs_L4[2][8][6][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[2][8][6][2]=l7r;
Curvature_Higgs_L4[2][8][6][5]=l5i/2.;
Curvature_Higgs_L4[2][8][8][0]=l7r;
Curvature_Higgs_L4[2][8][8][2]=l2;
Curvature_Higgs_L4[2][8][8][5]=l7i;
Curvature_Higgs_L4[2][9][0][9]=l10r;
Curvature_Higgs_L4[2][9][2][9]=l9;
Curvature_Higgs_L4[2][9][5][9]=l10i;
Curvature_Higgs_L4[2][9][9][0]=l10r;
Curvature_Higgs_L4[2][9][9][2]=l9;
Curvature_Higgs_L4[2][9][9][5]=l10i;
Curvature_Higgs_L4[3][0][0][1]=l6r;
Curvature_Higgs_L4[3][0][0][3]=l3;
Curvature_Higgs_L4[3][0][0][6]=l6i;
Curvature_Higgs_L4[3][0][1][0]=l6r;
Curvature_Higgs_L4[3][0][1][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][0][1][7]=-0.5*l5i;
Curvature_Higgs_L4[3][0][2][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][0][2][3]=l7r;
Curvature_Higgs_L4[3][0][2][6]=l5i/2.;
Curvature_Higgs_L4[3][0][3][0]=l3;
Curvature_Higgs_L4[3][0][3][2]=l7r;
Curvature_Higgs_L4[3][0][3][7]=-l7i;
Curvature_Higgs_L4[3][0][6][0]=l6i;
Curvature_Higgs_L4[3][0][6][2]=l5i/2.;
Curvature_Higgs_L4[3][0][6][7]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][0][7][1]=-0.5*l5i;
Curvature_Higgs_L4[3][0][7][3]=-l7i;
Curvature_Higgs_L4[3][0][7][6]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][1][0][0]=l6r;
Curvature_Higgs_L4[3][1][0][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][1][0][7]=-0.5*l5i;
Curvature_Higgs_L4[3][1][1][1]=3*l6r;
Curvature_Higgs_L4[3][1][1][3]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[3][1][1][6]=l6i;
Curvature_Higgs_L4[3][1][1][8]=-l5i;
Curvature_Higgs_L4[3][1][2][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][1][2][2]=l7r;
Curvature_Higgs_L4[3][1][2][5]=l5i/2.;
Curvature_Higgs_L4[3][1][3][1]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[3][1][3][3]=3*l7r;
Curvature_Higgs_L4[3][1][3][6]=l5i;
Curvature_Higgs_L4[3][1][3][8]=-l7i;
Curvature_Higgs_L4[3][1][4][4]=l10r;
Curvature_Higgs_L4[3][1][5][2]=l5i/2.;
Curvature_Higgs_L4[3][1][5][5]=l6r;
Curvature_Higgs_L4[3][1][5][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][1][6][1]=l6i;
Curvature_Higgs_L4[3][1][6][3]=l5i;
Curvature_Higgs_L4[3][1][6][6]=l6r;
Curvature_Higgs_L4[3][1][6][8]=l5r;
Curvature_Higgs_L4[3][1][7][0]=-0.5*l5i;
Curvature_Higgs_L4[3][1][7][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][1][7][7]=l7r;
Curvature_Higgs_L4[3][1][8][1]=-l5i;
Curvature_Higgs_L4[3][1][8][3]=-l7i;
Curvature_Higgs_L4[3][1][8][6]=l5r;
Curvature_Higgs_L4[3][1][8][8]=l7r;
Curvature_Higgs_L4[3][1][9][9]=l10r;
Curvature_Higgs_L4[3][2][0][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][2][0][3]=l7r;
Curvature_Higgs_L4[3][2][0][6]=l5i/2.;
Curvature_Higgs_L4[3][2][1][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][2][1][2]=l7r;
Curvature_Higgs_L4[3][2][1][5]=l5i/2.;
Curvature_Higgs_L4[3][2][2][1]=l7r;
Curvature_Higgs_L4[3][2][2][3]=l2;
Curvature_Higgs_L4[3][2][2][6]=l7i;
Curvature_Higgs_L4[3][2][3][0]=l7r;
Curvature_Higgs_L4[3][2][3][2]=l2;
Curvature_Higgs_L4[3][2][3][5]=l7i;
Curvature_Higgs_L4[3][2][5][1]=l5i/2.;
Curvature_Higgs_L4[3][2][5][3]=l7i;
Curvature_Higgs_L4[3][2][5][6]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[3][2][6][0]=l5i/2.;
Curvature_Higgs_L4[3][2][6][2]=l7i;
Curvature_Higgs_L4[3][2][6][5]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[3][3][0][0]=l3;
Curvature_Higgs_L4[3][3][0][2]=l7r;
Curvature_Higgs_L4[3][3][0][7]=-l7i;
Curvature_Higgs_L4[3][3][1][1]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[3][3][1][3]=3*l7r;
Curvature_Higgs_L4[3][3][1][6]=l5i;
Curvature_Higgs_L4[3][3][1][8]=-l7i;
Curvature_Higgs_L4[3][3][2][0]=l7r;
Curvature_Higgs_L4[3][3][2][2]=l2;
Curvature_Higgs_L4[3][3][2][5]=l7i;
Curvature_Higgs_L4[3][3][3][1]=3*l7r;
Curvature_Higgs_L4[3][3][3][3]=3*l2;
Curvature_Higgs_L4[3][3][3][6]=3*l7i;
Curvature_Higgs_L4[3][3][4][4]=l9;
Curvature_Higgs_L4[3][3][5][2]=l7i;
Curvature_Higgs_L4[3][3][5][5]=l3;
Curvature_Higgs_L4[3][3][5][7]=l7r;
Curvature_Higgs_L4[3][3][6][1]=l5i;
Curvature_Higgs_L4[3][3][6][3]=3*l7i;
Curvature_Higgs_L4[3][3][6][6]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[3][3][6][8]=l7r;
Curvature_Higgs_L4[3][3][7][0]=-l7i;
Curvature_Higgs_L4[3][3][7][5]=l7r;
Curvature_Higgs_L4[3][3][7][7]=l2;
Curvature_Higgs_L4[3][3][8][1]=-l7i;
Curvature_Higgs_L4[3][3][8][6]=l7r;
Curvature_Higgs_L4[3][3][8][8]=l2;
Curvature_Higgs_L4[3][3][9][9]=l9;
Curvature_Higgs_L4[3][4][1][4]=l10r;
Curvature_Higgs_L4[3][4][3][4]=l9;
Curvature_Higgs_L4[3][4][4][1]=l10r;
Curvature_Higgs_L4[3][4][4][3]=l9;
Curvature_Higgs_L4[3][4][4][6]=l10i;
Curvature_Higgs_L4[3][4][6][4]=l10i;
Curvature_Higgs_L4[3][5][1][2]=l5i/2.;
Curvature_Higgs_L4[3][5][1][5]=l6r;
Curvature_Higgs_L4[3][5][1][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][5][2][1]=l5i/2.;
Curvature_Higgs_L4[3][5][2][3]=l7i;
Curvature_Higgs_L4[3][5][2][6]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[3][5][3][2]=l7i;
Curvature_Higgs_L4[3][5][3][5]=l3;
Curvature_Higgs_L4[3][5][3][7]=l7r;
Curvature_Higgs_L4[3][5][5][1]=l6r;
Curvature_Higgs_L4[3][5][5][3]=l3;
Curvature_Higgs_L4[3][5][5][6]=l6i;
Curvature_Higgs_L4[3][5][6][2]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[3][5][6][5]=l6i;
Curvature_Higgs_L4[3][5][6][7]=l5i/2.;
Curvature_Higgs_L4[3][5][7][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][5][7][3]=l7r;
Curvature_Higgs_L4[3][5][7][6]=l5i/2.;
Curvature_Higgs_L4[3][6][0][0]=l6i;
Curvature_Higgs_L4[3][6][0][2]=l5i/2.;
Curvature_Higgs_L4[3][6][0][7]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][6][1][1]=l6i;
Curvature_Higgs_L4[3][6][1][3]=l5i;
Curvature_Higgs_L4[3][6][1][6]=l6r;
Curvature_Higgs_L4[3][6][1][8]=l5r;
Curvature_Higgs_L4[3][6][2][0]=l5i/2.;
Curvature_Higgs_L4[3][6][2][2]=l7i;
Curvature_Higgs_L4[3][6][2][5]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[3][6][3][1]=l5i;
Curvature_Higgs_L4[3][6][3][3]=3*l7i;
Curvature_Higgs_L4[3][6][3][6]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[3][6][3][8]=l7r;
Curvature_Higgs_L4[3][6][4][4]=l10i;
Curvature_Higgs_L4[3][6][5][2]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[3][6][5][5]=l6i;
Curvature_Higgs_L4[3][6][5][7]=l5i/2.;
Curvature_Higgs_L4[3][6][6][1]=l6r;
Curvature_Higgs_L4[3][6][6][3]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[3][6][6][6]=3*l6i;
Curvature_Higgs_L4[3][6][6][8]=l5i;
Curvature_Higgs_L4[3][6][7][0]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][6][7][5]=l5i/2.;
Curvature_Higgs_L4[3][6][7][7]=l7i;
Curvature_Higgs_L4[3][6][8][1]=l5r;
Curvature_Higgs_L4[3][6][8][3]=l7r;
Curvature_Higgs_L4[3][6][8][6]=l5i;
Curvature_Higgs_L4[3][6][8][8]=l7i;
Curvature_Higgs_L4[3][6][9][9]=l10i;
Curvature_Higgs_L4[3][7][0][1]=-0.5*l5i;
Curvature_Higgs_L4[3][7][0][3]=-l7i;
Curvature_Higgs_L4[3][7][0][6]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][7][1][0]=-0.5*l5i;
Curvature_Higgs_L4[3][7][1][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][7][1][7]=l7r;
Curvature_Higgs_L4[3][7][3][0]=-l7i;
Curvature_Higgs_L4[3][7][3][5]=l7r;
Curvature_Higgs_L4[3][7][3][7]=l2;
Curvature_Higgs_L4[3][7][5][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][7][5][3]=l7r;
Curvature_Higgs_L4[3][7][5][6]=l5i/2.;
Curvature_Higgs_L4[3][7][6][0]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[3][7][6][5]=l5i/2.;
Curvature_Higgs_L4[3][7][6][7]=l7i;
Curvature_Higgs_L4[3][7][7][1]=l7r;
Curvature_Higgs_L4[3][7][7][3]=l2;
Curvature_Higgs_L4[3][7][7][6]=l7i;
Curvature_Higgs_L4[3][8][1][1]=-l5i;
Curvature_Higgs_L4[3][8][1][3]=-l7i;
Curvature_Higgs_L4[3][8][1][6]=l5r;
Curvature_Higgs_L4[3][8][1][8]=l7r;
Curvature_Higgs_L4[3][8][3][1]=-l7i;
Curvature_Higgs_L4[3][8][3][6]=l7r;
Curvature_Higgs_L4[3][8][3][8]=l2;
Curvature_Higgs_L4[3][8][6][1]=l5r;
Curvature_Higgs_L4[3][8][6][3]=l7r;
Curvature_Higgs_L4[3][8][6][6]=l5i;
Curvature_Higgs_L4[3][8][6][8]=l7i;
Curvature_Higgs_L4[3][8][8][1]=l7r;
Curvature_Higgs_L4[3][8][8][3]=l2;
Curvature_Higgs_L4[3][8][8][6]=l7i;
Curvature_Higgs_L4[3][9][1][9]=l10r;
Curvature_Higgs_L4[3][9][3][9]=l9;
Curvature_Higgs_L4[3][9][6][9]=l10i;
Curvature_Higgs_L4[3][9][9][1]=l10r;
Curvature_Higgs_L4[3][9][9][3]=l9;
Curvature_Higgs_L4[3][9][9][6]=l10i;
Curvature_Higgs_L4[4][0][0][4]=l8;
Curvature_Higgs_L4[4][0][2][4]=l10r;
Curvature_Higgs_L4[4][0][4][0]=l8;
Curvature_Higgs_L4[4][0][4][2]=l10r;
Curvature_Higgs_L4[4][0][4][7]=-l10i;
Curvature_Higgs_L4[4][0][7][4]=-l10i;
Curvature_Higgs_L4[4][1][1][4]=l8;
Curvature_Higgs_L4[4][1][3][4]=l10r;
Curvature_Higgs_L4[4][1][4][1]=l8;
Curvature_Higgs_L4[4][1][4][3]=l10r;
Curvature_Higgs_L4[4][1][4][8]=-l10i;
Curvature_Higgs_L4[4][1][8][4]=-l10i;
Curvature_Higgs_L4[4][2][0][4]=l10r;
Curvature_Higgs_L4[4][2][2][4]=l9;
Curvature_Higgs_L4[4][2][4][0]=l10r;
Curvature_Higgs_L4[4][2][4][2]=l9;
Curvature_Higgs_L4[4][2][4][5]=l10i;
Curvature_Higgs_L4[4][2][5][4]=l10i;
Curvature_Higgs_L4[4][3][1][4]=l10r;
Curvature_Higgs_L4[4][3][3][4]=l9;
Curvature_Higgs_L4[4][3][4][1]=l10r;
Curvature_Higgs_L4[4][3][4][3]=l9;
Curvature_Higgs_L4[4][3][4][6]=l10i;
Curvature_Higgs_L4[4][3][6][4]=l10i;
Curvature_Higgs_L4[4][4][0][0]=l8;
Curvature_Higgs_L4[4][4][0][2]=l10r;
Curvature_Higgs_L4[4][4][0][7]=-l10i;
Curvature_Higgs_L4[4][4][1][1]=l8;
Curvature_Higgs_L4[4][4][1][3]=l10r;
Curvature_Higgs_L4[4][4][1][8]=-l10i;
Curvature_Higgs_L4[4][4][2][0]=l10r;
Curvature_Higgs_L4[4][4][2][2]=l9;
Curvature_Higgs_L4[4][4][2][5]=l10i;
Curvature_Higgs_L4[4][4][3][1]=l10r;
Curvature_Higgs_L4[4][4][3][3]=l9;
Curvature_Higgs_L4[4][4][3][6]=l10i;
Curvature_Higgs_L4[4][4][4][4]=6*lh;
Curvature_Higgs_L4[4][4][5][2]=l10i;
Curvature_Higgs_L4[4][4][5][5]=l8;
Curvature_Higgs_L4[4][4][5][7]=l10r;
Curvature_Higgs_L4[4][4][6][3]=l10i;
Curvature_Higgs_L4[4][4][6][6]=l8;
Curvature_Higgs_L4[4][4][6][8]=l10r;
Curvature_Higgs_L4[4][4][7][0]=-l10i;
Curvature_Higgs_L4[4][4][7][5]=l10r;
Curvature_Higgs_L4[4][4][7][7]=l9;
Curvature_Higgs_L4[4][4][8][1]=-l10i;
Curvature_Higgs_L4[4][4][8][6]=l10r;
Curvature_Higgs_L4[4][4][8][8]=l9;
Curvature_Higgs_L4[4][4][9][9]=2*lh;
Curvature_Higgs_L4[4][5][2][4]=l10i;
Curvature_Higgs_L4[4][5][4][2]=l10i;
Curvature_Higgs_L4[4][5][4][5]=l8;
Curvature_Higgs_L4[4][5][4][7]=l10r;
Curvature_Higgs_L4[4][5][5][4]=l8;
Curvature_Higgs_L4[4][5][7][4]=l10r;
Curvature_Higgs_L4[4][6][3][4]=l10i;
Curvature_Higgs_L4[4][6][4][3]=l10i;
Curvature_Higgs_L4[4][6][4][6]=l8;
Curvature_Higgs_L4[4][6][4][8]=l10r;
Curvature_Higgs_L4[4][6][6][4]=l8;
Curvature_Higgs_L4[4][6][8][4]=l10r;
Curvature_Higgs_L4[4][7][0][4]=-l10i;
Curvature_Higgs_L4[4][7][4][0]=-l10i;
Curvature_Higgs_L4[4][7][4][5]=l10r;
Curvature_Higgs_L4[4][7][4][7]=l9;
Curvature_Higgs_L4[4][7][5][4]=l10r;
Curvature_Higgs_L4[4][7][7][4]=l9;
Curvature_Higgs_L4[4][8][1][4]=-l10i;
Curvature_Higgs_L4[4][8][4][1]=-l10i;
Curvature_Higgs_L4[4][8][4][6]=l10r;
Curvature_Higgs_L4[4][8][4][8]=l9;
Curvature_Higgs_L4[4][8][6][4]=l10r;
Curvature_Higgs_L4[4][8][8][4]=l9;
Curvature_Higgs_L4[4][9][4][9]=2*lh;
Curvature_Higgs_L4[4][9][9][4]=2*lh;
Curvature_Higgs_L4[5][0][0][2]=l6i;
Curvature_Higgs_L4[5][0][0][5]=l1;
Curvature_Higgs_L4[5][0][0][7]=l6r;
Curvature_Higgs_L4[5][0][2][0]=l6i;
Curvature_Higgs_L4[5][0][2][2]=l5i;
Curvature_Higgs_L4[5][0][2][5]=l6r;
Curvature_Higgs_L4[5][0][2][7]=l5r;
Curvature_Higgs_L4[5][0][5][0]=l1;
Curvature_Higgs_L4[5][0][5][2]=l6r;
Curvature_Higgs_L4[5][0][5][7]=-l6i;
Curvature_Higgs_L4[5][0][7][0]=l6r;
Curvature_Higgs_L4[5][0][7][2]=l5r;
Curvature_Higgs_L4[5][0][7][5]=-l6i;
Curvature_Higgs_L4[5][0][7][7]=-l5i;
Curvature_Higgs_L4[5][1][1][2]=l6i;
Curvature_Higgs_L4[5][1][1][5]=l1;
Curvature_Higgs_L4[5][1][1][7]=l6r;
Curvature_Higgs_L4[5][1][2][1]=l6i;
Curvature_Higgs_L4[5][1][2][3]=l5i/2.;
Curvature_Higgs_L4[5][1][2][8]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][1][3][2]=l5i/2.;
Curvature_Higgs_L4[5][1][3][5]=l6r;
Curvature_Higgs_L4[5][1][3][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][1][5][1]=l1;
Curvature_Higgs_L4[5][1][5][3]=l6r;
Curvature_Higgs_L4[5][1][5][8]=-l6i;
Curvature_Higgs_L4[5][1][7][1]=l6r;
Curvature_Higgs_L4[5][1][7][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][1][7][8]=-0.5*l5i;
Curvature_Higgs_L4[5][1][8][2]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][1][8][5]=-l6i;
Curvature_Higgs_L4[5][1][8][7]=-0.5*l5i;
Curvature_Higgs_L4[5][2][0][0]=l6i;
Curvature_Higgs_L4[5][2][0][2]=l5i;
Curvature_Higgs_L4[5][2][0][5]=l6r;
Curvature_Higgs_L4[5][2][0][7]=l5r;
Curvature_Higgs_L4[5][2][1][1]=l6i;
Curvature_Higgs_L4[5][2][1][3]=l5i/2.;
Curvature_Higgs_L4[5][2][1][8]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][2][2][0]=l5i;
Curvature_Higgs_L4[5][2][2][2]=3*l7i;
Curvature_Higgs_L4[5][2][2][5]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[5][2][2][7]=l7r;
Curvature_Higgs_L4[5][2][3][1]=l5i/2.;
Curvature_Higgs_L4[5][2][3][3]=l7i;
Curvature_Higgs_L4[5][2][3][6]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[5][2][4][4]=l10i;
Curvature_Higgs_L4[5][2][5][0]=l6r;
Curvature_Higgs_L4[5][2][5][2]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[5][2][5][5]=3*l6i;
Curvature_Higgs_L4[5][2][5][7]=l5i;
Curvature_Higgs_L4[5][2][6][3]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[5][2][6][6]=l6i;
Curvature_Higgs_L4[5][2][6][8]=l5i/2.;
Curvature_Higgs_L4[5][2][7][0]=l5r;
Curvature_Higgs_L4[5][2][7][2]=l7r;
Curvature_Higgs_L4[5][2][7][5]=l5i;
Curvature_Higgs_L4[5][2][7][7]=l7i;
Curvature_Higgs_L4[5][2][8][1]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][2][8][6]=l5i/2.;
Curvature_Higgs_L4[5][2][8][8]=l7i;
Curvature_Higgs_L4[5][2][9][9]=l10i;
Curvature_Higgs_L4[5][3][1][2]=l5i/2.;
Curvature_Higgs_L4[5][3][1][5]=l6r;
Curvature_Higgs_L4[5][3][1][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][3][2][1]=l5i/2.;
Curvature_Higgs_L4[5][3][2][3]=l7i;
Curvature_Higgs_L4[5][3][2][6]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[5][3][3][2]=l7i;
Curvature_Higgs_L4[5][3][3][5]=l3;
Curvature_Higgs_L4[5][3][3][7]=l7r;
Curvature_Higgs_L4[5][3][5][1]=l6r;
Curvature_Higgs_L4[5][3][5][3]=l3;
Curvature_Higgs_L4[5][3][5][6]=l6i;
Curvature_Higgs_L4[5][3][6][2]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[5][3][6][5]=l6i;
Curvature_Higgs_L4[5][3][6][7]=l5i/2.;
Curvature_Higgs_L4[5][3][7][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][3][7][3]=l7r;
Curvature_Higgs_L4[5][3][7][6]=l5i/2.;
Curvature_Higgs_L4[5][4][2][4]=l10i;
Curvature_Higgs_L4[5][4][4][2]=l10i;
Curvature_Higgs_L4[5][4][4][5]=l8;
Curvature_Higgs_L4[5][4][4][7]=l10r;
Curvature_Higgs_L4[5][4][5][4]=l8;
Curvature_Higgs_L4[5][4][7][4]=l10r;
Curvature_Higgs_L4[5][5][0][0]=l1;
Curvature_Higgs_L4[5][5][0][2]=l6r;
Curvature_Higgs_L4[5][5][0][7]=-l6i;
Curvature_Higgs_L4[5][5][1][1]=l1;
Curvature_Higgs_L4[5][5][1][3]=l6r;
Curvature_Higgs_L4[5][5][1][8]=-l6i;
Curvature_Higgs_L4[5][5][2][0]=l6r;
Curvature_Higgs_L4[5][5][2][2]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[5][5][2][5]=3*l6i;
Curvature_Higgs_L4[5][5][2][7]=l5i;
Curvature_Higgs_L4[5][5][3][1]=l6r;
Curvature_Higgs_L4[5][5][3][3]=l3;
Curvature_Higgs_L4[5][5][3][6]=l6i;
Curvature_Higgs_L4[5][5][4][4]=l8;
Curvature_Higgs_L4[5][5][5][2]=3*l6i;
Curvature_Higgs_L4[5][5][5][5]=3*l1;
Curvature_Higgs_L4[5][5][5][7]=3*l6r;
Curvature_Higgs_L4[5][5][6][3]=l6i;
Curvature_Higgs_L4[5][5][6][6]=l1;
Curvature_Higgs_L4[5][5][6][8]=l6r;
Curvature_Higgs_L4[5][5][7][0]=-l6i;
Curvature_Higgs_L4[5][5][7][2]=l5i;
Curvature_Higgs_L4[5][5][7][5]=3*l6r;
Curvature_Higgs_L4[5][5][7][7]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[5][5][8][1]=-l6i;
Curvature_Higgs_L4[5][5][8][6]=l6r;
Curvature_Higgs_L4[5][5][8][8]=l3;
Curvature_Higgs_L4[5][5][9][9]=l8;
Curvature_Higgs_L4[5][6][2][3]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[5][6][2][6]=l6i;
Curvature_Higgs_L4[5][6][2][8]=l5i/2.;
Curvature_Higgs_L4[5][6][3][2]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[5][6][3][5]=l6i;
Curvature_Higgs_L4[5][6][3][7]=l5i/2.;
Curvature_Higgs_L4[5][6][5][3]=l6i;
Curvature_Higgs_L4[5][6][5][6]=l1;
Curvature_Higgs_L4[5][6][5][8]=l6r;
Curvature_Higgs_L4[5][6][6][2]=l6i;
Curvature_Higgs_L4[5][6][6][5]=l1;
Curvature_Higgs_L4[5][6][6][7]=l6r;
Curvature_Higgs_L4[5][6][7][3]=l5i/2.;
Curvature_Higgs_L4[5][6][7][6]=l6r;
Curvature_Higgs_L4[5][6][7][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][6][8][2]=l5i/2.;
Curvature_Higgs_L4[5][6][8][5]=l6r;
Curvature_Higgs_L4[5][6][8][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][7][0][0]=l6r;
Curvature_Higgs_L4[5][7][0][2]=l5r;
Curvature_Higgs_L4[5][7][0][5]=-l6i;
Curvature_Higgs_L4[5][7][0][7]=-l5i;
Curvature_Higgs_L4[5][7][1][1]=l6r;
Curvature_Higgs_L4[5][7][1][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][7][1][8]=-0.5*l5i;
Curvature_Higgs_L4[5][7][2][0]=l5r;
Curvature_Higgs_L4[5][7][2][2]=l7r;
Curvature_Higgs_L4[5][7][2][5]=l5i;
Curvature_Higgs_L4[5][7][2][7]=l7i;
Curvature_Higgs_L4[5][7][3][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][7][3][3]=l7r;
Curvature_Higgs_L4[5][7][3][6]=l5i/2.;
Curvature_Higgs_L4[5][7][4][4]=l10r;
Curvature_Higgs_L4[5][7][5][0]=-l6i;
Curvature_Higgs_L4[5][7][5][2]=l5i;
Curvature_Higgs_L4[5][7][5][5]=3*l6r;
Curvature_Higgs_L4[5][7][5][7]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[5][7][6][3]=l5i/2.;
Curvature_Higgs_L4[5][7][6][6]=l6r;
Curvature_Higgs_L4[5][7][6][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][7][7][0]=-l5i;
Curvature_Higgs_L4[5][7][7][2]=l7i;
Curvature_Higgs_L4[5][7][7][5]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[5][7][7][7]=3*l7r;
Curvature_Higgs_L4[5][7][8][1]=-0.5*l5i;
Curvature_Higgs_L4[5][7][8][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][7][8][8]=l7r;
Curvature_Higgs_L4[5][7][9][9]=l10r;
Curvature_Higgs_L4[5][8][1][2]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][8][1][5]=-l6i;
Curvature_Higgs_L4[5][8][1][7]=-0.5*l5i;
Curvature_Higgs_L4[5][8][2][1]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][8][2][6]=l5i/2.;
Curvature_Higgs_L4[5][8][2][8]=l7i;
Curvature_Higgs_L4[5][8][5][1]=-l6i;
Curvature_Higgs_L4[5][8][5][6]=l6r;
Curvature_Higgs_L4[5][8][5][8]=l3;
Curvature_Higgs_L4[5][8][6][2]=l5i/2.;
Curvature_Higgs_L4[5][8][6][5]=l6r;
Curvature_Higgs_L4[5][8][6][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][8][7][1]=-0.5*l5i;
Curvature_Higgs_L4[5][8][7][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[5][8][7][8]=l7r;
Curvature_Higgs_L4[5][8][8][2]=l7i;
Curvature_Higgs_L4[5][8][8][5]=l3;
Curvature_Higgs_L4[5][8][8][7]=l7r;
Curvature_Higgs_L4[5][9][2][9]=l10i;
Curvature_Higgs_L4[5][9][5][9]=l8;
Curvature_Higgs_L4[5][9][7][9]=l10r;
Curvature_Higgs_L4[5][9][9][2]=l10i;
Curvature_Higgs_L4[5][9][9][5]=l8;
Curvature_Higgs_L4[5][9][9][7]=l10r;
Curvature_Higgs_L4[6][0][0][3]=l6i;
Curvature_Higgs_L4[6][0][0][6]=l1;
Curvature_Higgs_L4[6][0][0][8]=l6r;
Curvature_Higgs_L4[6][0][2][3]=l5i/2.;
Curvature_Higgs_L4[6][0][2][6]=l6r;
Curvature_Higgs_L4[6][0][2][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][0][3][0]=l6i;
Curvature_Higgs_L4[6][0][3][2]=l5i/2.;
Curvature_Higgs_L4[6][0][3][7]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][0][6][0]=l1;
Curvature_Higgs_L4[6][0][6][2]=l6r;
Curvature_Higgs_L4[6][0][6][7]=-l6i;
Curvature_Higgs_L4[6][0][7][3]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][0][7][6]=-l6i;
Curvature_Higgs_L4[6][0][7][8]=-0.5*l5i;
Curvature_Higgs_L4[6][0][8][0]=l6r;
Curvature_Higgs_L4[6][0][8][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][0][8][7]=-0.5*l5i;
Curvature_Higgs_L4[6][1][1][3]=l6i;
Curvature_Higgs_L4[6][1][1][6]=l1;
Curvature_Higgs_L4[6][1][1][8]=l6r;
Curvature_Higgs_L4[6][1][3][1]=l6i;
Curvature_Higgs_L4[6][1][3][3]=l5i;
Curvature_Higgs_L4[6][1][3][6]=l6r;
Curvature_Higgs_L4[6][1][3][8]=l5r;
Curvature_Higgs_L4[6][1][6][1]=l1;
Curvature_Higgs_L4[6][1][6][3]=l6r;
Curvature_Higgs_L4[6][1][6][8]=-l6i;
Curvature_Higgs_L4[6][1][8][1]=l6r;
Curvature_Higgs_L4[6][1][8][3]=l5r;
Curvature_Higgs_L4[6][1][8][6]=-l6i;
Curvature_Higgs_L4[6][1][8][8]=-l5i;
Curvature_Higgs_L4[6][2][0][3]=l5i/2.;
Curvature_Higgs_L4[6][2][0][6]=l6r;
Curvature_Higgs_L4[6][2][0][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][2][2][3]=l7i;
Curvature_Higgs_L4[6][2][2][6]=l3;
Curvature_Higgs_L4[6][2][2][8]=l7r;
Curvature_Higgs_L4[6][2][3][0]=l5i/2.;
Curvature_Higgs_L4[6][2][3][2]=l7i;
Curvature_Higgs_L4[6][2][3][5]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[6][2][5][3]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[6][2][5][6]=l6i;
Curvature_Higgs_L4[6][2][5][8]=l5i/2.;
Curvature_Higgs_L4[6][2][6][0]=l6r;
Curvature_Higgs_L4[6][2][6][2]=l3;
Curvature_Higgs_L4[6][2][6][5]=l6i;
Curvature_Higgs_L4[6][2][8][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][2][8][2]=l7r;
Curvature_Higgs_L4[6][2][8][5]=l5i/2.;
Curvature_Higgs_L4[6][3][0][0]=l6i;
Curvature_Higgs_L4[6][3][0][2]=l5i/2.;
Curvature_Higgs_L4[6][3][0][7]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][3][1][1]=l6i;
Curvature_Higgs_L4[6][3][1][3]=l5i;
Curvature_Higgs_L4[6][3][1][6]=l6r;
Curvature_Higgs_L4[6][3][1][8]=l5r;
Curvature_Higgs_L4[6][3][2][0]=l5i/2.;
Curvature_Higgs_L4[6][3][2][2]=l7i;
Curvature_Higgs_L4[6][3][2][5]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[6][3][3][1]=l5i;
Curvature_Higgs_L4[6][3][3][3]=3*l7i;
Curvature_Higgs_L4[6][3][3][6]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[6][3][3][8]=l7r;
Curvature_Higgs_L4[6][3][4][4]=l10i;
Curvature_Higgs_L4[6][3][5][2]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[6][3][5][5]=l6i;
Curvature_Higgs_L4[6][3][5][7]=l5i/2.;
Curvature_Higgs_L4[6][3][6][1]=l6r;
Curvature_Higgs_L4[6][3][6][3]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[6][3][6][6]=3*l6i;
Curvature_Higgs_L4[6][3][6][8]=l5i;
Curvature_Higgs_L4[6][3][7][0]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][3][7][5]=l5i/2.;
Curvature_Higgs_L4[6][3][7][7]=l7i;
Curvature_Higgs_L4[6][3][8][1]=l5r;
Curvature_Higgs_L4[6][3][8][3]=l7r;
Curvature_Higgs_L4[6][3][8][6]=l5i;
Curvature_Higgs_L4[6][3][8][8]=l7i;
Curvature_Higgs_L4[6][3][9][9]=l10i;
Curvature_Higgs_L4[6][4][3][4]=l10i;
Curvature_Higgs_L4[6][4][4][3]=l10i;
Curvature_Higgs_L4[6][4][4][6]=l8;
Curvature_Higgs_L4[6][4][4][8]=l10r;
Curvature_Higgs_L4[6][4][6][4]=l8;
Curvature_Higgs_L4[6][4][8][4]=l10r;
Curvature_Higgs_L4[6][5][2][3]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[6][5][2][6]=l6i;
Curvature_Higgs_L4[6][5][2][8]=l5i/2.;
Curvature_Higgs_L4[6][5][3][2]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[6][5][3][5]=l6i;
Curvature_Higgs_L4[6][5][3][7]=l5i/2.;
Curvature_Higgs_L4[6][5][5][3]=l6i;
Curvature_Higgs_L4[6][5][5][6]=l1;
Curvature_Higgs_L4[6][5][5][8]=l6r;
Curvature_Higgs_L4[6][5][6][2]=l6i;
Curvature_Higgs_L4[6][5][6][5]=l1;
Curvature_Higgs_L4[6][5][6][7]=l6r;
Curvature_Higgs_L4[6][5][7][3]=l5i/2.;
Curvature_Higgs_L4[6][5][7][6]=l6r;
Curvature_Higgs_L4[6][5][7][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][5][8][2]=l5i/2.;
Curvature_Higgs_L4[6][5][8][5]=l6r;
Curvature_Higgs_L4[6][5][8][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][6][0][0]=l1;
Curvature_Higgs_L4[6][6][0][2]=l6r;
Curvature_Higgs_L4[6][6][0][7]=-l6i;
Curvature_Higgs_L4[6][6][1][1]=l1;
Curvature_Higgs_L4[6][6][1][3]=l6r;
Curvature_Higgs_L4[6][6][1][8]=-l6i;
Curvature_Higgs_L4[6][6][2][0]=l6r;
Curvature_Higgs_L4[6][6][2][2]=l3;
Curvature_Higgs_L4[6][6][2][5]=l6i;
Curvature_Higgs_L4[6][6][3][1]=l6r;
Curvature_Higgs_L4[6][6][3][3]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[6][6][3][6]=3*l6i;
Curvature_Higgs_L4[6][6][3][8]=l5i;
Curvature_Higgs_L4[6][6][4][4]=l8;
Curvature_Higgs_L4[6][6][5][2]=l6i;
Curvature_Higgs_L4[6][6][5][5]=l1;
Curvature_Higgs_L4[6][6][5][7]=l6r;
Curvature_Higgs_L4[6][6][6][3]=3*l6i;
Curvature_Higgs_L4[6][6][6][6]=3*l1;
Curvature_Higgs_L4[6][6][6][8]=3*l6r;
Curvature_Higgs_L4[6][6][7][0]=-l6i;
Curvature_Higgs_L4[6][6][7][5]=l6r;
Curvature_Higgs_L4[6][6][7][7]=l3;
Curvature_Higgs_L4[6][6][8][1]=-l6i;
Curvature_Higgs_L4[6][6][8][3]=l5i;
Curvature_Higgs_L4[6][6][8][6]=3*l6r;
Curvature_Higgs_L4[6][6][8][8]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[6][6][9][9]=l8;
Curvature_Higgs_L4[6][7][0][3]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][7][0][6]=-l6i;
Curvature_Higgs_L4[6][7][0][8]=-0.5*l5i;
Curvature_Higgs_L4[6][7][3][0]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][7][3][5]=l5i/2.;
Curvature_Higgs_L4[6][7][3][7]=l7i;
Curvature_Higgs_L4[6][7][5][3]=l5i/2.;
Curvature_Higgs_L4[6][7][5][6]=l6r;
Curvature_Higgs_L4[6][7][5][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][7][6][0]=-l6i;
Curvature_Higgs_L4[6][7][6][5]=l6r;
Curvature_Higgs_L4[6][7][6][7]=l3;
Curvature_Higgs_L4[6][7][7][3]=l7i;
Curvature_Higgs_L4[6][7][7][6]=l3;
Curvature_Higgs_L4[6][7][7][8]=l7r;
Curvature_Higgs_L4[6][7][8][0]=-0.5*l5i;
Curvature_Higgs_L4[6][7][8][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][7][8][7]=l7r;
Curvature_Higgs_L4[6][8][0][0]=l6r;
Curvature_Higgs_L4[6][8][0][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][8][0][7]=-0.5*l5i;
Curvature_Higgs_L4[6][8][1][1]=l6r;
Curvature_Higgs_L4[6][8][1][3]=l5r;
Curvature_Higgs_L4[6][8][1][6]=-l6i;
Curvature_Higgs_L4[6][8][1][8]=-l5i;
Curvature_Higgs_L4[6][8][2][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][8][2][2]=l7r;
Curvature_Higgs_L4[6][8][2][5]=l5i/2.;
Curvature_Higgs_L4[6][8][3][1]=l5r;
Curvature_Higgs_L4[6][8][3][3]=l7r;
Curvature_Higgs_L4[6][8][3][6]=l5i;
Curvature_Higgs_L4[6][8][3][8]=l7i;
Curvature_Higgs_L4[6][8][4][4]=l10r;
Curvature_Higgs_L4[6][8][5][2]=l5i/2.;
Curvature_Higgs_L4[6][8][5][5]=l6r;
Curvature_Higgs_L4[6][8][5][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][8][6][1]=-l6i;
Curvature_Higgs_L4[6][8][6][3]=l5i;
Curvature_Higgs_L4[6][8][6][6]=3*l6r;
Curvature_Higgs_L4[6][8][6][8]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[6][8][7][0]=-0.5*l5i;
Curvature_Higgs_L4[6][8][7][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[6][8][7][7]=l7r;
Curvature_Higgs_L4[6][8][8][1]=-l5i;
Curvature_Higgs_L4[6][8][8][3]=l7i;
Curvature_Higgs_L4[6][8][8][6]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[6][8][8][8]=3*l7r;
Curvature_Higgs_L4[6][8][9][9]=l10r;
Curvature_Higgs_L4[6][9][3][9]=l10i;
Curvature_Higgs_L4[6][9][6][9]=l8;
Curvature_Higgs_L4[6][9][8][9]=l10r;
Curvature_Higgs_L4[6][9][9][3]=l10i;
Curvature_Higgs_L4[6][9][9][6]=l8;
Curvature_Higgs_L4[6][9][9][8]=l10r;
Curvature_Higgs_L4[7][0][0][0]=-3*l6i;
Curvature_Higgs_L4[7][0][0][2]=-l5i;
Curvature_Higgs_L4[7][0][0][5]=l6r;
Curvature_Higgs_L4[7][0][0][7]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[7][0][1][1]=-l6i;
Curvature_Higgs_L4[7][0][1][3]=-0.5*l5i;
Curvature_Higgs_L4[7][0][1][8]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[7][0][2][0]=-l5i;
Curvature_Higgs_L4[7][0][2][2]=-l7i;
Curvature_Higgs_L4[7][0][2][5]=l5r;
Curvature_Higgs_L4[7][0][2][7]=l7r;
Curvature_Higgs_L4[7][0][3][1]=-0.5*l5i;
Curvature_Higgs_L4[7][0][3][3]=-l7i;
Curvature_Higgs_L4[7][0][3][6]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][0][4][4]=-l10i;
Curvature_Higgs_L4[7][0][5][0]=l6r;
Curvature_Higgs_L4[7][0][5][2]=l5r;
Curvature_Higgs_L4[7][0][5][5]=-l6i;
Curvature_Higgs_L4[7][0][5][7]=-l5i;
Curvature_Higgs_L4[7][0][6][3]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][0][6][6]=-l6i;
Curvature_Higgs_L4[7][0][6][8]=-0.5*l5i;
Curvature_Higgs_L4[7][0][7][0]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[7][0][7][2]=l7r;
Curvature_Higgs_L4[7][0][7][5]=-l5i;
Curvature_Higgs_L4[7][0][7][7]=-3*l7i;
Curvature_Higgs_L4[7][0][8][1]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[7][0][8][6]=-0.5*l5i;
Curvature_Higgs_L4[7][0][8][8]=-l7i;
Curvature_Higgs_L4[7][0][9][9]=-l10i;
Curvature_Higgs_L4[7][1][0][1]=-l6i;
Curvature_Higgs_L4[7][1][0][3]=-0.5*l5i;
Curvature_Higgs_L4[7][1][0][8]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[7][1][1][0]=-l6i;
Curvature_Higgs_L4[7][1][1][5]=l6r;
Curvature_Higgs_L4[7][1][1][7]=l3;
Curvature_Higgs_L4[7][1][3][0]=-0.5*l5i;
Curvature_Higgs_L4[7][1][3][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][1][3][7]=l7r;
Curvature_Higgs_L4[7][1][5][1]=l6r;
Curvature_Higgs_L4[7][1][5][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][1][5][8]=-0.5*l5i;
Curvature_Higgs_L4[7][1][7][1]=l3;
Curvature_Higgs_L4[7][1][7][3]=l7r;
Curvature_Higgs_L4[7][1][7][8]=-l7i;
Curvature_Higgs_L4[7][1][8][0]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[7][1][8][5]=-0.5*l5i;
Curvature_Higgs_L4[7][1][8][7]=-l7i;
Curvature_Higgs_L4[7][2][0][0]=-l5i;
Curvature_Higgs_L4[7][2][0][2]=-l7i;
Curvature_Higgs_L4[7][2][0][5]=l5r;
Curvature_Higgs_L4[7][2][0][7]=l7r;
Curvature_Higgs_L4[7][2][2][0]=-l7i;
Curvature_Higgs_L4[7][2][2][5]=l7r;
Curvature_Higgs_L4[7][2][2][7]=l2;
Curvature_Higgs_L4[7][2][5][0]=l5r;
Curvature_Higgs_L4[7][2][5][2]=l7r;
Curvature_Higgs_L4[7][2][5][5]=l5i;
Curvature_Higgs_L4[7][2][5][7]=l7i;
Curvature_Higgs_L4[7][2][7][0]=l7r;
Curvature_Higgs_L4[7][2][7][2]=l2;
Curvature_Higgs_L4[7][2][7][5]=l7i;
Curvature_Higgs_L4[7][3][0][1]=-0.5*l5i;
Curvature_Higgs_L4[7][3][0][3]=-l7i;
Curvature_Higgs_L4[7][3][0][6]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][3][1][0]=-0.5*l5i;
Curvature_Higgs_L4[7][3][1][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][3][1][7]=l7r;
Curvature_Higgs_L4[7][3][3][0]=-l7i;
Curvature_Higgs_L4[7][3][3][5]=l7r;
Curvature_Higgs_L4[7][3][3][7]=l2;
Curvature_Higgs_L4[7][3][5][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][3][5][3]=l7r;
Curvature_Higgs_L4[7][3][5][6]=l5i/2.;
Curvature_Higgs_L4[7][3][6][0]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][3][6][5]=l5i/2.;
Curvature_Higgs_L4[7][3][6][7]=l7i;
Curvature_Higgs_L4[7][3][7][1]=l7r;
Curvature_Higgs_L4[7][3][7][3]=l2;
Curvature_Higgs_L4[7][3][7][6]=l7i;
Curvature_Higgs_L4[7][4][0][4]=-l10i;
Curvature_Higgs_L4[7][4][4][0]=-l10i;
Curvature_Higgs_L4[7][4][4][5]=l10r;
Curvature_Higgs_L4[7][4][4][7]=l9;
Curvature_Higgs_L4[7][4][5][4]=l10r;
Curvature_Higgs_L4[7][4][7][4]=l9;
Curvature_Higgs_L4[7][5][0][0]=l6r;
Curvature_Higgs_L4[7][5][0][2]=l5r;
Curvature_Higgs_L4[7][5][0][5]=-l6i;
Curvature_Higgs_L4[7][5][0][7]=-l5i;
Curvature_Higgs_L4[7][5][1][1]=l6r;
Curvature_Higgs_L4[7][5][1][3]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][5][1][8]=-0.5*l5i;
Curvature_Higgs_L4[7][5][2][0]=l5r;
Curvature_Higgs_L4[7][5][2][2]=l7r;
Curvature_Higgs_L4[7][5][2][5]=l5i;
Curvature_Higgs_L4[7][5][2][7]=l7i;
Curvature_Higgs_L4[7][5][3][1]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][5][3][3]=l7r;
Curvature_Higgs_L4[7][5][3][6]=l5i/2.;
Curvature_Higgs_L4[7][5][4][4]=l10r;
Curvature_Higgs_L4[7][5][5][0]=-l6i;
Curvature_Higgs_L4[7][5][5][2]=l5i;
Curvature_Higgs_L4[7][5][5][5]=3*l6r;
Curvature_Higgs_L4[7][5][5][7]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[7][5][6][3]=l5i/2.;
Curvature_Higgs_L4[7][5][6][6]=l6r;
Curvature_Higgs_L4[7][5][6][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][5][7][0]=-l5i;
Curvature_Higgs_L4[7][5][7][2]=l7i;
Curvature_Higgs_L4[7][5][7][5]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[7][5][7][7]=3*l7r;
Curvature_Higgs_L4[7][5][8][1]=-0.5*l5i;
Curvature_Higgs_L4[7][5][8][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][5][8][8]=l7r;
Curvature_Higgs_L4[7][5][9][9]=l10r;
Curvature_Higgs_L4[7][6][0][3]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][6][0][6]=-l6i;
Curvature_Higgs_L4[7][6][0][8]=-0.5*l5i;
Curvature_Higgs_L4[7][6][3][0]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][6][3][5]=l5i/2.;
Curvature_Higgs_L4[7][6][3][7]=l7i;
Curvature_Higgs_L4[7][6][5][3]=l5i/2.;
Curvature_Higgs_L4[7][6][5][6]=l6r;
Curvature_Higgs_L4[7][6][5][8]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][6][6][0]=-l6i;
Curvature_Higgs_L4[7][6][6][5]=l6r;
Curvature_Higgs_L4[7][6][6][7]=l3;
Curvature_Higgs_L4[7][6][7][3]=l7i;
Curvature_Higgs_L4[7][6][7][6]=l3;
Curvature_Higgs_L4[7][6][7][8]=l7r;
Curvature_Higgs_L4[7][6][8][0]=-0.5*l5i;
Curvature_Higgs_L4[7][6][8][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][6][8][7]=l7r;
Curvature_Higgs_L4[7][7][0][0]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[7][7][0][2]=l7r;
Curvature_Higgs_L4[7][7][0][5]=-l5i;
Curvature_Higgs_L4[7][7][0][7]=-3*l7i;
Curvature_Higgs_L4[7][7][1][1]=l3;
Curvature_Higgs_L4[7][7][1][3]=l7r;
Curvature_Higgs_L4[7][7][1][8]=-l7i;
Curvature_Higgs_L4[7][7][2][0]=l7r;
Curvature_Higgs_L4[7][7][2][2]=l2;
Curvature_Higgs_L4[7][7][2][5]=l7i;
Curvature_Higgs_L4[7][7][3][1]=l7r;
Curvature_Higgs_L4[7][7][3][3]=l2;
Curvature_Higgs_L4[7][7][3][6]=l7i;
Curvature_Higgs_L4[7][7][4][4]=l9;
Curvature_Higgs_L4[7][7][5][0]=-l5i;
Curvature_Higgs_L4[7][7][5][2]=l7i;
Curvature_Higgs_L4[7][7][5][5]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[7][7][5][7]=3*l7r;
Curvature_Higgs_L4[7][7][6][3]=l7i;
Curvature_Higgs_L4[7][7][6][6]=l3;
Curvature_Higgs_L4[7][7][6][8]=l7r;
Curvature_Higgs_L4[7][7][7][0]=-3*l7i;
Curvature_Higgs_L4[7][7][7][5]=3*l7r;
Curvature_Higgs_L4[7][7][7][7]=3*l2;
Curvature_Higgs_L4[7][7][8][1]=-l7i;
Curvature_Higgs_L4[7][7][8][6]=l7r;
Curvature_Higgs_L4[7][7][8][8]=l2;
Curvature_Higgs_L4[7][7][9][9]=l9;
Curvature_Higgs_L4[7][8][0][1]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[7][8][0][6]=-0.5*l5i;
Curvature_Higgs_L4[7][8][0][8]=-l7i;
Curvature_Higgs_L4[7][8][1][0]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[7][8][1][5]=-0.5*l5i;
Curvature_Higgs_L4[7][8][1][7]=-l7i;
Curvature_Higgs_L4[7][8][5][1]=-0.5*l5i;
Curvature_Higgs_L4[7][8][5][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][8][5][8]=l7r;
Curvature_Higgs_L4[7][8][6][0]=-0.5*l5i;
Curvature_Higgs_L4[7][8][6][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[7][8][6][7]=l7r;
Curvature_Higgs_L4[7][8][7][1]=-l7i;
Curvature_Higgs_L4[7][8][7][6]=l7r;
Curvature_Higgs_L4[7][8][7][8]=l2;
Curvature_Higgs_L4[7][8][8][0]=-l7i;
Curvature_Higgs_L4[7][8][8][5]=l7r;
Curvature_Higgs_L4[7][8][8][7]=l2;
Curvature_Higgs_L4[7][9][0][9]=-l10i;
Curvature_Higgs_L4[7][9][5][9]=l10r;
Curvature_Higgs_L4[7][9][7][9]=l9;
Curvature_Higgs_L4[7][9][9][0]=-l10i;
Curvature_Higgs_L4[7][9][9][5]=l10r;
Curvature_Higgs_L4[7][9][9][7]=l9;
Curvature_Higgs_L4[8][0][0][1]=-l6i;
Curvature_Higgs_L4[8][0][0][6]=l6r;
Curvature_Higgs_L4[8][0][0][8]=l3;
Curvature_Higgs_L4[8][0][1][0]=-l6i;
Curvature_Higgs_L4[8][0][1][2]=-0.5*l5i;
Curvature_Higgs_L4[8][0][1][7]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[8][0][2][1]=-0.5*l5i;
Curvature_Higgs_L4[8][0][2][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][0][2][8]=l7r;
Curvature_Higgs_L4[8][0][6][0]=l6r;
Curvature_Higgs_L4[8][0][6][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][0][6][7]=-0.5*l5i;
Curvature_Higgs_L4[8][0][7][1]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[8][0][7][6]=-0.5*l5i;
Curvature_Higgs_L4[8][0][7][8]=-l7i;
Curvature_Higgs_L4[8][0][8][0]=l3;
Curvature_Higgs_L4[8][0][8][2]=l7r;
Curvature_Higgs_L4[8][0][8][7]=-l7i;
Curvature_Higgs_L4[8][1][0][0]=-l6i;
Curvature_Higgs_L4[8][1][0][2]=-0.5*l5i;
Curvature_Higgs_L4[8][1][0][7]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[8][1][1][1]=-3*l6i;
Curvature_Higgs_L4[8][1][1][3]=-l5i;
Curvature_Higgs_L4[8][1][1][6]=l6r;
Curvature_Higgs_L4[8][1][1][8]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[8][1][2][0]=-0.5*l5i;
Curvature_Higgs_L4[8][1][2][2]=-l7i;
Curvature_Higgs_L4[8][1][2][5]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][1][3][1]=-l5i;
Curvature_Higgs_L4[8][1][3][3]=-l7i;
Curvature_Higgs_L4[8][1][3][6]=l5r;
Curvature_Higgs_L4[8][1][3][8]=l7r;
Curvature_Higgs_L4[8][1][4][4]=-l10i;
Curvature_Higgs_L4[8][1][5][2]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][1][5][5]=-l6i;
Curvature_Higgs_L4[8][1][5][7]=-0.5*l5i;
Curvature_Higgs_L4[8][1][6][1]=l6r;
Curvature_Higgs_L4[8][1][6][3]=l5r;
Curvature_Higgs_L4[8][1][6][6]=-l6i;
Curvature_Higgs_L4[8][1][6][8]=-l5i;
Curvature_Higgs_L4[8][1][7][0]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[8][1][7][5]=-0.5*l5i;
Curvature_Higgs_L4[8][1][7][7]=-l7i;
Curvature_Higgs_L4[8][1][8][1]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[8][1][8][3]=l7r;
Curvature_Higgs_L4[8][1][8][6]=-l5i;
Curvature_Higgs_L4[8][1][8][8]=-3*l7i;
Curvature_Higgs_L4[8][1][9][9]=-l10i;
Curvature_Higgs_L4[8][2][0][1]=-0.5*l5i;
Curvature_Higgs_L4[8][2][0][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][2][0][8]=l7r;
Curvature_Higgs_L4[8][2][1][0]=-0.5*l5i;
Curvature_Higgs_L4[8][2][1][2]=-l7i;
Curvature_Higgs_L4[8][2][1][5]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][2][2][1]=-l7i;
Curvature_Higgs_L4[8][2][2][6]=l7r;
Curvature_Higgs_L4[8][2][2][8]=l2;
Curvature_Higgs_L4[8][2][5][1]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][2][5][6]=l5i/2.;
Curvature_Higgs_L4[8][2][5][8]=l7i;
Curvature_Higgs_L4[8][2][6][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][2][6][2]=l7r;
Curvature_Higgs_L4[8][2][6][5]=l5i/2.;
Curvature_Higgs_L4[8][2][8][0]=l7r;
Curvature_Higgs_L4[8][2][8][2]=l2;
Curvature_Higgs_L4[8][2][8][5]=l7i;
Curvature_Higgs_L4[8][3][1][1]=-l5i;
Curvature_Higgs_L4[8][3][1][3]=-l7i;
Curvature_Higgs_L4[8][3][1][6]=l5r;
Curvature_Higgs_L4[8][3][1][8]=l7r;
Curvature_Higgs_L4[8][3][3][1]=-l7i;
Curvature_Higgs_L4[8][3][3][6]=l7r;
Curvature_Higgs_L4[8][3][3][8]=l2;
Curvature_Higgs_L4[8][3][6][1]=l5r;
Curvature_Higgs_L4[8][3][6][3]=l7r;
Curvature_Higgs_L4[8][3][6][6]=l5i;
Curvature_Higgs_L4[8][3][6][8]=l7i;
Curvature_Higgs_L4[8][3][8][1]=l7r;
Curvature_Higgs_L4[8][3][8][3]=l2;
Curvature_Higgs_L4[8][3][8][6]=l7i;
Curvature_Higgs_L4[8][4][1][4]=-l10i;
Curvature_Higgs_L4[8][4][4][1]=-l10i;
Curvature_Higgs_L4[8][4][4][6]=l10r;
Curvature_Higgs_L4[8][4][4][8]=l9;
Curvature_Higgs_L4[8][4][6][4]=l10r;
Curvature_Higgs_L4[8][4][8][4]=l9;
Curvature_Higgs_L4[8][5][1][2]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][5][1][5]=-l6i;
Curvature_Higgs_L4[8][5][1][7]=-0.5*l5i;
Curvature_Higgs_L4[8][5][2][1]=(-4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][5][2][6]=l5i/2.;
Curvature_Higgs_L4[8][5][2][8]=l7i;
Curvature_Higgs_L4[8][5][5][1]=-l6i;
Curvature_Higgs_L4[8][5][5][6]=l6r;
Curvature_Higgs_L4[8][5][5][8]=l3;
Curvature_Higgs_L4[8][5][6][2]=l5i/2.;
Curvature_Higgs_L4[8][5][6][5]=l6r;
Curvature_Higgs_L4[8][5][6][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][5][7][1]=-0.5*l5i;
Curvature_Higgs_L4[8][5][7][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][5][7][8]=l7r;
Curvature_Higgs_L4[8][5][8][2]=l7i;
Curvature_Higgs_L4[8][5][8][5]=l3;
Curvature_Higgs_L4[8][5][8][7]=l7r;
Curvature_Higgs_L4[8][6][0][0]=l6r;
Curvature_Higgs_L4[8][6][0][2]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][6][0][7]=-0.5*l5i;
Curvature_Higgs_L4[8][6][1][1]=l6r;
Curvature_Higgs_L4[8][6][1][3]=l5r;
Curvature_Higgs_L4[8][6][1][6]=-l6i;
Curvature_Higgs_L4[8][6][1][8]=-l5i;
Curvature_Higgs_L4[8][6][2][0]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][6][2][2]=l7r;
Curvature_Higgs_L4[8][6][2][5]=l5i/2.;
Curvature_Higgs_L4[8][6][3][1]=l5r;
Curvature_Higgs_L4[8][6][3][3]=l7r;
Curvature_Higgs_L4[8][6][3][6]=l5i;
Curvature_Higgs_L4[8][6][3][8]=l7i;
Curvature_Higgs_L4[8][6][4][4]=l10r;
Curvature_Higgs_L4[8][6][5][2]=l5i/2.;
Curvature_Higgs_L4[8][6][5][5]=l6r;
Curvature_Higgs_L4[8][6][5][7]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][6][6][1]=-l6i;
Curvature_Higgs_L4[8][6][6][3]=l5i;
Curvature_Higgs_L4[8][6][6][6]=3*l6r;
Curvature_Higgs_L4[8][6][6][8]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[8][6][7][0]=-0.5*l5i;
Curvature_Higgs_L4[8][6][7][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][6][7][7]=l7r;
Curvature_Higgs_L4[8][6][8][1]=-l5i;
Curvature_Higgs_L4[8][6][8][3]=l7i;
Curvature_Higgs_L4[8][6][8][6]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[8][6][8][8]=3*l7r;
Curvature_Higgs_L4[8][6][9][9]=l10r;
Curvature_Higgs_L4[8][7][0][1]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[8][7][0][6]=-0.5*l5i;
Curvature_Higgs_L4[8][7][0][8]=-l7i;
Curvature_Higgs_L4[8][7][1][0]=(4*l4 - 4*l5r)/8.;
Curvature_Higgs_L4[8][7][1][5]=-0.5*l5i;
Curvature_Higgs_L4[8][7][1][7]=-l7i;
Curvature_Higgs_L4[8][7][5][1]=-0.5*l5i;
Curvature_Higgs_L4[8][7][5][6]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][7][5][8]=l7r;
Curvature_Higgs_L4[8][7][6][0]=-0.5*l5i;
Curvature_Higgs_L4[8][7][6][5]=(4*l4 + 4*l5r)/8.;
Curvature_Higgs_L4[8][7][6][7]=l7r;
Curvature_Higgs_L4[8][7][7][1]=-l7i;
Curvature_Higgs_L4[8][7][7][6]=l7r;
Curvature_Higgs_L4[8][7][7][8]=l2;
Curvature_Higgs_L4[8][7][8][0]=-l7i;
Curvature_Higgs_L4[8][7][8][5]=l7r;
Curvature_Higgs_L4[8][7][8][7]=l2;
Curvature_Higgs_L4[8][8][0][0]=l3;
Curvature_Higgs_L4[8][8][0][2]=l7r;
Curvature_Higgs_L4[8][8][0][7]=-l7i;
Curvature_Higgs_L4[8][8][1][1]=(8*l3 + 8*l4 - 8*l5r)/8.;
Curvature_Higgs_L4[8][8][1][3]=l7r;
Curvature_Higgs_L4[8][8][1][6]=-l5i;
Curvature_Higgs_L4[8][8][1][8]=-3*l7i;
Curvature_Higgs_L4[8][8][2][0]=l7r;
Curvature_Higgs_L4[8][8][2][2]=l2;
Curvature_Higgs_L4[8][8][2][5]=l7i;
Curvature_Higgs_L4[8][8][3][1]=l7r;
Curvature_Higgs_L4[8][8][3][3]=l2;
Curvature_Higgs_L4[8][8][3][6]=l7i;
Curvature_Higgs_L4[8][8][4][4]=l9;
Curvature_Higgs_L4[8][8][5][2]=l7i;
Curvature_Higgs_L4[8][8][5][5]=l3;
Curvature_Higgs_L4[8][8][5][7]=l7r;
Curvature_Higgs_L4[8][8][6][1]=-l5i;
Curvature_Higgs_L4[8][8][6][3]=l7i;
Curvature_Higgs_L4[8][8][6][6]=(8*l3 + 8*l4 + 8*l5r)/8.;
Curvature_Higgs_L4[8][8][6][8]=3*l7r;
Curvature_Higgs_L4[8][8][7][0]=-l7i;
Curvature_Higgs_L4[8][8][7][5]=l7r;
Curvature_Higgs_L4[8][8][7][7]=l2;
Curvature_Higgs_L4[8][8][8][1]=-3*l7i;
Curvature_Higgs_L4[8][8][8][6]=3*l7r;
Curvature_Higgs_L4[8][8][8][8]=3*l2;
Curvature_Higgs_L4[8][8][9][9]=l9;
Curvature_Higgs_L4[8][9][1][9]=-l10i;
Curvature_Higgs_L4[8][9][6][9]=l10r;
Curvature_Higgs_L4[8][9][8][9]=l9;
Curvature_Higgs_L4[8][9][9][1]=-l10i;
Curvature_Higgs_L4[8][9][9][6]=l10r;
Curvature_Higgs_L4[8][9][9][8]=l9;
Curvature_Higgs_L4[9][0][0][9]=l8;
Curvature_Higgs_L4[9][0][2][9]=l10r;
Curvature_Higgs_L4[9][0][7][9]=-l10i;
Curvature_Higgs_L4[9][0][9][0]=l8;
Curvature_Higgs_L4[9][0][9][2]=l10r;
Curvature_Higgs_L4[9][0][9][7]=-l10i;
Curvature_Higgs_L4[9][1][1][9]=l8;
Curvature_Higgs_L4[9][1][3][9]=l10r;
Curvature_Higgs_L4[9][1][8][9]=-l10i;
Curvature_Higgs_L4[9][1][9][1]=l8;
Curvature_Higgs_L4[9][1][9][3]=l10r;
Curvature_Higgs_L4[9][1][9][8]=-l10i;
Curvature_Higgs_L4[9][2][0][9]=l10r;
Curvature_Higgs_L4[9][2][2][9]=l9;
Curvature_Higgs_L4[9][2][5][9]=l10i;
Curvature_Higgs_L4[9][2][9][0]=l10r;
Curvature_Higgs_L4[9][2][9][2]=l9;
Curvature_Higgs_L4[9][2][9][5]=l10i;
Curvature_Higgs_L4[9][3][1][9]=l10r;
Curvature_Higgs_L4[9][3][3][9]=l9;
Curvature_Higgs_L4[9][3][6][9]=l10i;
Curvature_Higgs_L4[9][3][9][1]=l10r;
Curvature_Higgs_L4[9][3][9][3]=l9;
Curvature_Higgs_L4[9][3][9][6]=l10i;
Curvature_Higgs_L4[9][4][4][9]=2*lh;
Curvature_Higgs_L4[9][4][9][4]=2*lh;
Curvature_Higgs_L4[9][5][2][9]=l10i;
Curvature_Higgs_L4[9][5][5][9]=l8;
Curvature_Higgs_L4[9][5][7][9]=l10r;
Curvature_Higgs_L4[9][5][9][2]=l10i;
Curvature_Higgs_L4[9][5][9][5]=l8;
Curvature_Higgs_L4[9][5][9][7]=l10r;
Curvature_Higgs_L4[9][6][3][9]=l10i;
Curvature_Higgs_L4[9][6][6][9]=l8;
Curvature_Higgs_L4[9][6][8][9]=l10r;
Curvature_Higgs_L4[9][6][9][3]=l10i;
Curvature_Higgs_L4[9][6][9][6]=l8;
Curvature_Higgs_L4[9][6][9][8]=l10r;
Curvature_Higgs_L4[9][7][0][9]=-l10i;
Curvature_Higgs_L4[9][7][5][9]=l10r;
Curvature_Higgs_L4[9][7][7][9]=l9;
Curvature_Higgs_L4[9][7][9][0]=-l10i;
Curvature_Higgs_L4[9][7][9][5]=l10r;
Curvature_Higgs_L4[9][7][9][7]=l9;
Curvature_Higgs_L4[9][8][1][9]=-l10i;
Curvature_Higgs_L4[9][8][6][9]=l10r;
Curvature_Higgs_L4[9][8][8][9]=l9;
Curvature_Higgs_L4[9][8][9][1]=-l10i;
Curvature_Higgs_L4[9][8][9][6]=l10r;
Curvature_Higgs_L4[9][8][9][8]=l9;
Curvature_Higgs_L4[9][9][0][0]=l8;
Curvature_Higgs_L4[9][9][0][2]=l10r;
Curvature_Higgs_L4[9][9][0][7]=-l10i;
Curvature_Higgs_L4[9][9][1][1]=l8;
Curvature_Higgs_L4[9][9][1][3]=l10r;
Curvature_Higgs_L4[9][9][1][8]=-l10i;
Curvature_Higgs_L4[9][9][2][0]=l10r;
Curvature_Higgs_L4[9][9][2][2]=l9;
Curvature_Higgs_L4[9][9][2][5]=l10i;
Curvature_Higgs_L4[9][9][3][1]=l10r;
Curvature_Higgs_L4[9][9][3][3]=l9;
Curvature_Higgs_L4[9][9][3][6]=l10i;
Curvature_Higgs_L4[9][9][4][4]=2*lh;
Curvature_Higgs_L4[9][9][5][2]=l10i;
Curvature_Higgs_L4[9][9][5][5]=l8;
Curvature_Higgs_L4[9][9][5][7]=l10r;
Curvature_Higgs_L4[9][9][6][3]=l10i;
Curvature_Higgs_L4[9][9][6][6]=l8;
Curvature_Higgs_L4[9][9][6][8]=l10r;
Curvature_Higgs_L4[9][9][7][0]=-l10i;
Curvature_Higgs_L4[9][9][7][5]=l10r;
Curvature_Higgs_L4[9][9][7][7]=l9;
Curvature_Higgs_L4[9][9][8][1]=-l10i;
Curvature_Higgs_L4[9][9][8][6]=l10r;
Curvature_Higgs_L4[9][9][8][8]=l9;
Curvature_Higgs_L4[9][9][9][9]=6*lh;

Curvature_Gauge_G2H2[0][0][0][0]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][1][1]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][2][2]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][3][3]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][4][4]=8*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][5][5]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][6][6]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][7][7]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][8][8]=2*std::pow(gp,2);
Curvature_Gauge_G2H2[0][0][9][9]=8*std::pow(gp,2);
Curvature_Gauge_G2H2[0][1][0][1]=2*g*gp;
Curvature_Gauge_G2H2[0][1][1][0]=2*g*gp;
Curvature_Gauge_G2H2[0][1][2][3]=2*g*gp;
Curvature_Gauge_G2H2[0][1][3][2]=2*g*gp;
Curvature_Gauge_G2H2[0][1][5][6]=2*g*gp;
Curvature_Gauge_G2H2[0][1][6][5]=2*g*gp;
Curvature_Gauge_G2H2[0][1][7][8]=2*g*gp;
Curvature_Gauge_G2H2[0][1][8][7]=2*g*gp;
Curvature_Gauge_G2H2[0][2][0][6]=2*g*gp;
Curvature_Gauge_G2H2[0][2][1][5]=-2*g*gp;
Curvature_Gauge_G2H2[0][2][2][8]=2*g*gp;
Curvature_Gauge_G2H2[0][2][3][7]=-2*g*gp;
Curvature_Gauge_G2H2[0][2][5][1]=-2*g*gp;
Curvature_Gauge_G2H2[0][2][6][0]=2*g*gp;
Curvature_Gauge_G2H2[0][2][7][3]=-2*g*gp;
Curvature_Gauge_G2H2[0][2][8][2]=2*g*gp;
Curvature_Gauge_G2H2[0][3][0][0]=2*g*gp;
Curvature_Gauge_G2H2[0][3][1][1]=-2*g*gp;
Curvature_Gauge_G2H2[0][3][2][2]=2*g*gp;
Curvature_Gauge_G2H2[0][3][3][3]=-2*g*gp;
Curvature_Gauge_G2H2[0][3][5][5]=2*g*gp;
Curvature_Gauge_G2H2[0][3][6][6]=-2*g*gp;
Curvature_Gauge_G2H2[0][3][7][7]=2*g*gp;
Curvature_Gauge_G2H2[0][3][8][8]=-2*g*gp;
Curvature_Gauge_G2H2[1][0][0][1]=2*g*gp;
Curvature_Gauge_G2H2[1][0][1][0]=2*g*gp;
Curvature_Gauge_G2H2[1][0][2][3]=2*g*gp;
Curvature_Gauge_G2H2[1][0][3][2]=2*g*gp;
Curvature_Gauge_G2H2[1][0][5][6]=2*g*gp;
Curvature_Gauge_G2H2[1][0][6][5]=2*g*gp;
Curvature_Gauge_G2H2[1][0][7][8]=2*g*gp;
Curvature_Gauge_G2H2[1][0][8][7]=2*g*gp;
Curvature_Gauge_G2H2[1][1][0][0]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][1][1]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][2][2]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][3][3]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][5][5]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][6][6]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][7][7]=2*std::pow(g,2);
Curvature_Gauge_G2H2[1][1][8][8]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][0][0][6]=2*g*gp;
Curvature_Gauge_G2H2[2][0][1][5]=-2*g*gp;
Curvature_Gauge_G2H2[2][0][2][8]=2*g*gp;
Curvature_Gauge_G2H2[2][0][3][7]=-2*g*gp;
Curvature_Gauge_G2H2[2][0][5][1]=-2*g*gp;
Curvature_Gauge_G2H2[2][0][6][0]=2*g*gp;
Curvature_Gauge_G2H2[2][0][7][3]=-2*g*gp;
Curvature_Gauge_G2H2[2][0][8][2]=2*g*gp;
Curvature_Gauge_G2H2[2][2][0][0]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][1][1]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][2][2]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][3][3]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][5][5]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][6][6]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][7][7]=2*std::pow(g,2);
Curvature_Gauge_G2H2[2][2][8][8]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][0][0][0]=2*g*gp;
Curvature_Gauge_G2H2[3][0][1][1]=-2*g*gp;
Curvature_Gauge_G2H2[3][0][2][2]=2*g*gp;
Curvature_Gauge_G2H2[3][0][3][3]=-2*g*gp;
Curvature_Gauge_G2H2[3][0][5][5]=2*g*gp;
Curvature_Gauge_G2H2[3][0][6][6]=-2*g*gp;
Curvature_Gauge_G2H2[3][0][7][7]=2*g*gp;
Curvature_Gauge_G2H2[3][0][8][8]=-2*g*gp;
Curvature_Gauge_G2H2[3][3][0][0]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][1][1]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][2][2]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][3][3]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][5][5]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][6][6]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][7][7]=2*std::pow(g,2);
Curvature_Gauge_G2H2[3][3][8][8]=2*std::pow(g,2);


}

bool Class_ZeeModel::CalculateDebyeSimplified()
{
// IMPLEMENT LATER
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_ZeeModel::CalculateDebyeGaugeSimplified()
{
// IMPLEMENT LATER
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double Class_ZeeModel::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
double v1ar = v[0];
double v1br = v[1];
double v2ar = v[2];
double v2br = v[3];
double vr = v[4];
double v1ai = v[5];
double v1bi = v[6];
double v2ai = v[7];
double v2bi = v[8];
double vi = v[9];
res=(4*mu1s*std::pow(v1ai,2) + l1*std::pow(v1ai,4) + 4*mu1s*std::pow(v1ar,2) + 2*l1*std::pow(v1ai,2)*std::pow(v1ar,2) + l1*std::pow(v1ar,4) + 4*mu1s*std::pow(v1bi,2) + 2*l1*std::pow(v1ai,2)*std::pow(v1bi,2) + 2*l1*std::pow(v1ar,2)*std::pow(v1bi,2) + l1*std::pow(v1bi,4) + 4*mu1s*std::pow(v1br,2) + 2*l1*std::pow(v1ai,2)*std::pow(v1br,2) + 2*l1*std::pow(v1ar,2)*std::pow(v1br,2) + 2*l1*std::pow(v1bi,2)*std::pow(v1br,2) + l1*std::pow(v1br,4) + 8*mu3s*v1ai*v2ai + 4*l6r*std::pow(v1ai,3)*v2ai - 4*l6i*std::pow(v1ai,2)*v1ar*v2ai + 4*l6r*v1ai*std::pow(v1ar,2)*v2ai - 4*l6i*std::pow(v1ar,3)*v2ai + 4*l6r*v1ai*std::pow(v1bi,2)*v2ai - 4*l6i*v1ar*std::pow(v1bi,2)*v2ai + 4*l6r*v1ai*std::pow(v1br,2)*v2ai - 4*l6i*v1ar*std::pow(v1br,2)*v2ai + 4*mu2s*std::pow(v2ai,2) + 2*l3*std::pow(v1ai,2)*std::pow(v2ai,2) + 2*l4*std::pow(v1ai,2)*std::pow(v2ai,2) + 2*l5r*std::pow(v1ai,2)*std::pow(v2ai,2) - 4*l5i*v1ai*v1ar*std::pow(v2ai,2) + 2*l3*std::pow(v1ar,2)*std::pow(v2ai,2) + 2*l4*std::pow(v1ar,2)*std::pow(v2ai,2) - 2*l5r*std::pow(v1ar,2)*std::pow(v2ai,2) + 2*l3*std::pow(v1bi,2)*std::pow(v2ai,2) + 2*l3*std::pow(v1br,2)*std::pow(v2ai,2) + 4*l7r*v1ai*std::pow(v2ai,3) - 4*l7i*v1ar*std::pow(v2ai,3) + l2*std::pow(v2ai,4) + 4*l6i*std::pow(v1ai,3)*v2ar + 8*mu3s*v1ar*v2ar + 4*l6r*std::pow(v1ai,2)*v1ar*v2ar + 4*l6i*v1ai*std::pow(v1ar,2)*v2ar + 4*l6r*std::pow(v1ar,3)*v2ar + 4*l6i*v1ai*std::pow(v1bi,2)*v2ar + 4*l6r*v1ar*std::pow(v1bi,2)*v2ar + 4*l6i*v1ai*std::pow(v1br,2)*v2ar + 4*l6r*v1ar*std::pow(v1br,2)*v2ar + 4*l5i*std::pow(v1ai,2)*v2ai*v2ar + 8*l5r*v1ai*v1ar*v2ai*v2ar - 4*l5i*std::pow(v1ar,2)*v2ai*v2ar + 4*l7i*v1ai*std::pow(v2ai,2)*v2ar + 4*l7r*v1ar*std::pow(v2ai,2)*v2ar + 4*mu2s*std::pow(v2ar,2) + 2*l3*std::pow(v1ai,2)*std::pow(v2ar,2) + 2*l4*std::pow(v1ai,2)*std::pow(v2ar,2) - 2*l5r*std::pow(v1ai,2)*std::pow(v2ar,2) + 4*l5i*v1ai*v1ar*std::pow(v2ar,2) + 2*l3*std::pow(v1ar,2)*std::pow(v2ar,2) + 2*l4*std::pow(v1ar,2)*std::pow(v2ar,2) + 2*l5r*std::pow(v1ar,2)*std::pow(v2ar,2) + 2*l3*std::pow(v1bi,2)*std::pow(v2ar,2) + 2*l3*std::pow(v1br,2)*std::pow(v2ar,2) + 4*l7r*v1ai*v2ai*std::pow(v2ar,2) - 4*l7i*v1ar*v2ai*std::pow(v2ar,2) + 2*l2*std::pow(v2ai,2)*std::pow(v2ar,2) + 4*l7i*v1ai*std::pow(v2ar,3) + 4*l7r*v1ar*std::pow(v2ar,3) + l2*std::pow(v2ar,4) + 8*mu3s*v1bi*v2bi + 4*l6r*std::pow(v1ai,2)*v1bi*v2bi + 4*l6r*std::pow(v1ar,2)*v1bi*v2bi + 4*l6r*std::pow(v1bi,3)*v2bi - 4*l6i*std::pow(v1ai,2)*v1br*v2bi - 4*l6i*std::pow(v1ar,2)*v1br*v2bi - 4*l6i*std::pow(v1bi,2)*v1br*v2bi + 4*l6r*v1bi*std::pow(v1br,2)*v2bi - 4*l6i*std::pow(v1br,3)*v2bi + 4*l4*v1ai*v1bi*v2ai*v2bi + 4*l5r*v1ai*v1bi*v2ai*v2bi - 4*l5i*v1ar*v1bi*v2ai*v2bi - 4*l5i*v1ai*v1br*v2ai*v2bi + 4*l4*v1ar*v1br*v2ai*v2bi - 4*l5r*v1ar*v1br*v2ai*v2bi + 4*l7r*v1bi*std::pow(v2ai,2)*v2bi - 4*l7i*v1br*std::pow(v2ai,2)*v2bi + 4*l5i*v1ai*v1bi*v2ar*v2bi + 4*l4*v1ar*v1bi*v2ar*v2bi + 4*l5r*v1ar*v1bi*v2ar*v2bi - 4*l4*v1ai*v1br*v2ar*v2bi + 4*l5r*v1ai*v1br*v2ar*v2bi - 4*l5i*v1ar*v1br*v2ar*v2bi + 4*l7r*v1bi*std::pow(v2ar,2)*v2bi - 4*l7i*v1br*std::pow(v2ar,2)*v2bi + 4*mu2s*std::pow(v2bi,2) + 2*l3*std::pow(v1ai,2)*std::pow(v2bi,2) + 2*l3*std::pow(v1ar,2)*std::pow(v2bi,2) + 2*l3*std::pow(v1bi,2)*std::pow(v2bi,2) + 2*l4*std::pow(v1bi,2)*std::pow(v2bi,2) + 2*l5r*std::pow(v1bi,2)*std::pow(v2bi,2) - 4*l5i*v1bi*v1br*std::pow(v2bi,2) + 2*l3*std::pow(v1br,2)*std::pow(v2bi,2) + 2*l4*std::pow(v1br,2)*std::pow(v2bi,2) - 2*l5r*std::pow(v1br,2)*std::pow(v2bi,2) + 4*l7r*v1ai*v2ai*std::pow(v2bi,2) - 4*l7i*v1ar*v2ai*std::pow(v2bi,2) + 2*l2*std::pow(v2ai,2)*std::pow(v2bi,2) + 4*l7i*v1ai*v2ar*std::pow(v2bi,2) + 4*l7r*v1ar*v2ar*std::pow(v2bi,2) + 2*l2*std::pow(v2ar,2)*std::pow(v2bi,2) + 4*l7r*v1bi*std::pow(v2bi,3) - 4*l7i*v1br*std::pow(v2bi,3) + l2*std::pow(v2bi,4) + 4*l6i*std::pow(v1ai,2)*v1bi*v2br + 4*l6i*std::pow(v1ar,2)*v1bi*v2br + 4*l6i*std::pow(v1bi,3)*v2br + 8*mu3s*v1br*v2br + 4*l6r*std::pow(v1ai,2)*v1br*v2br + 4*l6r*std::pow(v1ar,2)*v1br*v2br + 4*l6r*std::pow(v1bi,2)*v1br*v2br + 4*l6i*v1bi*std::pow(v1br,2)*v2br + 4*l6r*std::pow(v1br,3)*v2br + 4*l5i*v1ai*v1bi*v2ai*v2br - 4*l4*v1ar*v1bi*v2ai*v2br + 4*l5r*v1ar*v1bi*v2ai*v2br + 4*l4*v1ai*v1br*v2ai*v2br + 4*l5r*v1ai*v1br*v2ai*v2br - 4*l5i*v1ar*v1br*v2ai*v2br + 4*l7i*v1bi*std::pow(v2ai,2)*v2br + 4*l7r*v1br*std::pow(v2ai,2)*v2br + 4*l4*v1ai*v1bi*v2ar*v2br - 4*l5r*v1ai*v1bi*v2ar*v2br + 4*l5i*v1ar*v1bi*v2ar*v2br + 4*l5i*v1ai*v1br*v2ar*v2br + 4*l4*v1ar*v1br*v2ar*v2br + 4*l5r*v1ar*v1br*v2ar*v2br + 4*l7i*v1bi*std::pow(v2ar,2)*v2br + 4*l7r*v1br*std::pow(v2ar,2)*v2br + 4*l5i*std::pow(v1bi,2)*v2bi*v2br + 8*l5r*v1bi*v1br*v2bi*v2br - 4*l5i*std::pow(v1br,2)*v2bi*v2br + 4*l7i*v1bi*std::pow(v2bi,2)*v2br + 4*l7r*v1br*std::pow(v2bi,2)*v2br + 4*mu2s*std::pow(v2br,2) + 2*l3*std::pow(v1ai,2)*std::pow(v2br,2) + 2*l3*std::pow(v1ar,2)*std::pow(v2br,2) + 2*l3*std::pow(v1bi,2)*std::pow(v2br,2) + 2*l4*std::pow(v1bi,2)*std::pow(v2br,2) - 2*l5r*std::pow(v1bi,2)*std::pow(v2br,2) + 4*l5i*v1bi*v1br*std::pow(v2br,2) + 2*l3*std::pow(v1br,2)*std::pow(v2br,2) + 2*l4*std::pow(v1br,2)*std::pow(v2br,2) + 2*l5r*std::pow(v1br,2)*std::pow(v2br,2) + 4*l7r*v1ai*v2ai*std::pow(v2br,2) - 4*l7i*v1ar*v2ai*std::pow(v2br,2) + 2*l2*std::pow(v2ai,2)*std::pow(v2br,2) + 4*l7i*v1ai*v2ar*std::pow(v2br,2) + 4*l7r*v1ar*v2ar*std::pow(v2br,2) + 2*l2*std::pow(v2ar,2)*std::pow(v2br,2) + 4*l7r*v1bi*v2bi*std::pow(v2br,2) - 4*l7i*v1br*v2bi*std::pow(v2br,2) + 2*l2*std::pow(v2bi,2)*std::pow(v2br,2) + 4*l7i*v1bi*std::pow(v2br,3) + 4*l7r*v1br*std::pow(v2br,3) + l2*std::pow(v2br,4) + 4*std::sqrt(2)*mu*(v1br*v2ai + v1bi*v2ar - v1ar*v2bi - v1ai*v2br)*vi + 2*lh*std::pow(vi,4) + 4*std::sqrt(2)*mu*(v1bi*v2ai - v1br*v2ar - v1ai*v2bi + v1ar*v2br)*vr + 2*(2*muhs + l8*(std::pow(v1ai,2) + std::pow(v1ar,2) + std::pow(v1bi,2) + std::pow(v1br,2)) + 2*l10r*v1ai*v2ai - 2*l10i*v1ar*v2ai + l9*std::pow(v2ai,2) + 2*l10i*v1ai*v2ar + 2*l10r*v1ar*v2ar + l9*std::pow(v2ar,2) + 2*l10r*v1bi*v2bi - 2*l10i*v1br*v2bi + l9*std::pow(v2bi,2) + 2*l10i*v1bi*v2br + 2*l10r*v1br*v2br + l9*std::pow(v2br,2))*std::pow(vr,2) + 2*lh*std::pow(vr,4) + 2*std::pow(vi,2)*(2*muhs + l8*std::pow(v1ai,2) + l8*std::pow(v1ar,2) + l8*std::pow(v1bi,2) + l8*std::pow(v1br,2) + 2*l10r*v1ai*v2ai - 2*l10i*v1ar*v2ai + l9*std::pow(v2ai,2) + 2*l10i*v1ai*v2ar + 2*l10r*v1ar*v2ar + l9*std::pow(v2ar,2) + 2*l10r*v1bi*v2bi - 2*l10i*v1br*v2bi + l9*std::pow(v2bi,2) + 2*l10i*v1bi*v2br + 2*l10r*v1br*v2br + l9*std::pow(v2br,2) + 2*lh*std::pow(vr,2)))/8.;


  return res;
}

double Class_ZeeModel::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0;
double v1ar = v[0];
double v1br = v[1];
double v2ar = v[2];
double v2br = v[3];
double vr = v[4];
double v1ai = v[5];
double v1bi = v[6];
double v2ai = v[7];
double v2bi = v[8];
double vi = v[9];
res=dT1ai*v1ai + dT1ar*v1ar + dT1bi*v1bi + dT1br*v1br + dT2ai*v2ai + dT2ar*v2ar + dT2bi*v2bi + dT2br*v2br + dTi*vi + dTr*vr;

  return res;
}

void Class_ZeeModel::Debugging(const std::vector<double> &input,
                               std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
