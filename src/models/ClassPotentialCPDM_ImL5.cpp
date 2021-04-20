/*
 * ClassTemplate.cpp
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

#include <ext/alloc_traits.h>     // for __alloc_traits<>::value_type
#include <stddef.h>               // for std::size_t
#include <algorithm>              // for max, copy
#include <iostream>               // for operator<<, endl, basic_o...
#include <memory>                 // for allocator_traits<>::value...
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"

#include <BSMPT/models/ClassPotentialCPDM_ImL5.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>
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
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
                Class_CPDM_ImL5::Class_CPDM_ImL5()
                {
                        Model = ModelID::ModelIDs::CPDM; // global int constant which will be used to tell the program which model is called
                        NNeutralHiggs = 5;               // number of neutral Higgs bosons at T = 0
                        NChargedHiggs = 4;               // number of charged Higgs bosons  at T = 0 (all d.o.f.)

                        nPar = 13;   // number of parameters in the tree-Level Lagrangian
                        nParCT = 19; // number of parameters in the counterterm potential

                        nVEV = 5; // number of VEVs to minimize the potential

                        NHiggs = NNeutralHiggs + NChargedHiggs;
                        NGauge = 4;
                        NLepton = 9;
                        NQuarks = 12;

                        VevOrder.resize(nVEV);
                        // Here you have to tell which scalar field gets which VEV.
                        VevOrder[0] = 2; //CB
                        VevOrder[1] = 4; //v1
                        VevOrder[2] = 6; //v2
                        VevOrder[3] = 7; //CP
                        VevOrder[4] = 8; //vS

                        // Set UseVTreeSimplified to use the tree-level potential defined in VTreeSimplified
                        UseVTreeSimplified = false;

                        // Set UseVCounterSimplified to use the counterterm potential defined in VCounterSimplified
                        UseVCounterSimplified = false;
                }

                Class_CPDM_ImL5::~Class_CPDM_ImL5()
                {
                }

                /**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
                std::vector<std::string> Class_CPDM_ImL5::addLegendCT() const
                {
                        std::vector<std::string> labels;
                        labels.push_back("CTL1");
                        labels.push_back("CTL2");
                        labels.push_back("CTL3");
                        labels.push_back("CTL4");
                        labels.push_back("CTL5");
                        labels.push_back("CTL6");
                        labels.push_back("CTL7");
                        labels.push_back("CTL8");
                        labels.push_back("CTReA");
                        labels.push_back("CTImA");
                        labels.push_back("CTm11sqrt");
                        labels.push_back("CTm22sqrt");
                        labels.push_back("CTmSsqrt");
                        labels.push_back("dTCB");
                        labels.push_back("dTCP");
                        labels.push_back("dT1");
                        labels.push_back("dT2");
                        labels.push_back("dTs");
                        labels.push_back("CTImL5");
                        return labels;
                }

                /**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
                std::vector<std::string> Class_CPDM_ImL5::addLegendTemp() const
                {
                        std::vector<std::string> labels;
                        labels.push_back("T_c");     // Label for the critical temperature
                        labels.push_back("v_c");     // Label for the critical vev
                        labels.push_back("v_c/T_c"); // Label for v_c/T_c, you could use xi_c also for example
                        //out += "Your VEV order"; // Now you have to put the label for your vevs
                        labels.push_back("wCB");
                        labels.push_back("w1");
                        labels.push_back("w2");
                        labels.push_back("wS");
                        labels.push_back("wCP");
                        return labels;
                }

                /**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
                std::vector<std::string> Class_CPDM_ImL5::addLegendTripleCouplings() const
                {
                        std::vector<std::string> labels;
                        std::vector<std::string> particles;
                        particles.resize(NHiggs);
                        //here you have to define the particle names in the vector particles

                        particles[0] = "H";

                        std::string out = "Tree_";
                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                for (std::size_t j = i; j < NHiggs; j++)
                                {
                                        for (std::size_t k = j; k < NHiggs; k++)
                                        {
                                                labels.push_back("Tree_" + particles.at(i) + particles.at(j) + particles.at(k));
                                                labels.push_back("CT_" + particles.at(i) + particles.at(j) + particles.at(k));
                                                labels.push_back("CW_" + particles.at(i) + particles.at(j) + particles.at(k));
                                        }
                                }
                        }

                        return labels;
                }

                /**
 * returns a string which tells the user the chronological order of the VEVs. Use this to
 * complement the legend of the given input file
 */
                std::vector<std::string> Class_CPDM_ImL5::addLegendVEV() const
                {
                        std::vector<std::string> labels;
                        //out = "Your VEV order";
                        labels.push_back("wCB");
                        labels.push_back("w1");
                        labels.push_back("w2");
                        labels.push_back("wS");
                        labels.push_back("wCP");
                        return labels;
                }

                /**
 * Reads the string linestr and sets the parameter point
 */
                void Class_CPDM_ImL5::ReadAndSet(const std::string &linestr, std::vector<double> &par)
                {
                        std::stringstream ss(linestr);
                        double temp;

                        // double temp_m11sqrt, temp_m22sqrt, temp_mSsqrt, temp_ReA, temp_ImA, temp_L1, temp_L2, temp_L3, temp_L4, temp_L5, temp_L6, temp_L7, temp_L8;

                        if (UseIndexCol)
                        {
                                ss >> temp;
                        }

                        for (int k = 1; k <= 30; k++)
                        {
                                ss >> temp;
                                //				if (k == 1)
                                //				{
                                //					if (temp == 0)
                                //					{
                                //						if (debug)
                                //							std::cout << "Free CT param set to = " << temp << std::endl;
                                //						tCT_free = temp;
                                //					}
                                //					if (temp > 0)
                                //					{
                                //						if (debug)
                                //							std::cout << "Free CT param is set to " << temp << std::endl;
                                //						tCT_free = temp;
                                //					}
                                //				}
                                //				if (k == 2)
                                //				{
                                //					if (temp == 0)
                                //					{
                                //						if (debug)
                                //							std::cout << "Dep CT param set to = " << temp << std::endl;
                                //						tCT_dep = temp;
                                //					}
                                //					if (temp > 0)
                                //					{
                                //						if (debug)
                                //							std::cout << "Dep CT param is set to " << temp << std::endl;
                                //						tCT_dep = temp;
                                //					}
                                //				}

                                if (k == 15)
                                {
                                        L1 = temp;
                                        // std::cout << "\tL1 = " << temp << std::endl;
                                }
                                if (k == 16)
                                {
                                        L2 = temp;
                                        // std::cout << "\tL2 = " << temp << std::endl;
                                }
                                if (k == 17)
                                {
                                        L3 = temp;
                                        // std::cout << "\tL3 = " << temp << std::endl;
                                }
                                if (k == 18)
                                {
                                        L4 = temp;
                                        // std::cout << "\tL4 = " << temp << std::endl;
                                }
                                if (k == 19)
                                {
                                        L5 = temp;
                                        // std::cout << "\tL5 = " << temp << std::endl;
                                }
                                if (k == 20)
                                {
                                        L6 = temp;
                                        // std::cout << "\tL6 = " << temp << std::endl;
                                }
                                if (k == 21)
                                {
                                        L7 = temp;
                                        // std::cout << "\tL7 = " << temp << std::endl;
                                }
                                if (k == 22 + 2)
                                {
                                        L8 = temp;
                                        // std::cout << "\tL8 = " << temp << std::endl;
                                }
                                if (k == 23)
                                {
                                        ReA = temp;
                                        // std::cout << "\tLReA = " << temp << std::endl;
                                }
                                if (k == 24)
                                {
                                        ImA = temp;
                                        // std::cout << "\tImA = " << temp << std::endl;
                                }
                                if (k == 25)
                                {
                                        m11sqrt = temp;
                                        // std::cout << "\tm11sqrt= " << temp << std::endl;
                                }
                                if (k == 26)
                                {
                                        m22sqrt = temp;
                                        // std::cout << "\tm22sqrt = " << temp << std::endl;
                                }
                                if (k == 27)
                                {
                                        mSsqrt = temp;
                                        // std::cout << "\tmSsqrt = " << temp << std::endl;
                                }
                                // if (k == 28)
                                // {
                                // 	// std::cout << "This should be the VEV: " << temp << std::endl;
                                // }
                        }
                        par[0] = L1;
                        par[1] = L2;
                        par[2] = L3;
                        par[3] = L4;
                        par[4] = L5;
                        par[5] = L6;
                        par[6] = L7;
                        par[7] = L8;
                        par[8] = ReA;
                        par[9] = ImA;
                        par[10] = m11sqrt;
                        par[11] = m22sqrt;
                        par[12] = mSsqrt;

                        set_gen(par); // This you have to call so that everything will be set
                        return;
                }

                /**
 * Set Class Object as well as the VEV configuration
 */
                void Class_CPDM_ImL5::set_gen(const std::vector<double> &par)
                {
                        std::cout << "Calling" << __func__ << std::endl;
                        L1 = par[0];
                        L2 = par[1];
                        L3 = par[2];
                        L4 = par[3];
                        L5 = par[4];
                        L6 = par[5];
                        L7 = par[6];
                        L8 = par[7];
                        ReA = par[8];
                        ImA = par[9];
                        m11sqrt = par[10];
                        m22sqrt = par[11];
                        mSsqrt = par[12];

                        scale = C_vev0; // Renormalisation scale is set to the SM VEV
                        vevTreeMin.resize(nVEV);
                        vevTree.resize(NHiggs);
                        vevTreeMin[0] = 0;      //CB
                        vevTreeMin[1] = C_vev0; //v1
                        vevTreeMin[2] = 0;      //v2
                        vevTreeMin[3] = 0;      //vS
                        vevTreeMin[4] = 0;      //CP

                        vevTree = MinimizeOrderVEV(vevTreeMin);

                        if (!SetCurvatureDone)
                                SetCurvatureArrays();
                }

                /**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
                void Class_CPDM_ImL5::set_CT_Pot_Par(const std::vector<double> &par)
                {

                        CTL1 = par[0];
                        CTL2 = par[1];
                        CTL3 = par[2];
                        CTL4 = par[3];
                        CTL5 = par[4];
                        CTL6 = par[5];
                        CTL7 = par[6];
                        CTL8 = par[7];
                        CTReA = par[8];
                        CTImA = par[9];
                        CTm11sqrt = par[10];
                        CTm22sqrt = par[11];
                        CTmSsqrt = par[12];
                        dTCB = par[13];
                        dTCP = par[14];
                        dT1 = par[15];
                        dT2 = par[16];
                        dTs = par[17];
                        CTImL5 = par[18];

                        Curvature_Higgs_CT_L1[2] = dTCB;
                        Curvature_Higgs_CT_L1[4] = dT1;
                        Curvature_Higgs_CT_L1[6] = dT2;
                        Curvature_Higgs_CT_L1[7] = dTCP;
                        Curvature_Higgs_CT_L1[8] = dTs;

                        Curvature_Higgs_CT_L2[0][0] = CTm11sqrt;
                        Curvature_Higgs_CT_L2[1][1] = CTm11sqrt;
                        Curvature_Higgs_CT_L2[2][2] = CTm22sqrt;
                        Curvature_Higgs_CT_L2[3][3] = CTm22sqrt;
                        Curvature_Higgs_CT_L2[4][4] = CTm11sqrt;
                        Curvature_Higgs_CT_L2[5][5] = CTm11sqrt;
                        Curvature_Higgs_CT_L2[6][6] = CTm22sqrt;
                        Curvature_Higgs_CT_L2[7][7] = CTm22sqrt;
                        Curvature_Higgs_CT_L2[8][8] = CTmSsqrt;

                        Curvature_Higgs_CT_L3[8][2][0] = CTReA;
                        Curvature_Higgs_CT_L3[8][3][0] = -CTImA;
                        Curvature_Higgs_CT_L3[2][8][0] = CTReA;
                        Curvature_Higgs_CT_L3[3][8][0] = -CTImA;
                        Curvature_Higgs_CT_L3[8][2][1] = CTImA;
                        Curvature_Higgs_CT_L3[8][3][1] = CTReA;
                        Curvature_Higgs_CT_L3[2][8][1] = CTImA;
                        Curvature_Higgs_CT_L3[3][8][1] = CTReA;
                        Curvature_Higgs_CT_L3[8][0][2] = CTReA;
                        Curvature_Higgs_CT_L3[8][1][2] = CTImA;
                        Curvature_Higgs_CT_L3[0][8][2] = CTReA;
                        Curvature_Higgs_CT_L3[1][8][2] = CTImA;
                        Curvature_Higgs_CT_L3[8][0][3] = -CTImA;
                        Curvature_Higgs_CT_L3[8][1][3] = CTReA;
                        Curvature_Higgs_CT_L3[0][8][3] = -CTImA;
                        Curvature_Higgs_CT_L3[1][8][3] = CTReA;
                        Curvature_Higgs_CT_L3[8][6][4] = CTReA;
                        Curvature_Higgs_CT_L3[8][7][4] = -CTImA;
                        Curvature_Higgs_CT_L3[6][8][4] = CTReA;
                        Curvature_Higgs_CT_L3[7][8][4] = -CTImA;
                        Curvature_Higgs_CT_L3[8][6][5] = CTImA;
                        Curvature_Higgs_CT_L3[8][7][5] = CTReA;
                        Curvature_Higgs_CT_L3[6][8][5] = CTImA;
                        Curvature_Higgs_CT_L3[7][8][5] = CTReA;
                        Curvature_Higgs_CT_L3[8][4][6] = CTReA;
                        Curvature_Higgs_CT_L3[8][5][6] = CTImA;
                        Curvature_Higgs_CT_L3[4][8][6] = CTReA;
                        Curvature_Higgs_CT_L3[5][8][6] = CTImA;
                        Curvature_Higgs_CT_L3[8][4][7] = -CTImA;
                        Curvature_Higgs_CT_L3[8][5][7] = CTReA;
                        Curvature_Higgs_CT_L3[4][8][7] = -CTImA;
                        Curvature_Higgs_CT_L3[5][8][7] = CTReA;
                        Curvature_Higgs_CT_L3[2][0][8] = CTReA;
                        Curvature_Higgs_CT_L3[3][0][8] = -CTImA;
                        Curvature_Higgs_CT_L3[2][1][8] = CTImA;
                        Curvature_Higgs_CT_L3[3][1][8] = CTReA;
                        Curvature_Higgs_CT_L3[0][2][8] = CTReA;
                        Curvature_Higgs_CT_L3[1][2][8] = CTImA;
                        Curvature_Higgs_CT_L3[0][3][8] = -CTImA;
                        Curvature_Higgs_CT_L3[1][3][8] = CTReA;
                        Curvature_Higgs_CT_L3[6][4][8] = CTReA;
                        Curvature_Higgs_CT_L3[7][4][8] = -CTImA;
                        Curvature_Higgs_CT_L3[6][5][8] = CTImA;
                        Curvature_Higgs_CT_L3[7][5][8] = CTReA;
                        Curvature_Higgs_CT_L3[4][6][8] = CTReA;
                        Curvature_Higgs_CT_L3[5][6][8] = CTImA;
                        Curvature_Higgs_CT_L3[4][7][8] = -CTImA;
                        Curvature_Higgs_CT_L3[5][7][8] = CTReA;

                        Curvature_Higgs_CT_L4[0][0][0][0] = 3 * CTL1;
                        Curvature_Higgs_CT_L4[1][1][0][0] = CTL1;
                        Curvature_Higgs_CT_L4[2][2][0][0] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[3][2][0][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][3][0][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[3][3][0][0] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[4][4][0][0] = CTL1;
                        Curvature_Higgs_CT_L4[5][5][0][0] = CTL1;
                        Curvature_Higgs_CT_L4[6][6][0][0] = CTL3;
                        Curvature_Higgs_CT_L4[7][7][0][0] = CTL3;
                        Curvature_Higgs_CT_L4[8][8][0][0] = CTL7;
                        Curvature_Higgs_CT_L4[1][0][1][0] = CTL1;
                        Curvature_Higgs_CT_L4[0][1][1][0] = CTL1;
                        Curvature_Higgs_CT_L4[2][2][1][0] = CTImL5;
                        Curvature_Higgs_CT_L4[3][2][1][0] = CTL5;
                        Curvature_Higgs_CT_L4[2][3][1][0] = CTL5;
                        Curvature_Higgs_CT_L4[3][3][1][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][0][2][0] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[3][0][2][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][1][2][0] = CTImL5;
                        Curvature_Higgs_CT_L4[3][1][2][0] = CTL5;
                        Curvature_Higgs_CT_L4[0][2][2][0] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[1][2][2][0] = CTImL5;
                        Curvature_Higgs_CT_L4[0][3][2][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][3][2][0] = CTL5;
                        Curvature_Higgs_CT_L4[6][4][2][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][4][2][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][5][2][0] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][5][2][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][6][2][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][6][2][0] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][7][2][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][7][2][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][0][3][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[3][0][3][0] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[2][1][3][0] = CTL5;
                        Curvature_Higgs_CT_L4[3][1][3][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[0][2][3][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][2][3][0] = CTL5;
                        Curvature_Higgs_CT_L4[0][3][3][0] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[1][3][3][0] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][4][3][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][4][3][0] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][5][3][0] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][5][3][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][6][3][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][6][3][0] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][7][3][0] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][7][3][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][0][4][0] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][4][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][2][4][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][3][4][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][3][4][0] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][4][4][0] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][4][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][6][4][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][7][4][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][7][4][0] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][0][5][0] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][5][0] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][2][5][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][3][5][0] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][3][5][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][5][5][0] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][5][0] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][6][5][0] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][7][5][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][7][5][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][0][6][0] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][6][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][2][6][0] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][3][6][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][3][6][0] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][4][6][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][4][6][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][5][6][0] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][5][6][0] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][6][6][0] = CTL3;
                        Curvature_Higgs_CT_L4[7][0][7][0] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][7][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][2][7][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][3][7][0] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][3][7][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][4][7][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][4][7][0] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][5][7][0] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][5][7][0] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][7][7][0] = CTL3;
                        Curvature_Higgs_CT_L4[8][0][8][0] = CTL7;
                        Curvature_Higgs_CT_L4[0][8][8][0] = CTL7;
                        Curvature_Higgs_CT_L4[1][0][0][1] = CTL1;
                        Curvature_Higgs_CT_L4[0][1][0][1] = CTL1;
                        Curvature_Higgs_CT_L4[2][2][0][1] = CTImL5;
                        Curvature_Higgs_CT_L4[3][2][0][1] = CTL5;
                        Curvature_Higgs_CT_L4[2][3][0][1] = CTL5;
                        Curvature_Higgs_CT_L4[3][3][0][1] = -CTImL5;
                        Curvature_Higgs_CT_L4[0][0][1][1] = CTL1;
                        Curvature_Higgs_CT_L4[1][1][1][1] = 3 * CTL1;
                        Curvature_Higgs_CT_L4[2][2][1][1] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[3][2][1][1] = CTImL5;
                        Curvature_Higgs_CT_L4[2][3][1][1] = CTImL5;
                        Curvature_Higgs_CT_L4[3][3][1][1] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[4][4][1][1] = CTL1;
                        Curvature_Higgs_CT_L4[5][5][1][1] = CTL1;
                        Curvature_Higgs_CT_L4[6][6][1][1] = CTL3;
                        Curvature_Higgs_CT_L4[7][7][1][1] = CTL3;
                        Curvature_Higgs_CT_L4[8][8][1][1] = CTL7;
                        Curvature_Higgs_CT_L4[2][0][2][1] = CTImL5;
                        Curvature_Higgs_CT_L4[3][0][2][1] = CTL5;
                        Curvature_Higgs_CT_L4[2][1][2][1] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[3][1][2][1] = CTImL5;
                        Curvature_Higgs_CT_L4[0][2][2][1] = CTImL5;
                        Curvature_Higgs_CT_L4[1][2][2][1] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[0][3][2][1] = CTL5;
                        Curvature_Higgs_CT_L4[1][3][2][1] = CTImL5;
                        Curvature_Higgs_CT_L4[6][4][2][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][4][2][1] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][5][2][1] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][5][2][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][6][2][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][6][2][1] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][7][2][1] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][7][2][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][0][3][1] = CTL5;
                        Curvature_Higgs_CT_L4[3][0][3][1] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][1][3][1] = CTImL5;
                        Curvature_Higgs_CT_L4[3][1][3][1] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[0][2][3][1] = CTL5;
                        Curvature_Higgs_CT_L4[1][2][3][1] = CTImL5;
                        Curvature_Higgs_CT_L4[0][3][3][1] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][3][3][1] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[6][4][3][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][4][3][1] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][5][3][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][5][3][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][6][3][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][6][3][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][7][3][1] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][7][3][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][1][4][1] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][4][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][2][4][1] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][3][4][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][3][4][1] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][4][4][1] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][4][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][6][4][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][7][4][1] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][7][4][1] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][1][5][1] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][5][1] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][2][5][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][3][5][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][3][5][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][5][5][1] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][5][1] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][6][5][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][7][5][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][7][5][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][1][6][1] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][6][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][2][6][1] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][3][6][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][3][6][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][4][6][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][4][6][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][5][6][1] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][5][6][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][6][6][1] = CTL3;
                        Curvature_Higgs_CT_L4[7][1][7][1] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][7][1] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][2][7][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][3][7][1] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][3][7][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][4][7][1] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][4][7][1] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][5][7][1] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][5][7][1] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][7][7][1] = CTL3;
                        Curvature_Higgs_CT_L4[8][1][8][1] = CTL7;
                        Curvature_Higgs_CT_L4[1][8][8][1] = CTL7;
                        Curvature_Higgs_CT_L4[2][0][0][2] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[3][0][0][2] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][1][0][2] = CTImL5;
                        Curvature_Higgs_CT_L4[3][1][0][2] = CTL5;
                        Curvature_Higgs_CT_L4[0][2][0][2] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[1][2][0][2] = CTImL5;
                        Curvature_Higgs_CT_L4[0][3][0][2] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][3][0][2] = CTL5;
                        Curvature_Higgs_CT_L4[6][4][0][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][4][0][2] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][5][0][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][5][0][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][6][0][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][6][0][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][7][0][2] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][7][0][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][0][1][2] = CTImL5;
                        Curvature_Higgs_CT_L4[3][0][1][2] = CTL5;
                        Curvature_Higgs_CT_L4[2][1][1][2] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[3][1][1][2] = CTImL5;
                        Curvature_Higgs_CT_L4[0][2][1][2] = CTImL5;
                        Curvature_Higgs_CT_L4[1][2][1][2] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[0][3][1][2] = CTL5;
                        Curvature_Higgs_CT_L4[1][3][1][2] = CTImL5;
                        Curvature_Higgs_CT_L4[6][4][1][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][4][1][2] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][5][1][2] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][5][1][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][6][1][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][6][1][2] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][7][1][2] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][7][1][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][0][2][2] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[1][0][2][2] = CTImL5;
                        Curvature_Higgs_CT_L4[0][1][2][2] = CTImL5;
                        Curvature_Higgs_CT_L4[1][1][2][2] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[2][2][2][2] = 3 * CTL2;
                        Curvature_Higgs_CT_L4[3][3][2][2] = CTL2;
                        Curvature_Higgs_CT_L4[4][4][2][2] = CTL3;
                        Curvature_Higgs_CT_L4[5][5][2][2] = CTL3;
                        Curvature_Higgs_CT_L4[6][6][2][2] = CTL2;
                        Curvature_Higgs_CT_L4[7][7][2][2] = CTL2;
                        Curvature_Higgs_CT_L4[8][8][2][2] = CTL8;
                        Curvature_Higgs_CT_L4[0][0][3][2] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][0][3][2] = CTL5;
                        Curvature_Higgs_CT_L4[0][1][3][2] = CTL5;
                        Curvature_Higgs_CT_L4[1][1][3][2] = CTImL5;
                        Curvature_Higgs_CT_L4[3][2][3][2] = CTL2;
                        Curvature_Higgs_CT_L4[2][3][3][2] = CTL2;
                        Curvature_Higgs_CT_L4[6][0][4][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][0][4][2] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][1][4][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][1][4][2] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][2][4][2] = CTL3;
                        Curvature_Higgs_CT_L4[2][4][4][2] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][4][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][6][4][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][7][4][2] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][7][4][2] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][0][5][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][0][5][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][1][5][2] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][1][5][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][2][5][2] = CTL3;
                        Curvature_Higgs_CT_L4[2][5][5][2] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][5][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][6][5][2] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][7][5][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][7][5][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][0][6][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][0][6][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][1][6][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][1][6][2] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][2][6][2] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][6][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][4][6][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][5][6][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][5][6][2] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][6][6][2] = CTL2;
                        Curvature_Higgs_CT_L4[4][0][7][2] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][0][7][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][1][7][2] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][1][7][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][2][7][2] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][7][2] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][4][7][2] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][5][7][2] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][5][7][2] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][7][7][2] = CTL2;
                        Curvature_Higgs_CT_L4[8][2][8][2] = CTL8;
                        Curvature_Higgs_CT_L4[2][8][8][2] = CTL8;
                        Curvature_Higgs_CT_L4[2][0][0][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[3][0][0][3] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[2][1][0][3] = CTL5;
                        Curvature_Higgs_CT_L4[3][1][0][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[0][2][0][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][2][0][3] = CTL5;
                        Curvature_Higgs_CT_L4[0][3][0][3] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[1][3][0][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][4][0][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][4][0][3] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][5][0][3] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][5][0][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][6][0][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][6][0][3] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][7][0][3] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][7][0][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][0][1][3] = CTL5;
                        Curvature_Higgs_CT_L4[3][0][1][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][1][1][3] = CTImL5;
                        Curvature_Higgs_CT_L4[3][1][1][3] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[0][2][1][3] = CTL5;
                        Curvature_Higgs_CT_L4[1][2][1][3] = CTImL5;
                        Curvature_Higgs_CT_L4[0][3][1][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][3][1][3] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[6][4][1][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][4][1][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][5][1][3] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][5][1][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][6][1][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][6][1][3] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][7][1][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][7][1][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][0][2][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][0][2][3] = CTL5;
                        Curvature_Higgs_CT_L4[0][1][2][3] = CTL5;
                        Curvature_Higgs_CT_L4[1][1][2][3] = CTImL5;
                        Curvature_Higgs_CT_L4[3][2][2][3] = CTL2;
                        Curvature_Higgs_CT_L4[2][3][2][3] = CTL2;
                        Curvature_Higgs_CT_L4[0][0][3][3] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[1][0][3][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[0][1][3][3] = -CTImL5;
                        Curvature_Higgs_CT_L4[1][1][3][3] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[2][2][3][3] = CTL2;
                        Curvature_Higgs_CT_L4[3][3][3][3] = 3 * CTL2;
                        Curvature_Higgs_CT_L4[4][4][3][3] = CTL3;
                        Curvature_Higgs_CT_L4[5][5][3][3] = CTL3;
                        Curvature_Higgs_CT_L4[6][6][3][3] = CTL2;
                        Curvature_Higgs_CT_L4[7][7][3][3] = CTL2;
                        Curvature_Higgs_CT_L4[8][8][3][3] = CTL8;
                        Curvature_Higgs_CT_L4[6][0][4][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][0][4][3] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][1][4][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][1][4][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][3][4][3] = CTL3;
                        Curvature_Higgs_CT_L4[3][4][4][3] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][4][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][6][4][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][7][4][3] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][7][4][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][0][5][3] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][0][5][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][1][5][3] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][1][5][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][3][5][3] = CTL3;
                        Curvature_Higgs_CT_L4[3][5][5][3] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][5][3] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][6][5][3] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][7][5][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][7][5][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][0][6][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][0][6][3] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][1][6][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][1][6][3] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][3][6][3] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][6][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][4][6][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][5][6][3] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][5][6][3] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][6][6][3] = CTL2;
                        Curvature_Higgs_CT_L4[4][0][7][3] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][0][7][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][1][7][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][1][7][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][3][7][3] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][7][3] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][4][7][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][5][7][3] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][5][7][3] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][7][7][3] = CTL2;
                        Curvature_Higgs_CT_L4[8][3][8][3] = CTL8;
                        Curvature_Higgs_CT_L4[3][8][8][3] = CTL8;
                        Curvature_Higgs_CT_L4[4][0][0][4] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][0][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][2][0][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][3][0][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][3][0][4] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][4][0][4] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][0][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][6][0][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][7][0][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][7][0][4] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][1][1][4] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][1][4] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][2][1][4] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][3][1][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][3][1][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][4][1][4] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][1][4] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][6][1][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][7][1][4] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][7][1][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][0][2][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][0][2][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][1][2][4] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][1][2][4] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][2][2][4] = CTL3;
                        Curvature_Higgs_CT_L4[2][4][2][4] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][2][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][6][2][4] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][7][2][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][7][2][4] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][0][3][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][0][3][4] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][1][3][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][1][3][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][3][3][4] = CTL3;
                        Curvature_Higgs_CT_L4[3][4][3][4] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][3][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][6][3][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][7][3][4] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][7][3][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][0][4][4] = CTL1;
                        Curvature_Higgs_CT_L4[1][1][4][4] = CTL1;
                        Curvature_Higgs_CT_L4[2][2][4][4] = CTL3;
                        Curvature_Higgs_CT_L4[3][3][4][4] = CTL3;
                        Curvature_Higgs_CT_L4[4][4][4][4] = 3 * CTL1;
                        Curvature_Higgs_CT_L4[5][5][4][4] = CTL1;
                        Curvature_Higgs_CT_L4[6][6][4][4] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[7][6][4][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][7][4][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[7][7][4][4] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[8][8][4][4] = CTL7;
                        Curvature_Higgs_CT_L4[5][4][5][4] = CTL1;
                        Curvature_Higgs_CT_L4[4][5][5][4] = CTL1;
                        Curvature_Higgs_CT_L4[6][6][5][4] = CTImL5;
                        Curvature_Higgs_CT_L4[7][6][5][4] = CTL5;
                        Curvature_Higgs_CT_L4[6][7][5][4] = CTL5;
                        Curvature_Higgs_CT_L4[7][7][5][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][0][6][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][0][6][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][1][6][4] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][1][6][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][2][6][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][2][6][4] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][3][6][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][3][6][4] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][4][6][4] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[7][4][6][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][5][6][4] = CTImL5;
                        Curvature_Higgs_CT_L4[7][5][6][4] = CTL5;
                        Curvature_Higgs_CT_L4[4][6][6][4] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[5][6][6][4] = CTImL5;
                        Curvature_Higgs_CT_L4[4][7][6][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][7][6][4] = CTL5;
                        Curvature_Higgs_CT_L4[2][0][7][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][0][7][4] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][1][7][4] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][1][7][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][2][7][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][2][7][4] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][3][7][4] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][3][7][4] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][4][7][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[7][4][7][4] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[6][5][7][4] = CTL5;
                        Curvature_Higgs_CT_L4[7][5][7][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[4][6][7][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][6][7][4] = CTL5;
                        Curvature_Higgs_CT_L4[4][7][7][4] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[5][7][7][4] = -CTImL5;
                        Curvature_Higgs_CT_L4[8][4][8][4] = CTL7;
                        Curvature_Higgs_CT_L4[4][8][8][4] = CTL7;
                        Curvature_Higgs_CT_L4[5][0][0][5] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][0][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][2][0][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][3][0][5] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][3][0][5] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][5][0][5] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][0][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][6][0][5] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][7][0][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][7][0][5] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][1][1][5] = CTL1;
                        Curvature_Higgs_CT_L4[6][2][1][5] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][2][1][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][3][1][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][3][1][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][5][1][5] = CTL1;
                        Curvature_Higgs_CT_L4[2][6][1][5] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][6][1][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][7][1][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][7][1][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][0][2][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][0][2][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][1][2][5] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][1][2][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][2][2][5] = CTL3;
                        Curvature_Higgs_CT_L4[2][5][2][5] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][2][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][6][2][5] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][7][2][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][7][2][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][0][3][5] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][0][3][5] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][1][3][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][1][3][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][3][3][5] = CTL3;
                        Curvature_Higgs_CT_L4[3][5][3][5] = CTL3;
                        Curvature_Higgs_CT_L4[0][6][3][5] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][6][3][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][7][3][5] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][7][3][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][4][4][5] = CTL1;
                        Curvature_Higgs_CT_L4[4][5][4][5] = CTL1;
                        Curvature_Higgs_CT_L4[6][6][4][5] = CTImL5;
                        Curvature_Higgs_CT_L4[7][6][4][5] = CTL5;
                        Curvature_Higgs_CT_L4[6][7][4][5] = CTL5;
                        Curvature_Higgs_CT_L4[7][7][4][5] = -CTImL5;
                        Curvature_Higgs_CT_L4[0][0][5][5] = CTL1;
                        Curvature_Higgs_CT_L4[1][1][5][5] = CTL1;
                        Curvature_Higgs_CT_L4[2][2][5][5] = CTL3;
                        Curvature_Higgs_CT_L4[3][3][5][5] = CTL3;
                        Curvature_Higgs_CT_L4[4][4][5][5] = CTL1;
                        Curvature_Higgs_CT_L4[5][5][5][5] = 3 * CTL1;
                        Curvature_Higgs_CT_L4[6][6][5][5] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[7][6][5][5] = CTImL5;
                        Curvature_Higgs_CT_L4[6][7][5][5] = CTImL5;
                        Curvature_Higgs_CT_L4[7][7][5][5] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[8][8][5][5] = CTL7;
                        Curvature_Higgs_CT_L4[2][0][6][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][0][6][5] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][1][6][5] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][1][6][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][2][6][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][2][6][5] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][3][6][5] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][3][6][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][4][6][5] = CTImL5;
                        Curvature_Higgs_CT_L4[7][4][6][5] = CTL5;
                        Curvature_Higgs_CT_L4[6][5][6][5] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[7][5][6][5] = CTImL5;
                        Curvature_Higgs_CT_L4[4][6][6][5] = CTImL5;
                        Curvature_Higgs_CT_L4[5][6][6][5] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[4][7][6][5] = CTL5;
                        Curvature_Higgs_CT_L4[5][7][6][5] = CTImL5;
                        Curvature_Higgs_CT_L4[2][0][7][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][0][7][5] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][1][7][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][1][7][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][2][7][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][2][7][5] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][3][7][5] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][3][7][5] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][4][7][5] = CTL5;
                        Curvature_Higgs_CT_L4[7][4][7][5] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][5][7][5] = CTImL5;
                        Curvature_Higgs_CT_L4[7][5][7][5] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[4][6][7][5] = CTL5;
                        Curvature_Higgs_CT_L4[5][6][7][5] = CTImL5;
                        Curvature_Higgs_CT_L4[4][7][7][5] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][7][7][5] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[8][5][8][5] = CTL7;
                        Curvature_Higgs_CT_L4[5][8][8][5] = CTL7;
                        Curvature_Higgs_CT_L4[6][0][0][6] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][0][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][2][0][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][3][0][6] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][3][0][6] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][4][0][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][4][0][6] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][5][0][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][5][0][6] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][6][0][6] = CTL3;
                        Curvature_Higgs_CT_L4[6][1][1][6] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][1][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][2][1][6] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][3][1][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][3][1][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][4][1][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][4][1][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][5][1][6] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][5][1][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][6][1][6] = CTL3;
                        Curvature_Higgs_CT_L4[4][0][2][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][0][2][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][1][2][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][1][2][6] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][2][2][6] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][2][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][4][2][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][5][2][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][5][2][6] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][6][2][6] = CTL2;
                        Curvature_Higgs_CT_L4[4][0][3][6] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][0][3][6] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][1][3][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][1][3][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][3][3][6] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][3][6] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][4][3][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][5][3][6] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][5][3][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][6][3][6] = CTL2;
                        Curvature_Higgs_CT_L4[2][0][4][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][0][4][6] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][1][4][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][1][4][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][2][4][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][2][4][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][3][4][6] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][3][4][6] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][4][4][6] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[7][4][4][6] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][5][4][6] = CTImL5;
                        Curvature_Higgs_CT_L4[7][5][4][6] = CTL5;
                        Curvature_Higgs_CT_L4[4][6][4][6] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[5][6][4][6] = CTImL5;
                        Curvature_Higgs_CT_L4[4][7][4][6] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][7][4][6] = CTL5;
                        Curvature_Higgs_CT_L4[2][0][5][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][0][5][6] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][1][5][6] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][1][5][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][2][5][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][2][5][6] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][3][5][6] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][3][5][6] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][4][5][6] = CTImL5;
                        Curvature_Higgs_CT_L4[7][4][5][6] = CTL5;
                        Curvature_Higgs_CT_L4[6][5][5][6] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[7][5][5][6] = CTImL5;
                        Curvature_Higgs_CT_L4[4][6][5][6] = CTImL5;
                        Curvature_Higgs_CT_L4[5][6][5][6] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[4][7][5][6] = CTL5;
                        Curvature_Higgs_CT_L4[5][7][5][6] = CTImL5;
                        Curvature_Higgs_CT_L4[0][0][6][6] = CTL3;
                        Curvature_Higgs_CT_L4[1][1][6][6] = CTL3;
                        Curvature_Higgs_CT_L4[2][2][6][6] = CTL2;
                        Curvature_Higgs_CT_L4[3][3][6][6] = CTL2;
                        Curvature_Higgs_CT_L4[4][4][6][6] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[5][4][6][6] = CTImL5;
                        Curvature_Higgs_CT_L4[4][5][6][6] = CTImL5;
                        Curvature_Higgs_CT_L4[5][5][6][6] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[6][6][6][6] = 3 * CTL2;
                        Curvature_Higgs_CT_L4[7][7][6][6] = CTL2;
                        Curvature_Higgs_CT_L4[8][8][6][6] = CTL8;
                        Curvature_Higgs_CT_L4[4][4][7][6] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][4][7][6] = CTL5;
                        Curvature_Higgs_CT_L4[4][5][7][6] = CTL5;
                        Curvature_Higgs_CT_L4[5][5][7][6] = CTImL5;
                        Curvature_Higgs_CT_L4[7][6][7][6] = CTL2;
                        Curvature_Higgs_CT_L4[6][7][7][6] = CTL2;
                        Curvature_Higgs_CT_L4[8][6][8][6] = CTL8;
                        Curvature_Higgs_CT_L4[6][8][8][6] = CTL8;
                        Curvature_Higgs_CT_L4[7][0][0][7] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][0][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][2][0][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][3][0][7] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][3][0][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][4][0][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][4][0][7] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][5][0][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][5][0][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][7][0][7] = CTL3;
                        Curvature_Higgs_CT_L4[7][1][1][7] = CTL3;
                        Curvature_Higgs_CT_L4[4][2][1][7] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][2][1][7] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][3][1][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][3][1][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][4][1][7] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][4][1][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][5][1][7] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][5][1][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][7][1][7] = CTL3;
                        Curvature_Higgs_CT_L4[4][0][2][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][0][2][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[4][1][2][7] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][1][2][7] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[7][2][2][7] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][2][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][4][2][7] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][5][2][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][5][2][7] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][7][2][7] = CTL2;
                        Curvature_Higgs_CT_L4[4][0][3][7] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[5][0][3][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[4][1][3][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[5][1][3][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[7][3][3][7] = CTL2;
                        Curvature_Higgs_CT_L4[0][4][3][7] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][4][3][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][5][3][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][5][3][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][7][3][7] = CTL2;
                        Curvature_Higgs_CT_L4[2][0][4][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][0][4][7] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[2][1][4][7] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][1][4][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][2][4][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][2][4][7] = (-CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][3][4][7] = (CTL4 - CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][3][4][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[6][4][4][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[7][4][4][7] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[6][5][4][7] = CTL5;
                        Curvature_Higgs_CT_L4[7][5][4][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[4][6][4][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][6][4][7] = CTL5;
                        Curvature_Higgs_CT_L4[4][7][4][7] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[5][7][4][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[2][0][5][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[3][0][5][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[2][1][5][7] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[3][1][5][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[0][2][5][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[1][2][5][7] = CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[0][3][5][7] = -CTImL5 / 2.;
                        Curvature_Higgs_CT_L4[1][3][5][7] = (CTL4 + CTL5) / 2.;
                        Curvature_Higgs_CT_L4[6][4][5][7] = CTL5;
                        Curvature_Higgs_CT_L4[7][4][5][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[6][5][5][7] = CTImL5;
                        Curvature_Higgs_CT_L4[7][5][5][7] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[4][6][5][7] = CTL5;
                        Curvature_Higgs_CT_L4[5][6][5][7] = CTImL5;
                        Curvature_Higgs_CT_L4[4][7][5][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][7][5][7] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[4][4][6][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][4][6][7] = CTL5;
                        Curvature_Higgs_CT_L4[4][5][6][7] = CTL5;
                        Curvature_Higgs_CT_L4[5][5][6][7] = CTImL5;
                        Curvature_Higgs_CT_L4[7][6][6][7] = CTL2;
                        Curvature_Higgs_CT_L4[6][7][6][7] = CTL2;
                        Curvature_Higgs_CT_L4[0][0][7][7] = CTL3;
                        Curvature_Higgs_CT_L4[1][1][7][7] = CTL3;
                        Curvature_Higgs_CT_L4[2][2][7][7] = CTL2;
                        Curvature_Higgs_CT_L4[3][3][7][7] = CTL2;
                        Curvature_Higgs_CT_L4[4][4][7][7] = CTL3 + CTL4 - CTL5;
                        Curvature_Higgs_CT_L4[5][4][7][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[4][5][7][7] = -CTImL5;
                        Curvature_Higgs_CT_L4[5][5][7][7] = CTL3 + CTL4 + CTL5;
                        Curvature_Higgs_CT_L4[6][6][7][7] = CTL2;
                        Curvature_Higgs_CT_L4[7][7][7][7] = 3 * CTL2;
                        Curvature_Higgs_CT_L4[8][8][7][7] = CTL8;
                        Curvature_Higgs_CT_L4[8][7][8][7] = CTL8;
                        Curvature_Higgs_CT_L4[7][8][8][7] = CTL8;
                        Curvature_Higgs_CT_L4[8][0][0][8] = CTL7;
                        Curvature_Higgs_CT_L4[0][8][0][8] = CTL7;
                        Curvature_Higgs_CT_L4[8][1][1][8] = CTL7;
                        Curvature_Higgs_CT_L4[1][8][1][8] = CTL7;
                        Curvature_Higgs_CT_L4[8][2][2][8] = CTL8;
                        Curvature_Higgs_CT_L4[2][8][2][8] = CTL8;
                        Curvature_Higgs_CT_L4[8][3][3][8] = CTL8;
                        Curvature_Higgs_CT_L4[3][8][3][8] = CTL8;
                        Curvature_Higgs_CT_L4[8][4][4][8] = CTL7;
                        Curvature_Higgs_CT_L4[4][8][4][8] = CTL7;
                        Curvature_Higgs_CT_L4[8][5][5][8] = CTL7;
                        Curvature_Higgs_CT_L4[5][8][5][8] = CTL7;
                        Curvature_Higgs_CT_L4[8][6][6][8] = CTL8;
                        Curvature_Higgs_CT_L4[6][8][6][8] = CTL8;
                        Curvature_Higgs_CT_L4[8][7][7][8] = CTL8;
                        Curvature_Higgs_CT_L4[7][8][7][8] = CTL8;
                        Curvature_Higgs_CT_L4[0][0][8][8] = CTL7;
                        Curvature_Higgs_CT_L4[1][1][8][8] = CTL7;
                        Curvature_Higgs_CT_L4[2][2][8][8] = CTL8;
                        Curvature_Higgs_CT_L4[3][3][8][8] = CTL8;
                        Curvature_Higgs_CT_L4[4][4][8][8] = CTL7;
                        Curvature_Higgs_CT_L4[5][5][8][8] = CTL7;
                        Curvature_Higgs_CT_L4[6][6][8][8] = CTL8;
                        Curvature_Higgs_CT_L4[7][7][8][8] = CTL8;
                        Curvature_Higgs_CT_L4[8][8][8][8] = 6 * CTL6;

                        return;
                }

                /**
 * console output of all Parameters
 */
                void Class_CPDM_ImL5::write() const
                {

                        std::cout << "Model = " << Model << std::endl;
                        std::cout << "Counterterms of the " << Model << ":" << std::endl;
                        std::cout << "\tCTL1 = " << CTL1 << std::endl;
                        std::cout << "\tCTL2 = " << CTL2 << std::endl;
                        std::cout << "\tCTL3 = " << CTL3 << std::endl;
                        std::cout << "\tCTL4 = " << CTL4 << std::endl;
                        std::cout << "\tCTL5 = " << CTL5 << std::endl;
                        std::cout << "\tCTImL5 = " << CTImL5 << std::endl;
                        std::cout << "\tCTL6 = " << CTL6 << std::endl;
                        std::cout << "\tCTL7 = " << CTL7 << std::endl;
                        std::cout << "\tCTL8 = " << CTL8 << std::endl;
                        std::cout << "\tCTm11sqrt = " << CTm11sqrt << std::endl;
                        std::cout << "\tCTm22sqrt = " << CTm22sqrt << std::endl;
                        std::cout << "\tCTmSsqrt = " << CTmSsqrt << std::endl;
                        std::cout << "\tCTReA = " << CTReA << std::endl;
                        std::cout << "\tCTImA = " << CTImA << std::endl;
                        std::cout << "\tCTdT1 = " << dT1 << std::endl;
                        std::cout << "\tCTdT2 = " << dT2 << std::endl;
                        std::cout << "\tCTdTCP = " << dTCP << std::endl;
                        std::cout << "\tCTdTCB = " << dTCB << std::endl;
                        std::cout << "\tt_free = " << tCT_free << std::endl;
                        std::cout << "\tt_dep  = " << tCT_dep << std::endl;
                }

                /**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
                std::vector<double> Class_CPDM_ImL5::calc_CT() const
                {
                        std::cout << "Calling" << __func__ << std::endl;
                        std::vector<double> parCT;
                        parCT.clear();

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

                        VectorXd NCW(NHiggs);
                        MatrixXd HCW(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                NCW[i] = WeinbergNabla[i];
                                for (std::size_t j = 0; j < NHiggs; j++)
                                        HCW(i, j) = WeinbergHesse.at(j * NHiggs + i);
                        }

                        // double dummyCTL2 = tCT_free;
                        // double dummyCTL3 = tCT_dep;
                        // double dummyCTL7 = tCT_dep;
                        // double dummyCTL6 = tCT_free;
                        // double dummyCTL8 = tCT_free;
                        double dummyCTL2 = 0;
                        double dummyCTL3 = 0;
                        double dummyCTL7 = 0;
                        double dummyCTL6 = 0;
                        double dummyCTL8 = 0;

                        double temp_CTm11sqrt = HCW(4, 4) / 2. - (3 * HCW(5, 5)) / 2.;
                        double temp_CTm22sqrt = -(dummyCTL3 * std::pow(C_vev0, 2)) / 2. - HCW(2, 2);
                        double temp_CTmSsqrt = -(dummyCTL7 * std::pow(C_vev0, 2)) / 2. - HCW(8, 8);
                        double temp_CTReA = -(HCW(6, 8) / C_vev0);
                        double temp_CTImA = HCW(7, 8) / C_vev0;
                        double temp_CTL1 = -(HCW(4, 4) / std::pow(C_vev0, 2)) + HCW(5, 5) / std::pow(C_vev0, 2);
                        double temp_CTL2 = dummyCTL2;
                        double temp_CTL3 = dummyCTL3;
                        double temp_CTL4 = (2 * HCW(2, 2)) / std::pow(C_vev0, 2) - HCW(6, 6) / std::pow(C_vev0, 2) - HCW(7, 7) / std::pow(C_vev0, 2);
                        double temp_CTL5 = -(HCW(6, 6) / std::pow(C_vev0, 2)) + HCW(7, 7) / std::pow(C_vev0, 2);
                        double temp_CTImL5 = (2 * HCW(6, 7)) / std::pow(C_vev0, 2);
                        double temp_CTL6 = dummyCTL6;
                        double temp_CTL7 = dummyCTL7;
                        double temp_CTL8 = dummyCTL8;
                        double temp_dTCB = -NCW(2);
                        double temp_dT1 = C_vev0 * HCW(5, 5) - NCW(4);
                        double temp_dT2 = -NCW(6);
                        double temp_dTCP = -NCW(7);
                        double temp_dTS = -NCW(8);

                        std::vector<int> NonConstrainedNableIndices{0, 1, 3, 5};
                        std::vector<std::vector<int>> NonConstrainedHesseIndices{{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, {2, 7}, {2, 8}, {3, 4}, {3, 5}, {3, 6}, {3, 7}, {3, 8}, {4, 5}, {4, 6}, {4, 7}, {4, 8}, {5, 6}, {5, 7}, {5, 8}};

                        std::cout << "\n\nCT-Solution Check:\n"
                                  << std::endl;
                        std::cout << "------------------------------------------------------" << std::endl;
                        for (int x : NonConstrainedNableIndices)
                        {
                                if (std::abs(NCW(x)) > 1e-8)
                                        std::cout << "\tLoop-induced coupling @ Nab[" << x << "] = " << NCW(x) << std::endl;
                        }
                        for (auto indexpair : NonConstrainedHesseIndices)
                        {
                                if (std::abs(HCW(indexpair.at(0), indexpair.at(1))) > 1e-8)
                                        std::cout << "\tLoop-induced coupling @ Hes[" << indexpair.at(0) << "][" << indexpair.at(1) << "] = " << HCW(indexpair.at(0), indexpair.at(1)) << std::endl;
                        }
                        std::vector<double> ConsistencyRelations{-HCW(0, 0) + HCW(5, 5), (-2 * HCW(2, 2)) / std::pow(C_vev0, 2) + (2 * HCW(3, 3)) / std::pow(C_vev0, 2), -HCW(1, 1) + HCW(5, 5)};
                        for (auto x : ConsistencyRelations)
                                if (std::abs(x) > 1e-8)
                                        std::cout << "Consistency Relation is not fulfilled: 0=!" << x << std::endl;
                        std::cout << "------------------------------------------------------" << std::endl;

                        parCT.push_back(temp_CTL1);
                        parCT.push_back(temp_CTL2);
                        parCT.push_back(temp_CTL3);
                        parCT.push_back(temp_CTL4);
                        parCT.push_back(temp_CTL5);
                        parCT.push_back(temp_CTL6);
                        parCT.push_back(temp_CTL7);
                        parCT.push_back(temp_CTL8);
                        parCT.push_back(temp_CTReA);
                        parCT.push_back(temp_CTImA);
                        parCT.push_back(temp_CTm11sqrt);
                        parCT.push_back(temp_CTm22sqrt);
                        parCT.push_back(temp_CTmSsqrt);
                        parCT.push_back(temp_dTCB);
                        parCT.push_back(temp_dTCP);
                        parCT.push_back(temp_dT1);
                        parCT.push_back(temp_dT2);
                        parCT.push_back(temp_dTS);
                        parCT.push_back(temp_CTImL5);

                        return parCT;
                }

                void Class_CPDM_ImL5::TripleHiggsCouplings()
                {
                        if (!SetCurvatureDone)
                                SetCurvatureArrays();
                        if (!CalcCouplingsdone)
                                CalculatePhysicalCouplings();

                        std::vector<double> HiggsOrder(NHiggs);
                        // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
                        // particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

                        // example for keeping the mass order
                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                HiggsOrder[i] = i;
                        }

                        std::vector<double> TripleDeriv;
                        TripleDeriv = WeinbergThirdDerivative();
                        std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                                          std::vector<double>(NHiggs)));
                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                for (std::size_t j = 0; j < NHiggs; j++)
                                {
                                        for (std::size_t k = 0; k < NHiggs; k++)
                                        {
                                                GaugeBasis[i][j][k] = TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
                                        }
                                }
                        }

                        MatrixXd HiggsRot(NHiggs, NHiggs);
                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                for (std::size_t j = 0; j < NHiggs; j++)
                                {
                                        HiggsRot(i, j) = HiggsRotationMatrix[i][j];
                                }
                        }

                        MatrixXd HiggsRotSort(NHiggs, NHiggs);

                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
                        }

                        TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
                        TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
                        TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
                                TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
                                TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
                                for (std::size_t j = 0; j < NHiggs; j++)
                                {
                                        TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
                                        TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
                                        TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
                                }
                        }

                        for (std::size_t i = 0; i < NHiggs; i++)
                        {
                                for (std::size_t j = 0; j < NHiggs; j++)
                                {
                                        for (std::size_t k = 0; k < NHiggs; k++)
                                        {
                                                TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
                                                TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
                                                TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
                                                for (std::size_t l = 0; l < NHiggs; l++)
                                                {
                                                        for (std::size_t m = 0; m < NHiggs; m++)
                                                        {
                                                                for (std::size_t n = 0; n < NHiggs; n++)
                                                                {
                                                                        double RotFac = HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
                                                                        TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac * GaugeBasis[l][m][n];
                                                                        TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac * LambdaHiggs_3[l][m][n];
                                                                        TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac * LambdaHiggs_3_CT[l][m][n];
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }

                void Class_CPDM_ImL5::SetCurvatureArrays()
                {
                        /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

                        initVectors();
                        for (std::size_t i = 0; i < NHiggs; i++)
                                HiggsVev[i] = vevTree[i];

                        Curvature_Higgs_L2[0][0] = m11sqrt;
                        Curvature_Higgs_L2[1][1] = m11sqrt;
                        Curvature_Higgs_L2[2][2] = m22sqrt;
                        Curvature_Higgs_L2[3][3] = m22sqrt;
                        Curvature_Higgs_L2[4][4] = m11sqrt;
                        Curvature_Higgs_L2[5][5] = m11sqrt;
                        Curvature_Higgs_L2[6][6] = m22sqrt;
                        Curvature_Higgs_L2[7][7] = m22sqrt;
                        Curvature_Higgs_L2[8][8] = mSsqrt;

                        Curvature_Higgs_L3[0][2][8] = ReA;
                        Curvature_Higgs_L3[0][3][8] = -ImA;
                        Curvature_Higgs_L3[0][8][2] = ReA;
                        Curvature_Higgs_L3[0][8][3] = -ImA;
                        Curvature_Higgs_L3[1][2][8] = ImA;
                        Curvature_Higgs_L3[1][3][8] = ReA;
                        Curvature_Higgs_L3[1][8][2] = ImA;
                        Curvature_Higgs_L3[1][8][3] = ReA;
                        Curvature_Higgs_L3[2][0][8] = ReA;
                        Curvature_Higgs_L3[2][1][8] = ImA;
                        Curvature_Higgs_L3[2][8][0] = ReA;
                        Curvature_Higgs_L3[2][8][1] = ImA;
                        Curvature_Higgs_L3[3][0][8] = -ImA;
                        Curvature_Higgs_L3[3][1][8] = ReA;
                        Curvature_Higgs_L3[3][8][0] = -ImA;
                        Curvature_Higgs_L3[3][8][1] = ReA;
                        Curvature_Higgs_L3[4][6][8] = ReA;
                        Curvature_Higgs_L3[4][7][8] = -ImA;
                        Curvature_Higgs_L3[4][8][6] = ReA;
                        Curvature_Higgs_L3[4][8][7] = -ImA;
                        Curvature_Higgs_L3[5][6][8] = ImA;
                        Curvature_Higgs_L3[5][7][8] = ReA;
                        Curvature_Higgs_L3[5][8][6] = ImA;
                        Curvature_Higgs_L3[5][8][7] = ReA;
                        Curvature_Higgs_L3[6][4][8] = ReA;
                        Curvature_Higgs_L3[6][5][8] = ImA;
                        Curvature_Higgs_L3[6][8][4] = ReA;
                        Curvature_Higgs_L3[6][8][5] = ImA;
                        Curvature_Higgs_L3[7][4][8] = -ImA;
                        Curvature_Higgs_L3[7][5][8] = ReA;
                        Curvature_Higgs_L3[7][8][4] = -ImA;
                        Curvature_Higgs_L3[7][8][5] = ReA;
                        Curvature_Higgs_L3[8][0][2] = ReA;
                        Curvature_Higgs_L3[8][0][3] = -ImA;
                        Curvature_Higgs_L3[8][1][2] = ImA;
                        Curvature_Higgs_L3[8][1][3] = ReA;
                        Curvature_Higgs_L3[8][2][0] = ReA;
                        Curvature_Higgs_L3[8][2][1] = ImA;
                        Curvature_Higgs_L3[8][3][0] = -ImA;
                        Curvature_Higgs_L3[8][3][1] = ReA;
                        Curvature_Higgs_L3[8][4][6] = ReA;
                        Curvature_Higgs_L3[8][4][7] = -ImA;
                        Curvature_Higgs_L3[8][5][6] = ImA;
                        Curvature_Higgs_L3[8][5][7] = ReA;
                        Curvature_Higgs_L3[8][6][4] = ReA;
                        Curvature_Higgs_L3[8][6][5] = ImA;
                        Curvature_Higgs_L3[8][7][4] = -ImA;
                        Curvature_Higgs_L3[8][7][5] = ReA;

                        Curvature_Higgs_L4[0][0][0][0] = 3 * L1;

                        Curvature_Higgs_L4[0][0][1][1] = L1;

                        Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + L5;

                        Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - L5;

                        Curvature_Higgs_L4[0][0][4][4] = L1;

                        Curvature_Higgs_L4[0][0][5][5] = L1;

                        Curvature_Higgs_L4[0][0][6][6] = L3;

                        Curvature_Higgs_L4[0][0][7][7] = L3;

                        Curvature_Higgs_L4[0][0][8][8] = L7;

                        Curvature_Higgs_L4[0][1][0][1] = L1;

                        Curvature_Higgs_L4[0][1][1][0] = L1;

                        Curvature_Higgs_L4[0][1][2][3] = L5;

                        Curvature_Higgs_L4[0][1][3][2] = L5;

                        Curvature_Higgs_L4[0][2][0][2] = L3 + L4 + L5;

                        Curvature_Higgs_L4[0][2][1][3] = L5;

                        Curvature_Higgs_L4[0][2][2][0] = L3 + L4 + L5;

                        Curvature_Higgs_L4[0][2][3][1] = L5;

                        Curvature_Higgs_L4[0][2][4][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][2][5][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][2][6][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][2][7][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][3][0][3] = L3 + L4 - L5;

                        Curvature_Higgs_L4[0][3][1][2] = L5;

                        Curvature_Higgs_L4[0][3][2][1] = L5;

                        Curvature_Higgs_L4[0][3][3][0] = L3 + L4 - L5;

                        Curvature_Higgs_L4[0][3][4][7] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[0][3][5][6] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][3][6][5] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][3][7][4] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[0][4][0][4] = L1;

                        Curvature_Higgs_L4[0][4][2][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][4][3][7] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[0][4][4][0] = L1;

                        Curvature_Higgs_L4[0][4][6][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][4][7][3] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[0][5][0][5] = L1;

                        Curvature_Higgs_L4[0][5][2][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][5][3][6] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][5][5][0] = L1;

                        Curvature_Higgs_L4[0][5][6][3] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][5][7][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][6][0][6] = L3;

                        Curvature_Higgs_L4[0][6][2][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][6][3][5] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][6][4][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][6][5][3] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][6][6][0] = L3;

                        Curvature_Higgs_L4[0][7][0][7] = L3;

                        Curvature_Higgs_L4[0][7][2][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][7][3][4] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[0][7][4][3] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[0][7][5][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[0][7][7][0] = L3;

                        Curvature_Higgs_L4[0][8][0][8] = L7;

                        Curvature_Higgs_L4[0][8][8][0] = L7;

                        Curvature_Higgs_L4[1][0][0][1] = L1;

                        Curvature_Higgs_L4[1][0][1][0] = L1;

                        Curvature_Higgs_L4[1][0][2][3] = L5;

                        Curvature_Higgs_L4[1][0][3][2] = L5;

                        Curvature_Higgs_L4[1][1][0][0] = L1;

                        Curvature_Higgs_L4[1][1][1][1] = 3 * L1;

                        Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - L5;

                        Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + L5;

                        Curvature_Higgs_L4[1][1][4][4] = L1;

                        Curvature_Higgs_L4[1][1][5][5] = L1;

                        Curvature_Higgs_L4[1][1][6][6] = L3;

                        Curvature_Higgs_L4[1][1][7][7] = L3;

                        Curvature_Higgs_L4[1][1][8][8] = L7;

                        Curvature_Higgs_L4[1][2][0][3] = L5;

                        Curvature_Higgs_L4[1][2][1][2] = L3 + L4 - L5;

                        Curvature_Higgs_L4[1][2][2][1] = L3 + L4 - L5;

                        Curvature_Higgs_L4[1][2][3][0] = L5;

                        Curvature_Higgs_L4[1][2][4][7] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][2][5][6] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[1][2][6][5] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[1][2][7][4] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][3][0][2] = L5;

                        Curvature_Higgs_L4[1][3][1][3] = L3 + L4 + L5;

                        Curvature_Higgs_L4[1][3][2][0] = L5;

                        Curvature_Higgs_L4[1][3][3][1] = L3 + L4 + L5;

                        Curvature_Higgs_L4[1][3][4][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][3][5][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][3][6][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][3][7][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][4][1][4] = L1;

                        Curvature_Higgs_L4[1][4][2][7] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][4][3][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][4][4][1] = L1;

                        Curvature_Higgs_L4[1][4][6][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][4][7][2] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][5][1][5] = L1;

                        Curvature_Higgs_L4[1][5][2][6] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[1][5][3][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][5][5][1] = L1;

                        Curvature_Higgs_L4[1][5][6][2] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[1][5][7][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][6][1][6] = L3;

                        Curvature_Higgs_L4[1][6][2][5] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[1][6][3][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][6][4][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][6][5][2] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[1][6][6][1] = L3;

                        Curvature_Higgs_L4[1][7][1][7] = L3;

                        Curvature_Higgs_L4[1][7][2][4] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][7][3][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][7][4][2] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][7][5][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[1][7][7][1] = L3;

                        Curvature_Higgs_L4[1][8][1][8] = L7;

                        Curvature_Higgs_L4[1][8][8][1] = L7;

                        Curvature_Higgs_L4[2][0][0][2] = L3 + L4 + L5;

                        Curvature_Higgs_L4[2][0][1][3] = L5;

                        Curvature_Higgs_L4[2][0][2][0] = L3 + L4 + L5;

                        Curvature_Higgs_L4[2][0][3][1] = L5;

                        Curvature_Higgs_L4[2][0][4][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][0][5][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][0][6][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][0][7][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][1][0][3] = L5;

                        Curvature_Higgs_L4[2][1][1][2] = L3 + L4 - L5;

                        Curvature_Higgs_L4[2][1][2][1] = L3 + L4 - L5;

                        Curvature_Higgs_L4[2][1][3][0] = L5;

                        Curvature_Higgs_L4[2][1][4][7] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][1][5][6] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[2][1][6][5] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[2][1][7][4] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][2][0][0] = L3 + L4 + L5;

                        Curvature_Higgs_L4[2][2][1][1] = L3 + L4 - L5;

                        Curvature_Higgs_L4[2][2][2][2] = 3 * L2;

                        Curvature_Higgs_L4[2][2][3][3] = L2;

                        Curvature_Higgs_L4[2][2][4][4] = L3;

                        Curvature_Higgs_L4[2][2][5][5] = L3;

                        Curvature_Higgs_L4[2][2][6][6] = L2;

                        Curvature_Higgs_L4[2][2][7][7] = L2;

                        Curvature_Higgs_L4[2][2][8][8] = L8;

                        Curvature_Higgs_L4[2][3][0][1] = L5;

                        Curvature_Higgs_L4[2][3][1][0] = L5;

                        Curvature_Higgs_L4[2][3][2][3] = L2;

                        Curvature_Higgs_L4[2][3][3][2] = L2;

                        Curvature_Higgs_L4[2][4][0][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][4][1][7] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][4][2][4] = L3;

                        Curvature_Higgs_L4[2][4][4][2] = L3;

                        Curvature_Higgs_L4[2][4][6][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][4][7][1] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][5][0][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][5][1][6] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[2][5][2][5] = L3;

                        Curvature_Higgs_L4[2][5][5][2] = L3;

                        Curvature_Higgs_L4[2][5][6][1] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[2][5][7][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][6][0][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][6][1][5] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[2][6][2][6] = L2;

                        Curvature_Higgs_L4[2][6][4][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][6][5][1] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[2][6][6][2] = L2;

                        Curvature_Higgs_L4[2][7][0][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][7][1][4] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][7][2][7] = L2;

                        Curvature_Higgs_L4[2][7][4][1] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][7][5][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[2][7][7][2] = L2;

                        Curvature_Higgs_L4[2][8][2][8] = L8;

                        Curvature_Higgs_L4[2][8][8][2] = L8;

                        Curvature_Higgs_L4[3][0][0][3] = L3 + L4 - L5;

                        Curvature_Higgs_L4[3][0][1][2] = L5;

                        Curvature_Higgs_L4[3][0][2][1] = L5;

                        Curvature_Higgs_L4[3][0][3][0] = L3 + L4 - L5;

                        Curvature_Higgs_L4[3][0][4][7] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[3][0][5][6] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][0][6][5] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][0][7][4] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[3][1][0][2] = L5;

                        Curvature_Higgs_L4[3][1][1][3] = L3 + L4 + L5;

                        Curvature_Higgs_L4[3][1][2][0] = L5;

                        Curvature_Higgs_L4[3][1][3][1] = L3 + L4 + L5;

                        Curvature_Higgs_L4[3][1][4][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][1][5][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][1][6][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][1][7][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][2][0][1] = L5;

                        Curvature_Higgs_L4[3][2][1][0] = L5;

                        Curvature_Higgs_L4[3][2][2][3] = L2;

                        Curvature_Higgs_L4[3][2][3][2] = L2;

                        Curvature_Higgs_L4[3][3][0][0] = L3 + L4 - L5;

                        Curvature_Higgs_L4[3][3][1][1] = L3 + L4 + L5;

                        Curvature_Higgs_L4[3][3][2][2] = L2;

                        Curvature_Higgs_L4[3][3][3][3] = 3 * L2;

                        Curvature_Higgs_L4[3][3][4][4] = L3;

                        Curvature_Higgs_L4[3][3][5][5] = L3;

                        Curvature_Higgs_L4[3][3][6][6] = L2;

                        Curvature_Higgs_L4[3][3][7][7] = L2;

                        Curvature_Higgs_L4[3][3][8][8] = L8;

                        Curvature_Higgs_L4[3][4][0][7] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[3][4][1][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][4][3][4] = L3;

                        Curvature_Higgs_L4[3][4][4][3] = L3;

                        Curvature_Higgs_L4[3][4][6][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][4][7][0] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[3][5][0][6] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][5][1][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][5][3][5] = L3;

                        Curvature_Higgs_L4[3][5][5][3] = L3;

                        Curvature_Higgs_L4[3][5][6][0] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][5][7][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][6][0][5] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][6][1][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][6][3][6] = L2;

                        Curvature_Higgs_L4[3][6][4][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][6][5][0] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][6][6][3] = L2;

                        Curvature_Higgs_L4[3][7][0][4] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[3][7][1][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][7][3][7] = L2;

                        Curvature_Higgs_L4[3][7][4][0] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[3][7][5][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[3][7][7][3] = L2;

                        Curvature_Higgs_L4[3][8][3][8] = L8;

                        Curvature_Higgs_L4[3][8][8][3] = L8;

                        Curvature_Higgs_L4[4][0][0][4] = L1;

                        Curvature_Higgs_L4[4][0][2][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][0][3][7] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[4][0][4][0] = L1;

                        Curvature_Higgs_L4[4][0][6][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][0][7][3] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[4][1][1][4] = L1;

                        Curvature_Higgs_L4[4][1][2][7] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][1][3][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][1][4][1] = L1;

                        Curvature_Higgs_L4[4][1][6][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][1][7][2] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][2][0][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][2][1][7] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][2][2][4] = L3;

                        Curvature_Higgs_L4[4][2][4][2] = L3;

                        Curvature_Higgs_L4[4][2][6][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][2][7][1] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][3][0][7] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[4][3][1][6] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][3][3][4] = L3;

                        Curvature_Higgs_L4[4][3][4][3] = L3;

                        Curvature_Higgs_L4[4][3][6][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][3][7][0] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[4][4][0][0] = L1;

                        Curvature_Higgs_L4[4][4][1][1] = L1;

                        Curvature_Higgs_L4[4][4][2][2] = L3;

                        Curvature_Higgs_L4[4][4][3][3] = L3;

                        Curvature_Higgs_L4[4][4][4][4] = 3 * L1;

                        Curvature_Higgs_L4[4][4][5][5] = L1;

                        Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + L5;

                        Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - L5;

                        Curvature_Higgs_L4[4][4][8][8] = L7;

                        Curvature_Higgs_L4[4][5][4][5] = L1;

                        Curvature_Higgs_L4[4][5][5][4] = L1;

                        Curvature_Higgs_L4[4][5][6][7] = L5;

                        Curvature_Higgs_L4[4][5][7][6] = L5;

                        Curvature_Higgs_L4[4][6][0][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][6][1][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][6][2][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][6][3][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][6][4][6] = L3 + L4 + L5;

                        Curvature_Higgs_L4[4][6][5][7] = L5;

                        Curvature_Higgs_L4[4][6][6][4] = L3 + L4 + L5;

                        Curvature_Higgs_L4[4][6][7][5] = L5;

                        Curvature_Higgs_L4[4][7][0][3] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[4][7][1][2] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][7][2][1] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[4][7][3][0] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[4][7][4][7] = L3 + L4 - L5;

                        Curvature_Higgs_L4[4][7][5][6] = L5;

                        Curvature_Higgs_L4[4][7][6][5] = L5;

                        Curvature_Higgs_L4[4][7][7][4] = L3 + L4 - L5;

                        Curvature_Higgs_L4[4][8][4][8] = L7;

                        Curvature_Higgs_L4[4][8][8][4] = L7;

                        Curvature_Higgs_L4[5][0][0][5] = L1;

                        Curvature_Higgs_L4[5][0][2][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][0][3][6] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][0][5][0] = L1;

                        Curvature_Higgs_L4[5][0][6][3] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][0][7][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][1][1][5] = L1;

                        Curvature_Higgs_L4[5][1][2][6] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[5][1][3][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][1][5][1] = L1;

                        Curvature_Higgs_L4[5][1][6][2] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[5][1][7][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][2][0][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][2][1][6] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[5][2][2][5] = L3;

                        Curvature_Higgs_L4[5][2][5][2] = L3;

                        Curvature_Higgs_L4[5][2][6][1] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[5][2][7][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][3][0][6] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][3][1][7] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][3][3][5] = L3;

                        Curvature_Higgs_L4[5][3][5][3] = L3;

                        Curvature_Higgs_L4[5][3][6][0] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][3][7][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][4][4][5] = L1;

                        Curvature_Higgs_L4[5][4][5][4] = L1;

                        Curvature_Higgs_L4[5][4][6][7] = L5;

                        Curvature_Higgs_L4[5][4][7][6] = L5;

                        Curvature_Higgs_L4[5][5][0][0] = L1;

                        Curvature_Higgs_L4[5][5][1][1] = L1;

                        Curvature_Higgs_L4[5][5][2][2] = L3;

                        Curvature_Higgs_L4[5][5][3][3] = L3;

                        Curvature_Higgs_L4[5][5][4][4] = L1;

                        Curvature_Higgs_L4[5][5][5][5] = 3 * L1;

                        Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - L5;

                        Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + L5;

                        Curvature_Higgs_L4[5][5][8][8] = L7;

                        Curvature_Higgs_L4[5][6][0][3] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][6][1][2] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[5][6][2][1] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[5][6][3][0] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][6][4][7] = L5;

                        Curvature_Higgs_L4[5][6][5][6] = L3 + L4 - L5;

                        Curvature_Higgs_L4[5][6][6][5] = L3 + L4 - L5;

                        Curvature_Higgs_L4[5][6][7][4] = L5;

                        Curvature_Higgs_L4[5][7][0][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][7][1][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][7][2][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][7][3][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[5][7][4][6] = L5;

                        Curvature_Higgs_L4[5][7][5][7] = L3 + L4 + L5;

                        Curvature_Higgs_L4[5][7][6][4] = L5;

                        Curvature_Higgs_L4[5][7][7][5] = L3 + L4 + L5;

                        Curvature_Higgs_L4[5][8][5][8] = L7;

                        Curvature_Higgs_L4[5][8][8][5] = L7;

                        Curvature_Higgs_L4[6][0][0][6] = L3;

                        Curvature_Higgs_L4[6][0][2][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][0][3][5] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][0][4][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][0][5][3] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][0][6][0] = L3;

                        Curvature_Higgs_L4[6][1][1][6] = L3;

                        Curvature_Higgs_L4[6][1][2][5] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[6][1][3][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][1][4][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][1][5][2] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[6][1][6][1] = L3;

                        Curvature_Higgs_L4[6][2][0][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][2][1][5] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[6][2][2][6] = L2;

                        Curvature_Higgs_L4[6][2][4][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][2][5][1] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[6][2][6][2] = L2;

                        Curvature_Higgs_L4[6][3][0][5] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][3][1][4] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][3][3][6] = L2;

                        Curvature_Higgs_L4[6][3][4][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][3][5][0] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][3][6][3] = L2;

                        Curvature_Higgs_L4[6][4][0][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][4][1][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][4][2][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][4][3][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][4][4][6] = L3 + L4 + L5;

                        Curvature_Higgs_L4[6][4][5][7] = L5;

                        Curvature_Higgs_L4[6][4][6][4] = L3 + L4 + L5;

                        Curvature_Higgs_L4[6][4][7][5] = L5;

                        Curvature_Higgs_L4[6][5][0][3] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][5][1][2] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[6][5][2][1] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[6][5][3][0] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[6][5][4][7] = L5;

                        Curvature_Higgs_L4[6][5][5][6] = L3 + L4 - L5;

                        Curvature_Higgs_L4[6][5][6][5] = L3 + L4 - L5;

                        Curvature_Higgs_L4[6][5][7][4] = L5;

                        Curvature_Higgs_L4[6][6][0][0] = L3;

                        Curvature_Higgs_L4[6][6][1][1] = L3;

                        Curvature_Higgs_L4[6][6][2][2] = L2;

                        Curvature_Higgs_L4[6][6][3][3] = L2;

                        Curvature_Higgs_L4[6][6][4][4] = L3 + L4 + L5;

                        Curvature_Higgs_L4[6][6][5][5] = L3 + L4 - L5;

                        Curvature_Higgs_L4[6][6][6][6] = 3 * L2;

                        Curvature_Higgs_L4[6][6][7][7] = L2;

                        Curvature_Higgs_L4[6][6][8][8] = L8;

                        Curvature_Higgs_L4[6][7][4][5] = L5;

                        Curvature_Higgs_L4[6][7][5][4] = L5;

                        Curvature_Higgs_L4[6][7][6][7] = L2;

                        Curvature_Higgs_L4[6][7][7][6] = L2;

                        Curvature_Higgs_L4[6][8][6][8] = L8;

                        Curvature_Higgs_L4[6][8][8][6] = L8;

                        Curvature_Higgs_L4[7][0][0][7] = L3;

                        Curvature_Higgs_L4[7][0][2][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][0][3][4] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[7][0][4][3] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[7][0][5][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][0][7][0] = L3;

                        Curvature_Higgs_L4[7][1][1][7] = L3;

                        Curvature_Higgs_L4[7][1][2][4] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][1][3][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][1][4][2] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][1][5][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][1][7][1] = L3;

                        Curvature_Higgs_L4[7][2][0][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][2][1][4] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][2][2][7] = L2;

                        Curvature_Higgs_L4[7][2][4][1] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][2][5][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][2][7][2] = L2;

                        Curvature_Higgs_L4[7][3][0][4] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[7][3][1][5] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][3][3][7] = L2;

                        Curvature_Higgs_L4[7][3][4][0] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[7][3][5][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][3][7][3] = L2;

                        Curvature_Higgs_L4[7][4][0][3] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[7][4][1][2] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][4][2][1] = -L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][4][3][0] = L4 / 2. - L5 / 2.;

                        Curvature_Higgs_L4[7][4][4][7] = L3 + L4 - L5;

                        Curvature_Higgs_L4[7][4][5][6] = L5;

                        Curvature_Higgs_L4[7][4][6][5] = L5;

                        Curvature_Higgs_L4[7][4][7][4] = L3 + L4 - L5;

                        Curvature_Higgs_L4[7][5][0][2] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][5][1][3] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][5][2][0] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][5][3][1] = L4 / 2. + L5 / 2.;

                        Curvature_Higgs_L4[7][5][4][6] = L5;

                        Curvature_Higgs_L4[7][5][5][7] = L3 + L4 + L5;

                        Curvature_Higgs_L4[7][5][6][4] = L5;

                        Curvature_Higgs_L4[7][5][7][5] = L3 + L4 + L5;

                        Curvature_Higgs_L4[7][6][4][5] = L5;

                        Curvature_Higgs_L4[7][6][5][4] = L5;

                        Curvature_Higgs_L4[7][6][6][7] = L2;

                        Curvature_Higgs_L4[7][6][7][6] = L2;

                        Curvature_Higgs_L4[7][7][0][0] = L3;

                        Curvature_Higgs_L4[7][7][1][1] = L3;

                        Curvature_Higgs_L4[7][7][2][2] = L2;

                        Curvature_Higgs_L4[7][7][3][3] = L2;

                        Curvature_Higgs_L4[7][7][4][4] = L3 + L4 - L5;

                        Curvature_Higgs_L4[7][7][5][5] = L3 + L4 + L5;

                        Curvature_Higgs_L4[7][7][6][6] = L2;

                        Curvature_Higgs_L4[7][7][7][7] = 3 * L2;

                        Curvature_Higgs_L4[7][7][8][8] = L8;

                        Curvature_Higgs_L4[7][8][7][8] = L8;

                        Curvature_Higgs_L4[7][8][8][7] = L8;

                        Curvature_Higgs_L4[8][0][0][8] = L7;

                        Curvature_Higgs_L4[8][0][8][0] = L7;

                        Curvature_Higgs_L4[8][1][1][8] = L7;

                        Curvature_Higgs_L4[8][1][8][1] = L7;

                        Curvature_Higgs_L4[8][2][2][8] = L8;

                        Curvature_Higgs_L4[8][2][8][2] = L8;

                        Curvature_Higgs_L4[8][3][3][8] = L8;

                        Curvature_Higgs_L4[8][3][8][3] = L8;

                        Curvature_Higgs_L4[8][4][4][8] = L7;

                        Curvature_Higgs_L4[8][4][8][4] = L7;

                        Curvature_Higgs_L4[8][5][5][8] = L7;

                        Curvature_Higgs_L4[8][5][8][5] = L7;

                        Curvature_Higgs_L4[8][6][6][8] = L8;

                        Curvature_Higgs_L4[8][6][8][6] = L8;

                        Curvature_Higgs_L4[8][7][7][8] = L8;

                        Curvature_Higgs_L4[8][7][8][7] = L8;

                        Curvature_Higgs_L4[8][8][0][0] = L7;

                        Curvature_Higgs_L4[8][8][1][1] = L7;

                        Curvature_Higgs_L4[8][8][2][2] = L8;

                        Curvature_Higgs_L4[8][8][3][3] = L8;

                        Curvature_Higgs_L4[8][8][4][4] = L7;

                        Curvature_Higgs_L4[8][8][5][5] = L7;

                        Curvature_Higgs_L4[8][8][6][6] = L8;

                        Curvature_Higgs_L4[8][8][7][7] = L8;

                        Curvature_Higgs_L4[8][8][8][8] = 6 * L6;

                        Curvature_Gauge_G2H2[0][0][0][0] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][1][1] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][2][2] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][3][3] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][4][4] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][5][5] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][6][6] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][0][7][7] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[0][3][0][4] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][1][5] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][2][6] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][3][7] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][4][0] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][5][1] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][6][2] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[0][3][7][3] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][1][0][0] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][1][1] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][2][2] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][3][3] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][4][4] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][5][5] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][6][6] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][1][7][7] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[1][3][0][5] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][1][4] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][2][7] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][3][6] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][4][1] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][5][0] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][6][3] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[1][3][7][2] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][2][0][0] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][1][1] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][2][2] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][3][3] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][4][4] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][5][5] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][6][6] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][2][7][7] = std::pow(C_g, 2) / 2.;
                        Curvature_Gauge_G2H2[2][3][0][0] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][1][1] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][2][2] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][3][3] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][4][4] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][5][5] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][6][6] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[2][3][7][7] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][0][4] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][1][5] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][2][6] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][3][7] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][4][0] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][5][1] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][6][2] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][0][7][3] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][0][5] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][1][4] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][2][7] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][3][6] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][4][1] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][5][0] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][6][3] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][1][7][2] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][0][0] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][1][1] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][2][2] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][3][3] = (C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][4][4] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][5][5] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][6][6] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][2][7][7] = -(C_gs * C_g) / 2.;
                        Curvature_Gauge_G2H2[3][3][0][0] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][1][1] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][2][2] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][3][3] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][4][4] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][5][5] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][6][6] = std::pow(C_gs, 2) / 2.;
                        Curvature_Gauge_G2H2[3][3][7][7] = std::pow(C_gs, 2) / 2.;

                        double v1 = C_vev0;

                        for (std::size_t i = 0; i < NQuarks; i++)
                        {
                                for (std::size_t j = 0; j < NQuarks; j++)
                                {
                                        for (std::size_t k = 0; k < NHiggs; k++)
                                        {
                                                Curvature_Quark_F2H1[i][j][k] = 0;
                                        }
                                }
                        }

                        std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
                        V11 = C_Vud;
                        V12 = C_Vus;
                        V13 = C_Vub;
                        V21 = C_Vcd;
                        V22 = C_Vcs;
                        V23 = C_Vcb;
                        V31 = C_Vtd;
                        V32 = C_Vts;
                        V33 = C_Vtb;

                        std::complex<double> II(0, 1); // define imaginary unit

                        Curvature_Quark_F2H1[0][1][4] = 0.1e1 / v1 * C_MassUp;
                        Curvature_Quark_F2H1[0][1][5] = -II / v1 * C_MassUp;
                        Curvature_Quark_F2H1[0][3][0] = 0.1e1 / v1 * C_MassDown * V11;
                        Curvature_Quark_F2H1[0][3][1] = II / v1 * C_MassDown * V11;
                        Curvature_Quark_F2H1[0][7][0] = V12 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[0][7][1] = II * V12 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[0][11][0] = V13 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[0][11][1] = II * V13 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[1][0][4] = 0.1e1 / v1 * C_MassUp;
                        Curvature_Quark_F2H1[1][0][5] = -II / v1 * C_MassUp;
                        Curvature_Quark_F2H1[1][2][0] = -0.1e1 / v1 * C_MassUp * std::conj(V11);
                        Curvature_Quark_F2H1[1][2][1] = II / v1 * C_MassUp * std::conj(V11);
                        Curvature_Quark_F2H1[1][6][0] = -0.1e1 / v1 * C_MassUp * std::conj(V12);
                        Curvature_Quark_F2H1[1][6][1] = II / v1 * C_MassUp * std::conj(V12);
                        Curvature_Quark_F2H1[1][10][0] = -0.1e1 / v1 * C_MassUp * std::conj(V13);
                        Curvature_Quark_F2H1[1][10][1] = II / v1 * C_MassUp * std::conj(V13);
                        Curvature_Quark_F2H1[2][1][0] = -0.1e1 / v1 * C_MassUp * std::conj(V11);
                        Curvature_Quark_F2H1[2][1][1] = II / v1 * C_MassUp * std::conj(V11);
                        Curvature_Quark_F2H1[2][3][4] = 0.1e1 / v1 * C_MassDown;
                        Curvature_Quark_F2H1[2][3][5] = II / v1 * C_MassDown;
                        Curvature_Quark_F2H1[2][5][0] = -0.1e1 / v1 * C_MassCharm * std::conj(V21);
                        Curvature_Quark_F2H1[2][5][1] = II / v1 * C_MassCharm * std::conj(V21);
                        Curvature_Quark_F2H1[2][9][0] = -0.1e1 / v1 * C_MassTop * std::conj(V31);
                        Curvature_Quark_F2H1[2][9][1] = II / v1 * C_MassTop * std::conj(V31);
                        Curvature_Quark_F2H1[3][0][0] = 0.1e1 / v1 * C_MassDown * V11;
                        Curvature_Quark_F2H1[3][0][1] = II / v1 * C_MassDown * V11;
                        Curvature_Quark_F2H1[3][2][4] = 0.1e1 / v1 * C_MassDown;
                        Curvature_Quark_F2H1[3][2][5] = II / v1 * C_MassDown;
                        Curvature_Quark_F2H1[3][4][0] = V21 / v1 * C_MassDown;
                        Curvature_Quark_F2H1[3][4][1] = II * V21 / v1 * C_MassDown;
                        Curvature_Quark_F2H1[3][8][0] = 0.1e1 / v1 * C_MassDown * V31;
                        Curvature_Quark_F2H1[3][8][1] = II / v1 * C_MassDown * V31;
                        Curvature_Quark_F2H1[4][3][0] = V21 / v1 * C_MassDown;
                        Curvature_Quark_F2H1[4][3][1] = II * V21 / v1 * C_MassDown;
                        Curvature_Quark_F2H1[4][5][4] = 0.1e1 / v1 * C_MassCharm;
                        Curvature_Quark_F2H1[4][5][5] = -II / v1 * C_MassCharm;
                        Curvature_Quark_F2H1[4][7][0] = V22 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[4][7][1] = II * V22 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[4][11][0] = V23 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[4][11][1] = II * V23 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[5][2][0] = -0.1e1 / v1 * C_MassCharm * std::conj(V21);
                        Curvature_Quark_F2H1[5][2][1] = II / v1 * C_MassCharm * std::conj(V21);
                        Curvature_Quark_F2H1[5][4][4] = 0.1e1 / v1 * C_MassCharm;
                        Curvature_Quark_F2H1[5][4][5] = -II / v1 * C_MassCharm;
                        Curvature_Quark_F2H1[5][6][0] = -0.1e1 / v1 * C_MassCharm * std::conj(V22);
                        Curvature_Quark_F2H1[5][6][1] = II / v1 * C_MassCharm * std::conj(V22);
                        Curvature_Quark_F2H1[5][10][0] = -0.1e1 / v1 * C_MassCharm * std::conj(V23);
                        Curvature_Quark_F2H1[5][10][1] = II / v1 * C_MassCharm * std::conj(V23);
                        Curvature_Quark_F2H1[6][1][0] = -0.1e1 / v1 * C_MassUp * std::conj(V12);
                        Curvature_Quark_F2H1[6][1][1] = II / v1 * C_MassUp * std::conj(V12);
                        Curvature_Quark_F2H1[6][5][0] = -0.1e1 / v1 * C_MassCharm * std::conj(V22);
                        Curvature_Quark_F2H1[6][5][1] = II / v1 * C_MassCharm * std::conj(V22);
                        Curvature_Quark_F2H1[6][7][4] = 0.1e1 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[6][7][5] = II / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[6][9][0] = -0.1e1 / v1 * C_MassTop * std::conj(V32);
                        Curvature_Quark_F2H1[6][9][1] = II / v1 * C_MassTop * std::conj(V32);
                        Curvature_Quark_F2H1[7][0][0] = V12 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][0][1] = II * V12 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][4][0] = V22 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][4][1] = II * V22 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][6][4] = 0.1e1 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][6][5] = II / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][8][0] = V32 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[7][8][1] = II * V32 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[8][3][0] = 0.1e1 / v1 * C_MassDown * V31;
                        Curvature_Quark_F2H1[8][3][1] = II / v1 * C_MassDown * V31;
                        Curvature_Quark_F2H1[8][7][0] = V32 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[8][7][1] = II * V32 / v1 * C_MassStrange;
                        Curvature_Quark_F2H1[8][9][4] = 0.1e1 / v1 * C_MassTop;
                        Curvature_Quark_F2H1[8][9][5] = -II / v1 * C_MassTop;
                        Curvature_Quark_F2H1[8][11][0] = V33 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[8][11][1] = II * V33 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[9][2][0] = -0.1e1 / v1 * C_MassTop * std::conj(V31);
                        Curvature_Quark_F2H1[9][2][1] = II / v1 * C_MassTop * std::conj(V31);
                        Curvature_Quark_F2H1[9][6][0] = -0.1e1 / v1 * C_MassTop * std::conj(V32);
                        Curvature_Quark_F2H1[9][6][1] = II / v1 * C_MassTop * std::conj(V32);
                        Curvature_Quark_F2H1[9][8][4] = 0.1e1 / v1 * C_MassTop;
                        Curvature_Quark_F2H1[9][8][5] = -II / v1 * C_MassTop;
                        Curvature_Quark_F2H1[9][10][0] = -0.1e1 / v1 * C_MassTop * std::conj(V33);
                        Curvature_Quark_F2H1[9][10][1] = II / v1 * C_MassTop * std::conj(V33);
                        Curvature_Quark_F2H1[10][1][0] = -0.1e1 / v1 * C_MassUp * std::conj(V13);
                        Curvature_Quark_F2H1[10][1][1] = II / v1 * C_MassUp * std::conj(V13);
                        Curvature_Quark_F2H1[10][5][0] = -0.1e1 / v1 * C_MassCharm * std::conj(V23);
                        Curvature_Quark_F2H1[10][5][1] = II / v1 * C_MassCharm * std::conj(V23);
                        Curvature_Quark_F2H1[10][9][0] = -0.1e1 / v1 * C_MassTop * std::conj(V33);
                        Curvature_Quark_F2H1[10][9][1] = II / v1 * C_MassTop * std::conj(V33);
                        Curvature_Quark_F2H1[10][11][4] = 0.1e1 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[10][11][5] = II / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][0][0] = V13 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][0][1] = II * V13 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][4][0] = V23 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][4][1] = II * V23 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][8][0] = V33 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][8][1] = II * V33 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][10][4] = 0.1e1 / v1 * C_MassBottom;
                        Curvature_Quark_F2H1[11][10][5] = II / v1 * C_MassBottom;

                        for (std::size_t i = 0; i < NLepton; i++)
                        {
                                for (std::size_t j = 0; j < NLepton; j++)
                                {
                                        for (std::size_t k = 0; k < NHiggs; k++)
                                        {
                                                Curvature_Lepton_F2H1[i][j][k] = 0;
                                        }
                                }
                        }

                        Curvature_Lepton_F2H1[0][1][4] = 0.1e1 / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[0][1][5] = II / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[1][0][4] = 0.1e1 / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[1][0][5] = II / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[1][6][1] = II / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[2][3][4] = 0.1e1 / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[2][3][5] = II / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[3][2][4] = 0.1e1 / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[3][2][5] = II / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[3][7][1] = II / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[4][5][4] = 0.1e1 / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[4][5][5] = II / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[5][4][4] = 0.1e1 / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[5][4][5] = II / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[5][8][1] = II / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[6][1][1] = II / v1 * C_MassElectron;
                        Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[7][3][1] = II / v1 * C_MassMu;
                        Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / v1 * C_MassTau;
                        Curvature_Lepton_F2H1[8][5][1] = II / v1 * C_MassTau;

                        // /*
                        //         Fermion Sector!
                        // */
                        // std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
                        // V11 = C_Vud;
                        // V12 = C_Vus;
                        // V13 = C_Vub;
                        // V21 = C_Vcd;
                        // V22 = C_Vcs;
                        // V23 = C_Vcb;
                        // V31 = C_Vtd;
                        // V32 = C_Vts;
                        // V33 = C_Vtb;

                        // std::complex<double> II(0, 1);
                        // double Yu = std::sqrt(2) * C_MassUp / C_vev0;
                        // double Yd = std::sqrt(2) * C_MassDown / C_vev0;
                        // double Ys = std::sqrt(2) * C_MassStrange / C_vev0;
                        // double Yc = std::sqrt(2) * C_MassCharm / C_vev0;
                        // double Yb = std::sqrt(2) * C_MassBottom / C_vev0;
                        // double Yt = std::sqrt(2) * C_MassTop / C_vev0;
                        // double Ye = std::sqrt(2) * C_MassElectron / C_vev0;
                        // double Ymu = std::sqrt(2) * C_MassMu / C_vev0;
                        // double Ytau = std::sqrt(2) * C_MassTau / C_vev0;

                        // // /* CKM = 1
                        // Curvature_Quark_F2H1[0][6][4] = -(Yu / std::sqrt(2));

                        // Curvature_Quark_F2H1[0][6][5] = (II * Yu) / std::sqrt(2);

                        // Curvature_Quark_F2H1[0][9][0] = Yu / std::sqrt(2);

                        // Curvature_Quark_F2H1[0][9][1] = -((II * Yu) / std::sqrt(2));

                        // Curvature_Quark_F2H1[1][7][4] = -(Yc / std::sqrt(2));

                        // Curvature_Quark_F2H1[1][7][5] = (II * Yc) / std::sqrt(2);

                        // Curvature_Quark_F2H1[1][10][0] = Yc / std::sqrt(2);

                        // Curvature_Quark_F2H1[1][10][1] = -((II * Yc) / std::sqrt(2));

                        // Curvature_Quark_F2H1[2][8][4] = -(Yt / std::sqrt(2));

                        // Curvature_Quark_F2H1[2][8][5] = (II * Yt) / std::sqrt(2);

                        // Curvature_Quark_F2H1[2][11][0] = Yt / std::sqrt(2);

                        // Curvature_Quark_F2H1[2][11][1] = -((II * Yt) / std::sqrt(2));

                        // Curvature_Quark_F2H1[3][6][0] = -(Yd / std::sqrt(2));

                        // Curvature_Quark_F2H1[3][6][1] = -((II * Yd) / std::sqrt(2));

                        // Curvature_Quark_F2H1[3][9][4] = -(Yd / std::sqrt(2));

                        // Curvature_Quark_F2H1[3][9][5] = -((II * Yd) / std::sqrt(2));

                        // Curvature_Quark_F2H1[4][7][0] = -(Ys / std::sqrt(2));

                        // Curvature_Quark_F2H1[4][7][1] = -((II * Ys) / std::sqrt(2));

                        // Curvature_Quark_F2H1[4][10][4] = -(Ys / std::sqrt(2));

                        // Curvature_Quark_F2H1[4][10][5] = -((II * Ys) / std::sqrt(2));

                        // Curvature_Quark_F2H1[5][8][0] = -(Yb / std::sqrt(2));

                        // Curvature_Quark_F2H1[5][8][1] = -((II * Yb) / std::sqrt(2));

                        // Curvature_Quark_F2H1[5][11][4] = -(Yb / std::sqrt(2));

                        // Curvature_Quark_F2H1[5][11][5] = -((II * Yb) / std::sqrt(2));

                        // Curvature_Quark_F2H1[6][0][4] = -(Yu / std::sqrt(2));

                        // Curvature_Quark_F2H1[6][0][5] = (II * Yu) / std::sqrt(2);

                        // Curvature_Quark_F2H1[6][3][0] = -(Yd / std::sqrt(2));

                        // Curvature_Quark_F2H1[6][3][1] = -((II * Yd) / std::sqrt(2));

                        // Curvature_Quark_F2H1[7][1][4] = -(Yc / std::sqrt(2));

                        // Curvature_Quark_F2H1[7][1][5] = (II * Yc) / std::sqrt(2);

                        // Curvature_Quark_F2H1[7][4][0] = -(Ys / std::sqrt(2));

                        // Curvature_Quark_F2H1[7][4][1] = -((II * Ys) / std::sqrt(2));

                        // Curvature_Quark_F2H1[8][2][4] = -(Yt / std::sqrt(2));

                        // Curvature_Quark_F2H1[8][2][5] = (II * Yt) / std::sqrt(2);

                        // Curvature_Quark_F2H1[8][5][0] = -(Yb / std::sqrt(2));

                        // Curvature_Quark_F2H1[8][5][1] = -((II * Yb) / std::sqrt(2));

                        // Curvature_Quark_F2H1[9][0][0] = Yu / std::sqrt(2);

                        // Curvature_Quark_F2H1[9][0][1] = -((II * Yu) / std::sqrt(2));

                        // Curvature_Quark_F2H1[9][3][4] = -(Yd / std::sqrt(2));

                        // Curvature_Quark_F2H1[9][3][5] = -((II * Yd) / std::sqrt(2));

                        // Curvature_Quark_F2H1[10][1][0] = Yc / std::sqrt(2);

                        // Curvature_Quark_F2H1[10][1][1] = -((II * Yc) / std::sqrt(2));

                        // Curvature_Quark_F2H1[10][4][4] = -(Ys / std::sqrt(2));

                        // Curvature_Quark_F2H1[10][4][5] = -((II * Ys) / std::sqrt(2));

                        // Curvature_Quark_F2H1[11][2][0] = Yt / std::sqrt(2);

                        // Curvature_Quark_F2H1[11][2][1] = -((II * Yt) / std::sqrt(2));

                        // Curvature_Quark_F2H1[11][5][4] = -(Yb / std::sqrt(2));

                        // Curvature_Quark_F2H1[11][5][5] = -((II * Yb) / std::sqrt(2));
                        // // */

                        // Curvature_Lepton_F2H1[0][1][4] = -(Ye / std::sqrt(2));
                        // Curvature_Lepton_F2H1[0][1][5] = -((II * Ye) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[1][0][4] = -(Ye / std::sqrt(2));
                        // Curvature_Lepton_F2H1[1][0][5] = -((II * Ye) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[1][6][0] = -(Ye / std::sqrt(2));
                        // Curvature_Lepton_F2H1[1][6][1] = -((II * Ye) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[2][3][4] = -(Ymu / std::sqrt(2));
                        // Curvature_Lepton_F2H1[2][3][5] = -((II * Ymu) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[3][2][4] = -(Ymu / std::sqrt(2));
                        // Curvature_Lepton_F2H1[3][2][5] = -((II * Ymu) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[3][7][0] = -(Ymu / std::sqrt(2));
                        // Curvature_Lepton_F2H1[3][7][1] = -((II * Ymu) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[4][5][4] = -(Ytau / std::sqrt(2));
                        // Curvature_Lepton_F2H1[4][5][5] = -((II * Ytau) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[5][4][4] = -(Ytau / std::sqrt(2));
                        // Curvature_Lepton_F2H1[5][4][5] = -((II * Ytau) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[5][8][0] = -(Ytau / std::sqrt(2));
                        // Curvature_Lepton_F2H1[5][8][1] = -((II * Ytau) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[6][1][0] = -(Ye / std::sqrt(2));
                        // Curvature_Lepton_F2H1[6][1][1] = -((II * Ye) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[7][3][0] = -(Ymu / std::sqrt(2));
                        // Curvature_Lepton_F2H1[7][3][1] = -((II * Ymu) / std::sqrt(2));
                        // Curvature_Lepton_F2H1[8][5][0] = -(Ytau / std::sqrt(2));
                        // Curvature_Lepton_F2H1[8][5][1] = -((II * Ytau) / std::sqrt(2));
                        SetCurvatureDone = true;
                }

                bool Class_CPDM_ImL5::CalculateDebyeSimplified()
                {
                        return false;
                        /*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
                }

                bool Class_CPDM_ImL5::CalculateDebyeGaugeSimplified()
                {
                        /*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */

                        return false;
                }
                double Class_CPDM_ImL5::VTreeSimplified(const std::vector<double> &v) const
                {
                        if (not UseVTreeSimplified)
                                return 0;
                        (void)v;
                        double res = 0;
                        return res;
                }

                double Class_CPDM_ImL5::VCounterSimplified(const std::vector<double> &v) const
                {
                        if (not UseVCounterSimplified)
                                return 0;
                        double res = 0;
                        (void)v;
                        return res;
                }

                void Class_CPDM_ImL5::Debugging(const std::vector<double> &input, std::vector<double> &output) const
                {
                        (void)input;
                        (void)output;
                }

        } // namespace Models
} // namespace BSMPT
