#include "Calculations.h"

#include<vector>
#include<math.h>
#include<cmath>
#include<algorithm>

using namespace std;

vector<double> weights;
vector<double> points;

//Single-point integration scheme
vector<double> weights1 = {2.0};
vector<double> points1 = {0.0};

//2-point integration scheme
vector<double> weights2 = {1.0, 1.0};
vector<double> points2 = {-(1.0 / sqrt(3.0)), 1.0 / sqrt(3.0)};

//3-point integration scheme
vector<double> weights3 = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
vector<double> points3 = {sqrt(3.0 / 5.0), 0.0, -sqrt(3.0 / 5.0)};

//4-point integration scheme
vector<double> weights4 = {(18.0 - sqrt((30.0))) / 36.0, (18.0 + sqrt((30.0))) / 36.0,
                        (18.0+sqrt((30.0)))/36.0, (18.0-sqrt((30.0)))/36.0};
vector<double> points4 = {-sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),
                          sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0)), sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0))};;

vector<int> BC;
vector<vector<double>> results_ksi;
vector<vector<double>> results_eta;
vector<vector<double>> jacobian_matrices;
vector<vector<vector<double>>> ElementsCalculations;
vector<vector<double>> determinantsJ;
vector<vector<double>> inverseDJ;
vector<vector<double>> vdndxpc;
vector<vector<double>> vdndypc;
vector<vector<vector<double>>> setwdndx;
vector<vector<vector<double>>> setwdndy;
vector<vector<vector<double>>> elementsMatricesH;
vector<vector<double>> compH;
vector<vector<vector<double>>> matricesHBC;
vector<vector<double>> vectorsP;
vector<vector<vector<double>>> finalH;
vector<vector<double>> aggregatedHmatrix;
vector<double> aggregatedPvector;
vector<vector<double>> results_N;
vector<vector<double>> results_C;
vector<vector<vector<double>>> results_CELEM;
vector<vector<double>> aggregatedCmatrix;
//----- FINAL RESULTS -----
vector<vector<double>> TemperaturesSet;
vector<vector<double>> MinMaxSet;

void integration_scheme(vector<double> chosen_weights, vector<double> chosen_points){
    for (int i = 0; i < chosen_weights.size(); ++i) {
        weights.push_back(chosen_weights[i]);
        points.push_back(chosen_points[i]);
    }
}

int BCs(int N){
    int BCs; int aux;
    aux  = N - (sqrt(N) * 2) - ((sqrt(N) - 2) * 2);
    BCs = N - aux;
    return BCs;
}

void integrate_ksi(vector<double> weight, vector<double> point){

    double N1 = 0.0;    double N2 = 0.0;    double N3 = 0.0;    double N4 = 0.0;
    results_ksi.clear();
    vector<double> intp;

    for (int i = 0; i < point.size(); i++) {
        for (int j = 0; j < point.size(); j++) {
            N1 = -0.25 * (1 - point[i]);
            intp.push_back(N1);
            N2 = 0.25 * (1 - point[i]);
            intp.push_back(N2);
            N3 = 0.25 * (1 + point[i]);
            intp.push_back(N3);
            N4 = -0.25 * (1 + point[i]);
            intp.push_back(N4);
            results_ksi.push_back(intp);
            intp.clear();
        }
    }
    intp.clear();
}

void integrate_eta(vector<double> weight, vector<double> point){

    double N1 = 0.0;    double N2 = 0.0;    double N3 = 0.0;    double N4 = 0.0;
    results_eta.clear();
    vector<double> intp;

    for (int i = 0; i < point.size(); i++) {
        for (int j = 0; j < point.size(); j++) {
            N1 = -0.25 * (1 - point[j]);
            intp.push_back(N1);
            N2 = -0.25 * (1 + point[j]);
            intp.push_back(N2);
            N3 = 0.25 * (1 + point[j]);
            intp.push_back(N3);
            N4 = 0.25 * (1 - point[j]);
            intp.push_back(N4);
            results_eta.push_back(intp);
            intp.clear();
        }
    }
    intp.clear();
}

void jacobian_matrix(Element vectors, GlobalData gd){
    //dx/dksi -> var1    dy/dksi -> var2      dx/deta -> var3       dy/deta -> var4
    double var1 = 0.0; double result1 = 0.0;
    double var2 = 0.0; double result2 = 0.0;
    double var3 = 0.0; double result3 = 0.0;
    double var4 = 0.0; double result4 = 0.0;
    int id_aux = 0;
    vector<double> jacobian;
    vector<double> jacobian1; //= {var1, var2, var3, var4};
    vector<double> jacobian2; //= {var1, var2, var3, var4};
    vector<double> jacobian3; //= {var1, var2, var3, var4};
    vector<double> jacobian4; //= {var1, var2, var3, var4}
    double aux = 0.0; double aux2 = 0.0;
    vector<double> detJfor4intp;
    vector<double> inverseDJfor4intp;
    vector<double> NewComponents;
    vector<vector<double>> NewJacobians;
    vector<vector<double>> NJ4intp;       //set of inverse Jacobian matrices for integration points
    vector<vector<vector<double>>> InverseMatrix;
    vector<vector<vector<double>>> NewMatrix;
    vector<double> vdndx;   //values of dn/dx
    vector<double> vdndy;   //values of dn/dy
    double resultH = 0.0;

    //jacobian_matrices
    for (int i = 0; i < vectors.ID.size(); ++i) {
        for (int j = 0; j < results_ksi.size(); ++j) {
            for (int k = 0; k < results_ksi[j].size(); ++k) {
                id_aux = vectors.ID[i][k];
                var1 += results_ksi[j][k] * vectors.x[id_aux - 1];
                var2 += results_ksi[j][k] * vectors.y[id_aux - 1];
                var3 += results_eta[j][k] * vectors.x[id_aux - 1];
                var4 += results_eta[j][k] * vectors.y[id_aux - 1];
                if (k==3){
                    result1 = var1;
                    result2 = var2;
                    result3 = var3;
                    result4 = var4;
                }
            }
            var1=0; var2=0; var3=0; var4=0;
            jacobian.push_back(result1);
            jacobian.push_back(result2);
            jacobian.push_back(result3);
            jacobian.push_back(result4);
            jacobian_matrices.push_back(jacobian);
            jacobian.clear();
        }
        ElementsCalculations.push_back(jacobian_matrices);
        jacobian_matrices.clear();
    }

    //determinants
    for (int i = 0; i < ElementsCalculations.size(); i++) {
        for (int j = 0; j < results_ksi.size(); ++j) {
            aux = (ElementsCalculations[i][j][0] * ElementsCalculations[i][j][3]) - (ElementsCalculations[i][j][1] * ElementsCalculations[i][j][2]);
            detJfor4intp.push_back(aux);
        }
        determinantsJ.push_back(detJfor4intp);
        detJfor4intp.clear();
    }

    //inverse jacobian
    for (int i = 0; i < determinantsJ.size(); i++) {
        for (int j = 0; j < determinantsJ[i].size(); ++j) {
            aux = 1 / determinantsJ[i][j];
            inverseDJfor4intp.push_back(aux);
        }
        inverseDJ.push_back(inverseDJfor4intp);
        inverseDJfor4intp.clear();
    }

    //inverse matrix!
    for (int i = 0; i < ElementsCalculations.size(); ++i) {  //elements
        for (int j = 0; j < ElementsCalculations[i].size(); ++j) {   //integration points
            var1 = ElementsCalculations[i][j][3];
            var2 = -ElementsCalculations[i][j][1];
            var3 = -ElementsCalculations[i][j][2];
            var4 = ElementsCalculations[i][j][0];
            jacobian.push_back(var1);
            jacobian.push_back(var2);
            jacobian.push_back(var3);
            jacobian.push_back(var4);
            NJ4intp.push_back(jacobian);
            jacobian.clear();
        }
        InverseMatrix.push_back(NJ4intp);
        NJ4intp.clear();
    }

    //new matrix - 1/detJ*inverted
    for (int a=0; a < ElementsCalculations.size(); a++){
        for (int i = 0; i < results_ksi.size(); i++) {
            for (int j = 0; j < results_ksi[i].size(); j++) {
                aux = InverseMatrix[a][i][j];
                aux2 = aux * inverseDJ[a][i];
                NewComponents.push_back(aux2);
            }
            NewJacobians.push_back(NewComponents);
            NewComponents.clear();
        }
        NewMatrix.push_back(NewJacobians);
        NewJacobians.clear();
    }

    //dN/dx
    for (int j = 0; j < NewMatrix.size(); ++j) {      //elems
        for (int k = 0; k < NewMatrix[j].size(); ++k) {       //intp
            for (int i = 0; i < NewMatrix[j][k].size(); ++i) {
                aux = (results_ksi[k][i] * NewMatrix[j][k][0]) + (results_eta[k][i] * NewMatrix[j][k][1]);
                vdndx.push_back(aux);
            }
            vdndxpc.push_back(vdndx);
            vdndx.clear();
        }
        setwdndx.push_back(vdndxpc);
        vdndxpc.clear();
    }

    //dN/dy
    for (int j = 0; j < NewMatrix.size(); ++j) {
        for (int k = 0; k < NewMatrix[j].size(); ++k) {
            for (int i = 0; i < NewMatrix[j][k].size(); ++i) {
                aux = (results_ksi[k][i] * NewMatrix[j][k][2]) + (results_eta[k][i] * NewMatrix[j][k][3]);
                vdndy.push_back(aux);
            }
            vdndypc.push_back(vdndy);
            vdndy.clear();
        }
        setwdndy.push_back(vdndypc);
        vdndypc.clear();
    }
}

void matrix_H(Element vectors, GlobalData gd){
    double aux = 0.0; double aux2 = 0.0;
    double aux3 = 0.0;  double aux4 = 0.0;
    vector<vector<double>> Hintpx;
    vector<vector<double>> Hintpy;
    vector<double> HintpComponents;
    vector<vector<vector<double>>> HintpxWhole;
    vector<vector<vector<double>>> HintpyWhole;
    vector<vector<vector<vector<double>>>> HintpxWholeELEM;
    vector<vector<vector<vector<double>>>> HintpyWholeELEM;

    //calculating Hintp (integral point)
    for (int i = 0; i < setwdndx.size(); ++i) {
        for (int j = 0; j < setwdndx[i].size(); ++j) {
            for (int k = 0; k < setwdndx[i][j].size(); ++k) {
                aux = setwdndx[i][j][k] * setwdndx[i][j][0];
                aux2 = setwdndx[i][j][k] * setwdndx[i][j][1];
                aux3 = setwdndx[i][j][k] * setwdndx[i][j][2];
                aux4 = setwdndx[i][j][k] * setwdndx[i][j][3];
                HintpComponents.push_back(aux);
                HintpComponents.push_back(aux2);
                HintpComponents.push_back(aux3);
                HintpComponents.push_back(aux4);
                Hintpx.push_back(HintpComponents);    //rows of one Hintp matrix from dx
                HintpComponents.clear();
                aux = aux2 = aux3 = aux4 = 0;
            }
            HintpxWhole.push_back(Hintpx);         //whole Hintp matrix, dx
            Hintpx.clear();
        }
        HintpxWholeELEM.push_back(HintpxWhole);
        HintpxWhole.clear();
    }

    for (int i = 0; i < setwdndy.size(); ++i) {
        for (int j = 0; j < setwdndy[i].size(); ++j) {
            for (int k = 0; k < setwdndy[i][j].size(); ++k) {
                aux = setwdndy[i][j][k] * setwdndy[i][j][0];
                aux2 = setwdndy[i][j][k] * setwdndy[i][j][1];
                aux3 = setwdndy[i][j][k] * setwdndy[i][j][2];
                aux4 = setwdndy[i][j][k] * setwdndy[i][j][3];
                HintpComponents.push_back(aux);
                HintpComponents.push_back(aux2);
                HintpComponents.push_back(aux3);
                HintpComponents.push_back(aux4);
                Hintpy.push_back(HintpComponents);    //rows of one Hintp matrix from dy
                HintpComponents.clear();
                aux = aux2 = aux3 = aux4 = 0;
            }
            HintpyWhole.push_back(Hintpy);         //whole Hintp matrix, dy
            Hintpy.clear();
        }
        HintpyWholeELEM.push_back(HintpyWhole);
        HintpyWhole.clear();
    }

    vector<double> elementsSum;
    vector<vector<double>> matricesSum;
    vector<vector<vector<double>>> allSums;
    vector<vector<vector<vector<double>>>> allSumsELEM;

    //sum of Hintpx and Hintpy matrices
    for (int a = 0; a < HintpxWholeELEM.size(); ++a) {
        for (int i = 0; i < HintpxWholeELEM[a].size(); ++i) {
            for (int j = 0; j < HintpxWholeELEM[a][i].size(); ++j) {
                for (int k = 0; k < HintpxWholeELEM[a][i][j].size(); ++k) {
                    aux= HintpxWholeELEM[a][i][j][k] + HintpyWholeELEM[a][i][j][k];
                    elementsSum.push_back(aux);
                    aux = 0;
                }
                matricesSum.push_back(elementsSum);
                elementsSum.clear();
            }
            allSums.push_back(matricesSum);
            matricesSum.clear();
        }
        allSumsELEM.push_back(allSums);
        allSums.clear();
    }

    vector<double> HComponents;
    vector<vector<double>> Hintp;
    vector<vector<vector<double>>> allHintp;
    vector<vector<vector<vector<double>>>> allHintpELEM;

    //H matrix
    for (int a = 0; a < allSumsELEM.size(); ++a) {
        for (int i = 0; i < allSumsELEM[a].size(); ++i) {
            for (int j = 0; j < allSumsELEM[a][i].size(); ++j) {
                for (int k = 0; k < allSumsELEM[a][i][j].size(); ++k) {
                    aux = allSumsELEM[a][i][j][k] * gd.Conductivity * determinantsJ[a][i];
                    HComponents.push_back(aux);
                    aux = 0;
                }
                Hintp.push_back(HComponents);
                HComponents.clear();
            }
            allHintp.push_back(Hintp);
            Hintp.clear();
        }
        allHintpELEM.push_back(allHintp);
        allHintp.clear();
    }

    aux = 0;
    //final matrix H
    vector<double> compHintp;
    for (int j = 0; j < allHintpELEM.size(); ++j) {
        for (int k = 0; k < 4; ++k) {
            for (int i = 0; i < 4; ++i) {
                for (int l = 0; l < allHintpELEM[j].size(); ++l) {
                    aux += allHintpELEM[j][l][k][i] * weights[l % (weights.size())] * weights[l / (weights.size())];
                }
                compHintp.push_back(aux);
                aux = 0;
            }
            compH.push_back(compHintp);
            compHintp.clear();
        }
        elementsMatricesH.push_back(compH);
        compH.clear();
    }
}

void matrix_Hbc_vector_P(Element vectors, GlobalData gd){
    double N11 = 0.0; double N12 = 0.0; double N13 = 0.0; double N14 = 0.0;
    double N21 = 0.0; double N22 = 0.0; double N23 = 0.0; double N24 = 0.0;
    double N31 = 0.0; double N32 = 0.0; double N33 = 0.0; double N34 = 0.0;
    double N41 = 0.0; double N42 = 0.0; double N43 = 0.0; double N44 = 0.0;
    double N1intp1 = 0; double N2intp1 = 0; double N3intp1 = 0; double N4intp1 = 0;
    double N1intp2 = 0; double N2intp2 = 0; double N3intp2 = 0; double N4intp2 = 0;
    double Naux = 0;
    vector<vector<double>> Hbc_side;
    int aux1=0; int aux2=0; double aux = 0;
    double ksi1=0;  double eta1=0;  double ksi2=0;  double eta2=0;  double ksi3=0;  double eta3=0;  double ksi4=0;  double eta4=0;
    double detJ=0;  double L=0;
    vector<double> vecaux1; vector<double> vecaux2;   vector<double> vecaux3; vector<double> vecaux4;
    vector<vector<double>> setvecaux;    //stores auxiliary vectors vecaux1, 2, ...
    vector<double> rows;                 //rows of matrices
    vector<double> test;                 //zeroing side's HBC
    vector<vector<double>> columns1;     //columns of HBC matrices for each side
    vector<vector<double>> columns2;
    vector<vector<double>> columns;
    vector<vector<double>> resultsHbcdetJalfa;    //final side's Hbc matrix
    vector<vector<vector<double>>> setHBCsides;
    vector<vector<double>> elementHBCmatrix;
    vector<double> vectorP;
    vector<double> vectorPmultiplied;
    vector<vector<double>> vectorPsides;
    vector<double> sidePvectorssum;

    vector<vector<vector<double>>> multipliedNs;

    elementHBCmatrix.clear();
    test.clear();

    //checking if there is a BC
    for (int i=0; i < vectors.ID.size(); i++){          //elements
        for(int j=0; j < vectors.ID[i].size(); j++){    //4 IDs in each element
            if (j==3){   //if(last ID)
                aux1 = vectors.ID[i][j];
                aux2 = vectors.ID[i][0];
            } else {     //other IDs
                aux1 = vectors.ID[i][j];
                aux2 = vectors.ID[i][j + 1];
            }
            //checking if there are 2 BCs next to each other
            if (vectors.BC_xy[aux1] == 1 && vectors.BC_xy[aux2] == 1){
                elementHBCmatrix.clear();
                resultsHbcdetJalfa.clear();
                if (weights.size() == 2){          //FOR 2-POINT INT SCHEME
                    //calculating shape functions and integrating
                    if(j==0){
                        ksi1=(-(1.0/sqrt(3.0)));
                        eta1=-1;
                        ksi2=(1.0/sqrt(3.0));
                        eta2=-1;
                    } else if (j==1){
                        ksi1=1;
                        eta1=(-(1.0/sqrt(3.0)));
                        ksi2=1;
                        eta2=(1.0/sqrt(3.0));
                    } else if (j==2){
                        ksi1=(1.0/sqrt(3.0));
                        eta1=1;
                        ksi2=-(1.0/sqrt(3.0));
                        eta2=1;
                    } else {
                        ksi1=-1;
                        eta1=(1.0/sqrt(3.0));
                        ksi2=-1;
                        eta2=(-(1.0/sqrt(3.0)));
                    }
                    //intp1
                    N11=0.25*(1-ksi1)*(1-eta1);
                    N12=0.25*(1+ksi1)*(1-eta1);
                    N13=0.25*(1+ksi1)*(1+eta1);
                    N14=0.25*(1-ksi1)*(1+eta1);
                    vecaux1.push_back(N11); vecaux1.push_back(N12); vecaux1.push_back(N13); vecaux1.push_back(N14);
                    setvecaux.push_back(vecaux1);
                    vecaux1.clear();
                    //intp2
                    N21=0.25*(1-ksi2)*(1-eta2);
                    N22=0.25*(1+ksi2)*(1-eta2);
                    N23=0.25*(1+ksi2)*(1+eta2);
                    N24=0.25*(1-ksi2)*(1+eta2);
                    vecaux2.push_back(N21); vecaux2.push_back(N22); vecaux2.push_back(N23); vecaux2.push_back(N24);
                    setvecaux.push_back(vecaux2);
                    vecaux2.clear();
                } else if (weights.size() == 3) {      //FOR 3-POINT INT SCHEME
                    //calculating shape functions and integrating
                    if(j==0){
                        ksi1=-sqrt(3.0/5.0);
                        eta1=-1;
                        ksi2=0;
                        eta2=-1;
                        ksi3=sqrt(3.0/5.0);
                        eta3=-1;
                    } else if (j==1){
                        ksi1=1;
                        eta1=-sqrt(3.0/5.0);
                        ksi2=1;
                        eta2=0;
                        ksi3=1;
                        eta3=sqrt(3.0/5.0);
                    } else if (j==2){
                        ksi1=sqrt(3.0/5.0);
                        eta1=1;
                        ksi2=0;
                        eta2=1;
                        ksi3=-sqrt(3.0/5.0);
                        eta3=1;
                    } else {
                        ksi1=-1;
                        eta1=sqrt(3.0/5.0);
                        ksi2=-1;
                        eta2=0;
                        ksi3=-1;
                        eta3=-sqrt(3.0/5.0);
                    }
                    //intp1
                    N11=0.25*(1-ksi1)*(1-eta1);
                    N12=0.25*(1+ksi1)*(1-eta1);
                    N13=0.25*(1+ksi1)*(1+eta1);
                    N14=0.25*(1-ksi1)*(1+eta1);
                    vecaux1.push_back(N11); vecaux1.push_back(N12); vecaux1.push_back(N13); vecaux1.push_back(N14);
                    setvecaux.push_back(vecaux1);
                    vecaux1.clear();
                    //intp2
                    N21=0.25*(1-ksi2)*(1-eta2);
                    N22=0.25*(1+ksi2)*(1-eta2);
                    N23=0.25*(1+ksi2)*(1+eta2);
                    N24=0.25*(1-ksi2)*(1+eta2);
                    vecaux2.push_back(N21); vecaux2.push_back(N22); vecaux2.push_back(N23); vecaux2.push_back(N24);
                    setvecaux.push_back(vecaux2);
                    vecaux2.clear();
                    //intp3
                    N31=0.25*(1-ksi3)*(1-eta3);
                    N32=0.25*(1+ksi3)*(1-eta3);
                    N33=0.25*(1+ksi3)*(1+eta3);
                    N34=0.25*(1-ksi3)*(1+eta3);
                    vecaux3.push_back(N31); vecaux3.push_back(N32); vecaux3.push_back(N33); vecaux3.push_back(N34);
                    setvecaux.push_back(vecaux3);
                    vecaux3.clear();
                } else if(weights.size() == 4) {       //FOR 4-POINT INT SCHEME
                    //calculating shape functions and integrating
                    if(j==0){
                        ksi1=-sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                        eta1=-1;
                        ksi2=-sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        eta2=-1;
                        ksi3=sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        eta3=-1;
                        ksi4=sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                        eta4=-1;
                    } else if (j==1){
                        ksi1=1;
                        eta1=-sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                        ksi2=1;
                        eta2=-sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        ksi3=1;
                        eta3=sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        ksi4=1;
                        eta4=sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                    } else if (j==2){
                        ksi1=sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                        eta1=1;
                        ksi2=sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        eta2=1;
                        ksi3=-sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        eta3=1;
                        ksi4=-sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                        eta4=1;
                    } else {
                        ksi1=-1;
                        eta1=sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                        ksi2=-1;
                        eta2=sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        ksi3=-1;
                        eta3=-sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0));
                        ksi4=-1;
                        eta4=-sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0));
                    }
                    //intp1
                    N11=0.25*(1-ksi1)*(1-eta1);
                    N12=0.25*(1+ksi1)*(1-eta1);
                    N13=0.25*(1+ksi1)*(1+eta1);
                    N14=0.25*(1-ksi1)*(1+eta1);
                    vecaux1.push_back(N11); vecaux1.push_back(N12); vecaux1.push_back(N13); vecaux1.push_back(N14);
                    setvecaux.push_back(vecaux1);
                    vecaux1.clear();
                    //intp2
                    N21=0.25*(1-ksi2)*(1-eta2);
                    N22=0.25*(1+ksi2)*(1-eta2);
                    N23=0.25*(1+ksi2)*(1+eta2);
                    N24=0.25*(1-ksi2)*(1+eta2);
                    vecaux2.push_back(N21); vecaux2.push_back(N22); vecaux2.push_back(N23); vecaux2.push_back(N24);
                    setvecaux.push_back(vecaux2);
                    vecaux2.clear();
                    //intp3
                    N31=0.25*(1-ksi3)*(1-eta3);
                    N32=0.25*(1+ksi3)*(1-eta3);
                    N33=0.25*(1+ksi3)*(1+eta3);
                    N34=0.25*(1-ksi3)*(1+eta3);
                    vecaux3.push_back(N31); vecaux3.push_back(N32); vecaux3.push_back(N33); vecaux3.push_back(N34);
                    setvecaux.push_back(vecaux3);
                    vecaux3.clear();
                    //intp4
                    N41=0.25*(1-ksi4)*(1-eta4);
                    N42=0.25*(1+ksi4)*(1-eta4);
                    N43=0.25*(1+ksi4)*(1+eta4);
                    N44=0.25*(1-ksi4)*(1+eta4);
                    vecaux4.push_back(N41); vecaux4.push_back(N42); vecaux4.push_back(N43); vecaux4.push_back(N44);
                    setvecaux.push_back(vecaux4);
                    vecaux4.clear();
                }
                //d=√((x_2-x_1)²+(y_2-y_1)²)
                L=sqrt(pow((vectors.x[aux2 - 1] - vectors.x[aux1 - 1]), 2) + pow((vectors.y[aux2 - 1] - vectors.y[aux1 - 1]), 2));
                detJ=L/2;

                //multiplying in loop to calculate N-vertical X N-horizontal matrix for each intp!
                for (int a = 0; a < setvecaux.size(); a++) {
                    for (int k = 0; k < 4; ++k) {
                        for (int l = 0; l < 4; ++l) {
                            Naux = setvecaux[a][k] * setvecaux[a][l] * weights[a];
                            rows.push_back(Naux);
                            Naux = 0;
                        }
                        columns.push_back(rows);
                        rows.clear();
                    }
                    multipliedNs.push_back(columns);
                    columns.clear();
                }

                for (int a = 0; a < multipliedNs.size(); a++) {
                    for (int k = 0; k < multipliedNs[a].size(); ++k) {
                        for (int l = 0; l < multipliedNs[a][k].size(); ++l) {
                        }
                    }
                }

                //vector P
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < setvecaux.size(); ++l) {
                        aux += (weights[l] * setvecaux[l][k] * gd.Tot * gd.Alfa * detJ);
                    }
                    vectorP.push_back(aux);
                    aux = 0;
                }

                vectorPsides.push_back(vectorP);
                vectorP.clear();
                vectorPmultiplied.clear();
                setvecaux.clear();

                //side matrix sum
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        for (int m = 0; m < multipliedNs.size(); ++m) {
                            Naux += multipliedNs[m][k][l];
                        }
                        rows.push_back(Naux);
                        Naux=0;
                    }
                    Hbc_side.push_back(rows);
                    rows.clear();
                }
                multipliedNs.clear();

                //25*detJ*matrixSum
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        Naux = gd.Alfa * Hbc_side[k][l] * detJ;
                        rows.push_back(Naux);
                        Naux=0;
                    }
                    resultsHbcdetJalfa.push_back(rows);
                    rows.clear();
                }
                Hbc_side.clear();
                setHBCsides.push_back(resultsHbcdetJalfa);
                resultsHbcdetJalfa.clear();

            } else {    //if nodes next each other don't have BCs
                for (int k = 0; k < 4; ++k) {
                    for (int l = 0; l < 4; ++l) {
                        rows.push_back(0.0);
                    }
                    resultsHbcdetJalfa.push_back(rows);
                    rows.clear();
                }
                setHBCsides.push_back(resultsHbcdetJalfa);
                resultsHbcdetJalfa.clear();

                //vector P
                for (int k = 0; k < 4; ++k) {
                    vectorP.push_back(0);
                }
                vectorPsides.push_back(vectorP);
                vectorP.clear();

            }
        }

        //sum of all side's HBC matrices
        for (int j = 0; j < 4; ++j) {
            vector<double> newRow;
            for (int k = 0; k < 4; ++k) {
                double sum = 0.0;
                for (int l = 0; l < 4; ++l) {
                    sum += setHBCsides[l][j][k];
                }
                newRow.push_back(sum);
            }
            elementHBCmatrix.push_back(newRow);
        }
        matricesHBC.push_back(elementHBCmatrix);
        elementHBCmatrix.clear();
        setHBCsides.clear();

        //sum of all side's P vectors
        double sum = 0;
        for (int j = 0; j < vectorPsides.size(); ++j) {
            for (int k = 0; k < vectorPsides[j].size(); ++k) {
                sum += vectorPsides[k][j];
            }
            sidePvectorssum.push_back(sum);
            sum = 0;
        }
        vectorsP.push_back(sidePvectorssum);
        sidePvectorssum.clear();
        vectorPsides.clear();
    }
}

void aggregation(Element vectors, GlobalData gd){

    double aux=0;
    vector<double> rows;
    vector<vector<double>> sumH;     //H+hbc of one element!
    vector<vector<double>>aggregatedHmatrixAUX(gd.NodesNumber, std::vector<double>(gd.NodesNumber, 0.0));
    vector<double> aggregatedPvectorAUX(gd.NodesNumber, 0.0);
    vector<vector<double>>aggregatedCmatrixAUX(gd.NodesNumber, std::vector<double>(gd.NodesNumber, 0.0));

    //H+HBC sum
    for (int i = 0; i < matricesHBC.size(); ++i) {
        for (int j = 0; j < matricesHBC[i].size(); ++j) {
            for (int k = 0; k < matricesHBC[i][j].size(); ++k) {
                aux = matricesHBC[i][j][k] + elementsMatricesH[i][j][k];
                rows.push_back(aux);
            }
            sumH.push_back(rows);
            rows.clear();
        }
        finalH.push_back(sumH);
        sumH.clear();
    }
    rows.clear();

    double var1 = 0; double var2 = 0;

    //H matrix aggregation
    for (int i = 0; i < vectors.ID.size(); ++i) {
            for (int j = 0; j < vectors.ID[i].size(); ++j) {
                for (int l = 0; l < vectors.ID[i].size(); ++l) {
                    var1 = vectors.ID[i][j] - 1;
                    var2 = vectors.ID[i][l] - 1;
                    aggregatedHmatrixAUX[var1][var2] += finalH[i][j][l];
                }
            }
    }

    //C matrix aggregation
    for (int i = 0; i < vectors.ID.size(); ++i) {
        for (int j = 0; j < vectors.ID[i].size(); ++j) {
            for (int l = 0; l < vectors.ID[i].size(); ++l) {
                var1 = vectors.ID[i][j] - 1;
                var2 = vectors.ID[i][l] - 1;
                aggregatedCmatrixAUX[var1][var2] += results_CELEM[i][j][l];
            }
        }
    }

    //copying H to another vector
    for (int i = 0; i < aggregatedHmatrixAUX.size(); ++i) {
        for (int j = 0; j < aggregatedHmatrixAUX[i].size(); ++j) {
                rows.push_back(aggregatedHmatrixAUX[i][j]);
        }
        aggregatedHmatrix.push_back(rows);
        rows.clear();
    }

    //copying C to another vector
    for (int i = 0; i < aggregatedCmatrixAUX.size(); ++i) {
        for (int j = 0; j < aggregatedCmatrixAUX[i].size(); ++j) {
            rows.push_back(aggregatedCmatrixAUX[i][j]);
        }
        aggregatedCmatrix.push_back(rows);
        rows.clear();
    }

    double var = 0;
    //P vector aggregation
    for (int i = 0; i < vectors.ID.size(); ++i) {
        for (int j = 0; j < vectors.ID[i].size(); ++j) {
            var = vectors.ID[i][j] - 1;
            aggregatedPvectorAUX[var] += vectorsP[i][j];
        }
    }
    //copying P to another vector
    for (int i = 0; i < aggregatedPvectorAUX.size(); ++i) {
        aggregatedPvector.push_back(aggregatedPvectorAUX[i]);
    }


}

void matrix_C(vector<double> weight, vector<double> point, GlobalData gd){

    double N1 = 0.0;    double N2 = 0.0;    double N3 = 0.0;    double N4 = 0.0;
    results_N.clear();
    vector<double> intp;
    int counter = 0;

    //shape functions
    for (int i = 0; i < point.size(); i++) {
        for (int j = 0; j < point.size(); j++) {
            N1 = 0.25 * (1 - point[j]) * (1 - point[i]);
            intp.push_back(N1);
            N2 = 0.25 * (1 + point[j]) * (1 - point[i]);
            intp.push_back(N2);
            N3 = 0.25 * (1 + point[j]) * (1 + point[i]);
            intp.push_back(N3);
            N4 = 0.25 * (1 - point[j]) * (1 + point[i]);
            intp.push_back(N4);
            results_N.push_back(intp);
            intp.clear();
        }
    }

    //N-horizontal X N-vertical
    vector<vector<vector<double>>> NintpMatrices;       //set of N matrices for every intp
    vector<vector<double>> Nmatrix;                     //after multiplying N-vertical X N-horizontal
    vector<double> rowsN;
    double aux = 0;

    for (int i = 0; i < results_N.size(); ++i) {
        for (int j = 0; j < results_N[i].size(); ++j) {
            for (int k = 0; k < results_N[i].size(); ++k) {
                aux = results_N[i][j] * results_N[i][k];
                rowsN.push_back(aux);
            }
            Nmatrix.push_back(rowsN);
            rowsN.clear();
        }
        NintpMatrices.push_back(Nmatrix);
        Nmatrix.clear();
    }

    //NUMBER OF ELEMENTS!!!
    vector<vector<vector<vector<double>>>> NintpMatricesELEM;
    for (int i = 0; i < gd.ElementsNumber; ++i) {
        NintpMatricesELEM.push_back(NintpMatrices);
    }

    //SpecificHeat*Density*N
    vector<vector<vector<double>>> sidesCmatrices;
    vector<vector<vector<vector<double>>>> sidesCmatricesELEMS;
    for (int i = 0; i < NintpMatricesELEM.size(); ++i) {                        //elems
        for (int j = 0; j < NintpMatricesELEM[i].size(); ++j) {                 //intps
            for (int k = 0; k < NintpMatricesELEM[i][j].size(); ++k) {          //row in Ns of each intp
                for (int l = 0; l < NintpMatricesELEM[i][j][k].size(); ++l) {   //Ns
                    aux = gd.SpecificHeat * gd.Density * NintpMatricesELEM[i][j][k][l] * determinantsJ[i][j];
                    rowsN.push_back(aux);
                }
                Nmatrix.push_back(rowsN);
                rowsN.clear();
            }
            sidesCmatrices.push_back(Nmatrix);
            Nmatrix.clear();
        }
        sidesCmatricesELEMS.push_back(sidesCmatrices);
        sidesCmatrices.clear();
    }

    aux = 0.0;
    vector<double> newRow;
    for (int a = 0; a < sidesCmatricesELEMS.size(); ++a) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < sidesCmatricesELEMS[a].size(); ++k) {
                    aux += sidesCmatricesELEMS[a][k][i][j] * weights[k % weights.size()] * weights[k / weights.size()];
                }
                newRow.push_back(aux);
                aux = 0;
            }
            results_C.push_back(newRow);
            newRow.clear();
        }
        results_CELEM.push_back(results_C);
        results_C.clear();
    }
}

void non_stationary_solution(GlobalData data){

    int N = data.SimulationTime / data.SimulationStepTime;    //number of iterations/step times

    vector<vector<double>> HplusC;      //matrix after adding H+C/deltatau
    vector<double> rows;
    double aux = 0;     double sum = 0;
    vector<double> Ct0;         //C matrix after multiplying by t0;
    vector<double> CplusP;      //C matrix after adding Ct0/deltatau+P
    vector<vector<double>> matrix;
    vector<vector<double>> HplusCAUX;
    vector<double> t1 (data.NodesNumber, 0);
    double min = 0;
    double max = 0;
    vector<double> minmax;

    //the addition of H+C/deltatau
    for (int i = 0; i < aggregatedHmatrix.size(); ++i) {
        for (int j = 0; j < aggregatedHmatrix[i].size(); ++j) {
            aux = aggregatedHmatrix[i][j] + aggregatedCmatrix[i][j] / data.SimulationStepTime;
            rows.push_back(aux);
        }
        HplusC.push_back(rows);
        rows.clear();
    }

    //vector full of t0;
    double t0 = data.InitialTemp;     //t0 = 100
    vector<double> tempVector;
    for (int i = 0; i < data.NodesNumber; i++){
        tempVector.push_back(t0);
    }

    TemperaturesSet.push_back(tempVector);      //T=0

    min = *min_element(tempVector.begin(), tempVector.end());
    minmax.push_back(min);
    max = *max_element(tempVector.begin(), tempVector.end());
    minmax.push_back(max);
    MinMaxSet.push_back(minmax);      //pair of min and max for T=0 added to the set
    minmax.clear();

    //CALCULATIONS LOOP
    for (int i = 0; i < N ; ++i) {
        //addition (C*t0)/deltatau + P
        for (int j = 0; j < aggregatedCmatrix.size(); ++j) {
            for (int k = 0; k < aggregatedCmatrix[i].size(); ++k) {
                sum += (aggregatedCmatrix[j][k] / data.SimulationStepTime) * tempVector[k];
            }
            Ct0.push_back(sum);
            sum = 0;
        }

        //c*t0 + P vector
        for (int j = 0; j < aggregatedPvector.size(); ++j) {
            aux = Ct0[j] + aggregatedPvector[j];
            CplusP.push_back(aux);
        }

        //copy of H+C matrix
        for (int i = 0; i < HplusC.size(); i++) {
            for (int j = 0; j < HplusC[i].size(); ++j) {
                rows.push_back(HplusC[i][j]);
            }
            matrix.push_back(rows);
            rows.clear();
        }

        //extension of H by P
        for (int i = 0; i < matrix.size(); i++) {
            matrix[i].push_back(CplusP[i]);
        }

        int n = matrix.size();

        for (int i = 0; i < n; i++) {
            //row normalisation
            double pivot = matrix[i][i];
            for (int j = i; j <= n; j++) {
                matrix[i][j] /= pivot;
            }

            //elimination of other rows
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = matrix[k][i];
                    for (int j = i; j <= n; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
        }

        //extraction of the result {t1}
        vector<double> tResult(n);
        for (int i = 0; i < n; i++) {
            tResult[i] = matrix[i][n];
        }

        //copying the result to t1
        for (int i = 0; i < n; i++) {
            t1[i] = tResult[i];
        }

        TemperaturesSet.push_back(t1);      //next step times added to the set

        min = *min_element(t1.begin(), t1.end());
        minmax.push_back(min);
        max = *max_element(t1.begin(), t1.end());
        minmax.push_back(max);
        MinMaxSet.push_back(minmax);      //pair of min and max for T
        minmax.clear();

        //copying t1 results to t0
        for (int j = 0; j < tempVector.size(); ++j) {
            tempVector[j] = t1[j];
        }

        Ct0.clear();
        CplusP.clear();
        matrix.clear();
    }
}