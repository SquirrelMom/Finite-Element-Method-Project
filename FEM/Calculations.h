#ifndef MES_CALCULATIONS_H
#define MES_CALCULATIONS_H

#include<vector>
using namespace std;

struct GlobalData {
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int NodesNumber;
    int ElementsNumber;
};

struct Element{
    vector<vector<double>> ID;
    vector<double> x;
    vector<double> y;
    vector<int> BC_xy;
};

extern vector<double> weights;
extern vector<double> points;

//Single-point integration scheme
extern vector<double> weights1;
extern vector<double> points1;

//2-point integration scheme
extern vector<double> weights2;
extern vector<double> points2;

//3-point integration scheme
extern vector<double> weights3;
extern vector<double> points3;

//4-point integration scheme
extern vector<double> weights4;
extern vector<double> points4;

//vectors storing results
extern vector<int> BC;                                          //vector that stores the IDs of nodes that have border condition
extern vector<vector<double>> results_ksi;                      //results of dN/dksi integration
extern vector<vector<double>> results_eta;                      //results of dN/dksi integration
extern vector<vector<double>> jacobian_matrices;                //the jacobian matrices
extern vector<vector<vector<double>>> ElementsCalculations;     //the jacobian_matrices of elements (jacobian matrix is 1 vector, 1 element is 4/9/16 jacobian_matrices for integration point)
extern vector<vector<double>> determinantsJ;                    //determinants of the jacobian matrices
extern vector<vector<double>> inverseDJ;                        //inverses of the jacobian matrices determinants
extern vector<vector<double>> vdndxpc;                          //values of dN/dx for the integration points
extern vector<vector<double>> vdndypc;                          //values of dN/dy for the integration points
extern vector<vector<vector<double>>> setwdndx;                 //values of dN/dx for the elements
extern vector<vector<vector<double>>> setwdndy;                 //values of dN/dy for the elements
extern vector<vector<double>> compH;                            //H matrices
extern vector<vector<vector<double>>> elementsMatricesH;        //set of H matrices
extern vector<vector<vector<double>>> matricesHBC;              //Hbc matrices of the elements
extern vector<vector<double>> vectorsP;                         //P vectors of the elements
extern vector<vector<vector<double>>> finalH;                   //H+Hbc for every element
extern vector<vector<double>> aggregatedHmatrix;                //aggregated H+Hbc matrix
extern vector<double> aggregatedPvector;                        //aggregated P vector
extern vector<vector<double>> results_N;                        //results of the shape function calculations
extern vector<vector<double>> results_C;                        //C matrices of the elements
extern vector<vector<vector<double>>> results_CELEM;            //set of C matrices
extern vector<vector<double>> aggregatedCmatrix;                //aggregated C matrix
extern vector<vector<double>> TemperaturesSet;                  //set of t1 vectors for every step time
extern vector<vector<double>> MinMaxSet;                        //min and max for a particular step time

//the function that assigns the chosen integration scheme
void integration_scheme(vector<double> chosen_weights, vector<double> chosen_points);

//the function that calculates the number of points with border condition
int BCs(int N);

//the function that integrates dN/dksi for every integration point
void integrate_ksi(vector<double> weight, vector<double> point);

//the function that integrates dN/deta for every integration point
void integrate_eta(vector<double> weight, vector<double> point);

//the function that calculates the jacobian matrix
void jacobian_matrix(Element vectors, GlobalData gd);

//the function calculates the H matrix for every element
void matrix_H (Element vectors, GlobalData gd);

//the function that calculates the Hbc matrix and P vector for every element
void matrix_Hbc_vector_P(Element vectors, GlobalData gd);

//the function that aggregates H+Hbc matrix, C matrix and P vector
void aggregation(Element vectors, GlobalData gd);

//the function that calculates the C matrix for every element
void matrix_C(vector<double> weight, vector<double> point, GlobalData gd);

//function that calculates the non-stationary solution of the problem
void non_stationary_solution(GlobalData data);

#endif //MES_CALCULATIONS_H