#include "Displaying.h"
#include "Calculations.h"

#include <iostream>
using namespace std;

void Display(GlobalData global_data, Element element, vector<double> points){

    cout << endl;
    cout << "---------------- DATA LOADED FROM THE FILE ----------------" << endl;
    cout << endl;
    cout << "SimulationTime: " << global_data.SimulationTime << endl;
    cout << "SimulationStepTime: " << global_data.SimulationStepTime << endl;
    cout << "Conductivity: " << global_data.Conductivity << endl;
    cout << "Alfa: " << global_data.Alfa << endl;
    cout << "Tot: " << global_data.Tot << endl;
    cout << "InitialTemp: " << global_data.InitialTemp << endl;
    cout << "Density: " << global_data.Density << endl;
    cout << "SpecificHeat: " << global_data.SpecificHeat << endl;
    cout << "Nodes Number: " << global_data.NodesNumber << endl;
    cout << "Elements Number: " << global_data.ElementsNumber << endl;

    cout << endl;
    cout << "Node:" << endl;
    for (int i = 0; i < global_data.NodesNumber; ++i) {
        cout << i+1 << ") " << "x=" << element.x[i] << ", y=" << element.y[i] << ", ";
        cout << "BC=" << element.BC_xy[i] << endl;
    }
    cout <<endl;
    cout << "Element, type=DC2D4:" << endl;
    for (int i = 0; i < element.ID.size(); ++i) {
        cout << "ID" << i+1 << " [ ";
        for (int j = 0; j < element.ID[i].size(); ++j) {
            cout << element.ID[i][j] << " ";
        }
        cout << "] " << endl;
    }
    cout <<endl;
    cout << "BC:" << endl;
    for (int i = 0; i < BC.size(); i++){
        cout << BC[i] << " ";
    }
    cout << endl;
    cout << endl;

    cout << "---------------- ORDER OF THE INTEGRATION POINTS ----------------" << endl;
    cout << endl;
    int count = 1;
    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            cout << "Point" << count << endl;
            cout << "ksi= " << points[j] << ", eta= " << points[i] << endl;
            count++;
        }
    }
    cout << endl;

    cout << "---------------- INTEGRATING THE SHAPE FUNCTIONS OVER KSI AND ETA ----------------" << endl;
    cout << endl;
    cout << "        N1/dksi  N2/dksi  N3/dksi  N4/dksi" << endl;
    for (int i = 0; i < results_ksi.size(); ++i) {
        cout << "intp" << i + 1 << ": ";
        for (int j = 0; j < results_ksi[i].size(); ++j) {
            cout << results_ksi[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    cout << "        N1/deta  N2/deta  N3/deta  N4/deta" << endl;
    for (int i = 0; i < results_eta.size(); ++i) {
        cout << "intp" << i + 1 << ": ";
        for (int j = 0; j < results_eta[i].size(); ++j) {
            cout << results_eta[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "---------------- COMPUTING THE JACOBIAN MATRIX, MATRICES H, HBC, AND VECTOR P ----------------" << endl;
    cout << endl;
    for (int a=0; a < ElementsCalculations.size(); a++){
        cout << ">>>> ELEMENT " << a + 1 << " <<<<" << endl;
        for (int i = 0; i < ElementsCalculations[a].size(); ++i) {
            cout << "Jacobian Matrix for the integration point " << i + 1 << ":" << endl;
            for (int j = 0; j < ElementsCalculations[a][i].size(); ++j) {
                cout << ElementsCalculations[a][i][j] << " ";
                if ((j + 1) % 2 == 0) {
                    cout << endl;
                }
            }
            cout << "The determinant of the Jacobian matrix for this point: det[J]=" << determinantsJ[a][i] << endl;
            cout << "The inverse of the Jacobian matrix for this point: 1/det[J]=" << inverseDJ[a][i] << endl;
            cout << endl;
            cout << "Values of dN/dx: " << endl;
            for (int j = 0; j < setwdndx[a][i].size(); ++j) {
                cout << setwdndx[a][i][j] << " ";
            }
            cout << endl;
            cout << "Values of dN/dy: " << endl;
            for (int j = 0; j < setwdndy[a][i].size(); ++j) {
                cout << setwdndy[a][i][j] << " ";
            }
            cout << endl; cout << endl;

        }
        cout << "H for the element " << a+1 << endl;
        for (int j = 0; j < elementsMatricesH[a].size(); ++j) {
            for (int k = 0; k < elementsMatricesH[a][j].size(); ++k) {
                cout << elementsMatricesH[a][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
        cout << "Hbc for the element " << a+1 << endl;
        for (int i = 0; i < matricesHBC[a].size(); ++i) {
            for (int j = 0; j < matricesHBC[i].size(); ++j) {
                cout << matricesHBC[a][i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        cout << "Vector P for the element " << a+1 << endl;
        for (int i = 0; i < vectorsP[a].size(); ++i) {
            cout << vectorsP[a][i] << " ";
        }
        cout << endl;
        cout << endl;
    }

    cout << "---------------- C MATRIX ----------------" << endl;
    cout << endl;

    cout << "C matrices for the elements" << endl;
    cout << endl;
    for (int i = 0; i < results_CELEM.size(); ++i) {
        cout << "ELEMENT " << i+1 << endl;
        for (int j = 0; j < results_CELEM[i].size(); ++j) {
            for (int k = 0; k < results_CELEM[i][j].size(); ++k) {
                cout << results_CELEM[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;

    cout << "---------------- AGGREGATION ----------------" << endl;
    cout << endl;
    cout << "Final H matrices" << endl;
    cout << endl;
    for (int i = 0; i < matricesHBC.size(); ++i) {
        cout << "ELEMENT " << i+1 << endl;
        for (int j = 0; j < matricesHBC[i].size(); ++j) {
            for (int k = 0; k < matricesHBC[i][j].size(); ++k) {
                cout << finalH[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    cout << "After the aggregation of the H matrix:" << endl;
    for (int i = 0; i < aggregatedHmatrix.size(); ++i) {
        for (int j = 0; j < aggregatedHmatrix[i].size(); ++j) {
            cout << aggregatedHmatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "After the aggregation of the P vector:" << endl;
    for (int i = 0; i < aggregatedPvector.size(); ++i) {
        cout << aggregatedPvector[i] << " ";
    }
    cout << endl;
    cout << endl;

    cout << "After the aggregation of the C matrix:" << endl;
    for (int i = 0; i < aggregatedCmatrix.size(); ++i) {
        for (int j = 0; j < aggregatedCmatrix[i].size(); ++j) {
            cout << aggregatedCmatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "---------------- SIMULATION RESULTS ----------------" << endl;
    cout << endl;

    int counter = 0;
    for (int i = 0; i < TemperaturesSet.size(); ++i) {
        cout << "T = " << counter << endl;
        for (int j = 0; j < TemperaturesSet[i].size(); ++j) {
            cout << TemperaturesSet[i][j] << " ";
        }
        cout << endl;
        cout << "MIN=" << MinMaxSet[i][0] << ", MAX=" << MinMaxSet[i][1] << endl;
        counter += global_data.SimulationStepTime;
        cout << endl;
    }

}

void Paraview(GlobalData global_data, Element element, int StepTime){

    cout << "# vtk DataFile Version 2.0" << endl;
    cout << "Unstructured Grid Example" << endl;
    cout << "ASCII" << endl;
    cout << "DATASET UNSTRUCTURED_GRID" << endl;
    cout << endl;

    cout << "POINTS " << global_data.NodesNumber << " float" << endl;
    for (int i = 0; i < global_data.NodesNumber; ++i) {
        cout << element.x[i] << " " << element.y[i] << " " << 0 << endl;
    }
    cout << endl;

    cout << "CELLS " << global_data.ElementsNumber << " " << (global_data.ElementsNumber*5) << endl;
    for (int i = 0; i < element.ID.size(); ++i) {
        cout << 4 << " ";
        for (int j = 0; j < element.ID[i].size(); ++j) {
            cout << element.ID[i][j]-1 << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "CELL_TYPES 9" << endl;
    for (int i = 0; i < 9; ++i) {
        cout << 9 << endl;
    }
    cout << endl;

    cout << "POINT_DATA " << global_data.NodesNumber << endl;
    cout << "SCALARS Temp float 1" << endl;
    cout << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < (global_data.NodesNumber); ++i) {
        cout << TemperaturesSet[StepTime][i] << endl;
    }
    cout << endl;

}

void Results(GlobalData global_data, Element element, vector<double> points){
    cout << "---------------- SIMULATION RESULTS ----------------" << endl;
    cout << endl;

    int counter = 0;
    for (int i = 0; i < TemperaturesSet.size(); ++i) {
        cout << "T = " << counter << endl;
        cout << "MIN=" << MinMaxSet[i][0] << ", MAX=" << MinMaxSet[i][1] << endl;
        counter += global_data.SimulationStepTime;
        cout << endl;
    }
}