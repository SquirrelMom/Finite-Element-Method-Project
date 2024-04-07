#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>

#include "Calculations.h"
#include "Displaying.h"

using namespace std;

int main() {

    GlobalData global_data{};
    Element element;
    ifstream file;
    string FileToOpen;    //grid file

    //----------------------------------- LOADING FILE -----------------------------------

    bool fileOK = false;

    while (!fileOK) {
        cout << "Type the file name: " << endl;
        cin >> FileToOpen;
        file.open(FileToOpen.c_str());
        if (file.good()) {
            fileOK = true;
            string name, name2;
            int value;
            for (int i = 1; i <= 10; i++) {
                file >> name >> name2;

                if (name2 == "number") {
                    file >> value;
                } else {
                    value = stoi(name2);
                }

                if (name == "SimulationTime") {
                    global_data.SimulationTime = value;
                } else if (name == "SimulationStepTime") {
                    global_data.SimulationStepTime = value;
                } else if (name == "Conductivity") {
                    global_data.Conductivity = value;
                } else if (name == "Alfa") {
                    global_data.Alfa = value;
                } else if (name == "Tot") {
                    global_data.Tot = value;
                } else if (name == "InitialTemp") {
                    global_data.InitialTemp = value;
                } else if (name == "Density") {
                    global_data.Density = value;
                } else if (name == "SpecificHeat") {
                    global_data.SpecificHeat = value;
                } else if (name == "Nodes") {
                    global_data.NodesNumber = value;
                } else if (name == "Elements") {
                    global_data.ElementsNumber = value;
                } else {
                    cout << "ERROR" << endl;
                }
            }

            string f;
            file >> f;  //processing not needed line in file

            for (int i = 0; i < global_data.NodesNumber; i++) {
                int number;
                double val1, val2;
                char comma;
                file >> number >> comma >> val1 >> comma >> val2;
                element.x.push_back(val1);
                element.y.push_back(val2);
            }

            file >> f;
            file >> f;

            vector<double> pom;

            for (int i = 0; i < global_data.ElementsNumber; ++i) {
                char comma;
                int id, v1, v2, v3, v4;
                file >> id >> comma >> v1 >> comma >> v2 >> comma >> v3 >> comma >> v4;
                pom.push_back(v1);
                pom.push_back(v2);
                pom.push_back(v3);
                pom.push_back(v4);
                element.ID.push_back(pom);
                pom.clear();
            }

            file >> f;

            int n_bc = BCs(global_data.NodesNumber);    //number of nodes with border condition
            for (int i = 0; i < n_bc; i++) {
                char comma;
                int id;
                file >> id >> comma;
                BC.push_back(id);
            }

            // Sorting BC vector
            sort(BC.begin(), BC.end());

            int counter = 0;
            for (int i = 0; i < global_data.NodesNumber + 1; i++) {
                if (BC[counter] == i) {
                    element.BC_xy.push_back(1);
                    counter++;
                } else {
                    element.BC_xy.push_back(0);
                }
            }
        } else {
            cout << "Invalid name!" << endl;
        }
    }

    //closing previously opened file
    file.close();

    //----------------------------------- INTEGRATION SCHEME SELECTION -----------------------------------

    int scheme;
    bool correctChoice = false;

    while(!correctChoice){
        cout << "Enter the number (2-4) of integration points: ";
        cin >> scheme;
        switch(scheme){
            case 2:
                integration_scheme(weights2, points2);
                correctChoice = true;
                break;
            case 3:
                integration_scheme(weights3, points3);
                correctChoice = true;
                break;
            case 4:
                integration_scheme(weights4, points4);
                correctChoice = true;
                break;
            default:
                cout << "INVALID NUMBER! Choose again:" << endl;
                break;
        }
    }


    //----------------------------------- CALCULATIONS -----------------------------------

    integrate_ksi(weights, points);
    integrate_eta(weights, points);
    jacobian_matrix(element, global_data);
    matrix_H(element, global_data);
    matrix_Hbc_vector_P(element, global_data);
    matrix_C(weights, points, global_data);
    aggregation(element, global_data);
    non_stationary_solution(global_data);

    //----------------------------------- DISPLAYING RESULTS -----------------------------------

    int version;
    int operation;
    correctChoice = false;

    while (!correctChoice) {
        cout << "What would you like to do?" << endl;
        cout << "1 - Show the calculations" << endl;
        cout << "2 - Show the results of calculations" << endl;
        cout << "3 - Generate the results in Paraview file format" << endl;
        cin >> operation;

        switch(operation){
            case 1:
                Display(global_data, element, points);
                correctChoice = true;
                break;
            case 2:
                Results(global_data, element, points);
                correctChoice = true;
                break;
            case 3:
                cout << "Enter which step time: ";
                cin >> version;
                cout << endl;
                Paraview(global_data, element, version);
                correctChoice = true;
                break;
            default:
                cout << "There is no such option! Choose again:" << endl;
                break;
        }
    }


}