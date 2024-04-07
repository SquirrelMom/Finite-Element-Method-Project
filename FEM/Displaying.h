#ifndef MES_DISPLAYING_H
#define MES_DISPLAYING_H

#include "Calculations.h"

//displaying the results
void Display(GlobalData global_data, Element element, vector<double> points);

//displaying the results in Paraview file standard
void Paraview(GlobalData global_data, Element element, int StepTime);

//displaying just the time simulation
void Results(GlobalData global_data, Element element, vector<double> points);

#endif //MES_DISPLAYING_H
