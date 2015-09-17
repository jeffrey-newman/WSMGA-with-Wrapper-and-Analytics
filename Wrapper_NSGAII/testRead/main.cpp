//
//  main.cpp
//  testRead
//
//  Created by a1091793 on 24/07/2015.
//  Copyright (c) 2015 University of Adelaide. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "decision_variables_2_ascii.h"

int main(int argc, const char * argv[]) {
    
    std::ofstream fs("testObj.txt");
    fs << "//Objective function values and constraints violations\n";
    fs << "2		/number of objectives (Jeff, feel free to skip this line if no_obj is input elsewhere)\n";
    fs << "1.111 	/objective 1 value\n";
    fs << "2.222 	/objective 2 value\n";
    fs << "7		/number of constraints (Jeff, feel free to skip this line if no_cons is input elsewhere)\n";
    fs << "1.0		/cons - value returned for constrained MUSIC outputs\n";
    fs << "1.0          \n";
    fs << "1.0\n";
    fs << "1.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0		/ overallcons - threshold for constrained MUSIC outputs - constraint met if = 0.0, constraint violated if > 0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs << "0.0\n";
    fs.close();
    
    Ascii2ObjectiveValues<double, double> myReader("testObj.txt", true, true);
    std::shared_ptr<ObjectivesAndConstraints<double, double> > obj_and_cnstr = myReader();
    
    if (obj_and_cnstr->objectives[0] == 1.111
        && obj_and_cnstr->objectives[1] == 2.222
        && obj_and_cnstr->constraints[0] == 1.0
        && obj_and_cnstr->constraints[1] == 1.0
        && obj_and_cnstr->constraints[2] == 1.0
        && obj_and_cnstr->constraints[3] == 1.0
        && obj_and_cnstr->constraints[4] == 0.0
        && obj_and_cnstr->constraints[5] == 0.0
        && obj_and_cnstr->constraints[6] == 0.0)
    {
        return EXIT_SUCCESS;
    }
    
    return EXIT_FAILURE;
}
