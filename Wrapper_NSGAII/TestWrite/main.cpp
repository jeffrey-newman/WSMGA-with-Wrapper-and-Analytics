//
//  main.cpp
//  TestWrite
//
//  Created by a1091793 on 24/07/2015.
//  Copyright (c) 2015 University of Adelaide. All rights reserved.
//

#include <iostream>
#include "decision_variables_2_ascii.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    DecisionVariable2Ascii<int> myWriter("dv.txt");
    std::vector<int> decision_variables { 5, 0, 5, 0, 0, 0, 1, 0, 5, 1 };
    myWriter(decision_variables);
    return EXIT_SUCCESS;
}
