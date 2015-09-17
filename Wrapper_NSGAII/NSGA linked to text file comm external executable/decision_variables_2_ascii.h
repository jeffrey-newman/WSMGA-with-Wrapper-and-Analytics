//
//  decision_variables_2_ascii.h
//  NSGA linked to text file comm external executable
//
//  Created by a1091793 on 24/07/2015.
//  Copyright (c) 2015 University of Adelaide. All rights reserved.
//

#ifndef __NSGA_linked_to_text_file_comm_external_executable__decision_variables_2_ascii__
#define __NSGA_linked_to_text_file_comm_external_executable__decision_variables_2_ascii__

#include "decision_variables_2_ascii.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


template<typename DecisionVariableType>
class DecisionVariable2Ascii
{
    
public:
    DecisionVariable2Ascii(std::string _filepath) :
    filepath(_filepath)
    {
        
    }
    
    void operator() (std::vector<DecisionVariableType> & _decision_variable)
    {
        std::ofstream dv_file(filepath.c_str(), std::ios_base::out | std::ios_base::trunc);
        for ( typename std::vector<DecisionVariableType>::const_iterator it = _decision_variable.begin();
             it != _decision_variable.end();
             ++it)
        {
            dv_file << *it << '\n';
        }
    }
    
private:
    std::string filepath;
};


bool
isComment (char _ch)
{
    return ((_ch == '/') ? true : false);
}

std::istringstream
getObjFunLine(std::ifstream & _fin)
{
    std::string line;
    std::getline(_fin, line);
    //the following line trims white space from the beginning of the string
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    
    //the following line trims out comments starting with '/'
    line.erase(std::find_if(line.begin(), line.end(), isComment), line.end());
    
    // Erase trailing spaces on the line.
    line.erase(std::find_if(line.rbegin(), line.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), line.end());
    
    if (line.size() == 0) return (getObjFunLine(_fin));
    
    return (std::istringstream(line));
    
}

template <typename ObjectivesTypes, typename ConstraintsTypes>
struct ObjectivesAndConstraints
{
    std::vector<ObjectivesTypes> objectives;
    std::vector<ConstraintsTypes> constraints;
    std::vector<ConstraintsTypes> constraint_thresholds;
};


template<typename ObjectivesType, typename ConstraintsType>
class Ascii2ObjectiveValues
{
public:
    Ascii2ObjectiveValues(std::string _filepath, bool _is_obj, bool _is_cnstr) :
    filepath(_filepath), is_obj(_is_obj), is_cnstr(_is_cnstr)
    {
        
    }
    
    std::shared_ptr<ObjectivesAndConstraints< ObjectivesType, ConstraintsType> > operator() ()
    {
        std::shared_ptr<ObjectivesAndConstraints< ObjectivesType, ConstraintsType> > objectives_and_constraints(new ObjectivesAndConstraints<ObjectivesType, ConstraintsType>);
        
        std::ifstream obj_file(filepath.c_str(), std::ios_base::in);
        
        if (is_obj)
        {
            int num_objectives;
            //        std::istringstream ss = getObjFunLine(obj_file);
            //        ss >> num_objectives;
            getObjFunLine(obj_file) >> num_objectives;
            
            ObjectivesType objective;
            for (int i = 0; i < num_objectives; ++i)
            {
                getObjFunLine(obj_file) >> objective;
                objectives_and_constraints->objectives.push_back(objective);
            }
        }
        
        if (is_cnstr)
        {
            int num_constraints;
            getObjFunLine(obj_file) >> num_constraints;
            
            ConstraintsType constraint;
            for (int i = 0; i < num_constraints; ++i)
            {
                getObjFunLine(obj_file) >> constraint;
                objectives_and_constraints->constraints.push_back(constraint);
            }
            
            for (int i = 0; i < num_constraints; ++i)
            {
                getObjFunLine(obj_file) >> constraint;
                objectives_and_constraints->constraint_thresholds.push_back(constraint);
            }
        }
        
        return (objectives_and_constraints);
    }
    
private:
    std::string filepath;
    bool is_obj;
    bool is_cnstr;
};





#endif /* defined(__NSGA_linked_to_text_file_comm_external_executable__decision_variables_2_ascii__) */
