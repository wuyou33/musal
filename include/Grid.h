//
//  Grid.h
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_Grid_h
#define MusAL_Grid_h

#include <iostream>
#include <algorithm>


class Grid {

public:
    
    /* Default constructor */
    Grid();
    
    /* Constructor from maximum length */
    Grid(int);
    
    /* Constructor from element dimension and maximum length */
    Grid(int,int);
    
    /* Copy constructor */
    Grid(const Grid&);
    
    /* Assignment operator */
    Grid & operator= (const Grid& obj){
        
        if (this != &obj){
            
            maxLen = obj.maxLen;
            numElems = obj.numElems;
            dim = obj.dim;
            
            double* c_time = new double[maxLen];
            std::copy(obj.time,obj.time+maxLen,c_time);
            
            double* c_elems = new double[maxLen*dim];
            std::copy(obj.elems,obj.elems+maxLen*dim,c_elems);
            
            if (time != NULL)
                delete [] time;
            
            if (elems != NULL)
                delete [] elems;
            
            time = c_time;
            elems = c_elems;
        }
        
        return *this;
    }
    
    /* Destructor */
    ~Grid();
    
    /* Sets */
    void setMaxLen(int);
    
    /* Gets */
    inline double* getTime(){
        return time;
    }
    
    inline double* getElems(){
        return elems;
    }
    
    inline int getNumElems(){
        return numElems;
    }
    
    inline int getMaxLen(){
        return maxLen;
    }
    
    inline int getDim(){
        return dim;
    }
    
private:
    
    /* Grid maximum length */
    int maxLen;
    
    /* Number of elements in grid */
    int numElems;
    
    /* Grid elements dimension */
    int dim;
    
    /* Vector of time values */
    double* time;
    
    /* Vector of elements */
    double* elems;
};


#endif
