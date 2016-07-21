//
//  Grid.cpp
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#include <assert.h>

#include "Grid.h"
#include "macros.h"


/* Default constructor */
Grid::Grid(){
    
    time = NULL;
    elems = NULL;
}


/* Constructor from max. length */
Grid::Grid(int maxLen){
    
    assert(maxLen>0);
    
    this->maxLen = maxLen;
    this->numElems = 0;
    
    time = new double[maxLen];
    
    elems = NULL;
}


/* Constructor from maximum length and element dimensions */
Grid::Grid(int maxLen,int dx){
    
    assert(maxLen>0);
    assert(dx>0);
    
    this->maxLen = maxLen;
    this->dim = dx;
    this->numElems = 0;
    
    time = new double[maxLen];
    elems = new double[dx*maxLen];
    
}


/* Copy constructor */
Grid::Grid(const Grid& obj){
    
    dim = obj.dim;
    maxLen = obj.maxLen;
    numElems = obj.numElems;
    
    time = new double[maxLen];
    std::copy(obj.time,obj.time+numElems,time);
    
    elems = new double[maxLen*dim];
    std::copy(obj.elems,obj.elems+numElems*dim,elems);
    
}


/* Default destructor */
Grid::~Grid(){
    
    if (time != NULL)
        delete [] time;
    
    if (elems != NULL)
        delete [] elems;
    
}


void Grid::setMaxLen(int maxLen){
    
    int minL;
    
    minL = MIN(this->maxLen,maxLen);
    this->maxLen = maxLen;
    
    double* ctime = new double[maxLen];
    double* celems = new double[maxLen*dim];
    
    // Re-allocate time vector
    if (time != NULL){
        std::copy(time,time+minL,ctime);
        delete [] time;
    }
    
    // Re-allocate elems vector
    if (elems != NULL){
        std::copy(elems,elems+minL*dim,celems);
        delete [] elems;
    }
    
    time = ctime;
    elems = celems;
}