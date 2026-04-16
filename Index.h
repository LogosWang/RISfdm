#pragma once
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <cmath>

inline int idx1d(int i, int nx){
    if (i<0 || i>=nx) {
        std::cerr<<"index out of bounds: "<<i<<std::endl;
    }
    return i;
}
inline int idx2d(int i, int j, int nx, int ny){
    if (i<0 || i>=nx || j<0 || j>=ny) {
        std::cerr<<"index out of bounds: "<<i<<", "<<j<<std::endl;
    }
    return i + j*nx;
}
inline int idx3d(int i, int j, int k, int nx, int ny, int nz){
    if (i<0 || i>=nx || j<0 || j>=ny || k<0 || k>=nz) {
        std::cerr<<"index out of bounds: "<<i<<", "<<j<<", "<<k<<std::endl;
    }
    return i + j*nx + k*nx*ny;
}