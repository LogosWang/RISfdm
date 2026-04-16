#pragma once
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <cmath>

enum class BCType { Dirichlet, Neumann };

struct Params
{
    // ---- geometry ----
    int dim = 1;
    int nvar = 1;

    int    nx = 10,  ny = 1,  nz = 1;
    double lx = 10.0, ly = 1.0, lz = 1.0;
    double dx = 0.0,  dy = 0.0, dz = 0.0;   // computed by finalize()

    // ---- time ----
    double dt = 100.0;
    double t_end = 1e8;
    int    n_steps = 0;
    int    output_interval = 100;

    // ---- per-variable diffusivities ----
    // Set D[v] for each variable before calling finalize().
    // If left empty, finalize() will fill with D_v (backward compat).
    std::vector<double> D;
    double D_v = 0.05;   // fallback single-variable diffusivity

    // ---- vacancy coupling (RIS) ----
    // iv: index of the vacancy variable in [0, nvar).
    // Set to -1 to disable RIS flux (fall back to independent log-diffusion).
    int iv = -1;

    // ---- source terms ----
    // f_v = src_const[v] + src_linear[v] * c_v   [1/s]
    // Filled with zeros by finalize() if empty.
    std::vector<double> src_const;
    std::vector<double> src_linear;

    // ---- initial concentrations ----
    // c_init[v]:  set before finalize().
    // If empty, finalize() fills with cv_init for all variables.
    std::vector<double> c_init;
    double cv_init = 1.0e-14;   // fallback

    // ---- boundary conditions ----
    struct BC { BCType type; double value; };
    std::vector<BC> BCset;
    // BCset[var * dim*2 + dim_idx*2 + side]
    //   side: 0 = left/bottom/front,  1 = right/top/back
    // Default: Neumann, value = 0  (zero flux)

    void finalize();
    void add_BC(int var, BCType type, double value, int dim_idx, int side);

private:
    void init_BC();
};

// ---- Convenience direction/side constants (include once in main.cpp) ----
static constexpr int Xdim = 0, Ydim = 1, Zdim = 2;
static constexpr int LEFT = 0, RIGHT = 1;
static constexpr int BOTTOM = 0, TOP = 1;
static constexpr int FRONT = 0, BACK = 1;
