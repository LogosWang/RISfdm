#pragma once
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <cmath>
enum class BCType
{
    Dirichlet,
    Neumann
};

struct Params
{
    // dimension
    int dim = 1;
    int nvar = 1;
    // grid size
    int nx = 10;
    int ny = 1;
    int nz = 1;

    // physical domain size
    double lx = 10.0;
    double ly = 1.0;
    double lz = 1.0;

    // grid spacing (computed)
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    // time control
    double dt = 100.0;
    double t_end = 1e8;
    int n_steps = 0;
    int output_interval = 100;

    // diffusion coefficients
    double D_v  = 0.05;
    // double D_cr = 0.05;
    // double D_ni = 0.002;

    // initial values
    double cv_init  = 1.0e-14;
    struct BC
    {
        BCType type;
        double value;
    };
    std::vector<BC> BCset;

    // double ccr_init = 0.18;
    // double cni_init = 0.82;

    // left boundary
    // BCType left_bc_type = BCType::Dirichlet;
    // double left_bc_cv   = 1.0e-14;
    // double left_bc_cr   = 0.18;
    // double left_bc_ni   = 0.10;

    // right boundary
   
    void finalize();
    void add_BC_condition(int var, BCType type, double value, int dimension, int side);
    void init_BC();
};