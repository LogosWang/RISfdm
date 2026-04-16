#include "Params.h"


void Params::finalize()
{
    if (dim < 1 || dim > 3)
        throw std::runtime_error("dim must be 1, 2, or 3.");
    if (nx < 2)
        throw std::runtime_error("nx must be >= 2.");
    if (dim >= 2 && ny < 2)
        throw std::runtime_error("ny must be >= 2 for 2D.");
    if (dim == 3 && nz < 2)
        throw std::runtime_error("nz must be >= 2 for 3D.");
    if (lx <= 0.0 || (dim >= 2 && ly <= 0.0) || (dim == 3 && lz <= 0.0))
        throw std::runtime_error("domain lengths must be positive.");

    dx = lx / (nx - 1);
    dy = (dim >= 2) ? ly / (ny - 1) : 1.0;
    dz = (dim == 3) ? lz / (nz - 1) : 1.0;

    if (dt <= 0.0 || t_end <= 0.0)
        throw std::runtime_error("dt and t_end must be positive.");
    n_steps = static_cast<int>(t_end / dt);
    if (n_steps < 1)
        throw std::runtime_error("n_steps < 1. Check dt and t_end.");
    if (output_interval < 1)
        throw std::runtime_error("output_interval must be >= 1.");

    // Fill per-variable vectors with defaults if the user left them empty.
    if (D.empty())
        D.assign(nvar, D_v);
    if (D.size() != static_cast<size_t>(nvar))
        throw std::runtime_error("D.size() != nvar");

    if (src_const.empty())  src_const.assign(nvar, 0.0);
    if (src_linear.empty()) src_linear.assign(nvar, 0.0);
    if (src_const.size()  != static_cast<size_t>(nvar) ||
        src_linear.size() != static_cast<size_t>(nvar))
        throw std::runtime_error("src_const/src_linear size != nvar");

    if (c_init.empty()) c_init.assign(nvar, cv_init);
    if (c_init.size() != static_cast<size_t>(nvar))
        throw std::runtime_error("c_init.size() != nvar");

    init_BC();   // sets all BCs to Neumann, value=0
}

void Params::add_BC(int var, BCType type, double value, int dim_idx, int side)
{
    int idx = var * dim * 2 + dim_idx * 2 + side;
    if (idx < 0 || idx >= static_cast<int>(BCset.size()))
        throw std::runtime_error("add_BC: index out of range");
    BCset[idx] = {type, value};
}

void Params::init_BC()
{
    int total = 2 * dim * nvar;
    BCset.assign(total, {BCType::Neumann, 0.0});
}
