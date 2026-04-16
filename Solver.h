#pragma once
#include <vector>
#include "Params.h"
#include "Field.h"

// Explicit (Forward Euler) finite difference solver for multi-component RIS diffusion.
//
// Flux (Manning form):
//   J_i = C_i * D_i * C_V * ( ∇ln(C_V) - ∇ln(C_i) )    for solute i
//   J_V = - Σ_{i≠V} J_i                                  vacancy (flux balance)
//
// When p.iv < 0, falls back to independent log-diffusion:
//   J_i = D_i * ∇ln(C_i)    (∂c/∂t = D∇²(ln c))
//
// Discretised in conservative (divergence) form using face-centred fluxes.
// Source term: f_v = src_const[v] + src_linear[v] * c_v
//
// Stability: D_v * dt / dx^2 <= 0.5  (1D)
//
// BC indexing: BCset[var * dim*2 + dim_idx*2 + side]
//   Neumann value  = prescribed normal flux J·n  at boundary face
//   Dirichlet value = prescribed concentration c  (overwrites after sweep)

class Solver
{
public:
    explicit Solver(const Params& params);

    // Set uniform initial condition for variable v from p.c_init,
    // or call directly with a value.
    void set_initial(int v, double value);
    void set_initial_from_params();

    // Run the full time loop.
    void run();

    // Current field — public for custom IC / post-processing.
    Field u;

private:
    const Params& p;
    Field u_new;
    double t;

    // ---- face-flux helpers ----
    // Compute fluxes at all internal faces and store into J[face * nvar + v].
    // face i is between nodes i and i+1, for i in [0, n-2].
    void compute_face_fluxes_1d(std::vector<double>& J) const;
    void compute_face_fluxes_2d(std::vector<double>& Jx,
                                std::vector<double>& Jy) const;
    void compute_face_fluxes_3d(std::vector<double>& Jx,
                                std::vector<double>& Jy,
                                std::vector<double>& Jz) const;

    // RIS face flux at one face, given left/right concentrations.
    // Returns J[nvar].  If iv < 0, uses simple log-grad flux instead.
    void face_flux(const std::vector<double>& CL,
                   const std::vector<double>& CR,
                   double h,
                   std::vector<double>& J) const;

    void step();
    void step1D();
    void step2D();
    void step3D();

    void apply_dirichlet_1d();
    void apply_dirichlet_2d();
    void apply_dirichlet_3d();

    Params::BC get_bc(int v, int dim_idx, int side) const;
    void output(int step_num) const;
    void check_cfl() const;
};
