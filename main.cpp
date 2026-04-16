#include <iostream>
#include "Params.h"
#include "Field.h"
#include "Solver.h"

// Variable indices  (Cr=0, Ni=1, V=2, Fe=3)
static constexpr int iCr = 0, iNi = 1, iV = 2, iFe = 3;

int main()
{
    // ------------------------------------------------------------------
    // 1-D RIS diffusion  [length: nm,  time: s]
    //
    // Domain: x in [0, 10 nm]
    //   Left  (x = 0):  symmetry — zero flux
    //   Right (x = L):  grain boundary sink — Dirichlet on C_V
    //
    // Grid: nx=51, dx = 10/(51-1) = 0.2 nm
    //
    // CFL stability:  D_eff_vac * dt / dx^2 <= 0.5
    //   D_eff_vac = 0.18*0.05 + 0.10*0.002 + 0.72*0.01 = 0.0164 nm^2/s
    //   dx^2 = 0.04 nm^2
    //   dt_max = 0.5 * 0.04 / 0.0164 = 1.22 s
    //   -> use dt = 20 s  ... wait, re-check with nx=51:
    //   dx = 10/50 = 0.2 nm,  dx^2 = 0.04 nm^2
    //   dt_max = 0.5 * 0.04 / 0.0164 = 1.22 s
    //   -> use dt = 1 s  (CFL = 0.41, safe)
    //
    // Source term stability: K_r * dt = 0.01 * 1 = 0.01 << 1  (safe)
    // ------------------------------------------------------------------

    Params p;
    p.dim  = 1;
    p.nvar = 4;

    p.nx = 51;
    p.lx = 10.0;    // nm

    p.dt    = 1.0;              // s
    p.t_end = 5.0e6;            // s  (~16 years)
    p.output_interval = 500000; // every 5e5 s (~6 days)

    // ---- diffusivities ----
    // D[iV] sets the CFL check for the vacancy equation via sum(C_i * D[i]).
    // It is NOT used in the Manning flux itself (J_V = -sum J_solutes).
    p.D.resize(4);
    p.D[iCr] = 0.05;
    p.D[iNi] = 0.002;
    p.D[iV]  = 0.0;    // vacancy does not appear in Manning flux as an explicit D
    p.D[iFe] = 0.01;
    // D_eff_vac = 0.18*0.05 + 0.10*0.002 + 0.72*0.01 = 0.009+0.0002+0.0072 = 0.0164

    // ---- vacancy variable index (enables RIS Manning flux) ----
    p.iv = iV;

    // ---- source terms: f_v = src_const[v] + src_linear[v]*c_v ----
    p.src_const.assign(4, 0.0);
    p.src_linear.assign(4, 0.0);
    p.src_const[iV]  =  1.0e-8;   // K_0: production rate
    p.src_linear[iV] = -1.0e-2;   // -K_r: recombination/sink  => C_V^eq = K_0/K_r = 1e-6

    // ---- initial concentrations ----
    p.c_init.resize(4);
    p.c_init[iCr] = 0.18;
    p.c_init[iNi] = 0.10;
    p.c_init[iV]  = 1.0e-6;       // near equilibrium: K0/Kr = 1e-8/1e-2
    p.c_init[iFe] = 1.0 - 0.18 - 0.10 - 1.0e-6;   // ~0.72

    p.finalize();

    // ---- boundary conditions ----
    // Right boundary: grain boundary sink — vacancy only.
    // C_V = thermal-equilibrium value (perfect vacancy sink).
    // Solutes keep Neumann zero-flux (GB is NOT a solute sink).
    // C_V^eq at GB ≈ thermal equilibrium at irradiation temperature (~500 °C for Fe-alloys).
    // Using 1e-8 rather than 1e-14: the ratio C_V_bulk/C_V_GB = 100, which is physically
    // reasonable and avoids the 1e8 span that makes the Manning flux drain all solutes from
    // the GB node → blocks vacancy transport → C_V piles up before sink.
    p.add_BC(iV, BCType::Dirichlet, 1.0e-8, Xdim, RIGHT);

    // ---- solver ----
    Solver solver(p);
    solver.set_initial_from_params();
    solver.u(iV, p.nx - 1) = 1.0e-14;   // enforce at t=0

    solver.run();
    return 0;
}
