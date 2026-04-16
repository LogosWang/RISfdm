#include "Solver.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <cmath>

// ---------------------------------------------------------------------------
// Construction / IC
// ---------------------------------------------------------------------------

Solver::Solver(const Params& params)
    : u(params), p(params), u_new(params), t(0.0)
{}

void Solver::set_initial(int v, double value)
{
    if (v < 0 || v >= p.nvar)
        throw std::runtime_error("var index out of range");
    u.fill(v, value);
}

void Solver::set_initial_from_params()
{
    for (int v = 0; v < p.nvar; ++v)
        u.fill(v, p.c_init[v]);
}

// ---------------------------------------------------------------------------
// BC helper
// ---------------------------------------------------------------------------

Params::BC Solver::get_bc(int v, int dim_idx, int side) const
{
    int idx = v * p.dim * 2 + dim_idx * 2 + side;
    if (idx < 0 || idx >= static_cast<int>(p.BCset.size()))
        throw std::runtime_error("BC index out of range");
    return p.BCset[idx];
}

// ---------------------------------------------------------------------------
// Core: face flux  (one face, given left/right nodal concentrations)
//
// RIS (Manning) flux when p.iv >= 0:
//   J_i = C_i_face * D_i * C_V_face * (ln C_V_R - ln C_V_L - ln C_i_R + ln C_i_L) / h
//   J_V = - sum_{i != V} J_i
//
// Independent log-diffusion fallback (p.iv < 0):
//   J_i = D_i * C_i_face * (ln C_i_R - ln C_i_L) / h
//
// h = grid spacing for this direction.
// ---------------------------------------------------------------------------

void Solver::face_flux(const std::vector<double>& CL,
                       const std::vector<double>& CR,
                       double h,
                       std::vector<double>& J) const
{
    const int iv = p.iv;

    if (iv >= 0) {
        // Manning flux — algebraically equivalent non-log form:
        //
        //   J_i = D_i * ( C_i * ∂C_V/∂x  -  C_V * ∂C_i/∂x )
        //
        // Derivation: C_i*C_V*(∇ln C_V - ∇ln C_i)
        //           = C_i * ∂C_V/∂x  -  C_V * ∂C_i/∂x      (exact identity)
        //
        // This form is numerically stable even when C_V spans many orders of
        // magnitude (e.g. bulk 1e-6 to GB sink 1e-14): the gradient is
        // (C_V_R - C_V_L)/h, bounded by the initial conditions, whereas the
        // log form uses C_V_face * log(C_V_R/C_V_L)/h which is ~9× larger
        // for a 1e-6-to-1e-14 step — enough to push C_V negative in one step.

        const double CV_face = 0.5 * (CL[iv] + CR[iv]);
        const double dCV     = (CR[iv] - CL[iv]) / h;   // ∂C_V/∂x at face

        double J_vac = 0.0;
        for (int v = 0; v < p.nvar; ++v) {
            if (v == iv) continue;
            const double Cv_face = 0.5 * (CL[v] + CR[v]);
            const double dCv     = (CR[v] - CL[v]) / h;
            const double Jv      = p.D[v] * (Cv_face * dCV - CV_face * dCv);
            J[v]    = Jv;
            J_vac  -= Jv;
        }
        J[iv] = J_vac;
    } else {
        // Independent Fickian diffusion (fallback when iv < 0)
        for (int v = 0; v < p.nvar; ++v)
            J[v] = p.D[v] * (CR[v] - CL[v]) / h;
    }
}

// ---------------------------------------------------------------------------
// 1D: compute face fluxes for all (nx-1) interior faces.
// J[face * nvar + v]  —  face i sits between nodes i and i+1.
// ---------------------------------------------------------------------------

void Solver::compute_face_fluxes_1d(std::vector<double>& J) const
{
    const int nx   = p.nx;
    const int nvar = p.nvar;
    J.assign((nx - 1) * nvar, 0.0);

    std::vector<double> CL(nvar), CR(nvar), Jf(nvar);
    for (int face = 0; face < nx - 1; ++face) {
        for (int v = 0; v < nvar; ++v) {
            CL[v] = u(v, face);
            CR[v] = u(v, face + 1);
        }
        face_flux(CL, CR, p.dx, Jf);
        for (int v = 0; v < nvar; ++v)
            J[face * nvar + v] = Jf[v];
    }
}

// ---------------------------------------------------------------------------
// 1D step
//
// Conservative update for node i:
//   c_new[v][i] = c[v][i] - dt/dx * (J_{i+1/2}[v] - J_{i-1/2}[v]) + dt*f[v][i]
//
// Boundary faces:
//   Neumann: prescribed flux at the boundary face (default 0).
//   Dirichlet: boundary face flux = 0 (node overwritten by apply_dirichlet_1d).
//
// For Neumann, the boundary flux is the normal flux J·n̂, so:
//   left  face: J_left  = -bc.value  (n̂ points into domain, i.e. +x; flux bc is outward)
//   right face: J_right =  bc.value
//
// For a zero-flux Neumann (bc.value = 0), both boundary fluxes are zero.
// ---------------------------------------------------------------------------

void Solver::step1D()
{
    const int nx   = p.nx;
    const int nvar = p.nvar;

    // Interior face fluxes
    std::vector<double> J;
    compute_face_fluxes_1d(J);

    // Boundary face fluxes: J_bnd[v]
    std::vector<double> J_left(nvar, 0.0), J_right(nvar, 0.0);
    for (int v = 0; v < nvar; ++v) {
        const Params::BC bcl = get_bc(v, 0, 0);
        const Params::BC bcr = get_bc(v, 0, 1);
        if (bcl.type == BCType::Neumann) J_left[v]  = -bcl.value; // inward = −outward
        if (bcr.type == BCType::Neumann) J_right[v] =  bcr.value;
        // Dirichlet: leave at 0; node will be overwritten below.
    }

    // Node updates
    for (int i = 0; i < nx; ++i) {
        const double* JL = (i == 0)    ? J_left.data()  : &J[(i-1) * nvar];
        const double* JR = (i == nx-1) ? J_right.data() : &J[i * nvar];

        // Lattice-site conservation correction (Manning mode only).
        // The vacancy source S_V = K0 - Kr*C_V violates Σ S_α = 0 unless we
        // distribute the compensating flux among solutes proportionally to their
        // concentrations: S_α = -S_V * C_α / C_solutes  (α ≠ vacancy).
        double src_vac_correction = 0.0;
        double C_solutes = 0.0;
        if (p.iv >= 0) {
            src_vac_correction = p.src_const[p.iv] + p.src_linear[p.iv] * u(p.iv, i);
            for (int v = 0; v < nvar; ++v)
                if (v != p.iv) C_solutes += u(v, i);
            if (C_solutes < 1.0e-30) C_solutes = 1.0e-30;
        }

        for (int v = 0; v < nvar; ++v) {
            double div_J = (JR[v] - JL[v]) / p.dx;
            double src;
            if (p.iv >= 0 && v != p.iv) {
                // Solute: intrinsic src (usually 0) plus lattice-conservation correction
                src = p.src_const[v] + p.src_linear[v] * u(v, i)
                      - src_vac_correction * u(v, i) / C_solutes;
            } else {
                src = p.src_const[v] + p.src_linear[v] * u(v, i);
            }
            u_new(v, i)  = u(v, i) - p.dt * div_J + p.dt * src;
        }
    }

    apply_dirichlet_1d();
    std::swap(u.data, u_new.data);
}

void Solver::apply_dirichlet_1d()
{
    for (int v = 0; v < p.nvar; ++v) {
        const Params::BC bcl = get_bc(v, 0, 0);
        const Params::BC bcr = get_bc(v, 0, 1);
        if (bcl.type == BCType::Dirichlet) u_new(v, 0)        = bcl.value;
        if (bcr.type == BCType::Dirichlet) u_new(v, p.nx - 1) = bcr.value;
    }
}

// ---------------------------------------------------------------------------
// 2D: face fluxes in x and y.
// Jx[face_x * nvar + v]: face_x = j*( nx-1 ) + i,  between (i,j) and (i+1,j)
// Jy[face_y * nvar + v]: face_y = j*(  nx  ) + i,  between (i,j) and (i,j+1)
// ---------------------------------------------------------------------------

void Solver::compute_face_fluxes_2d(std::vector<double>& Jx,
                                    std::vector<double>& Jy) const
{
    const int nx = p.nx, ny = p.ny, nvar = p.nvar;
    Jx.assign((nx-1) * ny      * nvar, 0.0);
    Jy.assign(nx     * (ny-1)  * nvar, 0.0);

    std::vector<double> CL(nvar), CR(nvar), Jf(nvar);

    // x-faces
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx-1; ++i) {
            for (int v = 0; v < nvar; ++v) { CL[v]=u(v,i,j); CR[v]=u(v,i+1,j); }
            face_flux(CL, CR, p.dx, Jf);
            int f = j*(nx-1) + i;
            for (int v = 0; v < nvar; ++v) Jx[f*nvar+v] = Jf[v];
        }

    // y-faces
    for (int j = 0; j < ny-1; ++j)
        for (int i = 0; i < nx; ++i) {
            for (int v = 0; v < nvar; ++v) { CL[v]=u(v,i,j); CR[v]=u(v,i,j+1); }
            face_flux(CL, CR, p.dy, Jf);
            int f = j*nx + i;
            for (int v = 0; v < nvar; ++v) Jy[f*nvar+v] = Jf[v];
        }
}

void Solver::step2D()
{
    const int nx = p.nx, ny = p.ny, nvar = p.nvar;
    std::vector<double> Jx, Jy;
    compute_face_fluxes_2d(Jx, Jy);

    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {

            // x-direction fluxes (left/right of node (i,j))
            auto jx_flux = [&](int ii, int jj, int side) -> double* {
                static std::vector<double> zero(64, 0.0);
                // side 0 = left face (between ii-1 and ii), side 1 = right face
                if (side == 0) {
                    if (ii == 0) return zero.data();  // boundary face (Neumann=0)
                    return &Jx[(jj*(nx-1) + (ii-1)) * nvar];
                } else {
                    if (ii == nx-1) return zero.data();
                    return &Jx[(jj*(nx-1) + ii) * nvar];
                }
            };
            auto jy_flux = [&](int ii, int jj, int side) -> double* {
                static std::vector<double> zero(64, 0.0);
                if (side == 0) {
                    if (jj == 0) return zero.data();
                    return &Jy[((jj-1)*nx + ii) * nvar];
                } else {
                    if (jj == ny-1) return zero.data();
                    return &Jy[(jj*nx + ii) * nvar];
                }
            };

            double* JxL = jx_flux(i, j, 0);
            double* JxR = jx_flux(i, j, 1);
            double* JyB = jy_flux(i, j, 0);
            double* JyT = jy_flux(i, j, 1);

            double src_vac_correction2d = 0.0, C_sol2d = 1.0e-30;
            if (p.iv >= 0) {
                src_vac_correction2d = p.src_const[p.iv] + p.src_linear[p.iv]*u(p.iv,i,j);
                C_sol2d = 0.0;
                for (int v=0;v<nvar;++v) if(v!=p.iv) C_sol2d+=u(v,i,j);
                if (C_sol2d < 1.0e-30) C_sol2d = 1.0e-30;
            }
            for (int v = 0; v < nvar; ++v) {
                double div_J = (JxR[v]-JxL[v])/p.dx + (JyT[v]-JyB[v])/p.dy;
                double src;
                if (p.iv >= 0 && v != p.iv)
                    src = p.src_const[v]+p.src_linear[v]*u(v,i,j) - src_vac_correction2d*u(v,i,j)/C_sol2d;
                else
                    src = p.src_const[v]+p.src_linear[v]*u(v,i,j);
                u_new(v,i,j) = u(v,i,j) - p.dt*div_J + p.dt*src;
            }
        }

    apply_dirichlet_2d();
    std::swap(u.data, u_new.data);
}

void Solver::apply_dirichlet_2d()
{
    const int nx = p.nx, ny = p.ny;
    for (int v = 0; v < p.nvar; ++v) {
        const Params::BC bcxl=get_bc(v,0,0), bcxr=get_bc(v,0,1);
        const Params::BC bcyb=get_bc(v,1,0), bcyt=get_bc(v,1,1);
        if (bcxl.type==BCType::Dirichlet) for(int j=0;j<ny;++j) u_new(v,0,   j)=bcxl.value;
        if (bcxr.type==BCType::Dirichlet) for(int j=0;j<ny;++j) u_new(v,nx-1,j)=bcxr.value;
        if (bcyb.type==BCType::Dirichlet) for(int i=0;i<nx;++i) u_new(v,i,   0)=bcyb.value;
        if (bcyt.type==BCType::Dirichlet) for(int i=0;i<nx;++i) u_new(v,i,ny-1)=bcyt.value;
    }
}

// ---------------------------------------------------------------------------
// 3D step
// ---------------------------------------------------------------------------

void Solver::compute_face_fluxes_3d(std::vector<double>& Jx,
                                    std::vector<double>& Jy,
                                    std::vector<double>& Jz) const
{
    const int nx=p.nx, ny=p.ny, nz=p.nz, nvar=p.nvar;
    Jx.assign((nx-1)*ny*nz*nvar, 0.0);
    Jy.assign(nx*(ny-1)*nz*nvar, 0.0);
    Jz.assign(nx*ny*(nz-1)*nvar, 0.0);

    std::vector<double> CL(nvar), CR(nvar), Jf(nvar);
    for (int k=0;k<nz;++k) for (int j=0;j<ny;++j) for (int i=0;i<nx-1;++i) {
        for (int v=0;v<nvar;++v){CL[v]=u(v,i,j,k);CR[v]=u(v,i+1,j,k);}
        face_flux(CL,CR,p.dx,Jf);
        int f=(k*ny+j)*(nx-1)+i;
        for(int v=0;v<nvar;++v) Jx[f*nvar+v]=Jf[v];
    }
    for (int k=0;k<nz;++k) for (int j=0;j<ny-1;++j) for (int i=0;i<nx;++i) {
        for (int v=0;v<nvar;++v){CL[v]=u(v,i,j,k);CR[v]=u(v,i,j+1,k);}
        face_flux(CL,CR,p.dy,Jf);
        int f=(k*(ny-1)+j)*nx+i;
        for(int v=0;v<nvar;++v) Jy[f*nvar+v]=Jf[v];
    }
    for (int k=0;k<nz-1;++k) for (int j=0;j<ny;++j) for (int i=0;i<nx;++i) {
        for (int v=0;v<nvar;++v){CL[v]=u(v,i,j,k);CR[v]=u(v,i,j,k+1);}
        face_flux(CL,CR,p.dz,Jf);
        int f=(k*ny+j)*nx+i;
        for(int v=0;v<nvar;++v) Jz[f*nvar+v]=Jf[v];
    }
}

void Solver::step3D()
{
    const int nx=p.nx,ny=p.ny,nz=p.nz,nvar=p.nvar;
    std::vector<double> Jx,Jy,Jz;
    compute_face_fluxes_3d(Jx,Jy,Jz);

    std::vector<double> zero(nvar, 0.0);
    for (int k=0;k<nz;++k)
    for (int j=0;j<ny;++j)
    for (int i=0;i<nx;++i) {
        auto gx=[&](int ii,int side)->const double*{
            if(side==0){if(ii==0)return zero.data();return &Jx[((k*ny+j)*(nx-1)+(ii-1))*nvar];}
            else       {if(ii==nx-1)return zero.data();return &Jx[((k*ny+j)*(nx-1)+ii)*nvar];}
        };
        auto gy=[&](int jj,int side)->const double*{
            if(side==0){if(jj==0)return zero.data();return &Jy[((k*(ny-1)+(jj-1))*nx+i)*nvar];}
            else       {if(jj==ny-1)return zero.data();return &Jy[((k*(ny-1)+jj)*nx+i)*nvar];}
        };
        auto gz=[&](int kk,int side)->const double*{
            if(side==0){if(kk==0)return zero.data();return &Jz[(((kk-1)*ny+j)*nx+i)*nvar];}
            else       {if(kk==nz-1)return zero.data();return &Jz[((kk*ny+j)*nx+i)*nvar];}
        };
        const double *JxL=gx(i,0),*JxR=gx(i,1);
        const double *JyB=gy(j,0),*JyT=gy(j,1);
        const double *JzF=gz(k,0),*JzK=gz(k,1);
        double svc3d=0.0, Cs3d=1.0e-30;
        if(p.iv>=0){
            svc3d=p.src_const[p.iv]+p.src_linear[p.iv]*u(p.iv,i,j,k);
            Cs3d=0.0; for(int v=0;v<nvar;++v) if(v!=p.iv) Cs3d+=u(v,i,j,k);
            if(Cs3d<1.0e-30) Cs3d=1.0e-30;
        }
        for(int v=0;v<nvar;++v){
            double dJ=(JxR[v]-JxL[v])/p.dx+(JyT[v]-JyB[v])/p.dy+(JzK[v]-JzF[v])/p.dz;
            double src;
            if(p.iv>=0&&v!=p.iv)
                src=p.src_const[v]+p.src_linear[v]*u(v,i,j,k)-svc3d*u(v,i,j,k)/Cs3d;
            else
                src=p.src_const[v]+p.src_linear[v]*u(v,i,j,k);
            u_new(v,i,j,k)=u(v,i,j,k)-p.dt*dJ+p.dt*src;
        }
    }
    apply_dirichlet_3d();
    std::swap(u.data,u_new.data);
}

void Solver::apply_dirichlet_3d()
{
    const int nx=p.nx,ny=p.ny,nz=p.nz;
    for(int v=0;v<p.nvar;++v){
        auto bcxl=get_bc(v,0,0),bcxr=get_bc(v,0,1);
        auto bcyb=get_bc(v,1,0),bcyt=get_bc(v,1,1);
        auto bczf=get_bc(v,2,0),bczk=get_bc(v,2,1);
        for(int k=0;k<nz;++k)for(int j=0;j<ny;++j){
            if(bcxl.type==BCType::Dirichlet) u_new(v,0,   j,k)=bcxl.value;
            if(bcxr.type==BCType::Dirichlet) u_new(v,nx-1,j,k)=bcxr.value;
        }
        for(int k=0;k<nz;++k)for(int i=0;i<nx;++i){
            if(bcyb.type==BCType::Dirichlet) u_new(v,i,0,   k)=bcyb.value;
            if(bcyt.type==BCType::Dirichlet) u_new(v,i,ny-1,k)=bcyt.value;
        }
        for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){
            if(bczf.type==BCType::Dirichlet) u_new(v,i,j,0)   =bczf.value;
            if(bczk.type==BCType::Dirichlet) u_new(v,i,j,nz-1)=bczk.value;
        }
    }
}

// ---------------------------------------------------------------------------
// Time advance + CFL
// ---------------------------------------------------------------------------

void Solver::step()
{
    if      (p.dim == 1) step1D();
    else if (p.dim == 2) step2D();
    else                 step3D();
    t += p.dt;
}

void Solver::check_cfl() const
{
    double sum_inv = 1.0/(p.dx*p.dx);
    if (p.dim >= 2) sum_inv += 1.0/(p.dy*p.dy);
    if (p.dim == 3) sum_inv += 1.0/(p.dz*p.dz);

    if (p.iv >= 0) {
        // Manning flux — effective diffusivities:
        //   solute v:   D_eff = D[v] * C_V_init
        //   vacancy iv: D_eff = sum_{i!=iv} C_i_init * D[i]
        double D_vac_eff = 0.0;
        for (int v = 0; v < p.nvar; ++v) {
            if (v == p.iv) continue;
            D_vac_eff += p.c_init[v] * p.D[v];
            double D_eff = p.D[v] * p.c_init[p.iv];
            double cfl   = D_eff * p.dt * sum_inv;
            if (cfl > 0.5)
                std::cerr << "CFL warning var " << v
                          << ": D_eff=" << D_eff << " CFL=" << cfl << "\n";
        }
        double cfl_vac = D_vac_eff * p.dt * sum_inv;
        if (cfl_vac > 0.5)
            std::cerr << "CFL warning vacancy (var " << p.iv << ")"
                      << ": D_eff=" << D_vac_eff << " CFL=" << cfl_vac << "\n";
        std::cout << "  CFL check (Manning): D_eff_vac = " << D_vac_eff
                  << "  CFL_vac = " << cfl_vac << "\n";
    } else {
        for (int v = 0; v < p.nvar; ++v) {
            double cfl = p.D[v] * p.dt * sum_inv;
            if (cfl > 0.5)
                std::cerr << "CFL warning var " << v
                          << ": D=" << p.D[v] << " CFL=" << cfl << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

void Solver::output(int step_num) const
{
    std::string filename = "output_" + std::to_string(step_num) + ".csv";
    std::ofstream ofs(filename);
    if (!ofs) { std::cerr << "Cannot open " << filename << "\n"; return; }

    ofs << std::scientific << std::setprecision(8);

    // Header row
    if (p.dim >= 1) ofs << "x";
    if (p.dim >= 2) ofs << ",y";
    if (p.dim == 3) ofs << ",z";
    for (int v = 0; v < p.nvar; ++v) ofs << ",c" << v;
    ofs << "\n";

    if (p.dim == 1) {
        for (int i = 0; i < p.nx; ++i) {
            ofs << i * p.dx;
            for (int v = 0; v < p.nvar; ++v) ofs << "," << u(v, i);
            ofs << "\n";
        }
    } else if (p.dim == 2) {
        for (int j = 0; j < p.ny; ++j)
            for (int i = 0; i < p.nx; ++i) {
                ofs << i*p.dx << "," << j*p.dy;
                for (int v = 0; v < p.nvar; ++v) ofs << "," << u(v, i, j);
                ofs << "\n";
            }
    } else {
        for (int k = 0; k < p.nz; ++k)
            for (int j = 0; j < p.ny; ++j)
                for (int i = 0; i < p.nx; ++i) {
                    ofs << i*p.dx << "," << j*p.dy << "," << k*p.dz;
                    for (int v = 0; v < p.nvar; ++v) ofs << "," << u(v, i, j, k);
                    ofs << "\n";
                }
    }
}

// ---------------------------------------------------------------------------
// Run
// ---------------------------------------------------------------------------

void Solver::run()
{
    std::cout << "Explicit FD  |  RIS Manning flux\n"
              << "  dim=" << p.dim << "  nvar=" << p.nvar
              << "  iv=" << p.iv
              << "  nx=" << p.nx;
    if (p.dim >= 2) std::cout << "  ny=" << p.ny;
    if (p.dim == 3) std::cout << "  nz=" << p.nz;
    std::cout << "  dx=" << p.dx << "\n"
              << "  dt=" << p.dt << "  t_end=" << p.t_end
              << "  n_steps=" << p.n_steps << "\n";

    check_cfl();
    output(0);

    for (int n = 1; n <= p.n_steps; ++n) {
        step();
        if (n % p.output_interval == 0) {
            std::cout << "  step " << std::setw(6) << n
                      << "  t=" << std::scientific << std::setprecision(4) << t << "\n";
            output(n);
        }
    }
    std::cout << "Done.  t=" << t << "\n";
}
