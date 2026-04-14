#include "Params.h"
const int V = 0;
const int Cr = 1;
const int Ni = 2;
const int Fe = 3;
const int Si = 4;
const int I = 5;
const int Xdim = 0;
const int Ydim = 1;
const int Zdim = 2;
const int LEFT = 0;
const int BOTTOM = 0;
const int FRONT = 0;
const int RIGHT = 1;
const int TOP = 1;
const int BACK = 1;


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

    if (dim >= 2)
        dy = ly / (ny - 1);
    else
        dy = 1.0;

    if (dim == 3)
        dz = lz / (nz - 1);
    else
        dz = 1.0;

    if (dt <= 0.0 || t_end <= 0.0)
        throw std::runtime_error("dt and t_end must be positive.");

    n_steps = static_cast<int>(t_end / dt);

    if (n_steps < 1)
        throw std::runtime_error("n_steps < 1. Check dt and t_end.");

    if (output_interval < 1)
        throw std::runtime_error("output_interval must be >= 1.");
    init_BC();
    add_BC_condition(V, BCType::Dirichlet,1e-14, Xdim, RIGHT);
}
void Params::add_BC_condition(int var, BCType type, double value, int dimension, int side) {
int total_num = 2*dim*nvar;
int BC_index = var*dim*2+dimension*2+side;
if (BC_index>=total_num) {
    std::cerr<<"no such BC"<<std::endl;
} 
BCset[BC_index].type = type;
BCset[BC_index].value = value;
}

void Params::init_BC() {
    int total_num = 2*dim*nvar;
    BCset.clear();
    BCset.resize(total_num);
    for (int i=0;i<total_num;++i){
        BC default_BC;
        default_BC.type = BCType::Neumann;
        default_BC.value = 0.0;
        BCset.push_back(default_BC);
    }
     
}