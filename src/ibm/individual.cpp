#include "individual.hpp"

Individual::Individual():
    genotype{g,g},
    epigenotype{w,w},
    modifier_w_z{{0.0,0.0},{0.0,0.0}},
    modifier_z_w{{0.0,0.0},{0.0,0.0}}
{
}

Individual::Individual(Individual const &other):
    genotype{other.genotype[0],other.genotype[1]},
    epigenotype{other.epigenotype[0],other.epigenotype[1]},
    modifier_w_z{
        {other.modifier_w_z[0][0],other.modifier_w_z[0][1]},
        {other.modifier_w_z[1][0],other.modifier_w_z[1][1]}},
    modifier_z_w{
        {other.modifier_z_w[0][0],other.modifier_z_w[0][1]},
        {other.modifier_z_w[1][0],other.modifier_z_w[1][1]}}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    genotype[0] = other.genotype[0];
    genotype[1] = other.genotype[1];

    epigenotype[0] = other.epigenotype[0];
    epigenotype[1] = other.epigenotype[1];

    modifier_w_z[0][0] = other.modifier_w_z[0][0];
    modifier_w_z[1][0] = other.modifier_w_z[1][0];
    modifier_w_z[0][1] = other.modifier_w_z[0][1];
    modifier_w_z[1][1] = other.modifier_w_z[1][1];

    modifier_z_w[0][0] = other.modifier_z_w[0][0];
    modifier_z_w[1][0] = other.modifier_z_w[1][0];
    modifier_z_w[0][1] = other.modifier_z_w[0][1];
    modifier_z_w[1][1] = other.modifier_z_w[1][1];
} // end void Individual::operator=()
