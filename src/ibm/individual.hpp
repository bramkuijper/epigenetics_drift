#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <vector>

// definition of individual properties

enum Allele {
    g = 0, // inactive
    G = 1 // active
};

enum EpiAllele {
    w = 0, // inactive
    z = 1 // active
};

class Individual
{
    public:

        // contains the two alleles governing adaptation
        Allele genotype[2];

        // contains w or z for each allele
        EpiAllele epigenotype[2];

        // modifier from inactive to active 
        // 2 x 2 because 2 environments (A, B) 
        // and 2 alleles
        double modifier_w_z[2][2];

        // modifier from active to inactive 
        double modifier_z_w[2][2];
        
        // default constructor
        Individual();

        // copy constructor
        Individual(Individual const &other);

        // assignment operator
        void operator=(Individual const &other);
};

#endif
