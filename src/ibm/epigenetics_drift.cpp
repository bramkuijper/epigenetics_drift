// Extending Leimar & McNamara's cue integration model
// to include cultural transmission
// Bram Kuijper 
// 2019
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>


// various functions, such as unique filename creation
#include "auxiliary.hpp"

// the individual class, which defines properties of each
// individual used here
#include "individual.hpp"

#define DEBUG

// standard namespace
using namespace std;

// C++ random number generation
int seed = get_nanoseconds();
mt19937 rng_r{static_cast<long unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);
// sample alleles from diploids with probability 0.5
bernoulli_distribution allele_sample(0.5);

// parameters & variables:
// a lot of the parameters declared here will be overridden
// in the init_arguments function which reads parameters from the command line

// number of individuals in population
const int NPatches = 400;

// make uniform distribution of patches so that we can 
// sample randomly amongst them
uniform_int_distribution<> random_patch(0, NPatches - 1);


// number of breeders of patch
int NBreeder = 100;

// number of generations
//int number_generations = 75000;
int number_generations = 100;

// environmental switch rate
//
// parameters below will be changed on the command line

// frequency of high state patches
double envt_switch[2] = {0.0,0.0};

// cost of being inactive in environment A
double sA = 0.5;
// cost of being active in environment B
double sB = 0.1;

// dominance coefficient
double h = 0.5;

// dispersal probability
double d = 0.1;

// genetic mutation rate of the genetic loci
double mu_g = 0.01;

// modifier mutation rate
double mu_mu = 0.01;
double sdmu_mu = 0.01;

//write out the data 
//every nth generation
int data_nth_generation = 10;


// the patch with its breeders
struct Patch
{
    // number of hermaphroditic breeder
    vector<Individual> breeders;

    // next generation breeders
    vector<Individual> breeders_t1;
    
    int n_breeders; // number of breeders currently in patch
    bool envt_state_is_A;  // high-state patch yes/no
};

Patch Pop[NPatches];

// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{


}

// write down all parameters to the file DataFile
void write_parameters(ofstream &DataFile)
{
    DataFile << endl << endl
        "mu_mu;" << mu_mu << endl
        "sdmu_mu;" << sdmu_mu << endl
        "mu_g;" << mu_g << endl
        "sigma_BA;" << envt_switch[0] << endl
        "sigma_AB;" << envt_switch[1] << endl
        "sA;" << sA << endl
        "sB;" << sB << endl
        "d;" << d << endl;
        "h;" << h << endl;
} // void write_parameters(ofstream &DataFile)


// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile << "generation;";
    
    for (int allele_1 = 0; allele_1 < 2; ++allele_1)
    {
        for (int allele_2 = 0; allele_2 < 2; ++allele_2)
        {
            for (int epi_allele_1 = 0; epi_allele_1 < 2; ++epi_allele_1)
            {
                for (int epi_allele_2 = 0; epi_allele_2 < 2; ++epi_allele_2)
                {
                    DataFile << "freq_" << (allele_1 == 0 ? "g" : "G") <<  (epi_allele_1 == 0 ? "w" : "z") << (allele_2 == 0 ? "g" : "G") << (epi_allele_2 == 0 ? "w" : "z") << ";";
                }
            }
        }
    }

    DataFile << "freq_w;freq_g;";

    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        DataFile 
            << "mean_mu_w_z_" << (envt_i == 0 ? "B" : "A") << ";"
            << "var_mu_w_z_" << (envt_i == 0 ? "B" : "A") << ";"
            << "mean_mu_z_w_" << (envt_i == 0 ? "B" : "A") << ";"
            << "var_mu_z_w_" << (envt_i == 0 ? "B" : "A") << ";" << endl;
    }

    DataFile << "freq_A;mean_surviving_breeders;" << endl;
} // void write_data_headers(ofstream &DataFile)



// write data both for winter and summer populations
void write_stats(ofstream &DataFile, int generation)
{
    // store frequencies of combinations of alleles and epialles
    int epigenotype_freq[2][2][2][2];

    // quadro-for-loop to initialize frequencies to 0
    for (int allele_1 = 0; allele_1 < 2; ++allele_1)
    {
        for (int allele_2 = 0; allele_2 < 2; ++allele_2)
        {
            for (int epi_allele_1 = 0; epi_allele_1 < 2; ++epi_allele_1)
            {
                for (int epi_allele_2 = 0; epi_allele_2 < 2; ++epi_allele_2)
                {
                    epigenotype_freq[allele_1][allele_2][epi_allele_1][epi_allele_2] = 0;
                }
            }
        }
    }

    double mean_modifier_w_z[2] = { 0.0, 0.0 };
    double mean_modifier_z_w[2] = { 0.0, 0.0 };
    double ss_modifier_w_z[2] = { 0.0, 0.0 };
    double ss_modifier_z_w[2] = { 0.0, 0.0 };

    double freq_A;

    int total_breeders = 0;

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        nbreeders = Pop[patch_i].breeders.size();

        total_breeders += nbreeders;

        freq_A += Pop[patch_i].envt_state_is_A;

        for (int breeder_i = 0; breeder_i < nbreeders; ++breeder_i)
        {
            ++epigenotype_freq[
                Pop[patch_i].breeders[breeder_i].genotype[0]
                ][
                Pop[patch_i].breeders[breeder_i].genotype[1]
                ][
                Pop[patch_i].breeders[breeder_i].epigenotype[0]
                ][
                Pop[patch_i].breeders[breeder_i].epigenotype[1]
                ];

            for (int envt_i = 0; envt_i < 2; ++envt_i)
            {
                val =  0.5 * (
                            Pop[patch_i].breeders[breeder_i].modifier_w_z[envt_i][0]
                            +
                            Pop[patch_i].breeders[breeder_i].modifier_w_z[envt_i][1]
                          );

                mean_modifier_w_z[envt_i] += val;
                ss_modifier_w_z[envt_i] += val * val;

                val =  0.5 * (
                            Pop[patch_i].breeders[breeder_i].modifier_z_w[envt_i][0]
                            +
                            Pop[patch_i].breeders[breeder_i].modifier_z_w[envt_i][1]
                          );

                mean_modifier_z_w[envt_i] += val;
                ss_modifier_z_w[envt_i] += val * val;
            }
        } //end for int breeder_i
    } // end for int patch_i


    DataFile << generation << ";";

    double freq_g = 0.0;
    double freq_w = 0.0;

    for (int allele_1 = 0; allele_1 < 2; ++allele_1)
    {
        for (int allele_2 = 0; allele_2 < 2; ++allele_2)
        {
            for (int epi_allele_1 = 0; epi_allele_1 < 2; ++epi_allele_1)
            {
                for (int epi_allele_2 = 0; epi_allele_2 < 2; ++epi_allele_2)
                {
                    mean_freq = (double) epigenotype_freq[allele_1][allele_2][epi_allele_1][epi_allele_2] / total_breeders;

                    DataFile << mean_freq << ";";

                    if (allele_1 == g)
                    {
                        freq_g += 0.5 * mean_freq; 
                    }

                    if (allele_2 == g)
                    {
                        freq_g += 0.5 * mean_freq; 
                    }
                    
                    if (epi_allele_1 == w)
                    {
                        freq_w += 0.5 * mean_freq; 
                    }

                    if (epi_allele_2 == w)
                    {
                        freq_w += 0.5 * mean_freq; 
                    }
                }
            }
        }
    }

    DataFile << freq_w << ";" << freq_g << ";";

    double mean, var;

    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        mean = mean_modifier_w_z[envt_i] / total_breeders;
        var =  ss_modifier_w_z[envt_i] / total_breeders - mean * mean;

        DataFile << mean << ";" << var << ";";
        
        mean = mean_modifier_z_w[envt_i] / total_breeders;
        var =  ss_modifier_z_w[envt_i] / total_breeders - mean * mean;

        DataFile << mean << ";" << var << ";";
    }

    DataFile << freq_A / NPatches << ";" << (double) total_breeders / NPatches << ";" << endl;
} // void write_stats(ofstream &DataFile, int generation, int timestep)

// initialize the population 
// at the start of the simulation
void init_population()
{
    // loop through all individuals 
    // and assign them values for the cue loci
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // patch in a high state or not
        Pop[patch_i].envt_state_is_A = uniform(rng_r) < 
            envt_switch[1] / (envt_switch[0] + envt_switch[1]);
      
        // initialize all the breeders within patch patch_i
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            // produce an individual 
            // initialize it with certain values
            // then add that to stack of breeders 
            Individual ind_init;

            // initialize every allele for each locus
            for (int allele_i = 0; allele_i < 2; ++allele_i)
            {
                // initialize the genotype for adaptation (g,G)
                ind_init.genotype[allele_i] = uniform(rng_r) < 0.5 ? g : G;

                // initialize the epigenotype that modifies the adaptation locus
                //
                // let's initialize this all as active alleles
                ind_init.epigenotype[allele_i] = z;

                // initialize the modifier loci that affect the epigenetic switch rate
                for (int envt_i = 0; envt_i < 2; ++envt_i)
                {
                    ind_init.modifier_w_z[envt_i][allele_i] = 0.0;
                    
                    ind_init.modifier_z_w[envt_i][allele_i] = 0.0;
                } // for (int envt_i =
            } // end for (int allele_i

            Pop[patch_i].breeders.append(ind_init);
        } // end for (int breeder_i

        // we should have initilized NBreeder individuals in this patch
        // assert that
        assert(Pop[patch_i].breeders.size() == NBreeder);

    } // end for (int patch_i
} // end void init_population()

// mutation of a certain allele with value val
// given mutation rate mu and mutational distribution stdev sdmu
double mutation(double val, double mu, double sdmu)
{
    if (uniform(rng_r) < mu)
    {
        normal_distribution<> mutational_effect(0.0, sdmu);
        val += mutational_effect(rng_r);
    }

    return(val);
}

// auxiliary function to bound values between [min,max]
void clamp(double &val, double min, double max)
{
    val = val > max ? max : val < min ? min : val;
}

// the survival probability of an individual
// with n active alleles either in envt A or B
double fitness_payoff(
        int const n_active_alleles
        , bool const envt_A)
{
    double val = 1.0;

    // G favored in A, g favored in B

    switch(n_active_alleles)
    {
        case 0: // 2 inactive alleles which fit envt B
            val = envt_A ? 1.0 - sA : 1.0;
            break;
        case 1: // 1 inactive allele
            val = envt_A ? 1.0 - h * sA : 1.0 - h * sB;
            break;
        case 2: // 2 active alleles, which fit envt A
            val = envt_A ? 1.0; 1.0 - sB;
            break;
        default:
            cout << "something goes quite horrifically wrong with the number of active alleles and the survival function. Abort." << endl;
            exit(1);
    }

    assert(val >= 0.0);
    assert(val <= 1.0);

    return(val);
} // end double survival_probability

// create a new offspring
void create_offspring(Individual &mother
        ,Individual &father
        ,Individual &offspring
        ,bool envt_is_A)
{
    // first inherit the modifiers
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {

    }


    Allele parental_allele;
    EpiAllele parental_epiallele;

    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        // inherit G alleles from mom or dad,
        // sample from random parental chromosome
        int parental_chromosome = allele_sample(rng_r);
        
        parental_allele = allele_i == 0 ? 
            mother.genotype[parental_chromosome]
            :
            father.genotype[parental_chromosome];

        // mutate parental allele
        if (uniform(rng_r) < mu_g)
        {
            parental_allele = !parental_allele;
        }

        offspring.genotype[allele_i] = parental_allele;

        // now the epialleles
        parental_epiallele = allele_i == 0 ? 
            mother.epigenotype[parental_chromosome]
            :
            father.epigenotype[parental_chromosome];

        if (parental_epiallele == w)
        {
            if (uniform(rng_r) < 0.5 * (
                        offspring.modifier_w_z[envt_is_A][0]
                        +
                        offspring.modifier_w_z[envt_is_A][1]))
            {
                parental_epiallele = !parental_epiallele;
            }
        } else
        {
            if (uniform(rng_r) < 0.5 * (
                        offspring.modifier_z_w[envt_is_A][0]
                        +
                        offspring.modifier_z_w[envt_is_A][1]))
            {
                parental_epiallele = !parental_epiallele;
            }
        }

    } // end for allele_i
} // end create_offspring()

void adult_survival()
{
    // 1. we need to loop through all the patches 
    // to get at the individual breeders
    //
    // 2. we need to look at the individual breeder's 
    // phenotype == genotype + epigenetic state
    //
    // 3. if phenotype is shit, kill individual
    //
    bool envt_is_A;
    
    
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        envt_is_A = Pop[patch_i].envt_state_is_A;

        for (int breeder_i = 0; breeder_i < Pop[patch_i].breeders[breeder_i].size(); ++breeder_i)
        {
            // count the number of alleles who are active
            int n_active_alleles = 0;

            for (int allele_i = 0; allele_i < 2; ++allele_i)
            {
                // look at the individual breeder's phenotype
                if (Pop[patch_i].breeders[breeder_i].genotype[allele_i] == G
                        && Pop[patch_i].breeders[breeder_i].epigenotype[allele_i] == z)
                {
                    ++n_active_alleles;
                }
            }

            assert(n_active_alleles >= 0);
            assert(n_active_alleles <= 2);


            // have individual die or not
            if (uniform(rng_r) < fitness_payoff(n_active_alleles, envt_is_A))
            {
                // remove breeder from the stack
                Pop[patch_i].breeders[breeder_i].erase(
                        Pop[patch_i].breeders[breeder_i].begin() + breeder_i);

                --breeder_i;
            }
        }
    }// end for (int patch_i = 0
} // end adult survival

 

// births of new offspring
// juvenile cue integration
// juvenile selection
// adult cue integration (aka horizontal social learning)
void replace()
{
    // auxiliary variable to keep track 
    // of floating value of clutch size
    // which will be rounded to a biologically 
    // relevant number of offspring later
    double clutch_d;
    bool envt_is_A;

    int clutch_i;

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int kid_i = 0; kid_i < NBreeders; ++kid_i)
        {

            // get lowest integer clutch size from the floating point
            // value
            clutch_i = floor(clutch_d);

            // of course, when rounding to the lowest integer, there is a 0<=remainder<=1
            // we draw a random number to see whether this remainder is an additional
            //
            if (uniform(rng_r) < clutch_d - clutch_i)
            {
                ++clutch_i;
            }

            // make offspring
            for (int egg_i = 0; egg_i < clutch_i; ++clutch_i)
            {
                // make new offspring
                Individual kid;

                // decide whether new offspring is immigrant
                // to the patch or not. If immigrant,
                // sample both parents from randomly chosen patch
                // (including the current patch with a probability 
                // of 1/NPatches)
                int patch_parents = uniform(rng_r) < d ? random_patch(rng_r) : patch_i;

                // generate a random sampler for this number of breeders
                uniform_int_distribution<> breeder_sampler(
                        0,
                        Pop[patch_parents].breeders.size() - 1);


                // inherit genes epigenes 
                // to offspring from parents,
                // where parents are randomly sampled
                // either locally or from a remote patch
                create_offspring(
                        Pop[patch_parents].breeders[breeder_sampler(rng_r)]
                        ,Pop[patch_parents].breeders[breeder_sampler(rng_r)]
                        ,kid
                        ,Pop[patch_i].envt_state_is_A);

                // add kid to the stack of future breeders
                Pop[patch_i].breeders_t1.append(kid);

            } // end for (int egg_i = 0;
        } // end for (int breeder_i = 0
    } // end for (int patch_i = 0

    // now replace breeders by breeders_t1 and change the environment
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        Pop[patch_i].breeders = Pop[patch_i].breeders_t1;

        // clear the next-gen vector as otherwise 
        // there'll be trouble
        Pop[patch_i].breeders_t1.clear();
    }
} // end replace()

// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    // create a file to output all the statistics about the population 
    string filename = "sim_epigenetics";
    // create unique filename
    create_filename(filename);
    ofstream DataFile(filename.c_str());  // output file 
    
    // prepare the output file by writing 
    // headers (which column is which
    write_data_headers(DataFile);
    
    // get command line arguments
    init_arguments(argc, argv);

    // initialize the population
    init_population();

    // auxiliary variable to store current generation
    int generation;

    for (generation = 0; generation < number_generations; ++generation)
    {
        // survival of adult breeders followed by reproduction
        adult_survival();

        // new breeder establishment, 
        // followed by horizontal learning
        // and finally enviromental change
        replace();

        if (generation % data_nth_generation == 0)
        {
            write_stats(DataFile, generation, 2);
        }
    }
            
    write_stats(DataFile, generation, 2);

    write_dist(DataFileDist);
    
    write_parameters(DataFile);
}
