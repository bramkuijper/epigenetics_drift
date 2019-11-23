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

// parameters & variables:
// a lot of the parameters declared here will be overridden
// in the init_arguments function which reads parameters from the command line

// number of individuals in population
const int NPatches = 400;

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


// genetic mutation rate of the genetic loci
double mu_g = 0.01;

// modifier mutation rate
double mu_mu = 0.01;
double sdmu_mu = 0.01;

// dispersal rates
double m = 0.0;

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

// survival function
double survival_probability(
        double const phen_ad, 
        bool const state_high)
{
    // eqns 3,4 Leimar & McNamara (2015) Amnat
    if (sigmoidal_survival)
    {
        double phen = state_high ? -phen_ad : phen_ad;
        return(1.0 / (1.0 + exp(survival_scalar[state_high] + 6.0 * phen)));
    }

    // eqns 1,2 Leimar & McNamara (2015) Amnat
    return(state_high ? 
            1.0 - survival_scalar[0] * (1.0 - phen_ad) * (1.0 - phen_ad)
            :
            1.0 - survival_scalar[0] * phen_ad * phen_ad
            );
}


// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{
    sigmoidal_survival = atoi(argv[1]);
    laplace = atoi(argv[2]);
    p = atof(argv[3]);
    survival_scalar[0] = atof(argv[4]);
    survival_scalar[1] = atof(argv[5]);
    qmat = atof(argv[6]);
    qjuv = atof(argv[7]);
    nloci_g = atoi(argv[8]);
    init_g = atof(argv[9]);
    init_amat = atof(argv[10]);
    init_ajuv = atof(argv[11]);
    init_agen = atof(argv[12]);
    init_asoc_horiz = atof(argv[13]);
    init_asoc_vert = atof(argv[14]);
    init_bmat_phen = atof(argv[15]);
    init_bmat_envt = atof(argv[16]);
    init_hp = atof(argv[17]);
    init_hc = atof(argv[18]);
    init_vp = atof(argv[19]);
    init_vc = atof(argv[20]);

    gmin = atof(argv[21]);
    gmax = atof(argv[22]);
    amin = atof(argv[23]);
    amax = atof(argv[24]);
    bmin = atof(argv[25]);
    bmax = atof(argv[26]);
    sdmat = atof(argv[27]);
    sdsoc_vert = atof(argv[28]);
    sdsoc_horiz = atof(argv[29]);

    mu_g = atof(argv[30]);
    mu_amat = atof(argv[31]);
    mu_ajuv = atof(argv[32]);
    mu_agen = atof(argv[33]);
    mu_asoc_horiz = atof(argv[34]);
    mu_asoc_vert = atof(argv[35]);
    mu_bmat_phen = atof(argv[36]);
    mu_bmat_envt = atof(argv[37]);
    mu_hp = atof(argv[38]);
    mu_hc = atof(argv[39]);
    mu_vp = atof(argv[40]);
    mu_vc = atof(argv[41]);
    sdmu_a = atof(argv[42]);
    sdmu_b = atof(argv[43]);
    sdmu_g = atof(argv[44]);
    m = atof(argv[45]);
    nph = atoi(argv[46]);
    nch = atoi(argv[47]);
    npv = atoi(argv[48]);
    ncv = atoi(argv[49]);
    juvenile_survival = atoi(argv[50]);
}

// write down all parameters to the file DataFile
void write_parameters(ofstream &DataFile)
{
    DataFile << endl << endl
        << "sigmoidal_survival;" << sigmoidal_survival << ";"<< endl
        << "laplace;" << laplace << ";"<< endl
        << "p;" << p << ";"<< endl
        << "qmat;" << qmat << ";"<< endl
        << "qjuv;" << qjuv << ";"<< endl
        << "nloci_g;" << nloci_g << ";"<< endl
        << "init_g;" << init_g << ";"<< endl
        << "init_amat;" << init_amat << ";"<< endl
        << "init_ajuv;" << init_ajuv << ";"<< endl
        << "init_agen;" << init_agen << ";"<< endl
        << "init_asoc_vert;" << init_asoc_vert << ";"<< endl
        << "init_asoc_horiz;" << init_asoc_horiz << ";"<< endl
        << "init_bmat_phen;" << init_bmat_phen << ";"<< endl
        << "init_bmat_envt;" << init_bmat_envt << ";"<< endl
        << "init_hp;" << init_hp << ";"<< endl
        << "init_hc;" << init_hc << ";"<< endl
        << "init_vp;" << init_vp << ";"<< endl
        << "init_vc;" << init_vc << ";"<< endl
        << "gmin;" << gmin << ";"<< endl
        << "gmax;" << gmax << ";"<< endl
        << "amin;" << amin << ";"<< endl
        << "amax;" << amax << ";"<< endl
        << "bmin;" << bmin << ";"<< endl
        << "bmax;" << bmax << ";"<< endl
        << "sdmat;" << sdmat << ";"<< endl
        << "sdsoc_vert;" << sdsoc_vert << ";"<< endl
        << "sdsoc_horiz;" << sdsoc_horiz << ";"<< endl
        << "mu_g;" << mu_g << ";"<< endl
        << "sdmu_g;" << sdmu_g << ";"<< endl
        << "mu_amat;" << mu_amat << ";"<< endl
        << "mu_ajuv;" << mu_ajuv << ";"<< endl
        << "mu_agen;" << mu_agen << ";"<< endl
        << "mu_asoc_horiz;" << mu_asoc_horiz << ";"<< endl
        << "juvenile_survival;" << juvenile_survival << ";"<< endl
        << "mu_asoc_vert;" << mu_asoc_vert << ";"<< endl
        << "mu_bmat_phen;" << mu_bmat_phen << ";"<< endl
        << "mu_bmat_envt;" << mu_bmat_envt << ";"<< endl
        << "mu_vc;" << mu_vc << ";"<< endl
        << "mu_vp;" << mu_vp << ";"<< endl
        << "mu_hc;" << mu_hc << ";"<< endl
        << "mu_hp;" << mu_hp << ";"<< endl
        << "sdmu_a;" << sdmu_a << ";"<< endl
        << "sdmu_b;" << sdmu_b << ";"<< endl
        << "m;" << m << ";"<< endl
        << "nph;" << nph << ";"<< endl
        << "nch;" << nch << ";"<< endl
        << "npv;" << npv << ";"<< endl
        << "ncv;" << ncv << ";"<< endl
        << "survival_scalar0;" << survival_scalar[0] << ";"<< endl
        << "survival_scalar1;" << survival_scalar[1] << ";"<< endl
        << "seed;" << seed << ";" << endl;

} // void write_parameters(ofstream &DataFile)


// write all properties of all individuals
// to the file DataFile (to obtain information about
// the distribution of traits)
void write_dist(ofstream &DataFile)
{
    double g; // auxiliary variable to temporarily 
                // store trait expression

    
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            DataFile << patch_i << ";" 
                << breeder_i << ";"
                << Pop[patch_i].breeders[breeder_i].phen_ad << ";"
                << Pop[patch_i].breeders[breeder_i].phen_mat << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_vert << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_horiz << ";"
                << Pop[patch_i].breeders[breeder_i].xmat << ";"
                << Pop[patch_i].breeders[breeder_i].xsoc_vert << ";"
                << Pop[patch_i].breeders[breeder_i].xsoc_horiz << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_vert << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_horiz << ";"

                // agen
                << 0.5 * (Pop[patch_i].breeders[breeder_i].agen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1]) << ";"

                // ajuv
                << 0.5 * (Pop[patch_i].breeders[breeder_i].ajuv[0]
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1]) << ";"

                // amat
                << 0.5 * (Pop[patch_i].breeders[breeder_i].amat[0]
                    +
                    Pop[patch_i].breeders[breeder_i].amat[1]) << ";"

                // asoc_vert 
                << 0.5 * (Pop[patch_i].breeders[breeder_i].asoc_vert[0]
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_vert[1]) << ";"
                
                // asoc_horiz 
                << 0.5 * (Pop[patch_i].breeders[breeder_i].asoc_horiz[0]
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[1]) << ";"

                // bmat_phen
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_phen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1]) << ";"

                // bmat_envt
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_envt[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1]) << ";"

                // hc
                << 0.5 * (Pop[patch_i].breeders[breeder_i].hc[0]
                    +
                    Pop[patch_i].breeders[breeder_i].hc[1]) << ";"

                // hp
                << 0.5 * (Pop[patch_i].breeders[breeder_i].hp[0]
                    +
                    Pop[patch_i].breeders[breeder_i].hp[1]) << ";"
                
                // vc
                << 0.5 * (Pop[patch_i].breeders[breeder_i].vc[0]
                    +
                    Pop[patch_i].breeders[breeder_i].vc[1]) << ";"

                // vp
                << 0.5 * (Pop[patch_i].breeders[breeder_i].vp[0]
                    +
                    Pop[patch_i].breeders[breeder_i].vp[1]) << ";";

            g = 0.0;

            // genetic cue
            for (int g_i = 0; g_i < nloci_g; ++g_i)
            {
                g += 0.5 * (
                    Pop[patch_i].breeders[breeder_i].g[0][g_i]
                    +
                    Pop[patch_i].breeders[breeder_i].g[1][g_i]);
            }

            DataFile << g << ";"
                << Pop[patch_i].envt_high << ";"
                << Pop[patch_i].breeders[breeder_i].cue_ad_envt_high << ";"
                << Pop[patch_i].breeders[breeder_i].cue_juv_envt_high << ";"
                << Pop[patch_i].breeders[breeder_i].maternal_cue << ";"
                << Pop[patch_i].breeders[breeder_i].mnoise << ";"
                << Pop[patch_i].breeders[breeder_i].svnoise << ";"
                << Pop[patch_i].breeders[breeder_i].shnoise << ";"
                << endl;
        } // end for (int breeder_i = 0; breeder_i < NPatches; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)
} // end  void write_dist(ofstream &DataFile)


// list of the data headers at the start of the file
// in which the distribution of evolved trait values
// and states is written at the end of the file
void write_data_headers_dist(ofstream &DataFile)
{
    DataFile 
        << "patch_id;" 
        << "id;" 
        << "phen_ad;" 
        << "phen_mat;" 
        << "phen_prestige_vert;" 
        << "phen_prestige_horiz;" 
        << "xmat;" 
        << "xsoc_vert;" 
        << "xsoc_horiz;" 
        << "xconformist_vert;" 
        << "xconformist_horiz;" 
        << "agen;" 
        << "ajuv;" 
        << "amat;" 
        << "asoc_vert;" 
        << "asoc_horiz;" 
        << "bmat_phen;" 
        << "bmat_envt;" 
        << "hc;" 
        << "hp;" 
        << "vc;" 
        << "vp;" 
        << "g;" 
        << "envt;" 
        << "cue_ad_envt_high;" 
        << "cue_juv_envt_high;" 
        << "maternal_cue;" 
        << "mnoise;" 
        << "svnoise;" 
        << "shnoise;" 
        << endl; 
}

// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile 
        << "generation;" 
        << "mean_phen_ad;" 
        << "mean_phen_juv;" 
        << "mean_phen_prestige_vert;" 
        << "mean_phen_prestige_horiz;" 
        << "mean_agen;" 
        << "mean_ajuv;" 
        << "mean_amat;" 
        << "mean_asoc_vert;" 
        << "mean_asoc_horiz;" 
        << "mean_bmat_phen;" 
        << "mean_bmat_envt;" 
        << "mean_hc;" 
        << "mean_hp;" 
        << "mean_vc;" 
        << "mean_vp;" 
        << "mean_g;" 
        << "var_phen_ad;" 
        << "var_phen_juv;" 
        << "var_phen_prestige_vert;" 
        << "var_phen_prestige_horiz;" 
        << "var_agen;" 
        << "var_ajuv;" 
        << "var_amat;" 
        << "var_asoc_vert;" 
        << "var_asoc_horiz;" 
        << "var_bmat_phen;" 
        << "var_bmat_envt;" 
        << "var_hc;" 
        << "var_hp;" 
        << "var_vc;" 
        << "var_vp;" 
        << "var_g;" 
        << "freq_high;" 
        << "mean_surv0;" 
        << "mean_surv1;" 
        << "var_surv0;" 
        << "var_surv1;" 
        << "var_component_gen;"
        << "var_component_ajuv;"
        << "var_component_amat;"
        << "var_component_amat_envt;"
        << "var_component_amat_phen;"
        << "cov_amat_phen_amat_envt;"
        << "cov_amat_ajuv;"
        << "cov_amat_envt_ajuv;"
        << "cov_amat_phen_ajuv;"
        << "var_component_asoc_vert;"
        << "var_component_asoc_vert_p;"
        << "var_component_asoc_vert_c;"
        << "var_component_asoc_horiz;"
        << "var_component_asoc_horiz_p;"
        << "var_component_asoc_horiz_c;"
        << "cov_amat_asoc_vert;"
        << "cov_amat_asoc_vert_c;"
        << "cov_amat_asoc_vert_p;"
        << "cov_amat_asoc_horiz;"
        << "cov_amat_asoc_horiz_c;"
        << "cov_amat_asoc_horiz_p;"
        << "cov_amat_envt_asoc_vert;"
        << "cov_amat_envt_asoc_vert_c;"
        << "cov_amat_envt_asoc_vert_p;"
        << "cov_amat_envt_asoc_horiz;"
        << "cov_amat_envt_asoc_horiz_c;"
        << "cov_amat_envt_asoc_horiz_p;"
        << "cov_amat_phen_asoc_vert;"
        << "cov_amat_phen_asoc_vert_c;"
        << "cov_amat_phen_asoc_vert_p;"
        << "cov_amat_phen_asoc_horiz;"
        << "cov_amat_phen_asoc_horiz_c;"
        << "cov_amat_phen_asoc_horiz_p;"

        << "cov_ajuv_asoc_vert;"
        << "cov_ajuv_asoc_vert_c;"
        << "cov_ajuv_asoc_vert_p;"
        
        << "cov_ajuv_asoc_horiz;"
        << "cov_ajuv_asoc_horiz_c;"
        << "cov_ajuv_asoc_horiz_p;"

        << "cov_agen_asoc_vert;"
        << "cov_agen_asoc_vert_c;"
        << "cov_agen_asoc_vert_p;"

        << "cov_agen_asoc_horiz;"
        << "cov_agen_asoc_horiz_c;"
        << "cov_agen_asoc_horiz_p;"

        << "cov_agen_ajuv;"
        << endl; 
}

// write data both for winter and summer populations
void write_stats(ofstream &DataFile, int generation, int timestep)
{
    // variables to store means and variances
    double mean_phen_ad = 0.0;
    double ss_phen_ad = 0.0;
    
    double mean_phen_juv = 0.0;
    double ss_phen_juv = 0.0;
    
    double mean_phen_prestige_vert = 0.0;
    double ss_phen_prestige_vert = 0.0;
    
    double mean_phen_prestige_horiz = 0.0;
    double ss_phen_prestige_horiz = 0.0;
   
    double mean_agen = 0.0;
    double ss_agen = 0.0;
   
    double mean_amat = 0.0;
    double ss_amat = 0.0;
    
    double mean_ajuv = 0.0;
    double ss_ajuv = 0.0;
    
    double mean_asoc_vert = 0.0;
    double ss_asoc_vert = 0.0;
    
    double mean_asoc_horiz = 0.0;
    double ss_asoc_horiz = 0.0;
    
    double mean_bmat_phen = 0.0;
    double ss_bmat_phen = 0.0;
    
    double mean_bmat_envt = 0.0;
    double ss_bmat_envt = 0.0;
    
    double mean_hc = 0.0;
    double ss_hc = 0.0;

    double mean_hp = 0.0;
    double ss_hp = 0.0;
    
    double mean_vc = 0.0;
    double ss_vc = 0.0;

    double mean_vp = 0.0;
    double ss_vp = 0.0;

    double mean_g = 0.0;
    double ss_g = 0.0;

    double freq_high = 0.0;

    // now some stats at the logit scale
    // maternal variance relative to total variance
    double ss1_gen_component = 0.0;
    double ss2_gen_component = 0.0;

    double ss1_ajuv_component = 0.0;
    double ss2_ajuv_component = 0.0;

    double ss1_amat_component = 0.0;
    double ss2_amat_component = 0.0;

    double ss1_amat_envt_component = 0.0;
    double ss2_amat_envt_component = 0.0;

    double ss1_amat_phen_component = 0.0;
    double ss2_amat_phen_component = 0.0;

    double ss1_asoc_vert_component = 0.0;
    double ss2_asoc_vert_component = 0.0;
    
    double ss1_asoc_horiz_component = 0.0;
    double ss2_asoc_horiz_component = 0.0;

    double ss1_asoc_vert_c_component = 0.0;
    double ss2_asoc_vert_c_component = 0.0;
    double ss1_asoc_vert_p_component = 0.0;
    double ss2_asoc_vert_p_component = 0.0;
    
    double ss1_asoc_horiz_c_component = 0.0;
    double ss2_asoc_horiz_c_component = 0.0;
    double ss1_asoc_horiz_p_component = 0.0;
    double ss2_asoc_horiz_p_component = 0.0;

    double ss2_cov_agen_ajuv = 0.0;
    double ss2_cov_amat_phen_amat_envt = 0.0;

    double ss2_cov_amat_ajuv = 0.0;
    double ss2_cov_amat_envt_ajuv = 0.0;
    double ss2_cov_amat_phen_ajuv = 0.0;

    double ss2_cov_agen_asoc_vert = 0.0;
    double ss2_cov_agen_asoc_vert_c = 0.0;
    double ss2_cov_agen_asoc_vert_p = 0.0;

    double ss2_cov_ajuv_asoc_vert = 0.0;
    double ss2_cov_ajuv_asoc_vert_c = 0.0;
    double ss2_cov_ajuv_asoc_vert_p = 0.0;
    
    double ss2_cov_ajuv_asoc_horiz = 0.0;
    double ss2_cov_ajuv_asoc_horiz_c = 0.0;
    double ss2_cov_ajuv_asoc_horiz_p = 0.0;

    double ss2_cov_amat_asoc_vert = 0.0;
    double ss2_cov_amat_asoc_vert_c = 0.0;
    double ss2_cov_amat_asoc_vert_p = 0.0;
    
    double ss2_cov_amat_asoc_horiz = 0.0;
    double ss2_cov_amat_asoc_horiz_c = 0.0;
    double ss2_cov_amat_asoc_horiz_p = 0.0;

    double ss2_cov_amat_phen_asoc_vert = 0.0;
    double ss2_cov_amat_phen_asoc_vert_c = 0.0;
    double ss2_cov_amat_phen_asoc_vert_p = 0.0;
    
    double ss2_cov_amat_envt_asoc_vert = 0.0;
    double ss2_cov_amat_envt_asoc_vert_c = 0.0;
    double ss2_cov_amat_envt_asoc_vert_p = 0.0;
    
    double ss2_cov_amat_phen_asoc_horiz = 0.0;
    double ss2_cov_amat_phen_asoc_horiz_c = 0.0;
    double ss2_cov_amat_phen_asoc_horiz_p = 0.0;
    
    double ss2_cov_amat_envt_asoc_horiz = 0.0;
    double ss2_cov_amat_envt_asoc_horiz_c = 0.0;
    double ss2_cov_amat_envt_asoc_horiz_p = 0.0;

    double ss2_cov_agen_asoc_horiz = 0.0;
    double ss2_cov_agen_asoc_horiz_c = 0.0;
    double ss2_cov_agen_asoc_horiz_p = 0.0;

    // auxiliary variables to calculate an individual's phenotype
    double g, z, agen, amat, ajuv, xmat_envt, xmat, xmat_phen, xsoc_vert, xsoc_vert_p, xsoc_vert_c, asoc_horiz, asoc_vert, cue_ad, xsoc_horiz_p, xsoc_horiz_c, xsoc_horiz;


    // summing means and sums of squares over all patches and breeders
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        freq_high += Pop[patch_i].envt_high;

        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            g = 0.0;
            for (int g_i = 0; g_i < nloci_g; ++g_i)
            {
                g +=  0.5 * (
                        Pop[patch_i].breeders[breeder_i].g[0][g_i]
                        +
                        Pop[patch_i].breeders[breeder_i].g[1][g_i]
                        );

            } // end for (int g_i = 0; g_i < nloci_g; ++g_i)

            mean_g += g;
            ss_g += g * g;

            // adult phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_ad;
            mean_phen_ad += z;
            ss_phen_ad += z * z;
            
            // juvenile phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_juv;
            mean_phen_juv += z;
            ss_phen_juv += z * z;
            
            // prestige phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_prestige_vert;
            mean_phen_prestige_vert += z;
            ss_phen_prestige_vert += z * z;
            
            // prestige phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_prestige_horiz;
            mean_phen_prestige_horiz += z;
            ss_phen_prestige_horiz += z * z;

            // sensitivity to genetic cues
            agen = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].agen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1] 
                    );

            mean_agen += agen;
            ss_agen += agen * agen;

            // variance components (liability scale)
            // sum of squares for the 
            // genetic variance component is
            // var(agen * sum(g)) = E[agen^2 * sum(g)^2] - E[agen * sum(g)]^2
            ss1_gen_component += agen * g;
            ss2_gen_component += agen * agen * g * g;

            // sensitivity to maternal cues
            amat = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].amat[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].amat[1] 
                    );

            mean_amat += amat;
            ss_amat += amat * amat;

            xmat_envt = Pop[patch_i].breeders[breeder_i].xmat_envt_only;
            xmat_phen  = Pop[patch_i].breeders[breeder_i].xmat_phen_only;
            xmat = Pop[patch_i].breeders[breeder_i].xmat;

            // variance components (liability scale)
            // sum of squares for the 
            // maternal environmental component
            ss1_amat_envt_component += amat * xmat_envt;
            ss2_amat_envt_component += amat * amat * xmat_envt * xmat_envt;

            // variance components (liability scale)
            // sum of squares for the 
            // maternal phenotypic component
            ss1_amat_phen_component += amat * xmat_phen;
            ss2_amat_phen_component += amat * amat * xmat_phen * xmat_phen;

            // variance components (liability scale)
            // sum of squares for the 
            // total maternal component
            ss1_amat_component += amat * xmat;
            ss2_amat_component += amat * amat * xmat * xmat;

            // covariance between different types of maternal effect
            ss2_cov_amat_phen_amat_envt += amat * xmat_envt * amat * xmat_phen;

            // sensitivity to juvenile cues
            ajuv = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].ajuv[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1] 
                    );
            
            mean_ajuv += ajuv;
            ss_ajuv += ajuv * ajuv;

            cue_ad = Pop[patch_i].breeders[breeder_i].cue_ad_envt_high;

            // variance components juvenile cue
            ss1_ajuv_component += ajuv * cue_ad;
            ss2_ajuv_component += ajuv * ajuv * cue_ad * cue_ad;

            // covariance between juvenile cue and maternal total effect
            // E[abcd] - E[ab]E[cd] (the latter product was what we already
            // calculated before)
            ss2_cov_amat_ajuv += amat * xmat * ajuv * cue_ad;
            ss2_cov_amat_envt_ajuv += amat * xmat_envt * ajuv * cue_ad;
            ss2_cov_amat_phen_ajuv += amat * xmat_phen * ajuv * cue_ad;


            ss2_cov_agen_ajuv += ajuv * cue_ad * agen * g;

            // sensitivity to socially learnt cues (vertically)
            asoc_vert = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].asoc_vert[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_vert[1] 
                    );

            mean_asoc_vert += asoc_vert;
            ss_asoc_vert += asoc_vert * asoc_vert;

            xsoc_vert = Pop[patch_i].breeders[breeder_i].xsoc_vert;
            xsoc_vert_c = Pop[patch_i].breeders[breeder_i].xsoc_vert_c;
            xsoc_vert_p = Pop[patch_i].breeders[breeder_i].xsoc_vert_p;

            // variance due to the total socially learnt component
            ss1_asoc_vert_component += asoc_vert * xsoc_vert;
            ss2_asoc_vert_component += asoc_vert * asoc_vert * xsoc_vert * xsoc_vert;

            // variance due to the conformity bias socially learnt component
            ss1_asoc_vert_c_component += asoc_vert * xsoc_vert_c;
            ss2_asoc_vert_c_component += asoc_vert * asoc_vert * xsoc_vert_c * xsoc_vert_c;
            
            // variance due to the prestige bias socially learnt component
            ss1_asoc_vert_p_component += asoc_vert * xsoc_vert_p;
            ss2_asoc_vert_p_component += asoc_vert * asoc_vert * xsoc_vert_p * xsoc_vert_p;
            
            // covariance between maternal and vertically learnt components
            ss2_cov_amat_asoc_vert += amat * xmat * asoc_vert * xsoc_vert;
            ss2_cov_amat_asoc_vert_c += amat * xmat * asoc_vert * xsoc_vert_c;
            ss2_cov_amat_asoc_vert_p += amat * xmat * asoc_vert * xsoc_vert_p;

            // covariance between maternal environment and vertically learnt components
            ss2_cov_amat_envt_asoc_vert += amat * xmat_envt * asoc_vert * xsoc_vert;
            ss2_cov_amat_envt_asoc_vert_c += amat * xmat_envt * asoc_vert * xsoc_vert_c;
            ss2_cov_amat_envt_asoc_vert_p += amat * xmat_envt * asoc_vert * xsoc_vert_p;

            // covariance between maternal phenotype and vertically learnt components
            ss2_cov_amat_phen_asoc_vert += amat * xmat_phen * asoc_vert * xsoc_vert;
            ss2_cov_amat_phen_asoc_vert_c += amat * xmat_phen * asoc_vert * xsoc_vert_c;
            ss2_cov_amat_phen_asoc_vert_p += amat * xmat_phen * asoc_vert * xsoc_vert_p;
            
            // covariance between juvenile cue and vertically learnt components
            ss2_cov_ajuv_asoc_vert += ajuv * cue_ad * asoc_vert * xsoc_vert;
            ss2_cov_ajuv_asoc_vert_c += ajuv * cue_ad * asoc_vert * xsoc_vert_c;
            ss2_cov_ajuv_asoc_vert_p += ajuv * cue_ad * asoc_vert * xsoc_vert_p;

            // covariance between genetic cue and vertically learnt components
            ss2_cov_agen_asoc_vert += agen * g * asoc_vert * xsoc_vert;
            ss2_cov_agen_asoc_vert_c += agen * g * asoc_vert * xsoc_vert_c;
            ss2_cov_agen_asoc_vert_p += agen * g * asoc_vert * xsoc_vert_p;

            // sensitivity to socially learnt cues (horizontally)
            asoc_horiz = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[1] 
                    );

            mean_asoc_horiz += asoc_horiz;
            ss_asoc_horiz += asoc_horiz * asoc_horiz;

            xsoc_horiz = Pop[patch_i].breeders[breeder_i].xsoc_horiz;
            xsoc_horiz_c = Pop[patch_i].breeders[breeder_i].xsoc_horiz_c;
            xsoc_horiz_p = Pop[patch_i].breeders[breeder_i].xsoc_horiz_p;

            ss1_asoc_horiz_component += asoc_horiz * xsoc_horiz;
            ss2_asoc_horiz_component += asoc_horiz * asoc_horiz * xsoc_horiz * xsoc_horiz;

            ss1_asoc_horiz_c_component += asoc_horiz * xsoc_horiz_c;
            ss2_asoc_horiz_c_component += asoc_horiz * asoc_horiz * xsoc_horiz_c * xsoc_horiz_c;

            ss1_asoc_horiz_p_component += asoc_horiz * xsoc_horiz_p;
            ss2_asoc_horiz_p_component += asoc_horiz * asoc_horiz * xsoc_horiz_p * xsoc_horiz_p;

            // covariance between maternal and vertically learnt components
            ss2_cov_amat_asoc_horiz += amat * xmat * asoc_horiz * xsoc_horiz;
            ss2_cov_amat_asoc_horiz_c += amat * xmat * asoc_horiz * xsoc_horiz_c;
            ss2_cov_amat_asoc_horiz_p += amat * xmat * asoc_horiz * xsoc_horiz_p;
            
            // covariance between maternal environmental cue and vertically learnt components
            ss2_cov_amat_envt_asoc_horiz += amat * xmat_envt * asoc_horiz * xsoc_horiz;
            ss2_cov_amat_envt_asoc_horiz_c += amat * xmat_envt * asoc_horiz * xsoc_horiz_c;
            ss2_cov_amat_envt_asoc_horiz_p += amat * xmat_envt * asoc_horiz * xsoc_horiz_p;
            
            // covariance between maternal phenotypic cue and vertically learnt components
            ss2_cov_amat_phen_asoc_horiz += amat * xmat_phen * asoc_horiz * xsoc_horiz;
            ss2_cov_amat_phen_asoc_horiz_c += amat * xmat_phen * asoc_horiz * xsoc_horiz_c;
            ss2_cov_amat_phen_asoc_horiz_p += amat * xmat_phen * asoc_horiz * xsoc_horiz_p;
            
            // covariance between juvenile cue and horizically learnt components
            ss2_cov_ajuv_asoc_horiz += ajuv * cue_ad * asoc_horiz * xsoc_horiz;
            ss2_cov_ajuv_asoc_horiz_c += ajuv * cue_ad * asoc_horiz * xsoc_horiz_c;
            ss2_cov_ajuv_asoc_horiz_p += ajuv * cue_ad * asoc_horiz * xsoc_horiz_p;

            // covariance between genetic cue and horizically learnt components
            ss2_cov_agen_asoc_horiz += agen * g * asoc_horiz * xsoc_horiz;
            ss2_cov_agen_asoc_horiz_c += agen * g * asoc_horiz * xsoc_horiz_c;
            ss2_cov_agen_asoc_horiz_p += agen * g * asoc_horiz * xsoc_horiz_p;
            
            // maternal sensitivity to phenotypic cues
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_phen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1] 
                    );

            mean_bmat_phen += z;
            ss_bmat_phen += z * z;
            
            // maternal sensitivity to environmental cues
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_envt[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1] 
                    );

            mean_bmat_envt += z;
            ss_bmat_envt += z * z;
            
            // sensitivity to performance-based cues when learning
            // horizontally
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hp[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].hp[1] 
                    );

            mean_hp += z;
            ss_hp += z * z;
            
            // sensitivity to confirmity-based cues when learning
            // horizontally
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hc[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].hc[1] 
                    );

            mean_hc += z;
            ss_hc += z * z;
            
            
            // sensitivity to performance-based cues when learning
            // vertically 
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].vp[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].vp[1] 
                    );

            mean_vp += z;
            ss_vp += z * z;
            
            // sensitivity to confirmity-based cues when learning
            // vertically
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].vc[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].vc[1] 
                    );

            mean_vc += z;
            ss_vc += z * z;
        } // end for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)

    mean_phen_ad /= NPatches * NBreeder;
    double var_phen_ad = ss_phen_ad / (NPatches * NBreeder) - mean_phen_ad * mean_phen_ad;
    
    mean_phen_juv /= NPatches * NBreeder;
    double var_phen_juv = ss_phen_juv / (NPatches * NBreeder) - mean_phen_juv * mean_phen_juv;
    
    mean_phen_prestige_vert /= NPatches * NBreeder;

    double var_phen_prestige_vert = ss_phen_prestige_vert / (NPatches * NBreeder) 
        - mean_phen_prestige_vert * mean_phen_prestige_vert;
   
    mean_phen_prestige_horiz /= NPatches * NBreeder;

    double var_phen_prestige_horiz = ss_phen_prestige_horiz / (NPatches * NBreeder) 
        - mean_phen_prestige_horiz * mean_phen_prestige_horiz;

    mean_agen /= NPatches * NBreeder;
    double var_agen = ss_agen / (NPatches * NBreeder) - mean_agen * mean_agen;
   
    mean_amat /= NPatches * NBreeder;
    double var_amat = NPatches * NBreeder - mean_amat * mean_amat;
    
    mean_ajuv /= NPatches * NBreeder;
    double var_ajuv = ss_ajuv / (NPatches * NBreeder) -  mean_ajuv * mean_ajuv;
    
    mean_asoc_vert /= NPatches * NBreeder;
    double var_asoc_vert = ss_asoc_vert / (NPatches * NBreeder) - mean_asoc_vert * mean_asoc_vert;
    
    mean_asoc_horiz /= NPatches * NBreeder;
    double var_asoc_horiz = ss_asoc_horiz / (NPatches * NBreeder) - mean_asoc_horiz * mean_asoc_horiz;
    
    mean_bmat_phen /= NPatches * NBreeder;
    double var_bmat_phen = ss_bmat_phen / (NPatches * NBreeder) - mean_bmat_phen * mean_bmat_phen;
    
    mean_bmat_envt /= NPatches * NBreeder;
    double var_bmat_envt = ss_bmat_envt / (NPatches * NBreeder) - mean_bmat_envt * mean_bmat_envt;
    
    mean_hp /= NPatches * NBreeder;
    double var_hp = ss_hp / (NPatches * NBreeder) - mean_hp * mean_hp;
    
    mean_hc /= NPatches * NBreeder;
    double var_hc = ss_hc / (NPatches * NBreeder) - mean_hc * mean_hc;
    
    mean_vp /= NPatches * NBreeder;
    double var_vp = ss_vp / (NPatches * NBreeder) - mean_vp * mean_vp;
    
    mean_vc /= NPatches * NBreeder;
    double var_vc = ss_vc / (NPatches * NBreeder) - mean_vc * mean_vc;

    mean_g /= NPatches * NBreeder;
    double var_g = ss_g / (NPatches * NBreeder) - mean_g * mean_g;



    // Variance components at the logit scale

    double var_component_gen = ss2_gen_component / (NPatches * NBreeder) - 
        pow(ss1_gen_component / (NPatches * NBreeder),2);
    
    double var_component_ajuv  = ss2_ajuv_component / (NPatches * NBreeder) - 
        pow(ss1_ajuv_component / (NPatches * NBreeder),2);

    double var_component_amat  = ss2_amat_component / (NPatches * NBreeder) - 
        pow(ss1_amat_component / (NPatches * NBreeder),2);

    double var_component_amat_envt  = ss2_amat_envt_component / (NPatches * NBreeder) - 
        pow(ss1_amat_envt_component / (NPatches * NBreeder),2);
    
    double var_component_amat_phen  = ss2_amat_phen_component / (NPatches * NBreeder) - 
        pow(ss1_amat_phen_component / (NPatches * NBreeder),2);
  
    // covariance between environmental and phenotypic maternal effects
    double cov_amat_phen_amat_envt = ss2_cov_amat_phen_amat_envt / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_amat_phen_component / (NPatches * NBreeder);

    // covariance between juvenile cue and maternal effects
    double cov_amat_ajuv = ss2_cov_amat_ajuv / (NPatches * NBreeder) - 
        ss1_amat_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);
    
    double cov_amat_envt_ajuv = ss2_cov_amat_envt_ajuv / (NPatches * NBreeder) - 
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);
    
    double cov_amat_phen_ajuv = ss2_cov_amat_phen_ajuv / (NPatches * NBreeder) - 
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);

    // variance of vertical social learning: total, conformism and prestige
    double var_component_asoc_vert = ss2_asoc_vert_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_vert_component / (NPatches * NBreeder),2);

    double var_component_asoc_vert_p = ss2_asoc_vert_p_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_vert_p_component / (NPatches * NBreeder),2);

    double var_component_asoc_vert_c = ss2_asoc_vert_c_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_vert_c_component / (NPatches * NBreeder),2);

    // variance of horizontal social learning: total, conformism and prestige
    double var_component_asoc_horiz = ss2_asoc_horiz_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_horiz_component / (NPatches * NBreeder),2);

    double var_component_asoc_horiz_p = ss2_asoc_horiz_p_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_horiz_p_component / (NPatches * NBreeder),2);

    double var_component_asoc_horiz_c = ss2_asoc_horiz_c_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_horiz_c_component / (NPatches * NBreeder),2);



    
    // covariance between maternal total cue and vertical social learning
    double cov_amat_asoc_vert = ss2_cov_amat_asoc_vert / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    double cov_amat_asoc_vert_c = ss2_cov_amat_asoc_vert_c / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);
    
    double cov_amat_asoc_vert_p = ss2_cov_amat_asoc_vert_p / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);


    
    // covariance between maternal total cue and vertical social learning
    double cov_amat_asoc_horiz = ss2_cov_amat_asoc_horiz / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    double cov_amat_asoc_horiz_c = ss2_cov_amat_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);
    
    double cov_amat_asoc_horiz_p = ss2_cov_amat_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);



    // covariance between maternal environmental cue and vertical social learning
    double cov_amat_envt_asoc_vert = ss2_cov_amat_envt_asoc_vert / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    double cov_amat_envt_asoc_vert_c = ss2_cov_amat_envt_asoc_vert_c / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);
    
    double cov_amat_envt_asoc_vert_p = ss2_cov_amat_envt_asoc_vert_p / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);


    // covariance between maternal environmental cue and horizontal social learning
    double cov_amat_envt_asoc_horiz = ss2_cov_amat_envt_asoc_horiz / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    double cov_amat_envt_asoc_horiz_c = ss2_cov_amat_envt_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);
    
    double cov_amat_envt_asoc_horiz_p = ss2_cov_amat_envt_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);


    // covariance between maternal phenotypic cue and vertical social learning
    double cov_amat_phen_asoc_vert = ss2_cov_amat_phen_asoc_vert / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    double cov_amat_phen_asoc_vert_c = ss2_cov_amat_phen_asoc_vert_c / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);
    
    double cov_amat_phen_asoc_vert_p = ss2_cov_amat_phen_asoc_vert_p / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);


    // covariance between maternal phenotypic cue and horizontal social learning
    double cov_amat_phen_asoc_horiz = ss2_cov_amat_phen_asoc_horiz / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    double cov_amat_phen_asoc_horiz_c = ss2_cov_amat_phen_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);
    
    double cov_amat_phen_asoc_horiz_p = ss2_cov_amat_phen_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);


    // covariance between vertical social learning and juvenile cues
    double cov_ajuv_asoc_vert = ss2_cov_ajuv_asoc_vert / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_vert_p = ss2_cov_ajuv_asoc_vert_p / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_vert_c = ss2_cov_ajuv_asoc_vert_c / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);

    
    // covariance between horizontal social learning and juvenile cues
    double cov_ajuv_asoc_horiz = ss2_cov_ajuv_asoc_horiz / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_horiz_p = ss2_cov_ajuv_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_horiz_c = ss2_cov_ajuv_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);


    // covariance between vertical social learning and genetic cues
    double cov_agen_asoc_vert = ss2_cov_agen_asoc_vert / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_vert_p = ss2_cov_agen_asoc_vert_p / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_vert_c = ss2_cov_agen_asoc_vert_c / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);

    // covariance between horizontal social learning and genetic cues
    double cov_agen_asoc_horiz = ss2_cov_agen_asoc_horiz / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_horiz_p = ss2_cov_agen_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_horiz_c = ss2_cov_agen_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);

    // covariance between genetic cue and juvenile cue 
    double cov_agen_ajuv = ss2_cov_agen_ajuv / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);

    freq_high /= NPatches;

    DataFile << generation << ";"
        << mean_phen_ad << ";"
        << mean_phen_juv << ";"
        << mean_phen_prestige_vert << ";"
        << mean_phen_prestige_horiz << ";"
        << mean_agen << ";"
        << mean_ajuv << ";"
        << mean_amat << ";"
        << mean_asoc_vert << ";"
        << mean_asoc_horiz << ";"
        << mean_bmat_phen << ";"
        << mean_bmat_envt << ";"
        << mean_hc << ";"
        << mean_hp << ";"
        << mean_vc << ";"
        << mean_vp << ";"
        << mean_g << ";"
        << var_phen_ad << ";"
        << var_phen_juv << ";"
        << var_phen_prestige_vert << ";"
        << var_phen_prestige_horiz << ";"
        << var_agen << ";"
        << var_ajuv << ";"
        << var_amat << ";"
        << var_asoc_vert << ";"
        << var_asoc_horiz << ";"
        << var_bmat_phen << ";"
        << var_bmat_envt << ";"
        << var_hc << ";"
        << var_hp << ";"
        << var_vc << ";"
        << var_vp << ";"
        << var_g << ";" 
        << freq_high << ";" 
        << mean_survival[0] << ";" 
        << mean_survival[1] << ";" 
        << var_survival[0] << ";" 
        << var_survival[1] << ";" 

        << var_component_gen << ";"
        << var_component_ajuv << ";"
        << var_component_amat << ";"
        << var_component_amat_envt << ";"
        << var_component_amat_phen << ";"
        << cov_amat_phen_amat_envt << ";"

        << cov_amat_ajuv << ";"
        << cov_amat_envt_ajuv << ";"
        << cov_amat_phen_ajuv << ";"

        << var_component_asoc_vert << ";"
        << var_component_asoc_vert_p << ";"
        << var_component_asoc_vert_c << ";"

        << var_component_asoc_horiz << ";"
        << var_component_asoc_horiz_p << ";"
        << var_component_asoc_horiz_c << ";"

        << cov_amat_asoc_vert << ";"
        << cov_amat_asoc_vert_c << ";"
        << cov_amat_asoc_vert_p << ";"

        << cov_amat_asoc_horiz << ";"
        << cov_amat_asoc_horiz_c << ";"
        << cov_amat_asoc_horiz_p << ";"

        << cov_amat_envt_asoc_vert << ";"
        << cov_amat_envt_asoc_vert_c << ";"
        << cov_amat_envt_asoc_vert_p << ";"
        
        << cov_amat_envt_asoc_horiz << ";"
        << cov_amat_envt_asoc_horiz_c << ";"
        << cov_amat_envt_asoc_horiz_p << ";"

        << cov_amat_phen_asoc_vert << ";"
        << cov_amat_phen_asoc_vert_c << ";"
        << cov_amat_phen_asoc_vert_p << ";"

        << cov_amat_phen_asoc_horiz << ";"
        << cov_amat_phen_asoc_horiz_c << ";"
        << cov_amat_phen_asoc_horiz_p << ";"

        << cov_ajuv_asoc_vert << ";"
        << cov_ajuv_asoc_vert_c << ";"
        << cov_ajuv_asoc_vert_p << ";"

        << cov_ajuv_asoc_horiz << ";"
        << cov_ajuv_asoc_horiz_c << ";"
        << cov_ajuv_asoc_horiz_p << ";"

        << cov_agen_asoc_vert << ";"
        << cov_agen_asoc_vert_c << ";"
        << cov_agen_asoc_vert_p << ";"
        
        << cov_agen_asoc_horiz << ";"
        << cov_agen_asoc_horiz_c << ";"
        << cov_agen_asoc_horiz_p << ";"

        << cov_agen_ajuv << ";"
        << endl;
}

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
                ind_init.epigenotype[allele_i] = w;

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

// create a new offspring
void create_offspring(Individual &mother
        ,Individual &father
        ,Individual &offspring
        ,bool const offspring_envt_high
        ,double const phen_prestige_vert
        ,double const xconformist_vert
        )
{
    // set up a bernoulli distribution that returns 0s or 1s
    // at equal probability to sample alleles from the first or
    // second set of chromosomes of diploid individuals
    bernoulli_distribution allele_sample(0.5);

    assert((int)mother.g[0].size() == nloci_g);
    
    // reset arrays for the genetic cue values
    if (offspring.g[0].size() > 0)
    {
        offspring.g[0].clear();
        offspring.g[1].clear();
    }

    // inherit genetic cue values
    // first auxiliary variables
    double sum_genes = 0.0;
    double allelic_val;

    // iterate over all gene loci and inherit
    for (int g_loc_i = 0; g_loc_i < nloci_g; ++g_loc_i)
    {
        // maternal values 
        allelic_val = mutation(
                mother.g[allele_sample(rng_r)][g_loc_i],
                mu_g,
                sdmu_g
                );

        clamp(allelic_val, gmin, gmax);

        offspring.g[0].push_back(allelic_val);
        
        // paternal values 
        allelic_val = mutation(
                father.g[allele_sample(rng_r)][g_loc_i],
                mu_g,
                sdmu_g);

        clamp(allelic_val, gmin, gmax);

        offspring.g[1].push_back(allelic_val);

        sum_genes += 0.5 * (offspring.g[0][g_loc_i] + offspring.g[1][g_loc_i]);
    }
    
    assert((int)offspring.g[0].size() == nloci_g);

    // inheritance of maternal cue values 
    offspring.amat[0] = mutation(
            mother.amat[allele_sample(rng_r)],
            mu_amat,
            sdmu_a);

    clamp(offspring.amat[0], amin, amax);

    offspring.amat[1] = mutation(
            father.amat[allele_sample(rng_r)],
            mu_amat,
            sdmu_a);

    clamp(offspring.amat[1], amin, amax);

    double amat_phen = 0.5 * (offspring.amat[0] + offspring.amat[1]);
   

    // inheritance of juvenile cue values 
    offspring.ajuv[0] = mutation(
            mother.ajuv[allele_sample(rng_r)],
            mu_ajuv,
            sdmu_a);

    clamp(offspring.ajuv[0], amin, amax);

    offspring.ajuv[1] = mutation(
            father.ajuv[allele_sample(rng_r)],
            mu_ajuv,
            sdmu_a);

    clamp(offspring.ajuv[1], amin, amax);

    double ajuv_phen = 0.5 * (offspring.ajuv[0] + offspring.ajuv[1]);
   


    // inheritance of genetic cue sensitivity values 
    offspring.agen[0] = mutation(
            mother.agen[allele_sample(rng_r)],
            mu_agen,
            sdmu_a);

    clamp(offspring.agen[0], amin, amax);

    offspring.agen[1] = mutation(
            father.agen[allele_sample(rng_r)],
            mu_agen,
            sdmu_a);

    clamp(offspring.agen[1], amin, amax);

    double agen_phen = 0.5 * (offspring.agen[0] + offspring.agen[1]);
    
    
    // inheritance of vertical social cue sensitivity values 
    offspring.asoc_vert[0] = mutation(
            mother.asoc_vert[allele_sample(rng_r)],
            mu_asoc_vert,
            sdmu_a);

    clamp(offspring.asoc_vert[0], amin, amax);

    offspring.asoc_vert[1] = mutation(
            father.asoc_vert[allele_sample(rng_r)],
            mu_asoc_vert,
            sdmu_a);

    clamp(offspring.asoc_vert[1], amin, amax);

    double asoc_vert_phen = 0.5 * (offspring.asoc_vert[0] + offspring.asoc_vert[1]);
    
    
    // inheritance of horizontal social cue sensitivity values 
    offspring.asoc_horiz[0] = mutation(
            mother.asoc_horiz[allele_sample(rng_r)],
            mu_asoc_horiz,
            sdmu_a);

    clamp(offspring.asoc_horiz[0], amin, amax);

    offspring.asoc_horiz[1] = mutation(
            father.asoc_horiz[allele_sample(rng_r)],
            mu_asoc_horiz,
            sdmu_a);

    clamp(offspring.asoc_horiz[1], amin, amax);


    // inheritance of maternal phenotypic cue values 
    offspring.bmat_phen[0] = mutation(
            mother.bmat_phen[allele_sample(rng_r)],
            mu_bmat_phen,
            sdmu_b);

    clamp(offspring.bmat_phen[0], bmin, bmax);

    offspring.bmat_phen[1] = mutation(
            father.bmat_phen[allele_sample(rng_r)],
            mu_bmat_phen,
            sdmu_b);

    clamp(offspring.bmat_phen[1], bmin, bmax);

    // inheritance of maternal environmental cue values 
    offspring.bmat_envt[0] = mutation(
            mother.bmat_envt[allele_sample(rng_r)],
            mu_bmat_envt,
            sdmu_b);

    clamp(offspring.bmat_envt[0], bmin, bmax);

    offspring.bmat_envt[1] = mutation(
            father.bmat_envt[allele_sample(rng_r)],
            mu_bmat_envt,
            sdmu_b);

    clamp(offspring.bmat_envt[1], bmin, bmax);


    
    // inheritance of rel. sensitivity to 
    // vertical prestige biases
    offspring.vp[0] = mutation(
            mother.vp[allele_sample(rng_r)],
            mu_vp,
            sdmu_b);

    clamp(offspring.vp[0], bmin, bmax);

    offspring.vp[1] = mutation(
            father.vp[allele_sample(rng_r)],
            mu_vp,
            sdmu_b);

    clamp(offspring.vp[1], bmin, bmax);

    double vp_phen = 0.5 * (offspring.vp[0] + offspring.vp[1]);
    
    // inheritance of rel. sensitivity to 
    // vertical conformity biases
    offspring.vc[0] = mutation(
            mother.vc[allele_sample(rng_r)],
            mu_vc,
            sdmu_b);

    clamp(offspring.vc[0], bmin, bmax);

    offspring.vc[1] = mutation(
            father.vc[allele_sample(rng_r)],
            mu_vc,
            sdmu_b);

    clamp(offspring.vc[1], bmin, bmax);

    double vc_phen = 0.5 * (offspring.vc[0] + offspring.vc[1]);
    
    
    // inheritance of rel. sensitivity to 
    // horizontal prestige biases
    offspring.hp[0] = mutation(
            mother.hp[allele_sample(rng_r)],
            mu_hp,
            sdmu_b);

    clamp(offspring.hp[0], bmin, bmax);

    offspring.hp[1] = mutation(
            father.hp[allele_sample(rng_r)],
            mu_hp,
            sdmu_b);

    clamp(offspring.hp[1], bmin, bmax);

    // inheritance of rel. sensitivity to 
    // horizontal conformity biases
    offspring.hc[0] = mutation(
            mother.hc[allele_sample(rng_r)],
            mu_hc,
            sdmu_b);

    clamp(offspring.hc[0], bmin, bmax);

    offspring.hc[1] = mutation(
            father.hc[allele_sample(rng_r)],
            mu_hc,
            sdmu_b);

    clamp(offspring.hc[1], bmin, bmax);

    // kid receives juvenile cue
    offspring.cue_juv_envt_high = uniform(rng_r) < qjuv ? 
        offspring_envt_high : !offspring_envt_high;

    // adult cue will be received after potential envt'al change
    //
    // has the mother observed a high cue or a low one?
    double dmat_weighting = mother.cue_ad_envt_high ? -1.0 : 1.0;

    // store the maternal phenotype for stats purposes
    offspring.phen_mat = mother.phen_ad;
    offspring.maternal_cue = mother.cue_ad_envt_high;

    // express sensitivity to maternal phenotype
    double b_phen = 0.5 * (offspring.bmat_phen[0] + offspring.bmat_phen[1]);
    assert(b_phen >= bmin);
    assert(b_phen <= bmax);

    // express sensitivity to maternal environment 
    double b_envt = 0.5 * (offspring.bmat_envt[0] + offspring.bmat_envt[1]);
    assert(b_envt >= bmin);
    assert(b_phen <= bmax);

    // generate maternal cue
    double xoff = 1.0 /
        (1.0 + exp(
                   -b_phen * (mother.phen_ad - 0.5) 
                   + 
                   dmat_weighting * b_envt));

    // generate maternal cue again 
    // but then holding the maternal environment constant
    // this is for stats purposes
    double xoff_phen_only = 1.0 /
        (1.0 + exp(
                   -b_phen * (mother.phen_ad - 0.5)));

    // generate maternal cue again 
    // but then holding the maternal phenotype constant
    // this is for stats purposes
    double xoff_envt_only = 1.0 /
        (1.0 + exp(dmat_weighting * b_envt));

    // noise in the maternal cue
    normal_distribution<> maternal_noise(0.0, sdmat);

    double mnoise = maternal_noise(rng_r);

    // store the noise deviate for stats purposes
    offspring.mnoise = mnoise;

    // calculate final value of maternal cue
    offspring.xmat = xoff + mnoise;
    offspring.xmat_phen_only = xoff_phen_only + mnoise;
    offspring.xmat_envt_only = xoff_envt_only + mnoise;

    clamp(offspring.xmat, 0.0, 1.0);
    clamp(offspring.xmat_phen_only, 0.0, 1.0);
    clamp(offspring.xmat_envt_only, 0.0, 1.0);

    // social learning
    offspring.xconformist_vert = xconformist_vert;
    offspring.phen_prestige_vert = phen_prestige_vert;


    // generate vertical socially learnt cue
    double xsoc_vert = 1.0 / (1.0 + exp(
                - vp_phen * (offspring.phen_prestige_vert - 0.5)
                - vc_phen * offspring.xconformist_vert));

    // also calculate vertical social cues for prestige only
    double xsoc_vert_c = 1.0 / (1.0 + exp(
                - vc_phen * offspring.xconformist_vert));

    double xsoc_vert_p = 1.0 / (1.0 + exp(
                - vp_phen * (offspring.phen_prestige_vert - 0.5)));

    normal_distribution<> social_noise(0.0, sdsoc_vert);

    double socnoise = social_noise(rng_r);

    // store the noise deviate for stats
    offspring.svnoise = socnoise;

    offspring.xsoc_vert = xsoc_vert + socnoise;
    offspring.xsoc_vert_c = xsoc_vert_c + socnoise;
    offspring.xsoc_vert_p = xsoc_vert_p + socnoise;

    clamp(offspring.xsoc_vert, 0.0, 1.0);
    clamp(offspring.xsoc_vert_c, 0.0, 1.0);
    clamp(offspring.xsoc_vert_p, 0.0, 1.0);

    // expressing a juvenile phenotype
    offspring.phen_juv = 1.0 / 
        (1.0 + exp(
                   -amat_phen * offspring.xmat +
                   -agen_phen * sum_genes +
                   -ajuv_phen * offspring.cue_juv_envt_high
                   -asoc_vert_phen * offspring.xsoc_vert
                   ));
    // 
    offspring.phen_ad = NAN;
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
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < NBreeders; ++breeder_i)
        {
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

        }
    }// for (int patch_i = 0

}

 

// births of new offspring
// juvenile cue integration
// juvenile selection
// adult cue integration (aka horizontal social learning)
void replace()
{
    // randomly chosen remote patch to obtain
    // individuals from
    int random_remote_patch;

    // auxiliary variables to calculate
    // socially learnt cues through prestige-based social learning
    // and through conformism-based social learning
    double prestige_phen, xconformist;

    // set up a random number generator to 
    // sample remote patches
    uniform_int_distribution<> patch_sampler(
            0,
            NPatches - 1);


    // count juveniles that have survived
    // and hence are recruited as local breeder
    int n_recruited_breeders;

    // auxiliary variable to keep track 
    // of the survival probability
    double surv;

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        n_recruited_breeders = 0;

        // now make NBreeder offspring and have them survive
        for (int offspring_i = 0; offspring_i < NBreeder; ++offspring_i)
        {
            Individual Kid;

            // offspring is born in local patch hence sample from local parents
            if (uniform(rng_r) < 1.0 - m 
                    &&
                    Pop[patch_i].n_breeders > 0)
            {
                // set up a random number generator to 
                // sample from the remaining breeders
                uniform_int_distribution<> random_local_breeder(
                        0,
                        Pop[patch_i].n_breeders - 1);

                // prepare cues in local patch from social learning
                social_learning(
                        patch_i
                        ,false
                        ,prestige_phen
                        ,xconformist);

                create_offspring(
                        Pop[patch_i].breeders[random_local_breeder(rng_r)]
                        ,Pop[patch_i].breeders[random_local_breeder(rng_r)]
                        ,Kid
                        ,Pop[patch_i].envt_high
                        ,prestige_phen
                        ,xconformist
                );
            }
            else // offspring born in remote patch
            {
                // sample a random remote pathc
                do {

                    random_remote_patch = patch_sampler(rng_r);

                }
                while(Pop[random_remote_patch].n_breeders < 1);
        
                uniform_int_distribution<> random_remote_breeder(
                0,
                Pop[random_remote_patch].n_breeders - 1);
               
                // prepare cues in remote patch (where parents are breeding)
                // from social learning
                social_learning(
                        random_remote_patch
                        ,false
                        ,prestige_phen
                        ,xconformist);

                create_offspring(
                        Pop[random_remote_patch].breeders[random_remote_breeder(rng_r)]
                        ,Pop[random_remote_patch].breeders[random_remote_breeder(rng_r)]
                        ,Kid
                        ,Pop[random_remote_patch].envt_high
                        ,prestige_phen
                        ,xconformist
                );
            }

            // in case juvenile selection acts,
            // and offspring dies, just continue on the next
            // iteration of this loop without adding the offspring
            // to the breeder population
            if (juvenile_survival)
            {
                surv = survival_probability(
                            Kid.phen_juv
                            ,Pop[patch_i].envt_high);

                cout << "sodeju" << endl;

                // offspring dies 
                if (uniform(rng_r) > surv)
                {
                    continue;
                }
            }

            // offspring survives
            Pop[patch_i].breeders_t1[n_recruited_breeders] = Kid;
            assert((int)Pop[patch_i].breeders_t1[n_recruited_breeders].g[0].size() == nloci_g);

            ++n_recruited_breeders;

            assert(n_recruited_breeders <= NBreeder);

        } // for (int offspring_i = 0; offspring_i < NBreeder; ++offspring_i)
    } // end for (int patch_i = 0

    bool cue_ad_envt_high;

    // auxiliary variables for horizontal social learning
    double hp, hc, asoc_horiz, xsoc_horiz, noise;

    // random number for errors in horizontal social learning
    normal_distribution<> social_noise(0.0, sdsoc_horiz);

    // all new breeders born etc, copy them over
    // and perform horiz social learning
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // calculate adult cue value supplied to mothers
        // the cue value is the same for all mothers
        cue_ad_envt_high = uniform(rng_r) < qmat ? 
            Pop[patch_i].envt_high 
            : 
            !Pop[patch_i].envt_high;
        
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            Pop[patch_i].breeders[breeder_i] = 
                Pop[patch_i].breeders_t1[breeder_i];
    
            assert((int)Pop[patch_i].breeders[breeder_i].g[0].size() == nloci_g);
            
            Pop[patch_i].n_breeders = NBreeder;
        
            // give breeder an environmental cue as adult
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_ad_envt_high;
        } // all juveniles are now assigned a breeding position

        // horizontal social learning 
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            // ad_phen should be NaN as it is yet to be set
            assert(abs(::isnan(Pop[patch_i].breeders[breeder_i].phen_ad)) > 0);

            social_learning(
                    patch_i
                    ,true
                    ,prestige_phen
                    ,xconformist);

            hp = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hp[0]
                    + 
                    Pop[patch_i].breeders[breeder_i].hp[1]);
            
            hc = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hc[0]
                    + 
                    Pop[patch_i].breeders[breeder_i].hc[1]);

            asoc_horiz = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[0]
                    + 
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[1]);

            Pop[patch_i].breeders[breeder_i].xconformist_horiz = xconformist;
            Pop[patch_i].breeders[breeder_i].phen_prestige_horiz = prestige_phen;
            noise = social_noise(rng_r);

            // generate the horizontally learnt cue and add error
            xsoc_horiz = 1.0 /
                (1.0 + exp(
                           -hp * (prestige_phen - 0.5)
                           -hc * xconformist)
                ) + noise;

            // store noise for stats purposes
            Pop[patch_i].breeders[breeder_i].shnoise = noise;

            clamp(xsoc_horiz, 0.0, 1.0);

            // assign the horizontal social learnt 
            // conformism bias cue to the breeding individual
            Pop[patch_i].breeders[breeder_i].xsoc_horiz = 
                xsoc_horiz;

            // add social learning to the adult phenotype
            Pop[patch_i].breeders[breeder_i].phen_ad = 1.0 / 
                (1.0 + exp(
                            log(1.0 / Pop[patch_i].breeders[breeder_i].phen_juv - 1.0) 
                            -asoc_horiz * xsoc_horiz));
        }

        // envtal change after breeder establishment
        if (uniform(rng_r) < 1.0 - p)
        {
            Pop[patch_i].envt_high = !Pop[patch_i].envt_high;
        }

    } // end for (int patch_i = 0
}

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
