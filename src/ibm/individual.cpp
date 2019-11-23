#include "individual.hpp"

Individual::Individual():
    phen_juv{0.0},
    phen_ad{0.0},
    phen_mat{0.0},
    phen_prestige_horiz{0.0},
    phen_prestige_vert{0.0},
    xmat{0.0},
    xmat_envt_only{0.0},
    xmat_phen_only{0.0},
    xsoc_horiz{0.0},
    xsoc_horiz_c{0.0},
    xsoc_horiz_p{0.0},
    xsoc_vert{0.0},
    xsoc_vert_c{0.0},
    xsoc_vert_p{0.0},
    xconformist_horiz{0.0},
    xconformist_vert{0.0},
    agen{0.0,0.0},
    amat{0.0,0.0},
    ajuv{0.0,0.0},
    asoc_vert{0.0,0.0},
    asoc_horiz{0.0,0.0},
    bmat_envt{0.0,0.0},
    bmat_phen{0.0,0.0},
    vp{0.0,0.0},
    vc{0.0,0.0},
    hp{0.0,0.0},
    hc{0.0,0.0},
    cue_ad_envt_high{false},
    cue_juv_envt_high{false},
    mnoise{0.0},
    svnoise{0.0},
    shnoise{0.0}
{
}

Individual::Individual(Individual const &other):
    phen_juv{other.phen_juv},
    phen_ad{other.phen_ad},
    phen_mat{other.phen_mat},
    phen_prestige_horiz{other.phen_prestige_horiz},
    phen_prestige_vert{other.phen_prestige_vert},
    g{other.g[0],other.g[1]},
    xmat{other.xmat},
    xmat_envt_only{other.xmat_envt_only},
    xmat_phen_only{other.xmat_phen_only},
    xsoc_horiz{other.xsoc_horiz},
    xsoc_horiz_c{other.xsoc_horiz_c},
    xsoc_horiz_p{other.xsoc_horiz_p},
    xsoc_vert{other.xsoc_vert},
    xsoc_vert_c{other.xsoc_vert_c},
    xsoc_vert_p{other.xsoc_vert_p},
    xconformist_horiz{other.xconformist_horiz},
    xconformist_vert{other.xconformist_vert},
    agen{other.agen[0],other.agen[1]},
    amat{other.amat[0],other.amat[1]},
    ajuv{other.ajuv[0],other.ajuv[1]},
    asoc_vert{other.asoc_vert[0],other.asoc_vert[1]},
    asoc_horiz{other.asoc_horiz[0],other.asoc_horiz[1]},
    bmat_envt{other.bmat_envt[0],other.bmat_envt[1]},
    bmat_phen{other.bmat_phen[0],other.bmat_phen[1]},
    vp{other.vp[0],other.vp[1]},
    vc{other.vc[0],other.vc[1]},
    hp{other.hp[0],other.hp[1]},
    hc{other.hc[0],other.hc[1]},
    cue_ad_envt_high{other.cue_ad_envt_high},
    cue_juv_envt_high{other.cue_juv_envt_high},
    mnoise{other.mnoise},
    svnoise{other.svnoise},
    shnoise{other.shnoise}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    phen_juv = other.phen_juv;
    phen_ad = other.phen_ad;
    phen_mat = other.phen_mat;
    phen_prestige_horiz = other.phen_prestige_horiz;
    phen_prestige_vert = other.phen_prestige_vert;
    xmat = other.xmat;
    xmat_envt_only = other.xmat_envt_only;
    xmat_phen_only = other.xmat_phen_only;
    xconformist_horiz = other.xconformist_horiz;
    xconformist_vert = other.xconformist_vert;
    xsoc_horiz = other.xsoc_horiz;
    xsoc_horiz_c = other.xsoc_horiz_c;
    xsoc_horiz_p = other.xsoc_horiz_p;
    xsoc_vert = other.xsoc_vert;
    xsoc_vert_c = other.xsoc_vert_c;
    xsoc_vert_p = other.xsoc_vert_p;

    agen[0] = other.agen[0];
    agen[1] = other.agen[1];

    amat[0] = other.amat[0];
    amat[1] = other.amat[1];

    ajuv[0] = other.ajuv[0];
    ajuv[1] = other.ajuv[1];
    
    asoc_vert[0] = other.asoc_vert[0];
    asoc_vert[1] = other.asoc_vert[1];
    
    asoc_horiz[0] = other.asoc_horiz[0];
    asoc_horiz[1] = other.asoc_horiz[1];
    
    bmat_phen[0] = other.bmat_phen[0];
    bmat_phen[1] = other.bmat_phen[1];
    
    bmat_envt[0] = other.bmat_envt[0];
    bmat_envt[1] = other.bmat_envt[1];

    hp[0] = other.hp[0];
    hp[1] = other.hp[1];
    
    hc[0] = other.hc[0];
    hc[1] = other.hc[1];
    
    vp[0] = other.vp[0];
    vp[1] = other.vp[1];
    
    vc[0] = other.vc[0];
    vc[1] = other.vc[1];

    cue_ad_envt_high = other.cue_ad_envt_high;
    cue_juv_envt_high = other.cue_juv_envt_high;

    g[0] = other.g[0];
    g[1] = other.g[1];

    mnoise = other.mnoise;
    svnoise = other.svnoise;
    shnoise = other.shnoise;
} // end void Individual::operator=()
