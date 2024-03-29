import("stdfaust.lib");
process = gate : envARE with { // * no.noise : fr <: _,_ with {
    a_ = 0.005;
    len = 8487;
    r_ = (len - (a_*ma.SR))/ma.SR;//0.1718;
    s_ = 0.1;
    f_ = 0.05;
    gate = (1-(1@(a_*ma.SR)));// + 0.5*(1@750-(1@1700));
    envASRFE = en.asrfe(a_,s_,r_,f_);
    envARFE = en.arfe(a_,r_,f_);
    envARE = en.are(a_,r_);
    fr = fi.iir(b,a) with {
        b = 1.0038,-0.16283,0.0062466,-0.10801,-0.24058, 0.029842,-0.121,-0.16796,-0.15775,-0.20561,0.0077204;
        a = -1.3267,0.61699,-0.75244,0.5751,-0.2797,0.497,-0.45368,0.3945,-0.22875,0.0441;
    };
};
