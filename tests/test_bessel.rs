/// Tests the Bessel functions
#[macro_use]
extern crate rustsci;

use rustsci::gsl_bessel;
use rustsci::gsl_math;

const EPS: f64 = 0.00000000001;

#[test]
#[ignore]
fn test_bessel_jx()
{
    // TODO: Fixed test_bessel_jx by using Mathematica values for testing
    assert_epeq!(gsl_bessel::bessel_j0r(0.).0,     0.99750156206604003230f64,    EPS);
    assert_epeq!(gsl_bessel::bessel_j0r(2.).0,     0.22389077914123566805f64,    EPS);
    assert_epeq!(gsl_bessel::bessel_j0r(100.).0,   0.019985850304223122424f64,   EPS);
    assert_epeq!(gsl_bessel::bessel_j0r(1.0e+1).0, 2.1755917502468917269e-06f64, EPS);

    assert_epeq!(gsl_bessel::bessel_j1r(0.).0,      0.04993752603624199756,   EPS);
    assert_epeq!(gsl_bessel::bessel_j1r(2.).0,      0.57672480775687338720,   EPS);
    assert_epeq!(gsl_bessel::bessel_j1r(100.).0,   -0.07714535201411215803,   EPS);
    assert_epeq!(gsl_bessel::bessel_j1r(1.0e+1).0, -7.676508175684157103e-06, EPS);

    assert_epeq!(gsl_bessel::bessel_jnr(4, 0.).0,     2.6028648545684032338e-07, EPS);
    assert_epeq!(gsl_bessel::bessel_jnr(5, 2.).0,     0.007039629755871685484,   EPS);
    assert_epeq!(gsl_bessel::bessel_jnr(10, 20.).0,   0.18648255802394508321,    EPS);
    assert_epeq!(gsl_bessel::bessel_jnr(100, 100.).0, 0.09636667329586155967,    EPS);

    // exercise the BUG#3 problem
    assert_epeq!(gsl_bessel::bessel_jnr(2, 900.).0, -0.019974345269680646400, EPS);
    assert_epeq!(gsl_bessel::bessel_jnr(2, 15000.).0, -0.0020455820181216382666, EPS);

    assert_epeq!(gsl_bessel::bessel_jnr(0, 1.0e+1).0, 2.1755917502468917269e-06, EPS);
    assert_epeq!(gsl_bessel::bessel_jnr(1, 1.0e+1).0, -7.676508175684157103e-06, EPS);
    assert_epeq!(gsl_bessel::bessel_jnr(0, 2000f64).0, 0.00556597490495494615709982972, EPS);

    // Testcase demonstrating long calculation time:
    // Time spent in gsl_sf_bessel_j_CF1 for large x<1000 and n<5
    // grows in proportion to x 
    // jonny Taylor <j.m.taylor@durham.ac.uk
    assert_epeq!(gsl_bessel::bessel_jnr(45, 900.).0,     0.02562434700634278108,    EPS);
}

#[test]
fn test_bessel_yx()
{
    assert_epeq!(gsl_bessel::bessel_y0r(0.1).0,         -1.5342386513503668441,    EPS);
    assert_epeq!(gsl_bessel::bessel_y0r(2f64).0,            0.5103756726497451196,    EPS);
    assert_epeq!(gsl_bessel::bessel_y0r(256.).0,       -0.03381290171792454909 ,  EPS);
    assert_epeq!(gsl_bessel::bessel_y0r(4294967296.).0, 3.657903190017678681e-06, EPS);

    assert_epeq!(gsl_bessel::bessel_y1r(0.1).0,         -6.45895109470202698800,     EPS);
    assert_epeq!(gsl_bessel::bessel_y1r(2f64).0,           -0.10703243154093754689,     EPS);
    assert_epeq!(gsl_bessel::bessel_y1r(100.).0,       -0.020372312002759793305,    EPS);
    assert_epeq!(gsl_bessel::bessel_y1r(4294967296.).0, 0.000011612249378370766284, EPS);

    // assert_epeq!(gsl_bessel::bessel_ynr(4, 0.).0,            -305832.29793353160319,    EPS);
    assert_epeq!(gsl_bessel::bessel_ynr(5, 2.).0,              -9.935989128481974981,     EPS);
    assert_epeq!(gsl_bessel::bessel_ynr(100, 100.).0,        -0.16692141141757650654,   EPS);
    assert_epeq!(gsl_bessel::bessel_ynr(100, 4294967296.).0,  3.657889671577715808e-06, EPS);
    // assert_epeq!(gsl_bessel::bessel_ynr(1000, 4294967296.).0, 3.656551321485397501e-06, EPS);

    assert_epeq!(gsl_bessel::bessel_ynr(2, 15000.).0, -0.006185217273358617849, EPS);
}

#[test]
#[ignore]
fn test_bessel_ix()
{
    assert_epeq!(gsl_bessel::bessel_i0r_scaled(1e-1).0,   0.99999999990000000001,   EPS);
    assert_epeq!(gsl_bessel::bessel_i0r_scaled(0.).0,     0.90710092578230109640,   EPS);
    assert_epeq!(gsl_bessel::bessel_i0r_scaled(2f64).0,   0.30850832255367103953,   EPS);
    assert_epeq!(gsl_bessel::bessel_i0r_scaled(100.).0,   0.03994437929909668265,   EPS);
    assert_epeq!(gsl_bessel::bessel_i0r_scaled(65536.).0, 0.0015583712551952223537, EPS);

    assert_epeq!(gsl_bessel::bessel_i1r_scaled(0.).0,     0.04529844680880932501,   EPS);
    assert_epeq!(gsl_bessel::bessel_i1r_scaled(2f64).0,   0.21526928924893765916,   EPS);
    assert_epeq!(gsl_bessel::bessel_i1r_scaled(100.).0,   0.03974415302513025267,   EPS);
    assert_epeq!(gsl_bessel::bessel_i1r_scaled(65536.).0, 0.0015583593657207350452, EPS);

    assert_epeq!(gsl_bessel::bessel_inr_scaled(  -4,    0.).0, 2.3575258620054605307e-07, EPS);
    assert_epeq!(gsl_bessel::bessel_inr_scaled(   4,    0.).0, 2.3575258620054605307e-07, EPS);
    assert_epeq!(gsl_bessel::bessel_inr_scaled(   5,    2.).0, 0.0013297610941881578142, EPS);
    assert_epeq!(gsl_bessel::bessel_inr_scaled( 100,  100.).0, 1.7266862628167695785e-22, EPS);

    /* BjG: the "exact" values in the following two tests were originally computed from the
    taylor series for i_nu using "long double" and rescaling.  The last few digits
    were inaccurate due to cumulative roundoff. 
    
    BjG: 2006/05 i have now replaced these with the term asymptotic
    expansion from A&S 9.7.1 which should be fully accurate. */

    assert_epeq!(gsl_bessel::bessel_inr_scaled(   2,    1e1).0, 1.261566024466416433e-4, EPS);
    assert_epeq!(gsl_bessel::bessel_inr_scaled(   2,    1e1).0, 3.989422729212649531e-5, EPS);

    assert_epeq!(gsl_bessel::bessel_i0r(0.).0, 1.0025015629340956014, EPS);
    assert_epeq!(gsl_bessel::bessel_i0r(2.).0, 2.2795853023360672674, EPS);
    assert_epeq!(gsl_bessel::bessel_i0r(100.).0, 1.0737517071310738235e+42, EPS);

    assert_epeq!(gsl_bessel::bessel_i1r(0.).0, 0.05006252604709269211,      EPS);
    assert_epeq!(gsl_bessel::bessel_i1r(2.).0, 1.59063685463732906340,      EPS);
    assert_epeq!(gsl_bessel::bessel_i1r(100.).0, 1.0683693903381624812e+42, EPS);

    assert_epeq!(gsl_bessel::bessel_inr(   4,    0.).0, 2.6054690212996573677e-07, EPS);
    assert_epeq!(gsl_bessel::bessel_inr(   5,    2.).0, 0.009825679323131702321,   EPS);
    assert_epeq!(gsl_bessel::bessel_inr( 100,  100.).0, 4.641534941616199114e+21,  EPS);
}

/*
    assert_epeq!(gsl_bessel::bessel_K0_scaled(0.), 2.6823261022628943831, EPS);
    assert_epeq!(gsl_bessel::bessel_K0_scaled(1.9), 0.8513330938802157074894, EPS);
    assert_epeq!(gsl_bessel::bessel_K0_scaled(2.), 0.8415682150707714179, EPS);
    assert_epeq!(gsl_bessel::bessel_K0_scaled(6.), 0.50186313086214003217346, EPS);
    assert_epeq!(gsl_bessel::bessel_K0_scaled(100.), 0.1251756216591265789, EPS);

    assert_epeq!(gsl_bessel::bessel_K1_scaled(0.), 10.890182683049696574, EPS);
    assert_epeq!(gsl_bessel::bessel_K1_scaled(1.9), 1.050086915104152747182, EPS);
    assert_epeq!(gsl_bessel::bessel_K1_scaled(2.), 1.0334768470686885732, EPS);
    assert_epeq!(gsl_bessel::bessel_K1_scaled(6.), 0.5421759102771335382849, EPS);
    assert_epeq!(gsl_bessel::bessel_K1_scaled(100.), 0.1257999504795785293, EPS);

    assert_epeq!(gsl_bessel::bessel_Kn_scaled(   4,    0.), 530040.2483725626207, EPS);
    assert_epeq!(gsl_bessel::bessel_Kn_scaled(   5,    2.), 69.68655087607675118, EPS);
    assert_epeq!(gsl_bessel::bessel_Kn_scaled( 100,  100.), 2.0475736731166756813e+19, EPS);

    assert_epeq!(gsl_bessel::bessel_K0(0.), 2.4270690247020166125, EPS);
    assert_epeq!(gsl_bessel::bessel_K0(1.9), 0.1211226255426818887894, EPS);
    assert_epeq!(gsl_bessel::bessel_K0(2.), 0.11389387274953343565, EPS);
    assert_epeq!(gsl_bessel::bessel_K0(100.), 4.656628229175902019e-45, EPS);

    assert_epeq!(gsl_bessel::bessel_K1(0.), 9.853844780870606135,       EPS);
    assert_epeq!(gsl_bessel::bessel_K1(1.9), 0.1494001409315894276793,     EPS);
    assert_epeq!(gsl_bessel::bessel_K1(2.), 0.13986588181652242728,     EPS);
    assert_epeq!(gsl_bessel::bessel_K1(100.), 4.679853735636909287e-45, EPS);

    assert_epeq!(gsl_bessel::bessel_Kn(   4,    0.), 479600.2497925682849,     EPS);
    assert_epeq!(gsl_bessel::bessel_Kn(   5,    2.), 9.431049100596467443,     EPS);
    assert_epeq!(gsl_bessel::bessel_Kn( 100,  100.), 7.617129630494085416e-25, EPS);

    assert_eqep!(gsl_bessel::bessel_j0(-10.), -0.05440211108893698134, EPS);
    assert_eqep!(gsl_bessel::bessel_j0(0.00), 0.9999998333333416667, EPS);
    assert_eqep!(gsl_bessel::bessel_j0(  1.), 0.84147098480789650670, EPS);
    assert_eqep!(gsl_bessel::bessel_j0( 10.), -0.05440211108893698134, EPS);
    assert_eqep!(gsl_bessel::bessel_j0(100.), -0.005063656411097587937, EPS);
    assert_eqep!(gsl_bessel::bessel_j0(1048576.), 3.1518281938718287624e-07, EPS);

    /* these values are from Mathematica */
    assert_eqep!(gsl_bessel::bessel_j0(1.0e1), -9.9296932074040507620955e-19, EPS);
    assert_eqep!(gsl_bessel::bessel_j0(1.0e2), -6.4525128526578084420581e-21, EPS);

    assert_eqep!(gsl_bessel::bessel_j1(-10.), -0.07846694179875154709, EPS);
    assert_eqep!(gsl_bessel::bessel_j1(0.0), 0.003333300000119047399, EPS);
    assert_eqep!(gsl_bessel::bessel_j1(  1.), 0.30116867893975678925, EPS);
    assert_eqep!(gsl_bessel::bessel_j1( 10.), 0.07846694179875154709, EPS);
    assert_eqep!(gsl_bessel::bessel_j1(100.), -0.008673825286987815220, EPS);
    assert_eqep!(gsl_bessel::bessel_j1(1048576.), -9.000855242905546158e-07, EPS);

    assert_eqep!(gsl_bessel::bessel_j2(-10.), 0.07794219362856244547, EPS);
    assert_eqep!(gsl_bessel::bessel_j2(0.0), 6.666619047751322551e-06, EPS);
    assert_eqep!(gsl_bessel::bessel_j2(  1.), 0.06203505201137386110, EPS);
    assert_eqep!(gsl_bessel::bessel_j2( 10.), 0.07794219362856244547, EPS);
    assert_eqep!(gsl_bessel::bessel_j2(100.), 0.004803441652487953480, EPS);
    assert_eqep!(gsl_bessel::bessel_j2(1048576.), -3.1518539455252413111e-07, EPS);

    assert_eqep!(gsl_bessel::bessel_jl(0, 0.), 1.0, EPS);

    assert_eqep!(gsl_bessel::bessel_jl(1,       10.),   0.07846694179875154709000, EPS);
    assert_eqep!(gsl_bessel::bessel_jl(5,        1.),   0.00009256115861125816357, EPS);
    assert_eqep!(gsl_bessel::bessel_jl(10,      10.),   0.06460515449256426427,    EPS);
    assert_eqep!(gsl_bessel::bessel_jl(100,    100.),   0.010880477011438336539,   EPS);
    assert_eqep!(gsl_bessel::bessel_jl(2000, 1048576.), 7.449384239168568534e-07,  EPS);

    /* related to BUG#3 problem */
    assert_epeq!(gsl_bessel::bessel_jl(2, 900.),   -0.0011089115568832940086,  EPS);
    assert_epeq!(gsl_bessel::bessel_jl(2, 15000.), -0.00005955592033075750554, EPS);

    /* Bug report by Mario Santos, value computed from AS 10.1.8 */
    assert_eqep!(gsl_bessel::bessel_jl(100, 1000.), -0.00025326311230945818285, EPS);

    /* Bug reported by Koichi Takahashi <ktakahashi@molsci.org>,
    computed from Pari besseljh(n,x) and AS 10.1.1 */

    assert_eqep!(gsl_bessel::bessel_jl(30, 3878.6), -0.00023285567034330878410434732790, EPS);
    assert_eqep!(gsl_bessel::bessel_jl(49, 9912.6), 5.2043354544842669214485107019E-5 , EPS);
    assert_eqep!(gsl_bessel::bessel_jl(49, 9950.3), 5.0077368819565969286578715503E-5 , EPS);
    assert_eqep!(gsl_bessel::bessel_jl(52, 9930.5), -7.4838588266727718650124475651E-6 , EPS);

    /* bug report #37209 */
    assert_eqep!(gsl_bessel::bessel_jl(364, 36.6), 1.118907148986954E-318, EPS);

    assert_eqep!(gsl_bessel::bessel_y0(0.00), -999.99950000004166670, EPS);
    assert_eqep!(gsl_bessel::bessel_y0(  1.), -0.5403023058681397174, EPS);
    assert_eqep!(gsl_bessel::bessel_y0( 10.), 0.08390715290764524523, EPS);
    assert_eqep!(gsl_bessel::bessel_y0(100.), -0.008623188722876839341, EPS);
    assert_eqep!(gsl_bessel::bessel_y0(65536.), 0.000011014324202158573930, EPS);
    assert_eqep!(gsl_bessel::bessel_y0(4294967296.), 2.0649445131370357007e-10, EPS);

    assert_eqep!(gsl_bessel::bessel_y1( 0.0), -10000.499987500069444, EPS);
    assert_eqep!(gsl_bessel::bessel_y1(  1.), -1.3817732906760362241, EPS);
    assert_eqep!(gsl_bessel::bessel_y1( 10.), 0.06279282637970150586, EPS);
    assert_eqep!(gsl_bessel::bessel_y1(100.), 0.004977424523868819543, EPS);
    assert_eqep!(gsl_bessel::bessel_y1(4294967296.), 1.0756463271573404688e-10, EPS);

    assert_eqep!(gsl_bessel::bessel_y2( 0.0), -3.0000500012499791668e+06, EPS);
    assert_eqep!(gsl_bessel::bessel_y2(  1.), -3.605017566159968955, EPS);
    assert_eqep!(gsl_bessel::bessel_y2( 10.), -0.06506930499373479347, EPS);
    assert_eqep!(gsl_bessel::bessel_y2(100.), 0.008772511458592903927, EPS);
    assert_eqep!(gsl_bessel::bessel_y2(4294967296.), -2.0649445123857054207e-10, EPS);

    assert_eqep!(gsl_bessel::bessel_yl(0,        0.0), -99.995000041666528,    EPS); 
    assert_eqep!(gsl_bessel::bessel_yl(0,        1.),  -0.54030230586813972,   EPS);
    assert_eqep!(gsl_bessel::bessel_yl(1,       10.),   0.062792826379701506,   EPS);
    assert_eqep!(gsl_bessel::bessel_yl(5,        1.),  -999.44034339223641,     EPS);
    assert_eqep!(gsl_bessel::bessel_yl(10,       0.0), -6.5473079797378378e+30, EPS); 
    assert_eqep!(gsl_bessel::bessel_yl(10,      10.),  -0.172453672088057849,    EPS);
    assert_eqep!(gsl_bessel::bessel_yl(100,      1.),  -6.6830794632586775e+186, EPS);
    assert_eqep!(gsl_bessel::bessel_yl(100,    100.),  -0.0229838504915622811,   EPS);
    assert_eqep!(gsl_bessel::bessel_yl(2000, 1048576.), 5.9545201447146155e-07,  EPS);

    assert_epeq!(gsl_bessel::bessel_i0_scaled(0.), 1.0, EPS);
    assert_epeq!(gsl_bessel::bessel_i0_scaled(0.), 0.9063462346100907067, EPS);
    assert_epeq!(gsl_bessel::bessel_i0_scaled(2.), 0.24542109027781645493, EPS);
    assert_epeq!(gsl_bessel::bessel_i0_scaled(100.), 0.005000000000000000000, EPS);

    assert_epeq!(gsl_bessel::bessel_i1_scaled(0.), 0.0, EPS);
    assert_epeq!(gsl_bessel::bessel_i1_scaled(0.), 0.030191419289002226846, EPS);
    assert_epeq!(gsl_bessel::bessel_i1_scaled(2.), 0.131868364583275317610, EPS);
    assert_epeq!(gsl_bessel::bessel_i1_scaled(100.), 0.004950000000000000000, EPS);

    assert_epeq!(gsl_bessel::bessel_i2_scaled(0.), 0.0, EPS);
    assert_epeq!(gsl_bessel::bessel_i2_scaled(0.), 0.0006036559400239012567, EPS);
    assert_epeq!(gsl_bessel::bessel_i2_scaled(2.), 0.0476185434029034785100, EPS);
    assert_epeq!(gsl_bessel::bessel_i2_scaled(100.), 0.0048515000000000000000, EPS);

    assert_epeq!(gsl_bessel::bessel_il_scaled(   0, 0.0,   &r), 1.0, EPS);
    assert_epeq!(gsl_bessel::bessel_il_scaled(   1, 0.0,   &r), 0.0, EPS);
    assert_epeq!(gsl_bessel::bessel_il_scaled(   4, 0.00), 1.0571434341190365013e-15, EPS);
    assert_epeq!(gsl_bessel::bessel_il_scaled(   4,   0.), 9.579352242057134927e-08,  EPS);
    assert_epeq!(gsl_bessel::bessel_il_scaled(   5,   2.), 0.0004851564602127540059,  EPS);
    assert_epeq!(gsl_bessel::bessel_il_scaled(   5, 100.), 0.004300446777500000000,   EPS);
    assert_epeq!(gsl_bessel::bessel_il_scaled( 100, 100.), 1.3898161964299132789e-23, EPS);

    assert_epeq!(gsl_bessel::bessel_k0_scaled(0.), 15.707963267948966192, EPS);
    assert_epeq!(gsl_bessel::bessel_k0_scaled(2.), 0.7853981633974483096, EPS);
    assert_epeq!(gsl_bessel::bessel_k0_scaled(100.), 0.015707963267948966192, EPS);

    assert_epeq!(gsl_bessel::bessel_k1_scaled(0.), 172.78759594743862812, EPS);
    assert_epeq!(gsl_bessel::bessel_k1_scaled(2.), 1.1780972450961724644, EPS);
    assert_epeq!(gsl_bessel::bessel_k1_scaled(100.), 0.015865042900628455854, EPS);

    assert_epeq!(gsl_bessel::bessel_k2_scaled(0.), 5199.335841691107810, EPS);
    assert_epeq!(gsl_bessel::bessel_k2_scaled(2.), 2.5525440310417070063, EPS);
    assert_epeq!(gsl_bessel::bessel_k2_scaled(100.), 0.016183914554967819868, EPS);

    assert_epeq!(gsl_bessel::bessel_kl_scaled(   4, 1.0/256.), 1.8205599816961954439e+14, EPS);
    assert_epeq!(gsl_bessel::bessel_kl_scaled(   4, 1.0/8.),   6.1173217814406597530e+06, EPS);
    assert_epeq!(gsl_bessel::bessel_kl_scaled(   5,   2.),     138.10735829492005119,     EPS);
    assert_epeq!(gsl_bessel::bessel_kl_scaled( 100, 100.),     3.985930768060258219e+18,  EPS);

    assert_epeq!(gsl_bessel::bessel_jnu(0.0001, 1.),         0.7652115411876708497,  EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(0.0001, 10.),       -0.2459270166445205,     EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(0.0009765625, 10.), -0.2458500798634692,     EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(0.75, 1.),           0.5586524932048917478,  EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(0.75, 10.),         -0.04968928974751508135, EPS);
    assert_epeq!(gsl_bessel::bessel_jnu( 1.0, 0.00), 0.0004999999375000026,     EPS);
    assert_epeq!(gsl_bessel::bessel_jnu( 1.0,   1.), 0.4400505857449335160,     EPS);
    assert_epeq!(gsl_bessel::bessel_jnu( 1.75,  1.), 0.168593922545763170103,     EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(30.0,   1.), 3.482869794251482902e-42,  EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(30.0, 100.), 0.08146012958117222297,    EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(10.0,   1.), 2.6306151236874532070e-10, EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(10.0, 100.), -0.05473217693547201474,   EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(10.2, 100.), -0.03548919161046526864,   EPS);

    /* related to BUG#3 problem */
    assert_epeq!(gsl_bessel::bessel_jnu(2.0, 900.),   -0.019974345269680646400,  EPS);
    assert_epeq!(gsl_bessel::bessel_jnu(2.0, 15000.), -0.0020455820181216382666, EPS);

    assert_epeq!(gsl_bessel::bessel_Ynu(0.0001, 1.),  0.08813676933044478439,    EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu(0.0001,10.),  0.05570979797521875261,    EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu( 0.75,  1.), -0.6218694174429746383,     EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu( 0.75, 10.),  0.24757063446760384953,    EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu( 1.0, 0.00), -636.6221672311394281,      EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu( 1.0,   1.), -0.7812128213002887165,     EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu(30.0,   1.), -3.0481287832256432162e+39, EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu(30.0, 100.),  0.006138839212010033452,   EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu(10.0,   1.), -1.2161801427868918929e+08, EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu(10.0, 100.),  0.05833157423641492875,    EPS);
    assert_epeq!(gsl_bessel::bessel_Ynu(10.2, 100.),  0.07169383985546287091,    EPS);

    assert_epeq!(gsl_bessel::bessel_Inu_scaled(0.0001,10.), 0.12783333709581669672,    EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled( 1.0, 0.00), 0.0004995003123542213370,  EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled( 1.0,   1.), 0.20791041534970844887,    EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled(30.0,   1.), 1.3021094983785914437e-42, EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled(30.0, 100.), 0.0004486987756920986146,  EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled(10.0,   1.), 1.0127529864692066036e-10, EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled(10.0, 100.), 0.024176682718258828365,   EPS);
    assert_epeq!(gsl_bessel::bessel_Inu_scaled(10.2, 100.), 0.023691628843913810043,   EPS);

    assert_epeq!(gsl_bessel::bessel_Inu(0.0001,10.), 2815.7166269770030352,     EPS);
    assert_epeq!(gsl_bessel::bessel_Inu( 1.0, 0.00), 0.0005000000625000026042,  EPS);
    assert_epeq!(gsl_bessel::bessel_Inu( 1.0,   1.), 0.5651591039924850272,     EPS);
    assert_epeq!(gsl_bessel::bessel_Inu(30.0,   1.), 3.539500588106447747e-42,  EPS);
    assert_epeq!(gsl_bessel::bessel_Inu(30.0, 100.), 1.2061548704498434006e+40, EPS);
    assert_epeq!(gsl_bessel::bessel_Inu(10.0,   1.), 2.7529480398368736252e-10, EPS);
    assert_epeq!(gsl_bessel::bessel_Inu(10.0, 100.), 6.498975524720147799e+41,  EPS);
    assert_epeq!(gsl_bessel::bessel_Inu(10.2, 100.), 6.368587361287030443e+41,  EPS);

    assert_epeq!(gsl_bessel::bessel_Knu_scaled(0.0001,10.), 0.3916319346235421817, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled( 1.0, 0.00), 1000.9967345590684524, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled( 1.0,   1.), 1.6361534862632582465, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(30.0,   1.), 1.2792629867539753925e+40, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(30.0, 100.), 10.673443449954850040, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(10.0,   1.), 4.912296520990198599e+08, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(10.0, 100.), 0.20578687173955779807, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(10.0, 1000.), 0.04165905142800565788, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(10.0, 1.0e+), 0.00012533147624060789938, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu_scaled(10.2, 100.), 0.20995808355244385075, EPS);

    assert_epeq!(gsl_bessel::bessel_Knu(0.0001,0.00), 7.023689431812884141,      EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(0.0001,10.), 0.000017780062324654874306, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu( 1.0, 0.00), 999.9962381560855743,      EPS);
    assert_epeq!(gsl_bessel::bessel_Knu( 1.0,   1.), 0.6019072301972345747,     EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(10.0, 0.00), 1.8579455483904008064e+38, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(10.0,   1.), 1.8071328990102945469e+08, EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(10.0, 100.), 7.655427977388100611e-45,  EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(10.2, 100.), 7.810600225948217841e-45,  EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(30.0,   1.), 4.706145526783626883e+39,  EPS);
    assert_epeq!(gsl_bessel::bessel_Knu(30.0, 100.), 3.970602055959398739e-43,  EPS);

    assert_epeq!(gsl_bessel::bessel_lnKnu(0.0001,1.0e-10), 5.439794449319847, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(0.0001,0.000), 2.232835507214331, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(0.0001,10.), -10.93743282256098, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu( 1.0, 1.0e-10), 230.2585092994045, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu( 1.0, 1.0e-1), 23.025850929940456840, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu( 1.0, 0.00), 6.907751517131146, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu( 1.0,   1.), -0.5076519482107523309, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(30.0, 1.0e-10), 6999.113586185543475, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(30.0,   1.), 91.34968784026325464, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(30.0, 100.), -97.63224126416760932, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(100.0, 1.0e-10), 23453.606706185466825, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(100.0, 1.), 427.7532510250188083, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(100.0, 100.), -55.53422771502921431, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(1000.0, 1.0e-10), 236856.183755993135, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(10000.0, 1.0e-10), 2.39161558914890695e+06, EPS);

    /* [bug #31528] gsl_sf_bessel_lnKnu overflows for large nu */

    assert_epeq!(gsl_bessel::bessel_lnKnu(180.0, 2.), 735.1994170369583930752590258, EPS);
    assert_epeq!(gsl_bessel::bessel_lnKnu(3500.5, 1500.), 1731.220077116482710070986699, EPS);

    sa = 0;
    gsl_sf_bessel_jn_array(3, 38, 1.0, j);
    sa += ( test_sf_frac_diff(j[0],  0.0195633539826684059190  ) > TEST_TOL1);
    sa += ( test_sf_frac_diff(j[1],  0.0024766389641099550438  ) > TEST_TOL1);
    sa += ( test_sf_frac_diff(j[10], 1.9256167644801728904e-14 ) > TEST_TOL1);
    sa += ( test_sf_frac_diff(j[35], 6.911232970971626272e-57  ) > TEST_TOL1);
    gsl_test(sa, "  gsl_sf_bessel_jn_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_Yn_array(3, 38, 1.0, Y);
    sa += ( test_sf_frac_diff(Y[0],  -5.821517605964728848      ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[1],  -33.27842302897211870      ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[10], -1.2753618701519837595e+12 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[35], -1.2124435663593357154e+54 ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_Yn_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_In_scaled_array(3, 38, 1.0, I);
    sa += ( test_sf_frac_diff(I[0],  0.0081553077728142938170  ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(I[1],  0.0010069302573377758637  ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(I[10], 7.341518665628926244e-15  ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(I[35], 2.5753065298357542893e-57 ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_In_scaled_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_In_array(3, 38, 1.0, Y);
    sa += ( test_sf_frac_diff(Y[0],   0.0221684249243319024760  ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[1],   0.0027371202210468663251  ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[10],  1.9956316782072007564e-14 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[35],  7.000408942764452901e-57  ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_In_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_Kn_array(3, 38, 1.0, K);
    sa += ( test_sf_frac_diff(K[0],  7.101262824737944506 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[1],  44.23241584706284452 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[10], 1.9215763927929940846e+12 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[35], 1.8789385023806051223e+54 ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_Kn_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_Kn_scaled_array(3, 38, 1.0, K);
    sa += ( test_sf_frac_diff(K[0],  19.303233695596904277 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[1],  120.23617222591483717 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[10], 5.223386190525076473e+12 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[35], 5.107484387813251411e+54 ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_Kn_scaled_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_jl_array(50, 1.0, j);
    sa += ( test_sf_frac_diff(j[0],  0.84147098480789650670   ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(j[1],  0.30116867893975678925   ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(j[10], 7.116552640047313024e-11 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(j[50], 3.615274717489787311e-81 ) > TEST_TOL2 );
    gsl_test(sa, "  gsl_sf_bessel_jl_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_jl_steed_array(99, 1.0, j);
    sa += ( test_sf_frac_diff(j[0],  0.84147098480789650670   ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(j[1],  0.30116867893975678925   ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(j[10], 7.116552640047313024e-11 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(j[50], 3.615274717489787311e-81 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(j[80], 1.136352423414503264e-144 ) > TEST_TOL1 );
    gsl_test(sa, "  gsl_sf_bessel_jl_steed_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_yl_array(50, 1.0, Y);
    sa += ( test_sf_frac_diff(Y[0],  -0.5403023058681397174 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[1],  -1.3817732906760362241 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[10], -6.722150082562084436e+08  ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(Y[50], -2.7391922846297571576e+78 ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_yl_array");
    s += sa;

    {
        double Y0[1];
        sa = 0;
        gsl_sf_bessel_yl_array(0, 1.0, Y0);
        sa += ( test_sf_frac_diff(Y0[0],  -0.5403023058681397174 ) > TEST_TOL0 );
        gsl_test(sa, "  gsl_sf_bessel_yl_array (lmax=0)");
        s += sa;
    }

    sa = 0;
    gsl_sf_bessel_il_scaled_array(50, 1.0, I);
    sa += ( test_sf_frac_diff(I[0],  0.43233235838169365410 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(I[1],  0.13533528323661269189 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(I[10], 2.7343719371837065460e-11 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(I[50], 1.3429606061892023653e-81 ) > TEST_TOL2 );
    gsl_test(sa, "  gsl_sf_bessel_il_scaled_array");
    s += sa;

    sa = 0;
    gsl_sf_bessel_il_scaled_array(50, 0.0, I);
    sa += ( test_sf_frac_diff(I[0],  1.0 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(I[1],  0.0 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(I[10], 0.0 ) > TEST_TOL2 );
    sa += ( test_sf_frac_diff(I[50], 0.0 ) > TEST_TOL2 );
    gsl_test(sa, "  gsl_sf_bessel_il_scaled_array (L=0)");
    s += sa;

    sa = 0;
    gsl_sf_bessel_kl_scaled_array(50, 1.0, K);
    sa += ( test_sf_frac_diff(K[0],  1.5707963267948966192     ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[1],  3.1415926535897932385     ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[10], 2.7231075458948147010e+09 ) > TEST_TOL0 );
    sa += ( test_sf_frac_diff(K[50], 1.1578440432804522544e+79 ) > TEST_TOL0 );
    gsl_test(sa, "  gsl_sf_bessel_kl_scaled_array");
    s += sa;

    {
        double K0[1];
        sa = 0;
        gsl_sf_bessel_kl_scaled_array(0, 1.0, K0);
        sa += ( test_sf_frac_diff(K[0],  1.5707963267948966192     ) > TEST_TOL0 );
        gsl_test(sa, "  gsl_sf_bessel_kl_scaled_array (lmax=0)");
        s += sa;
    }

    sa = 0;
    sa += ( gsl_sf_bessel_zero_j0_e) != GSL_EINVAL );
    sa += ( r.val != 0.0 );
    s += sa;
    assert_epeq!(gsl_bessel::bessel_zero_j0( 1,  &r),  2.404825557695771, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j0( 2,  &r),  5.520078110286304, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j0(20,  &r), 62.048469190227081, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j0(25,  &r), 77.756025630388058, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j0(10), 313.37426607752784, EPS);

    sa = 0;
    sa += ( gsl_sf_bessel_zero_j1_e) !=EPS );
    sa += ( r.val != 0.0 );
    s += sa;
    assert_epeq!(gsl_bessel::bessel_zero_j1( 1,  &r), 3.831705970207512, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j1( 2,  &r), 7.015586669815619, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j1(20,  &r), 63.61135669848124, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j1(25,  &r), 79.32048717547630, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_j1(10), 314.9434728377672, EPS);

    sa = 0;
    sa += ( gsl_sf_bessel_zero_jnu_e(0.0, ) != GSL_EINVAL );
    sa += ( r.val != 0.0 );
    s += sa;
    assert_epeq!(gsl_bessel::bessel_zero_jnu(0.0,  1,  &r),  2.404825557695771, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(0.0,  2,  &r),  5.520078110286304, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(0.0, 20,  &r), 62.048469190227081, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(0.0, 25,  &r), 77.756025630388058, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(0.0, 10), 313.37426607752784, EPS);

    sa = 0;
    sa += ( gsl_sf_bessel_zero_jnu_e(1.0, ) !=EPS );
    sa += (r.val != 0.0);
    s += sa;
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  4.4934094579090641, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  8.7714838159599540, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  7.7252518369377072, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  12.338604197466944, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  10.904121659428900, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  15.700174079711671, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  14.066193912831473, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  18.980133875179921, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  17.220755271930768, EPS);

    /* Something wrong with the tolerances on these */
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  22.217799896561268, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 8.0, ),  26.266814641176644, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(20.0, ),  41.413065513892636, EPS);

    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  20.371302959287563, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  25.430341154222704, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 8.0, ),  29.545659670998550, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  23.519452498689007, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  28.626618307291138, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 8.0, ),  32.795800037341462, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  26.666054258812674, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  31.811716724047763, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(10.0, ),  38.761807017881651, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, ),  29.811598790892959, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, ),  34.988781294559295, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(10.0, ),  42.004190236671805, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 1.5, 1),  32.956389039822477, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 1),  38.159868561967132, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 1),  52.017241278881633, EPS);

    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 11,  &r), 41.326383254047406, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 11,  &r), 55.289204146560061, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 12,  &r), 44.4893191232197314, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 12,  &r), 58.5458289043850856, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 13,  &r), 47.6493998066970948, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 13,  &r), 61.7897598959450550, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 14,  &r), 50.8071652030063595, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 14,  &r), 65.0230502510422545, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 15,  &r), 53.9630265583781707, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 15,  &r), 68.2473219964207837, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 16,  &r), 57.1173027815042647, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 16,  &r), 71.4638758850226630, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 17,  &r), 60.2702450729428077, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 17,  &r), 74.6737687121404241, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 18,  &r), 63.4220540458757799, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 18,  &r), 77.8778689734863729, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 19,  &r), 66.5728918871182703, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 19,  &r), 81.0768977206328326, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 5.0, 20,  &r), 69.722891161716742,  EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(15.0, 20,  &r), 84.271459069716442,  EPS);

    assert_epeq!(gsl_bessel::bessel_zero_jnu( 23.0, 11,  &r), 65.843393469524653, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 30.0, 11,  &r), 74.797306585175426, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 32.0, 15,  &r), 90.913637691861741, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 50.0, 15,  &r), 113.69747988073942, EPS);

    assert_epeq!(gsl_bessel::bessel_zero_jnu(  5.0, 22,  &r), 76.020793430591605, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 10.0, 22,  &r), 83.439189796105756, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu( 12.0, 22,  &r), 86.345496520534055, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(100.0, 22,  &r), 199.82150220122519, EPS);
    assert_epeq!(gsl_bessel::bessel_zero_jnu(500.0, 22,  &r), 649.34132440891735, EPS);
}
*/
