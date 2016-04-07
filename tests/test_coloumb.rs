//////////////////////////////////////////
// Special Functions: Coulomb Functions //
//////////////////////////////////////////
#[macro_use]
extern crate rustsci;

use rustsci::gsl_coulomb;
use rustsci::gsl_math;

const EPS: f64 = 0.00000000001;

#[test]
pub fn test_hydrogenic_r1()
{
    assert_epeq!(gsl_coulomb::hydrogenic_r1(3.0, 2.0).0,  0.025759948256148471036f64,  EPS);
    assert_epeq!(gsl_coulomb::hydrogenic_r1(3.0, 10.0).0, 9.724727052062819704e-13f64, EPS);
}

#[test]
pub fn test_hydrogenic_rn()
{
    assert_epeq!(gsl_coulomb::hydrogenic_rn(4, 1, 3.0, 0.0).0,  0f64,  EPS);
    assert_epeq!(gsl_coulomb::hydrogenic_rn(4, 0, 3.0, 2.0).0, -0.03623182256981820062f64,  EPS);
    assert_epeq!(gsl_coulomb::hydrogenic_rn(4, 1, 3.0, 2.0).0, -0.028065049083129581005f64, EPS);
    assert_epeq!(gsl_coulomb::hydrogenic_rn(4, 2, 3.0, 2.0).0,  0.14583027278668431009f64,  EPS);

    assert_epeq!(gsl_coulomb::hydrogenic_rn(100,  0, 3.0, 2.0).0, -0.00007938950980052281367f64, EPS);
    assert_epeq!(gsl_coulomb::hydrogenic_rn(100, 10, 3.0, 2.0).0,  7.112823375353605977e-12f64,  EPS);
    assert_epeq!(gsl_coulomb::hydrogenic_rn(100, 90, 3.0, 2.0).0,  5.845231751418131548e-245f64, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_1()
{
    let lam_f = 0f64;
    let k_g: usize   = 0;
    let eta = 1f64;
    let x = 5f64;  
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(f,  0.6849374120059439677, EPS);
    assert_epeq!(fp, -0.7236423862556063963, EPS);
    assert_epeq!(g, -0.8984143590920205487, EPS);
    assert_epeq!(gp, -0.5108047585190350106, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_2()
{
    let lam_f = 10f64;
    let k_g: usize = 2;
    let eta = 1f64;
    let x = 5f64;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  0.0006423773354915823698, EPS);
    assert_epeq!( fp,  0.0013299570958719702545, EPS);
    assert_epeq!(  g,  33.27615734455096130,     EPS);
    assert_epeq!( gp, -45.49180102261540580,     EPS);
}

#[test]
pub fn test_coloumb_wave_fg_3()
{
    let lam_f = 4f64;
    let k_g: usize = 2;
    let eta = 50f64;
    let x = 120f64;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    assert_epeq!(  f,  0.0735194711823798495, EPS);
    assert_epeq!( fp,  0.6368149124126783325, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_4()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = -1000.0;
    let x = 1.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  9.68222518991341e-02, EPS);
    assert_epeq!( fp,  5.12063396274631e+00, EPS);
    assert_epeq!(  g,  1.13936784379472e-01, EPS);
    assert_epeq!( gp, -4.30243486522438e+00, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_5()
{
    let lam_min = 0.0;
    let k_g: usize = 0;
    let eta = -50.0;
    let x = 5.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_min, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  1.52236975714236e-01, EPS);
    assert_epeq!( fp,  2.03091041166137e+00, EPS);
    assert_epeq!(  g,  4.41680690236251e-01, EPS);
    assert_epeq!( gp, -6.76485374766869e-01, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_6()
{
    let lam_min = 0.0;
    let eta = -50.0;
    let x = 1000.0;
    let k_g: usize = 0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_min, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f, -0.2267212182760888523, EPS);
    assert_epeq!( fp, -0.9961306810018401525, EPS);
    assert_epeq!(  g, -0.9497684438900352186, EPS);
    assert_epeq!( gp,  0.2377656295411961399, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_7()
{
    let lam_f = 10.0;
    let k_g: usize = 0;
    let eta = -50.0;
    let x = 5.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f, -3.681143602184922e-01, EPS);
    assert_epeq!( fp,  1.338467510317215e+00, EPS);
    assert_epeq!(  g,  3.315883246109351e-01, EPS);
    assert_epeq!( gp,  1.510888628136180e+00, EPS);
}

#[test]
pub fn test_coloumb_wave_fg_8()
{
    let lam_f = 0.0;
    let k_g: usize = 0;
    let eta = -4.0;
    let x = 5.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  4.078627230056172e-01, EPS);
    assert_epeq!( fp,  1.098212336357310e+00, EPS);
    assert_epeq!(  g,  6.743270353832442e-01, EPS);
    assert_epeq!( gp, -6.361104272804447e-01, EPS);
}

#[test]
pub fn test_coulomb_wave_fg_9()
{
    let lam_f = 3.0;
    let k_g: usize = 0;
    let eta = -4.0;
    let x = 5.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f, -2.568630935581323e-01, EPS);
    assert_epeq!( fp,  1.143229422014827e+00, EPS);
    assert_epeq!(  g,  7.879899223927996e-01, EPS);
    assert_epeq!( gp,  3.859905878106713e-01, EPS);
}

#[test]
pub fn test_coulomb_wave_fg_10()
{
    let lam_f = 0.0;
    let k_g: usize = 0;
    let eta = 1.0;
    let x = 2.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  6.61781613832681e-01, EPS);
    assert_epeq!( fp,  4.81557455709949e-01, EPS);
    assert_epeq!(  g,  1.27577878476828e+00, EPS);
    assert_epeq!( gp, -5.82728813097184e-01, EPS);
}

#[test]
pub fn test_coulomb_wave_fg_11()
{
    let lam_f = 0.0;
    let k_g: usize = 0;
    let eta = 1.0;
    let x = 0.5;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  0.08315404535022023302, EPS);
    assert_epeq!( fp,  0.22693874616222787568, EPS);
    assert_epeq!(  g,  3.1060069279548875140,  EPS);
    assert_epeq!( gp, -3.549156038719924236,   EPS);
}

#[test]
pub fn test_coulomb_wave_fg_12()
{
    let lam_f = 0.5;
    let k_g: usize = 0;
    let eta = 1.0;
    let x = 0.5;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  0.04049078073829290935, EPS);
    assert_epeq!( fp,  0.14194939168094778795, EPS);
    assert_epeq!(  g,  4.720553853049677897,   EPS);
    assert_epeq!( gp, -8.148033852319180005,   EPS);
}

#[test]
pub fn test_coulomb_wave_fg_13()
{
    let lam_f = 0.1;
    let k_g: usize = 0;
    let eta = 1.0;
    let x = 0.5;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  0.07365466672379703418, EPS);
    assert_epeq!( fp,  0.21147121807178518647, EPS);
    assert_epeq!(  g,  3.306705446241024890, EPS);
    assert_epeq!( gp, -4.082931670935696644, EPS);
}

#[test]
pub fn test_coulomb_wave_fg_14()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 8.0;
    let x = 1.05;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  9.882706082810274357e-09, EPS);
    assert_epeq!( fp,  4.005167028235547770e-08, EPS);
    assert_epeq!(  g,  1.333127992006686320e+07, 10000f64*EPS);
    assert_epeq!( gp, -4.715914530842402330e+07, 10000f64*EPS);
}

#[test]
pub fn test_coulomb_wave_fg_15()
{
    let lam_f = 0.1;
    let k_g = 0;
    let eta = 8.0;
    let x = 1.05;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  9.611416736061987761e-09, EPS);
    assert_epeq!( fp,  3.909628126126824140e-08, EPS);
    assert_epeq!(  g,  1.365928464219262581e+07, 100000f64*EPS);
    assert_epeq!( gp, -4.848117385783386850e+07, 100000f64*EPS);
}

#[test]
#[ignore]
pub fn test_coulomb_wave_fg_16()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 50.0;
    let x = 0.1;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  2.807788027954216071e-67, EPS);
    assert_epeq!( fp,  9.677600748751576606e-66, EPS);
    assert_epeq!(  g,  5.579810686998358766e+64, EPS);
    assert_epeq!( gp, -1.638329512756321424e+66, EPS);
}

#[test]
pub fn test_coulomb_wave_fg_17()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 10.0;
    let x = 5.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  1.7207454091787930614e-06, 100000000f64*EPS);
    assert_epeq!( fp,  3.0975994706405458046e-06, 100000000f64*EPS);
    assert_epeq!(  g,  167637.56609459967623, 100000000f64*EPS);
    assert_epeq!( gp, -279370.76655361803075, 100000000f64*EPS);
}

#[test]
pub fn test_coulomb_wave_fg_18()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 25.0;
    let x = 10.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  1.5451274501076114315e-16, 100000000f64*EPS);
    assert_epeq!( fp,  3.1390869393378630928e-16, 100000000f64*EPS);
    assert_epeq!(  g,  1.6177129008336318136e+15, 100000000f64*EPS);
    assert_epeq!( gp, -3.1854062013149740860e+15, 100000000f64*EPS);
}

#[test]
pub fn test_coulomb_wave_fg_19()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 1.0;
    let x = 9.2;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f, -0.25632012319757955655, EPS);
    assert_epeq!( fp,  0.91518792286724220370, EPS);
    assert_epeq!(  g,  1.03120585918973466110, EPS);
    assert_epeq!( gp,  0.21946326717491250193, EPS);
}

#[test]
pub fn test_coulomb_wave_fg_20()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 10.0;
    let x = 10.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  0.0016262711250135878249, 10000000f64*EPS);
    assert_epeq!( fp,  0.0017060476320792806014, 100000000f64*EPS);
    assert_epeq!(  g,  307.87321661090837987, 10000000f64*EPS);
    assert_epeq!( gp, -291.92772380826822871, 10000000f64*EPS);
}

#[test]
pub fn test_coulomb_wave_fg_21()
{
    let lam_f = 0.0;
    let k_g = 0;
    let eta = 100.0;
    let x = 1.0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  8.999367996930662705e-126, 1e8f64*EPS);
    assert_epeq!( fp,  1.292746745757069321e-124, 1e8f64*EPS);
    assert_epeq!(  g,  3.936654148133683610e+123, 1e8f64*EPS);
    assert_epeq!( gp, -5.456942268061526371e+124, 1e8f64*EPS);
}

#[test]
pub fn test_coulomb_wave_fg_22()
{
    let lam_f = 1.0;
    let eta = 0.0;
    let x = 3.25;
    let k_g = 0;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  x.sin()/x - x.cos(), EPS);
    assert_epeq!( fp,  -x.sin()/(x*x) + x.cos()/x + x.sin(), EPS);
    assert_epeq!(  g,  x.cos()/x + x.sin(), EPS);
    assert_epeq!( gp,  -x.cos()/(x*x) - x.sin()/x + x.cos(), EPS);
}

#[test]
pub fn test_coulomb_wave_fg_23()
{
    let lam_f = 1.0;
    let eta = 0.0;
    let x = 3.25;
    let k_g = 1;
    let (vals, _) = gsl_coulomb::coulomb_wave_fg(eta, x, lam_f, k_g);
    let f = vals[0];
    let fp = vals[1];
    let g = vals[2];
    let gp = vals[3];
    assert_epeq!(  f,  x.sin()/x - x.cos(), EPS);
    assert_epeq!( fp,  -x.sin()/(x*x) + x.cos()/x +x.sin(), EPS);
    assert_epeq!(  g,  x.cos(), EPS);
    assert_epeq!( gp,  -x.sin(), EPS);
}
