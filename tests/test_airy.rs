///////////////////////////////////////
// Special Functions: Airy Functions //
///////////////////////////////////////
#[macro_use]
extern crate rustsci;

use rustsci::gsl_airy;
use rustsci::gsl_math;

const EPS: f64 = 0.00000000001;

#[test]
fn test_gslairy_ai()
{
    let val = gsl_airy::airy_ai(1f64).0;
    assert!(gsl_math::gslmath_fcmp(val, 0.13529241631288141f64, 0.001f64));
}

#[test]
fn test_airy_ai()
{
    assert_epeq!(gsl_airy::airy_ai(-500.0).0,              0.0725901201040411396, EPS);
    assert_epeq!(gsl_airy::airy_ai(-5.0).0,                0.3507610090241142, EPS);
    assert_epeq!(gsl_airy::airy_ai(-0.3000000000000094).0, 0.4309030952855831, EPS);
    assert_epeq!(gsl_airy::airy_ai(0.6999999999999907).0,  0.1891624003981519, EPS);

    assert_epeq!(gsl_airy::airy_ai(1.649999999999991).0,   0.0583105861872088521, EPS);

    assert_epeq!(gsl_airy::airy_ai(2.54999999999999).0,    0.01446149513295428, EPS);
    assert_epeq!(gsl_airy::airy_ai(3.499999999999987).0,   0.002584098786989702, EPS);
    assert_epeq!(gsl_airy::airy_ai(5.39999999999998).0,    4.272986169411866e-05, EPS);
}

#[test]
fn test_airy_ai_scaled()
{
    assert_epeq!(gsl_airy::airy_ai_scaled(-5.0).0,                  0.3507610090241142, EPS);
    assert_epeq!(gsl_airy::airy_ai_scaled(0.6999999999999907).0, 0.2795125667681217, EPS);
    assert_epeq!(gsl_airy::airy_ai_scaled(1.649999999999991).0,  0.2395493001442741, EPS);
    assert_epeq!(gsl_airy::airy_ai_scaled(2.54999999999999).0,   0.2183658595899388, EPS);
    assert_epeq!(gsl_airy::airy_ai_scaled(3.499999999999987).0,  0.2032920808163519, EPS);
    assert_epeq!(gsl_airy::airy_ai_scaled(5.39999999999998).0,   0.1836050093282229, EPS);
}

#[test]
fn test_airy_bi()
{
    assert_epeq!(gsl_airy::airy_bi(-500.0).0,             -0.094688570132991028, EPS);
    assert_epeq!(gsl_airy::airy_bi(-5.0).0,               -0.1383691349016005,   EPS);
    assert_epeq!(gsl_airy::airy_bi(0.6999999999999907).0,  0.9733286558781599,   EPS);
    assert_epeq!(gsl_airy::airy_bi(1.649999999999991).0,   2.196407956850028,    EPS);
    assert_epeq!(gsl_airy::airy_bi(2.54999999999999).0,    6.973628612493443,    EPS);
    assert_epeq!(gsl_airy::airy_bi(3.499999999999987).0,   33.05550675461069,    EPS);
    assert_epeq!(gsl_airy::airy_bi(5.39999999999998).0,    1604.476078241272,    EPS);
}

#[test]
fn test_airy_bi_scaled()
{
    assert_epeq!(gsl_airy::airy_bi_scaled(-5.0).0,                -0.1383691349016005, EPS);
    assert_epeq!(gsl_airy::airy_bi_scaled(0.6999999999999907).0,  0.6587080754582302, EPS);
    assert_epeq!(gsl_airy::airy_bi_scaled(1.649999999999991).0,   0.5346449995597539, EPS);
    assert_epeq!(gsl_airy::airy_bi_scaled(2.54999999999999).0,    0.461835455542297,  EPS);
    assert_epeq!(gsl_airy::airy_bi_scaled(3.499999999999987).0,   0.4201771882353061, EPS);
    assert_epeq!(gsl_airy::airy_bi_scaled(5.39999999999998).0,    0.3734050675720473, EPS);
}

// derivatives

#[test]
fn test_airy_ai_deriv()
{
    assert_epeq!(gsl_airy::airy_ai_deriv(-5.0).0,                 0.3271928185544435,    EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv(-0.5500000000000094).0, -0.1914604987143629,    EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv(0.4999999999999906).0,  -0.2249105326646850,    EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv(1.899999999999992).0,   -0.06043678178575718,   EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv(3.249999999999988).0,   -0.007792687926790889,  EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv(5.199999999999981).0,   -0.0001589434526459543, EPS);
}

#[test]
fn test_airy_ai_deriv_scaled()
{
    assert_epeq!(gsl_airy::airy_ai_deriv_scaled(-5.0).0,                0.3271928185544435, EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv_scaled(0.5499999999999906).0, -0.2874057279170166, EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv_scaled(1.499999999999991).0,  -0.3314199796863637, EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv_scaled(2.49999999999999).0,   -0.3661089384751620, EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv_scaled(3.649999999999986).0,  -0.3974033831453963, EPS);
    assert_epeq!(gsl_airy::airy_ai_deriv_scaled(6.299999999999977).0,  -0.4508799189585947, EPS);
}

#[test]
fn test_airy_bi_deriv()
{
    assert_epeq!(gsl_airy::airy_bi_deriv(-5.0).0,                0.778411773001899,  EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv(-0.5500000000000094).0, 0.5155785358765014, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv(0.4999999999999906).0,  0.5445725641405883, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv(1.899999999999992).0,   3.495165862891568,  EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv(3.249999999999988).0,   36.55485149250338,  EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv(5.199999999999981).0,   2279.748293583233,  EPS);
}

#[test]
fn test_airy_bi_deriv_scaled()
{
    assert_epeq!(gsl_airy::airy_bi_deriv_scaled(-5.0).0,               0.778411773001899, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv_scaled(0.5499999999999906).0, 0.4322811281817566, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv_scaled(1.499999999999991).0,  0.5542307563918037, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv_scaled(2.49999999999999).0,   0.6755384441644985, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv_scaled(3.649999999999986).0,  0.7613959373000228, EPS);
    assert_epeq!(gsl_airy::airy_bi_deriv_scaled(6.299999999999977).0,  0.8852064139737571, EPS);
}

#[test]
fn test_airy_zero_ai()
{
    assert_epeq!(gsl_airy::airy_zero_ai(2).0,  -4.087949444130970617, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai(50).0, -38.02100867725525443, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai(100).0,  -60.45555727411669871, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai(110).0,  -64.43135670991324811, EPS);
}

#[test]
fn test_airy_zero_bi()
{
    assert_epeq!(gsl_airy::airy_zero_bi(2).0, -3.271093302836352716, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi(50).0, -37.76583438165180116, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi(100).0, -60.25336482580837088, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi(110).0, -64.2355167606561537, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi(111).0, -64.6268994819519378, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi(200).0, -95.88699147356682665, EPS);
}

#[test]
fn test_airy_zero_ai_deriv()
{
    assert_epeq!(gsl_airy::airy_zero_ai_deriv(2).0, -3.248197582179836561, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai_deriv(50).0, -37.76565910053887108, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai_deriv(100).0, -60.25329596442479317, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai_deriv(110).0, -64.23545617243546956, EPS);
    assert_epeq!(gsl_airy::airy_zero_ai_deriv(1000).0, -280.9378080358935071, EPS);
}

#[test]
fn test_airy_zero_bi_deriv()
{
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(2).0, -4.073155089071828216, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(50).0, -38.02083574095788210, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(100).0, -60.45548887257140819, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(110).0, -64.43129648944845060, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(111).0, -64.82208737584206093, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(200).0, -96.04731050310324450, EPS);
    assert_epeq!(gsl_airy::airy_zero_bi_deriv(1000).0, -281.0315164471118527, EPS);
}
