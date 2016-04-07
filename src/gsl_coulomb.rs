/// Coulomb functions from the GSL
use libc::{c_int, c_double};
use gsl_sf;
use array::Array;
use array::Order;

#[link(name = "gsl")]
extern
{
    /// Calculate the Clausen integral:
    ///   Cl_2(x) := Integrate[-Log[2 Sin[t/2]], {t,0,x}]
    /// Relation to dilogarithm:
    ///   Cl_2(theta) = Im[ Li_2(e^(i theta)) ]
    fn gsl_sf_clausen_e(x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Normalized hydrogenic bound states, radial dependence.
    /// R_1 := 2Z sqrt(Z) exp(-Z r)
    fn gsl_sf_hydrogenicR_1_e(
        Z: c_double,
        r: c_double,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
    /// normalization such that psi(n,l,r) = R_n Y_{lm}
    fn gsl_sf_hydrogenicR_e(
        n: c_int,
        l: c_int,
        Z: c_double,
        r: c_double,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Coulomb wave functions F_{lam_F}(eta,x), G_{lam_G}(eta,x)
    /// and their derivatives; lam_G := lam_F - k_lam_G
    /// lam_F, lam_G > -0.5
    /// x > 0.0
    /// Conventions of Abramowitz+Stegun.
    /// Because there can be a large dynamic range of values,
    /// overflows are handled gracefully. If an overflow occurs,
    /// GSL_EOVRFLW is signalled and exponent(s) are returned
    /// through exp_F, exp_G. These are such that
    ///   F_L(eta,x)  =  fc[k_L] /// exp(exp_F)
    ///   G_L(eta,x)  =  gc[k_L] /// exp(exp_G)
    ///   F_L'(eta,x) = fcp[k_L] /// exp(exp_F)
    ///   G_L'(eta,x) = gcp[k_L] /// exp(exp_G)
    fn gsl_sf_coulomb_wave_FG_e(
        eta: c_double,
        x: c_double,
        lam_F: c_double,
        k_lam_G: c_int,
        F: *mut gsl_sf::gsl_sf_result,
        Fp: *mut gsl_sf::gsl_sf_result,
        G: *mut gsl_sf::gsl_sf_result,
        Gp: *mut gsl_sf::gsl_sf_result,
        exp_F: *mut c_double,
        exp_G: *mut c_double) -> c_int;


    /// F_L(eta,x) as array
    fn gsl_sf_coulomb_wave_F_array(
        lam_min: c_double,
        kmax: c_int,
        eta: c_double,
        x: c_double,
        fc_array: *mut c_double,
        F_exponent: *mut c_double) -> c_int;

    /// F_L(eta,x), G_L(eta,x) as arrays
    fn gsl_sf_coulomb_wave_FG_array(
        lam_min: c_double,
        kmax: c_int,
        eta: c_double,
        x: c_double,
        fc_array: *mut c_double,
        gc_array: *mut c_double,
        F_exponent: *mut c_double,
        G_exponent: *mut c_double) -> c_int;

    /// F_L(eta,x), G_L(eta,x), F'_L(eta,x), G'_L(eta,x) as arrays
    fn gsl_sf_coulomb_wave_FGp_array(
        lam_min: c_double,
        kmax: c_int,
        eta: c_double,
        x: c_double,
        fc_array: *mut c_double,
        fcp_array: *mut c_double,
        gc_array: *mut c_double,
        gcp_array: *mut c_double,
        F_exponent: *mut c_double,
        G_exponent: *mut c_double) -> c_int;

    //// Coulomb wave function divided by the argument,
    /// F(eta, x)/x. This is the function which reduces to
    /// spherical Bessel functions in the limit eta->0.
    ////
    fn gsl_sf_coulomb_wave_sphF_array(
        lam_min: c_double,
        kmax: c_int,
        eta: c_double,
        x: c_double,
        fc_array: *mut c_double,
        F_exponent: *mut c_double) -> c_int;


    //// Coulomb wave function normalization constant.
    /// [Abramowitz+Stegun 14.1.8, 14.1.9]
    fn gsl_sf_coulomb_CL_e(
        L: c_double,
        eta: c_double,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    //// Coulomb wave function normalization constant, array version
    /// [Abramowitz+Stegun 14.1.8, 14.1.9]    
    fn gsl_sf_coulomb_CL_array(
        Lmin: c_double,
        kmax: c_int,
        eta: c_double,
        cl: *mut c_double) -> c_int;
}

/// Calculate the Clausen integral
pub fn clausen(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_clausen_e(x, &mut s) != 0
        {
            panic!("Clausen calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Normalized hydrogenic bound states, radial dependence.
pub fn hydrogenic_r1(z: f64, r: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_hydrogenicR_1_e(z, r, &mut s) != 0
        {
            panic!("Coloumb calculation failed");
        }
    }
    (s.val, s.err)    
}

/// R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
/// normalization such that psi(n,l,r) = R_n Y_{lm}
pub fn hydrogenic_rn(n: i32, l: i32, z: f64, r: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_hydrogenicR_e(n, l, z, r, &mut s) != 0
        {
            panic!("Coloumb calculation failed");
        }
    }
    (s.val, s.err)
}

/// Coulomb wave functions F_{lam_F}(eta,x), G_{lam_G}(eta,x)
/// and their derivatives; lam_g := lam_f - k_lam_g
pub fn coulomb_wave_fg(
    eta: f64,
    x: f64,
    lam_f: f64,
    k: usize) -> ([f64; 4], [f64; 4])
{
    if x < 0f64
    {
        panic!("x cannot be less than zero");
    }
    if lam_f - (k as f64) < -0.5f64
    {
        panic!("Invalid value for k");
    }
    let mut f = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    let mut fp = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    let mut g = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    let mut gp = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    let mut exp_f = 0.0f64;
    let mut exp_g = 0.0f64;
    unsafe {
        if gsl_sf_coulomb_wave_FG_e(
            eta,
            x,
            lam_f,
            k as i32,
            &mut f,
            &mut fp,
            &mut g,
            &mut gp,
            &mut exp_f as *mut c_double,
            &mut exp_g as *mut c_double) != 0
        {
            panic!("Coloumb function failed");
        }
    }
    ([f.val, fp.val, g.val, gp.val], [f.err, fp.err, g.err, gp.err])
}


/// F_L(eta,x) as array
pub fn coulomb_wave_f_array(
    lam_min: f64,
    kmax: usize,
    eta: f64,
    x: f64) -> Array<f64>
{
    let mut fc_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut f_exponent: c_double = 0.0;
    unsafe {
         if 0 != gsl_sf_coulomb_wave_F_array(
            lam_min,
            kmax as c_int,
            eta,
            x,
            fc_array.as_mut_ptr() as *mut c_double,
            &mut f_exponent as *mut c_double)
        {
            panic!("Coloumb function failed");
        }
    }
    fc_array
}

/// F_L(eta,x), G_L(eta,x) as arrays
pub fn coulomb_wave_fg_array(
    lam_min: f64,
    kmax: usize,
    eta: f64,
    x: f64) -> (Array<f64>, Array<f64>)
{
    let mut fc_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut gc_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut f_exponent: c_double = 0f64;
    let mut g_exponent: c_double = 0f64;
    unsafe {
        if 0 != gsl_sf_coulomb_wave_FG_array(
            lam_min,
            kmax as i32,
            eta,
            x,
            fc_array.as_mut_ptr() as *mut c_double,
            gc_array.as_mut_ptr() as *mut c_double,
            &mut f_exponent as *mut c_double,
            &mut g_exponent as *mut c_double)
        {
            panic!("Coloumb function failed");
        }
    }
    (fc_array, gc_array)
}

/// F_L(eta,x), G_L(eta,x), F'_L(eta,x), G'_L(eta,x) as arrays
pub fn coulomb_wave_fgp_array(
    lam_min: f64,
    kmax: usize,
    eta: f64,
    x: f64) -> (Array<f64>, Array<f64>, Array<f64>, Array<f64>)
{
    let mut fc_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut gc_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut fcp_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut gcp_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut f_exponent: c_double = 0f64;
    let mut g_exponent: c_double = 0f64;
    unsafe {
        if 0 != gsl_sf_coulomb_wave_FGp_array(
            lam_min,
            kmax as i32,
            eta,
            x,
            fc_array.as_mut_ptr() as *mut c_double,
            fcp_array.as_mut_ptr() as *mut c_double,
            gc_array.as_mut_ptr() as *mut c_double,
            gcp_array.as_mut_ptr() as *mut c_double,
            &mut f_exponent as *mut c_double,
            &mut g_exponent as *mut c_double)
        {
            panic!("Coloumb function failed");
        }
    }
    (fc_array, gc_array, fcp_array, gcp_array)
}

//// Coulomb wave function divided by the argument,
/// F(eta, x)/x. This is the function which reduces to
/// spherical Bessel functions in the limit eta->0.
pub fn coulomb_wave_sphf_array(
    lam_min: f64,
    kmax: usize,
    eta: f64,
    x: f64) -> Array<f64>
{
    let mut fc_array = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    let mut f_exponent: c_double = 0.0;
    unsafe {
        if 0 != gsl_sf_coulomb_wave_sphF_array(
            lam_min,
            kmax as i32,
            eta,
            x,
            fc_array.as_mut_ptr() as *mut c_double,
            &mut f_exponent as *mut c_double)
        {
            panic!("Coloumb function failed");
        }
    }
    fc_array
}


//// Coulomb wave function normalization constant.
/// [Abramowitz+Stegun 14.1.8, 14.1.9]
pub fn coulomb_cl_e(l: f64, eta: f64) -> (f64, f64)
{
    let mut result = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if 0 != gsl_sf_coulomb_CL_e(l, eta, &mut result)
        {
            panic!("Coloumb function failed");
        }
    }
    (result.val, result.err)
}

//// Coulomb wave function normalization constant, array version
/// [Abramowitz+Stegun 14.1.8, 14.1.9]    
pub fn coulomb_cl_array(lmin: f64, kmax: usize, eta: f64) -> Array<f64>
{
    let mut cl: Array<f64> = Array::<f64>::new_filled(0f64, kmax, Order::Row);
    unsafe {
        if 0 != gsl_sf_coulomb_CL_array(
            lmin,
            kmax as i32,
            eta,
            cl.as_mut_ptr() as *mut c_double)
        {
            panic!("Coloumb function failed");
        }
    }
    cl
}
