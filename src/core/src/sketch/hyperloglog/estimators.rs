use std::cmp;

pub type CounterType = u8;

pub fn counts(registers: &[CounterType], q: usize) -> Vec<u16> {
    let mut counts = vec![0; q + 2];

    for k in registers {
        counts[*k as usize] += 1;
    }

    counts
}

#[allow(clippy::many_single_char_names)]
pub fn mle(counts: &[u16], p: usize, q: usize, relerr: f64) -> f64 {
    let m = 1 << p;
    if counts[q + 1] == m {
        return std::f64::INFINITY;
    }

    let (k_min, _) = counts.iter().enumerate().find(|(_, v)| **v != 0).unwrap();
    let k_min_prime = cmp::max(1, k_min);

    let (k_max, _) = counts
        .iter()
        .enumerate()
        .rev()
        .find(|(_, v)| **v != 0)
        .unwrap();
    let k_max_prime = cmp::min(q, k_max as usize);

    let mut z = 0.;
    for i in num_iter::range_step_inclusive(k_max_prime as i32, k_min_prime as i32, -1) {
        z = 0.5 * z + (counts[i as usize] as f64);
    }

    // ldexp(x, i) = x * (2 ** i)
    z *= 2f64.powi(-(k_min_prime as i32));

    let mut c_prime = counts[q + 1];
    if q >= 1 {
        c_prime += counts[k_max_prime];
    }

    let mut g_prev = 0.;
    let a = z + (counts[0] as f64);
    let b = z + (counts[q + 1] as f64) * 2f64.powi(-(q as i32));
    let m_prime = (m - counts[0]) as f64;

    let mut x = if b <= 1.5 * a {
        // weak lower bound (47)
        m_prime / (0.5 * b + a)
    } else {
        // strong lower bound (46)
        m_prime / (b * (1. + b / a).ln())
    };

    let mut delta_x = x;
    let del = relerr / (m as f64).sqrt();
    while delta_x > x * del {
        // secant method iteration

        let kappa: usize = az::saturating_cast(2. + x.log2().floor());

        // x_prime in [0, 0.25]
        let mut x_prime = x * 2f64.powi(-(cmp::max(k_max_prime, kappa) as i32) - 1);
        let x_pp = x_prime * x_prime;

        // Taylor approximation (58)
        let mut h = x_prime - (x_pp / 3.) + (x_pp * x_pp) * (1. / 45. - x_pp / 472.5);

        // Calculate h(x/2^k), see (56), at this point x_prime = x / (2^(k+2))
        for _k in num_iter::range_step_inclusive(kappa as i32 - 1, k_max_prime as i32, -1) {
            let h_prime = 1. - h;
            h = (x_prime + h * h_prime) / (x_prime + h_prime);
            x_prime += x_prime;
        }

        // compare (53)
        let mut g = c_prime as f64 * h;

        for k in num_iter::range_step_inclusive(k_max_prime as i32 - 1, k_min_prime as i32, -1) {
            let h_prime = 1. - h;
            // Calculate h(x/2^k), see (56), at this point x_prime = x / (2^(k+2))
            h = (x_prime + h * h_prime) / (x_prime + h_prime);
            g += counts[k as usize] as f64 * h;
            x_prime += x_prime;
        }

        g += x * a;
        delta_x = if (g > g_prev) | (m_prime >= g) {
            // see (54)
            delta_x * (m_prime - g) / (g - g_prev)
        } else {
            0.
        };

        x += delta_x;
        g_prev = g
    }

    m as f64 * x
}

/// Calculate the joint maximum likelihood of A and B.
///
/// Returns a tuple (only in A, only in B, intersection)
pub fn joint_mle(
    k1: &[CounterType],
    k2: &[CounterType],
    p: usize,
    q: usize,
) -> (usize, usize, usize) {
    let mut c1 = vec![0; q + 2];
    let mut c2 = vec![0; q + 2];
    let mut cu = vec![0; q + 2];
    let mut cg1 = vec![0; q + 2];
    let mut cg2 = vec![0; q + 2];
    let mut ceq = vec![0; q + 2];

    for (k1_, k2_) in k1.iter().zip(k2.iter()) {
        match k1_.cmp(&k2_) {
            cmp::Ordering::Less => {
                c1[*k1_ as usize] += 1;
                cg2[*k2_ as usize] += 1;
            }
            cmp::Ordering::Greater => {
                cg1[*k1_ as usize] += 1;
                c2[*k2_ as usize] += 1;
            }
            cmp::Ordering::Equal => {
                ceq[*k1_ as usize] += 1;
            }
        }
        cu[*cmp::max(k1_, k2_) as usize] += 1;
    }

    for (i, (v, u)) in cg1.iter().zip(ceq.iter()).enumerate() {
        c1[i] += v + u;
    }

    for (i, (v, u)) in cg2.iter().zip(ceq.iter()).enumerate() {
        c2[i] += v + u;
    }

    let c_ax = mle(&c1, p, q, 0.01);
    let c_bx = mle(&c2, p, q, 0.01);
    let c_abx = mle(&cu, p, q, 0.01);

    let mut counts_axb_half = vec![0u16; q + 2];
    let mut counts_bxa_half = vec![0u16; q + 2];

    counts_axb_half[q] = k1.len() as u16;
    counts_bxa_half[q] = k2.len() as u16;

    for _q in 0..q {
        counts_axb_half[_q] = cg1[_q] + ceq[_q] + cg2[_q + 1];
        debug_assert!(counts_axb_half[q] >= counts_axb_half[_q]);
        counts_axb_half[q] -= counts_axb_half[_q];

        counts_bxa_half[_q] = cg2[_q] + ceq[_q] + cg1[_q + 1];
        debug_assert!(counts_bxa_half[q] >= counts_bxa_half[_q]);
        counts_bxa_half[q] -= counts_bxa_half[_q];
    }

    let c_axb_half = mle(&counts_axb_half, p, q - 1, 0.01);
    let c_bxa_half = mle(&counts_bxa_half, p, q - 1, 0.01);

    let cx1 = 1.5 * c_bx + 1.5 * c_ax - c_bxa_half - c_axb_half;
    let cx2 = 2. * (c_bxa_half + c_axb_half) - 3. * c_abx;

    (
        (c_abx - c_bx) as usize,
        (c_abx - c_ax) as usize,
        cmp::max(0, (0.5 * (cx1 + cx2)) as usize),
    )
}
