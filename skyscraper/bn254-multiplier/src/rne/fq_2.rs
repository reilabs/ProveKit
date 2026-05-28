use super::simd_mul_fq;

/// Fused two-term sum-of-products over Fq: returns `(a0*b0 + a1*b1) mod q` in
/// Montgomery form.
///
/// # Preconditions
///
/// - All four inputs must be < 2^255 (5×51-bit limb precondition).
/// - Inputs in Fq-Montgomery form (R = 2^256).
#[inline(always)]
pub fn sum_of_products_2_fq(a0: [u64; 4], b0: [u64; 4], a1: [u64; 4], b1: [u64; 4]) -> [u64; 4] {
    let (p, q) = simd_mul_fq(a0, b0, a1, b1);

    // Each of p, q lives in [0, 2*Fq) per the mul invariant, so p + q < 4*Fq.
    // Since 4*Fq < 2^256 (Fq is 254 bits), the 4-limb add can't overflow.
    let s = add_4limb(p, q);

    // Bring back into [0, 2*Fq): one conditional subtract of 2*Fq suffices.
    let two_q = U64_2P_FQ;
    if cmp_geq(&s, &two_q) {
        sub_4limb(s, two_q)
    } else {
        s
    }
}

#[inline(always)]
fn add_4limb(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let mut r = [0u64; 4];
    let mut carry: u64 = 0;
    for i in 0..4 {
        let (s1, c1) = a[i].overflowing_add(b[i]);
        let (s2, c2) = s1.overflowing_add(carry);
        r[i] = s2;
        carry = (c1 as u64) + (c2 as u64);
    }
    debug_assert_eq!(carry, 0, "4P + 4P should not overflow 2^256 for BN254 Fq");
    r
}

#[inline(always)]
fn sub_4limb(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let mut r = [0u64; 4];
    let mut borrow: u64 = 0;
    for i in 0..4 {
        let (d1, b1) = a[i].overflowing_sub(b[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        r[i] = d2;
        borrow = (b1 as u64) + (b2 as u64);
    }
    debug_assert_eq!(borrow, 0, "subtrahend was larger than minuend");
    r
}

#[inline(always)]
fn cmp_geq(a: &[u64; 4], b: &[u64; 4]) -> bool {
    for i in (0..4).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    true
}

/// 2 * Fq in 64-bit limbs — used as the subtractive bound for
/// [`sum_of_products_2_fq`]. Lifted from `U64_P_FQ_MULTIPLES[2]` so the
/// constants stay in one place; kept here as an inline copy to avoid pulling
/// the whole multiples table into the hot path.
const U64_2P_FQ: [u64; 4] = [
    0x7841182db0f9fa8e,
    0x2f02d522d0e3951a,
    0x70a08b6d0302b0bb,
    0x60c89ce5c2634053,
];
