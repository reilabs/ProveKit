//! Constants for RNE Montgomery multiplication over the BN254 scalar field.

use crate::pow_2;

/// Field specific constants
pub trait RNEMulParams {
    const P: [u64; 5];
    const NP0: u64;
    const RHO_1: [u64; 5];
    const RHO_2: [u64; 5];
    const RHO_3: [u64; 5];
    const RHO_4: [u64; 5];
}

pub struct FRParams;

impl RNEMulParams for FRParams {
    const P: [u64; 5] = U51_P_FR;

    const NP0: u64 = U51_NP0_FR;

    const RHO_1: [u64; 5] = RHO_1_FR;

    const RHO_2: [u64; 5] = RHO_2_FR;

    const RHO_3: [u64; 5] = RHO_3_FR;

    const RHO_4: [u64; 5] = RHO_4_FR;
}

pub struct FQParams;

impl RNEMulParams for FQParams {
    const P: [u64; 5] = U51_P_FQ;

    const NP0: u64 = U51_NP0_FQ;

    const RHO_1: [u64; 5] = RHO_1_FQ;

    const RHO_2: [u64; 5] = RHO_2_FQ;

    const RHO_3: [u64; 5] = RHO_3_FQ;

    const RHO_4: [u64; 5] = RHO_4_FQ;
}

/// Montgomery reduction constant: `-p⁻¹ mod 2⁵¹`
pub const U51_NP0_FR: u64 = 0x1f593efffffff;

/// BN254 scalar field prime
pub const U51_P_FR: [u64; 5] = U51_P_MULTIPLES[1];

/// Lookup table: `U51_P_MULTIPLES[k]` = `k * P` for k in 0..=5, in 51-bit
/// limbs.
pub const U51_P_MULTIPLES: [[u64; 5]; 6] = [
    [
        0x0000000000000,
        0x0000000000000,
        0x0000000000000,
        0x0000000000000,
        0x0000000000000,
    ], // 0P
    [
        0x1f593f0000001,
        0x10f372e12287c,
        0x6056174a0cfa1,
        0x014dc2822db40,
        0x30644e72e131a,
    ], // 1P
    [
        0x3eb27e0000002,
        0x21e6e5c2450f8,
        0x40ac2e9419f42,
        0x029b85045b681,
        0x60c89ce5c2634,
    ], // 2P
    [
        0x5e0bbd0000003,
        0x32da58a367974,
        0x210245de26ee3,
        0x03e94786891c2,
        0x112ceb58a394e,
    ], // 3P
    [
        0x7d64fc0000004,
        0x43cdcb848a1f0,
        0x01585d2833e84,
        0x05370a08b6d03,
        0x419139cb84c68,
    ], // 4P
    [
        0x1cbe3b0000005,
        0x54c13e65aca6d,
        0x61ae747240e25,
        0x0684cc8ae4843,
        0x71f5883e65f82,
    ], // 5P
];

/// Bit mask for 51-bit limbs.
pub const MASK51: u64 = 2_u64.pow(51) - 1;

/// Reduction constants: `RHO_i = 2^{-51*i} mod p` in 51-bit limbs.
pub const RHO_1_FR: [u64; 5] = [
    0x05cc89dc987a4,
    0x64e24f262c77a,
    0x237f02685263f,
    0x70aad55e2a1fd,
    0x0bda088fbd071,
];

pub const RHO_2_FR: [u64; 5] = [
    0x3459f4a69e5e7,
    0x25faeea4c9ca7,
    0x3e771def3ca40,
    0x46003708f7bc8,
    0x088b040ada652,
];

pub const RHO_3_FR: [u64; 5] = [
    0x76fe2f2b3ebb4,
    0x6d028b8f2441f,
    0x461c7904ae683,
    0x71824d0dd38b7,
    0x18c6b0be26ceb,
];

pub const RHO_4_FR: [u64; 5] = [
    0x30bf04e2f27cc,
    0x039b11bea2ed3,
    0x2fb7665568cc8,
    0x0cc99c143d8f0,
    0x0523513296c10,
];

/// Montgomery reduction constant for Fq: `-q⁻¹ mod 2⁵¹`.
pub const U51_NP0_FQ: u64 = 0x20782e4866389;

/// BN254 base field prime in 5×51-bit limbs.
pub const U51_P_FQ: [u64; 5] = [
    0x08c16d87cfd47,
    0x22d0e3951a784,
    0x60561765e05aa,
    0x014dc2822db40,
    0x30644e72e131a,
];

/// Reduction constants for Fq: `RHO_i = 2^{-51*i} mod q` in 51-bit limbs.
pub const RHO_1_FQ: [u64; 5] = [
    0x12747bdc458aa,
    0x057d3a8941ddf,
    0x3f85a4f3e6b1f,
    0x079e650a08953,
    0x0c468325af424,
];

pub const RHO_2_FQ: [u64; 5] = [
    0x2db617d8e0c0d,
    0x33c4fdbde44d9,
    0x3a7be6e2bf9de,
    0x6d4f15f591fa0,
    0x1c11456c86287,
];

pub const RHO_3_FQ: [u64; 5] = [
    0x76866e318c8ac,
    0x580fdd937314d,
    0x68dad1d169bf2,
    0x4b27493600721,
    0x04e6a4adeff09,
];

pub const RHO_4_FQ: [u64; 5] = [
    0x5b604a24bca48,
    0x712b927f808f0,
    0x3fc758b4db1ae,
    0x08e0587815152,
    0x259169009dbe2,
];

pub const C1: f64 = pow_2(103);
pub const C2: f64 = pow_2(103) + pow_2(52) + pow_2(51);
pub const C3: f64 = pow_2(52) + pow_2(51);
