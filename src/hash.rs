use std::hash::{BuildHasher, Hasher};
use rustc_hash::FxHasher;

#[inline]
pub fn compute_fingerprint(mut kmer: u64) -> u16 {
    kmer ^= kmer >> 33;
    kmer = kmer.wrapping_mul(0xff51afd7ed558ccd);
    kmer ^= kmer >> 33;
    (kmer & 0xFFFF) as u16
}

#[inline]
pub fn hash_slice(slice: &[u32]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = FxHasher::default();
    slice.hash(&mut hasher);
    hasher.finish()
}

#[inline(always)]
    pub fn bloom_hash_pair(item: u64) -> (u64, u64) {
        let mut h1 = item;
        h1 ^= h1 >> 33;
        h1 = h1.wrapping_mul(0xff51afd7ed558ccd);
        h1 ^= h1 >> 33;
        h1 = h1.wrapping_mul(0xc4ceb9fe1a85ec53);
        h1 ^= h1 >> 33;

        let mut h2 = item;
        h2 ^= h2 >> 31;
        h2 = h2.wrapping_mul(0x85ebca6b);
        h2 ^= h2 >> 13;
        h2 = h2.wrapping_mul(0xc2b2ae35);
        h2 ^= h2 >> 16;
        h2 |= 1;

        (h1, h2)
    }

// This is for HLL++
#[derive(Default)]
pub struct Fmix64Hasher(u64);

impl Hasher for Fmix64Hasher {
    #[inline(always)]
    fn write_u64(&mut self, i: u64) {
        let mut k = i;
        k ^= k >> 33;
        k = k.wrapping_mul(0xff51afd7ed558ccd);
        k ^= k >> 33;
        k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
        k ^= k >> 33;
        self.0 = k;
    }

    #[inline(always)]
    fn write(&mut self, bytes: &[u8]) {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(&bytes[..8]);
        self.write_u64(u64::from_ne_bytes(arr));
    }

    #[inline(always)]
    fn finish(&self) -> u64 {
        self.0
    }
}

#[derive(Default, Clone)]
pub struct BuildFmix64Hasher;

impl BuildHasher for BuildFmix64Hasher {
    type Hasher = Fmix64Hasher;
    fn build_hasher(&self) -> Self::Hasher {
        Fmix64Hasher::default()
    }
}