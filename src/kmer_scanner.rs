//! # K-mer Scanner Module
//!
//! This module provides high-performance k-mer extraction and manipulation functions.
//! It's designed to be memory-efficient and blazingly fast for genomic sequence processing.
//!
//! ## Key Features
//! - Zero-copy k-mer extraction where possible
//! - Efficient canonical k-mer representation (minimum of k-mer and reverse complement)
//! - Bit-packed k-mer storage (2 bits per nucleotide)
//! - Spaced seed support for sensitive alignment
//! - SIMD-friendly operations where applicable

use std::fmt;

/// Type alias for k-mer representation
/// We use u64 which can store up to 32-mers (32 * 2 bits = 64 bits)
pub type Kmer = u64;

/// Maximum k-mer size that fits in a u64 (32 nucleotides)
pub const MAX_KMER_SIZE: usize = 32;

/// Lookup table for converting ASCII nucleotides to 2-bit encoding
/// A=0, C=1, G=2, T=3
/// This is a const array initialized at compile time for zero runtime cost
const NUCLEOTIDE_TO_CODE: [u8; 256] = {
    let mut table = [4u8; 256]; // 4 means invalid
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// Lookup table for reverse complement
/// A(0) <-> T(3), C(1) <-> G(2)
const COMPLEMENT_TABLE: [u8; 4] = [3, 2, 1, 0];

/// Error types for k-mer operations
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KmerError {
    /// K-mer size exceeds maximum allowed (32)
    KmerTooLarge(usize),
    /// Sequence too short to extract k-mer
    SequenceTooShort { seq_len: usize, k: usize },
    /// Invalid nucleotide character found
    InvalidNucleotide(char),
}

impl fmt::Display for KmerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KmerError::KmerTooLarge(k) => {
                write!(f, "K-mer size {} exceeds maximum of {}", k, MAX_KMER_SIZE)
            }
            KmerError::SequenceTooShort { seq_len, k } => {
                write!(f, "Sequence length {} is shorter than k={}", seq_len, k)
            }
            KmerError::InvalidNucleotide(c) => {
                write!(f, "Invalid nucleotide character: '{}'", c)
            }
        }
    }
}

impl std::error::Error for KmerError {}

/// K-mer scanner that efficiently extracts k-mers from sequences
///
/// # Design Philosophy
/// - Store k as part of the scanner to avoid passing it everywhere
/// - Pre-compute masks for fast bit operations
/// - Use bit manipulation instead of string operations for speed
///
/// # Memory Layout
/// Each k-mer uses exactly 2k bits, packed into a u64
/// Example for "ACGT": 00 01 10 11 = 0b00011011 = 27
pub struct KmerScanner {
    /// K-mer size (number of nucleotides)
    k: usize,
    /// Bit mask to keep only the rightmost k nucleotides (2k bits)
    /// For k=31: mask = (1 << 62) - 1 = 0x3FFFFFFFFFFFFFFF
    kmer_mask: Kmer,
}

impl KmerScanner {
    /// Create a new k-mer scanner
    ///
    /// # Arguments
    /// * `k` - K-mer size (must be 1-32)
    ///
    /// # Examples
    /// ```
    /// use kraken_rs::kmer_scanner::KmerScanner;
    /// let scanner = KmerScanner::new(31).unwrap();
    /// ```
    pub fn new(k: usize) -> Result<Self, KmerError> {
        if k == 0 || k > MAX_KMER_SIZE {
            return Err(KmerError::KmerTooLarge(k));
        }

        // Create mask: for k=31, we need 62 bits set (2 bits per base)
        // Bit shift: (1u64 << 62) - 1 = 0x3FFFFFFFFFFFFFFF
        let kmer_mask = if k == MAX_KMER_SIZE {
            // Special case: all 64 bits
            u64::MAX
        } else {
            (1u64 << (k * 2)) - 1
        };

        Ok(Self { k, kmer_mask })
    }

    /// Get the k-mer size
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Encode a single nucleotide to 2-bit representation
    ///
    /// # Returns
    /// - `Some(code)` where code is 0-3 for valid nucleotides
    /// - `None` for invalid characters
    ///
    /// # Performance
    /// This is a single array lookup - extremely fast O(1) operation
    #[inline(always)] // Always inline for performance
    pub fn encode_nucleotide(nucleotide: u8) -> Option<u8> {
        let code = NUCLEOTIDE_TO_CODE[nucleotide as usize];
        if code < 4 {
            Some(code)
        } else {
            None
        }
    }

    /// Encode a sequence of nucleotides into a k-mer
    ///
    /// # Arguments
    /// * `sequence` - Byte slice of nucleotides (must be exactly k long)
    ///
    /// # Returns
    /// * `Ok(kmer)` - The encoded k-mer
    /// * `Err(KmerError)` - If sequence is wrong length or contains invalid nucleotides
    ///
    /// # Algorithm
    /// We build the k-mer left-to-right by:
    /// 1. Shift existing bits left by 2 (make room for new nucleotide)
    /// 2. OR in the new 2-bit code
    /// 3. Mask to keep only k nucleotides
    ///
    /// # Example
    /// ```text
    /// Sequence: "ACG"
    /// A=0, C=1, G=2
    ///
    /// Step 1: kmer = 0, add A(0):  0b000000 | 0b00 = 0b000000
    /// Step 2: shift left: 0b000000 << 2 = 0b000000, add C(1): 0b000001
    /// Step 3: shift left: 0b000001 << 2 = 0b000100, add G(2): 0b000110
    /// Final: 0b000110 = 6
    /// ```
    pub fn encode_sequence(&self, sequence: &[u8]) -> Result<Kmer, KmerError> {
        if sequence.len() != self.k {
            return Err(KmerError::SequenceTooShort {
                seq_len: sequence.len(),
                k: self.k,
            });
        }

        let mut kmer: Kmer = 0;

        // Process each nucleotide
        for &nucleotide in sequence {
            let code = Self::encode_nucleotide(nucleotide)
                .ok_or_else(|| KmerError::InvalidNucleotide(nucleotide as char))?;

            // Shift left by 2 bits and add new nucleotide
            kmer = ((kmer << 2) | (code as Kmer)) & self.kmer_mask;
        }

        Ok(kmer)
    }

    /// Compute reverse complement of a k-mer
    ///
    /// # Algorithm
    /// 1. Reverse the order of nucleotides (bit reversal)
    /// 2. Complement each nucleotide (A<->T, C<->G)
    ///
    /// # Bit Manipulation Magic
    /// We use parallel bit manipulation to reverse and complement in O(1) time:
    /// - Swap pairs of bits
    /// - Swap nibbles (4 bits)
    /// - Swap bytes
    /// - Finally complement by XOR
    ///
    /// # Example
    /// ```text
    /// ACG (k=3) = 0b000110
    /// Reverse: GCA = 0b100100
    /// Complement: CGT = 0b011011
    /// ```
    pub fn reverse_complement(&self, mut kmer: Kmer) -> Kmer {
        // Step 1: Reverse pairs of bits (reverse each nucleotide's position)
        // Pattern: 0xCCCC... = 11001100... (alternating pairs)
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);

        // Step 2: Reverse nibbles (4-bit groups)
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);

        // Step 3: Reverse bytes
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);

        // Step 4: Reverse 16-bit groups
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);

        // Step 5: Reverse 32-bit groups
        kmer = (kmer >> 32) | (kmer << 32);

        // Step 6: Shift right to align (we reversed all 64 bits, but only need 2k bits)
        kmer >>= 64 - (self.k * 2);

        // Step 7: Complement each nucleotide (XOR with 0b11 to flip bits)
        // For each 2-bit nucleotide: A(00)->T(11), C(01)->G(10), G(10)->C(01), T(11)->A(00)
        (!kmer) & self.kmer_mask
    }

    /// Get canonical k-mer (lexicographically smaller of k-mer and its reverse complement)
    ///
    /// # Why Canonical?
    /// DNA is double-stranded, so "ACGT" and "ACGT" (reverse complement) represent
    /// the same biological sequence. We always use the smaller one to save memory
    /// and enable efficient lookup.
    ///
    /// # Performance
    /// This is just one comparison after computing reverse complement
    #[inline]
    pub fn canonical(&self, kmer: Kmer) -> Kmer {
        let revcomp = self.reverse_complement(kmer);
        if kmer <= revcomp {
            kmer
        } else {
            revcomp
        }
    }

    /// Extract all k-mers from a sequence (sliding window)
    ///
    /// # Arguments
    /// * `sequence` - Byte slice of nucleotide sequence
    ///
    /// # Returns
    /// Vector of k-mers. If sequence contains invalid nucleotides,
    /// the k-mer stream is reset at that position.
    ///
    /// # Algorithm - Rolling Hash
    /// Instead of encoding each k-mer from scratch, we use a rolling window:
    /// 1. Encode the first k-mer completely
    /// 2. For each subsequent position:
    ///    - Shift left by 2 bits (remove leftmost nucleotide)
    ///    - OR in new rightmost nucleotide
    ///    - Mask to keep only k nucleotides
    ///
    /// This is O(n) instead of O(nk) - MUCH faster!
    ///
    /// # Example
    /// ```text
    /// Sequence: "ACGTAC" (k=3)
    /// 
    /// Position 0: ACG = encode(A,C,G) = 0b000110
    /// Position 1: CGT = (ACG << 2) | T = (0b000110 << 2) | 0b11 = 0b011011
    /// Position 2: GTA = (CGT << 2) | A = (0b011011 << 2) | 0b00 = 0b101100
    /// Position 3: TAC = (GTA << 2) | C = (0b101100 << 2) | 0b01 = 0b110001
    /// ```
    pub fn extract_kmers(&self, sequence: &[u8]) -> Vec<Kmer> {
        if sequence.len() < self.k {
            return Vec::new();
        }

        // Pre-allocate exact size needed
        let num_kmers = sequence.len() - self.k + 1;
        let mut kmers = Vec::with_capacity(num_kmers);

        let mut current_kmer: Kmer = 0;
        let mut valid_bases = 0;

        for (i, &nucleotide) in sequence.iter().enumerate() {
            // Try to encode the nucleotide
            if let Some(code) = Self::encode_nucleotide(nucleotide) {
                // Shift and add new nucleotide
                current_kmer = ((current_kmer << 2) | (code as Kmer)) & self.kmer_mask;
                valid_bases += 1;

                // Once we have k valid bases, start emitting k-mers
                if valid_bases >= self.k {
                    kmers.push(current_kmer);
                }
            } else {
                // Invalid nucleotide - reset the window
                current_kmer = 0;
                valid_bases = 0;
            }
        }

        kmers
    }

    /// Extract canonical k-mers from a sequence
    ///
    /// This is the most commonly used function - it extracts all k-mers
    /// and returns their canonical forms.
    pub fn extract_canonical_kmers(&self, sequence: &[u8]) -> Vec<Kmer> {
        self.extract_kmers(sequence)
            .into_iter()
            .map(|kmer| self.canonical(kmer))
            .collect()
    }

    /// Apply a spaced seed mask to a k-mer
    ///
    /// # What is a Spaced Seed?
    /// Instead of using all positions in a k-mer, we use only certain positions.
    /// This can increase sensitivity for divergent sequences.
    ///
    /// # Example
    /// ```text
    /// Pattern: "111011" (positions 0,1,2,4,5 are used, position 3 is ignored)
    /// Mask for nucleotides:
    ///   Position: 5 4 3 2 1 0
    ///   Pattern:  1 1 0 1 1 1
    ///   Mask:     11 11 00 11 11 11 (2 bits per position)
    ///   Binary:   0b111100111111
    ///
    /// K-mer ACGTAT = 0b00011011000011
    /// Masked:        0b00011000000011 (position 3's "TA" is zeroed out)
    /// ```
    #[inline]
    pub fn apply_spaced_seed(&self, kmer: Kmer, mask: Kmer) -> Kmer {
        kmer & mask
    }

    /// Decode a k-mer back to a DNA string (for debugging/display)
    ///
    /// # Note
    /// This is relatively slow and should only be used for debugging,
    /// not in performance-critical code.
    pub fn decode_kmer(&self, mut kmer: Kmer) -> String {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        let mut result = String::with_capacity(self.k);

        // Extract nucleotides from right to left
        for _ in 0..self.k {
            let code = (kmer & 0b11) as usize;
            result.push(BASES[code]);
            kmer >>= 2;
        }

        // Reverse because we extracted right-to-left
        result.chars().rev().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_nucleotide() {
        assert_eq!(KmerScanner::encode_nucleotide(b'A'), Some(0));
        assert_eq!(KmerScanner::encode_nucleotide(b'C'), Some(1));
        assert_eq!(KmerScanner::encode_nucleotide(b'G'), Some(2));
        assert_eq!(KmerScanner::encode_nucleotide(b'T'), Some(3));
        assert_eq!(KmerScanner::encode_nucleotide(b'N'), None);
    }

    #[test]
    fn test_encode_sequence() {
        let scanner = KmerScanner::new(4).unwrap();
        let kmer = scanner.encode_sequence(b"ACGT").unwrap();
        // A=00, C=01, G=10, T=11 -> 0b00011011 = 27
        assert_eq!(kmer, 0b00011011);
    }

    #[test]
    fn test_reverse_complement() {
        let scanner = KmerScanner::new(4).unwrap();
        let kmer = scanner.encode_sequence(b"ACGT").unwrap();
        let revcomp = scanner.reverse_complement(kmer);
        let expected = scanner.encode_sequence(b"ACGT").unwrap(); // ACGT -> ACGT
        assert_eq!(revcomp, expected);
    }

    #[test]
    fn test_extract_kmers() {
        let scanner = KmerScanner::new(3).unwrap();
        let kmers = scanner.extract_kmers(b"ACGTAC");
        assert_eq!(kmers.len(), 4); // ACG, CGT, GTA, TAC
    }

    #[test]
    fn test_decode_kmer() {
        let scanner = KmerScanner::new(4).unwrap();
        let kmer = scanner.encode_sequence(b"ACGT").unwrap();
        let decoded = scanner.decode_kmer(kmer);
        assert_eq!(decoded, "ACGT");
    }

    #[test]
    fn test_invalid_nucleotide_handling() {
        let scanner = KmerScanner::new(3).unwrap();
        let kmers = scanner.extract_kmers(b"ACNGT");
        // After N at position 2, we have "GT" left which is < k=3
        // So we should get 0 k-mers
        println!("Extracted {} k-mers from ACNGT", kmers.len());
        for (i, kmer) in kmers.iter().enumerate() {
            println!("  K-mer {}: {}", i, scanner.decode_kmer(*kmer));
        }
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_canonical() {
        let scanner = KmerScanner::new(3).unwrap();
        let acg = scanner.encode_sequence(b"ACG").unwrap();
        let cgt = scanner.encode_sequence(b"CGT").unwrap(); // reverse complement of ACG
        let canonical_acg = scanner.canonical(acg);
        let canonical_cgt = scanner.canonical(cgt);
        assert_eq!(canonical_acg, canonical_cgt);
    }
}