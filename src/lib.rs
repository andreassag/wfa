pub mod seq;
pub mod arena;
pub mod cigar;
pub mod align;
pub mod dp;
pub mod biwfa;

pub use seq::{PaddedSeq, read_fasta};
pub use cigar::{CigarBuffer, CigarOp};
pub use align::{wfa_align, lcp_extend, AlignmentResult, WfaParams};
pub use dp::dp_align;
pub use biwfa::biwfa_align;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lcp_branchless_extend() {
        let q = b"ACGTACGTTAGCTAGCTTTTTTTT";
        let t = b"ACGTACGTTAGCTAGCXXXXXXXX";
        let (eq, et) = lcp_extend(q, t, 0, 0);
        assert_eq!(eq, 16);
        assert_eq!(et, 16);
    }

    #[test]
    fn test_exact_match() {
        let q = PaddedSeq::new(b"ACGTACGTACGT", b'#');
        let t = PaddedSeq::new(b"ACGTACGTACGT", b'$');
        let params = WfaParams::default();

        let res = wfa_align(&q, &t, &params);
        assert_eq!(res.score, 0);
        assert_eq!(res.cigar.to_sam_string(), "12M");
    }

    #[test]
    fn test_single_mismatch() {
        let q = PaddedSeq::new(b"ACGTACGT", b'#');
        let t = PaddedSeq::new(b"ACGTTCGT", b'$');
        let params = WfaParams::default();

        let res_wfa = wfa_align(&q, &t, &params);
        let res_dp = dp_align(&q, &t, &params);

        assert_eq!(res_wfa.score, params.x);
        assert_eq!(res_wfa.score, res_dp.score);
    }

    #[test]
    fn test_single_insertion() {
        let q = PaddedSeq::new(b"ACGTACGT", b'#');
        let t = PaddedSeq::new(b"ACGTCGT", b'$');
        let params = WfaParams::default();

        let res_wfa = wfa_align(&q, &t, &params);
        let res_dp = dp_align(&q, &t, &params);

        assert_eq!(res_wfa.score, params.o + params.e);
        assert_eq!(res_wfa.score, res_dp.score);
    }

    #[test]
    fn test_complex_alignment_dp_equivalence() {
        let q = PaddedSeq::new(b"ACGTGTCGTAGCTAGCTAACGT", b'#');
        let t = PaddedSeq::new(b"ACGTCTAGCTACGTAGCTAACGT", b'$');
        let params = WfaParams { x: 4, o: 6, e: 2 };

        let res_wfa = wfa_align(&q, &t, &params);
        let res_dp = dp_align(&q, &t, &params);

        assert_eq!(res_wfa.score, res_dp.score);
        assert_eq!(res_wfa.cigar.score(params.x, params.o, params.e), res_dp.score);
    }
}
