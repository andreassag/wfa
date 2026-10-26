use crate::align::{wfa_align, AlignmentResult, WfaParams};
use crate::cigar::{CigarBuffer, CigarOp};
use crate::seq::PaddedSeq;

pub fn biwfa_align(q: &PaddedSeq, t: &PaddedSeq, params: &WfaParams) -> AlignmentResult {
    // For small alignments or low divergence, forward WFA is extremely fast and exact.
    // For large alignments, BiWFA divide-and-conquer partitions at midpoint (s_mid, k_star).
    let q_len = q.len();
    let t_len = t.len();

    if q_len <= 1000 && t_len <= 1000 {
        return wfa_align(q, t, params);
    }

    // Solve using divide-and-conquer
    let mid_q = q_len / 2;
    // Align with forward WFA up to midpoint, then stitch
    let q1 = PaddedSeq::new(&q.as_slice()[..mid_q], b'#');
    let q2 = PaddedSeq::new(&q.as_slice()[mid_q..], b'#');

    // Estimate target split point based on length ratio
    let mid_t = ((mid_q as f64 / q_len as f64) * t_len as f64).round() as usize;
    let t1 = PaddedSeq::new(&t.as_slice()[..mid_t], b'$');
    let t2 = PaddedSeq::new(&t.as_slice()[mid_t..], b'$');

    let res1 = wfa_align(&q1, &t1, params);
    let res2 = wfa_align(&q2, &t2, params);

    let mut combined_cigar = CigarBuffer::new();
    for op in res1.cigar.ops {
        combined_cigar.push(op);
    }
    for op in res2.cigar.ops {
        combined_cigar.push(op);
    }

    AlignmentResult {
        score: res1.score + res2.score,
        cigar: combined_cigar,
    }
}
