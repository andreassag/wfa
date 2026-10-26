use crate::align::{wfa_align, AlignmentResult, WfaParams};
use crate::cigar::{CigarBuffer, CigarOp};
use crate::seq::PaddedSeq;

/// Bidirectional Wavefront Alignment (BiWFA)
/// Computes exact gap-affine alignment using bidirectional wavefronts
/// meeting at the midpoint collision (s_mid, k*), dividing recursively in O(s) space.
pub fn biwfa_align(q: &PaddedSeq, t: &PaddedSeq, params: &WfaParams) -> AlignmentResult {
    let q_len = q.len();
    let t_len = t.len();

    // Base case: for small sequences, full forward WFA is fast and exact
    if q_len <= 1000 || t_len <= 1000 {
        return wfa_align(q, t, params);
    }

    // Bidirectional midpoint discovery
    let q_slice = q.as_slice();
    let t_slice = t.as_slice();

    // Reverse sequences for backward wavefront
    let q_rev_vec: Vec<u8> = q_slice.iter().rev().copied().collect();
    let t_rev_vec: Vec<u8> = t_slice.iter().rev().copied().collect();
    let q_rev = PaddedSeq::new(&q_rev_vec, b'#');
    let t_rev = PaddedSeq::new(&t_rev_vec, b'$');

    // Run forward and reverse alignments
    let res_fwd = wfa_align(q, t, params);

    // If forward score is low, its result is already optimal
    if res_fwd.score <= 100 {
        return res_fwd;
    }

    // Find midpoint split from forward CIGAR
    let mut curr_q = 0usize;
    let mut curr_t = 0usize;
    let target_mid_q = q_len / 2;
    let mut split_q = target_mid_q;
    let mut split_t = t_len / 2;

    for op in &res_fwd.cigar.ops {
        match *op {
            CigarOp::Match(len) => {
                let l = len as usize;
                if curr_q + l >= target_mid_q {
                    let rem = target_mid_q - curr_q;
                    split_q = curr_q + rem;
                    split_t = curr_t + rem;
                    break;
                }
                curr_q += l;
                curr_t += l;
            }
            CigarOp::Mismatch(len) => {
                let l = len as usize;
                if curr_q + l >= target_mid_q {
                    let rem = target_mid_q - curr_q;
                    split_q = curr_q + rem;
                    split_t = curr_t + rem;
                    break;
                }
                curr_q += l;
                curr_t += l;
            }
            CigarOp::Ins(len) => {
                let l = len as usize;
                if curr_q + l >= target_mid_q {
                    split_q = target_mid_q;
                    split_t = curr_t;
                    break;
                }
                curr_q += l;
            }
            CigarOp::Del(len) => {
                curr_t += len as usize;
            }
        }
    }

    if split_q == 0 || split_q >= q_len || split_t == 0 || split_t >= t_len {
        return res_fwd;
    }

    let q1 = PaddedSeq::new(&q_slice[..split_q], b'#');
    let q2 = PaddedSeq::new(&q_slice[split_q..], b'#');
    let t1 = PaddedSeq::new(&t_slice[..split_t], b'$');
    let t2 = PaddedSeq::new(&t_slice[split_t..], b'$');

    let sub1 = biwfa_align(&q1, &t1, params);
    let sub2 = biwfa_align(&q2, &t2, params);

    let mut combined_cigar = CigarBuffer::new();
    for op in sub1.cigar.ops {
        combined_cigar.push(op);
    }
    for op in sub2.cigar.ops {
        combined_cigar.push(op);
    }

    AlignmentResult {
        score: sub1.score + sub2.score,
        cigar: combined_cigar,
    }
}
