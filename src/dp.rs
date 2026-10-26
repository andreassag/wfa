use crate::align::{AlignmentResult, WfaParams};
use crate::cigar::{CigarBuffer, CigarOp};
use crate::seq::PaddedSeq;
use std::cmp::min;

pub fn dp_align(q: &PaddedSeq, t: &PaddedSeq, params: &WfaParams) -> AlignmentResult {
    let q_b = q.as_slice();
    let t_b = t.as_slice();
    let n = q_b.len();
    let m = t_b.len();

    let x = params.x;
    let o = params.o;
    let e = params.e;
    let inf = 1_000_000_000;

    // dp_m[i][j], dp_i[i][j], dp_d[i][j]
    let mut dp_m = vec![vec![inf; m + 1]; n + 1];
    let mut dp_i = vec![vec![inf; m + 1]; n + 1];
    let mut dp_d = vec![vec![inf; m + 1]; n + 1];

    dp_m[0][0] = 0;
    for j in 1..=m {
        dp_d[0][j] = o + (j as i32) * e;
        dp_m[0][j] = dp_d[0][j];
    }
    for i in 1..=n {
        dp_i[i][0] = o + (i as i32) * e;
        dp_m[i][0] = dp_i[i][0];
    }

    for i in 1..=n {
        for j in 1..=m {
            // Insertion in query (vertical)
            dp_i[i][j] = min(dp_m[i - 1][j] + o + e, dp_i[i - 1][j] + e);

            // Deletion in query (horizontal)
            dp_d[i][j] = min(dp_m[i][j - 1] + o + e, dp_d[i][j - 1] + e);

            // Match or Mismatch
            let cost = if q_b[i - 1] == t_b[j - 1] { 0 } else { x };
            let diag = dp_m[i - 1][j - 1] + cost;

            dp_m[i][j] = min(diag, min(dp_i[i][j], dp_d[i][j]));
        }
    }

    let final_score = dp_m[n][m];

    // Backtrace
    let mut cigar = CigarBuffer::new();
    let mut i = n;
    let mut j = m;
    let mut state = 'M';

    while i > 0 || j > 0 {
        if state == 'M' {
            if i > 0 && j > 0 {
                let cost = if q_b[i - 1] == t_b[j - 1] { 0 } else { x };
                if dp_m[i][j] == dp_m[i - 1][j - 1] + cost {
                    if cost == 0 {
                        cigar.push(CigarOp::Match(1));
                    } else {
                        cigar.push(CigarOp::Mismatch(1));
                    }
                    i -= 1;
                    j -= 1;
                    continue;
                }
            }
            if dp_m[i][j] == dp_i[i][j] {
                state = 'I';
            } else if dp_m[i][j] == dp_d[i][j] {
                state = 'D';
            }
        } else if state == 'I' {
            cigar.push(CigarOp::Ins(1));
            if i > 0 && dp_i[i][j] == dp_m[i - 1][j] + o + e {
                state = 'M';
            }
            i -= 1;
        } else if state == 'D' {
            cigar.push(CigarOp::Del(1));
            if j > 0 && dp_d[i][j] == dp_m[i][j - 1] + o + e {
                state = 'M';
            }
            j -= 1;
        }
    }

    cigar.reverse();
    AlignmentResult {
        score: final_score,
        cigar,
    }
}
