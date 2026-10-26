use crate::arena::{Wavefront, WavefrontArena};
use crate::cigar::{CigarBuffer, CigarOp};
use crate::seq::PaddedSeq;
use std::cmp::max;

#[derive(Clone, Copy, Debug)]
pub struct WfaParams {
    pub x: i32, // Mismatch penalty
    pub o: i32, // Gap open penalty
    pub e: i32, // Gap extend penalty
}

impl Default for WfaParams {
    fn default() -> Self {
        Self { x: 4, o: 6, e: 2 }
    }
}

#[inline(always)]
pub fn lcp_extend(q: &[u8], t: &[u8], mut i: usize, mut j: usize) -> (usize, usize) {
    while i + 8 <= q.len() && j + 8 <= t.len() {
        let qw = u64::from_le_bytes(q[i..i + 8].try_into().unwrap());
        let tw = u64::from_le_bytes(t[j..j + 8].try_into().unwrap());
        let xor = qw ^ tw;
        if xor != 0 {
            let match_bytes = (xor.trailing_zeros() >> 3) as usize;
            return (i + match_bytes, j + match_bytes);
        }
        i += 8;
        j += 8;
    }
    while i < q.len() && j < t.len() && q[i] == t[j] {
        i += 1;
        j += 1;
    }
    (i, j)
}

#[derive(Debug)]
pub struct AlignmentResult {
    pub score: i32,
    pub cigar: CigarBuffer,
}

pub fn wfa_align(q: &PaddedSeq, t: &PaddedSeq, params: &WfaParams) -> AlignmentResult {
    let q_len = q.len();
    let t_len = t.len();
    let target_k = (t_len as i32) - (q_len as i32);

    let mut arena = WavefrontArena::new();
    let mut m_wf: Vec<Wavefront> = Vec::new();
    let mut i_wf: Vec<Wavefront> = Vec::new();
    let mut d_wf: Vec<Wavefront> = Vec::new();

    // Score s = 0: diagonal 0 extended from (0, 0)
    let (init_q, _) = lcp_extend(q.as_slice(), t.as_slice(), 0, 0);
    let wf0 = arena.allocate(0, 0, init_q as i32);
    m_wf.push(wf0);
    i_wf.push(Wavefront::NULL);
    d_wf.push(Wavefront::NULL);

    if target_k == 0 && init_q == q_len {
        let mut cigar = CigarBuffer::new();
        if q_len > 0 {
            cigar.push(CigarOp::Match(q_len as u32));
        }
        return AlignmentResult { score: 0, cigar };
    }

    let x = params.x as usize;
    let o = params.o as usize;
    let e = params.e as usize;
    let oe = o + e;

    let mut final_score = 0;

    for s in 1.. {
        // Compute lo and hi bounds for diagonal k
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;

        if s >= oe && !m_wf[s - oe].is_null() {
            lo = lo.min(m_wf[s - oe].lo - 1);
            hi = hi.max(m_wf[s - oe].hi + 1);
        }
        if s >= e && !i_wf[s - e].is_null() {
            lo = lo.min(i_wf[s - e].lo - 1);
            hi = hi.max(i_wf[s - e].hi - 1);
        }
        if s >= e && !d_wf[s - e].is_null() {
            lo = lo.min(d_wf[s - e].lo + 1);
            hi = hi.max(d_wf[s - e].hi + 1);
        }
        if s >= x && !m_wf[s - x].is_null() {
            lo = lo.min(m_wf[s - x].lo);
            hi = hi.max(m_wf[s - x].hi);
        }

        if lo > hi {
            m_wf.push(Wavefront::NULL);
            i_wf.push(Wavefront::NULL);
            d_wf.push(Wavefront::NULL);
            continue;
        }

        let cur_i = arena.allocate(lo, hi, i32::MIN / 2);
        let cur_d = arena.allocate(lo, hi, i32::MIN / 2);
        let cur_m = arena.allocate(lo, hi, i32::MIN / 2);

        // Compute I_s[k] and D_s[k]
        for k in lo..=hi {
            let mut i_val = i32::MIN / 2;
            if s >= oe && !m_wf[s - oe].is_null() {
                let prev_m = arena.get(&m_wf[s - oe], k + 1);
                if prev_m >= 0 { i_val = max(i_val, prev_m + 1); }
            }
            if s >= e && !i_wf[s - e].is_null() {
                let prev_i = arena.get(&i_wf[s - e], k + 1);
                if prev_i >= 0 { i_val = max(i_val, prev_i + 1); }
            }
            arena.set(&cur_i, k, i_val);

            let mut d_val = i32::MIN / 2;
            if s >= oe && !m_wf[s - oe].is_null() {
                let prev_m = arena.get(&m_wf[s - oe], k - 1);
                if prev_m >= 0 { d_val = max(d_val, prev_m); }
            }
            if s >= e && !d_wf[s - e].is_null() {
                let prev_d = arena.get(&d_wf[s - e], k - 1);
                if prev_d >= 0 { d_val = max(d_val, prev_d); }
            }
            arena.set(&cur_d, k, d_val);
        }

        // Compute M_s[k] and LCP extend
        let mut reached_end = false;

        for k in lo..=hi {
            let mut m_val = i32::MIN / 2;
            if s >= x && !m_wf[s - x].is_null() {
                let prev_m = arena.get(&m_wf[s - x], k);
                if prev_m >= 0 { m_val = max(m_val, prev_m + 1); }
            }
            let i_val = arena.get(&cur_i, k);
            let d_val = arena.get(&cur_d, k);
            m_val = max(m_val, max(i_val, d_val));

            if m_val >= 0 {
                let q_pos = m_val as usize;
                let t_pos_i64 = (m_val as i64) + (k as i64);
                if t_pos_i64 >= 0 && q_pos <= q_len && (t_pos_i64 as usize) <= t_len {
                    let (ext_q, _) = lcp_extend(q.as_slice(), t.as_slice(), q_pos, t_pos_i64 as usize);
                    arena.set(&cur_m, k, ext_q as i32);

                    if k == target_k && ext_q >= q_len {
                        reached_end = true;
                    }
                } else {
                    arena.set(&cur_m, k, m_val);
                }
            } else {
                arena.set(&cur_m, k, i32::MIN / 2);
            }
        }

        m_wf.push(cur_m);
        i_wf.push(cur_i);
        d_wf.push(cur_d);

        if reached_end {
            final_score = s as i32;
            break;
        }
    }

    // Backtrace to construct CIGAR
    let mut cigar = CigarBuffer::new();
    let mut curr_s = final_score as usize;
    let mut curr_k = target_k;
    let mut curr_state = 'M'; // 'M', 'I', 'D'

    while curr_s > 0 || curr_k != 0 {
        if curr_state == 'M' {
            let off = arena.get(&m_wf[curr_s], curr_k);
            if curr_s == 0 {
                if off > 0 {
                    cigar.push(CigarOp::Match(off as u32));
                }
                break;
            }

            // Find how many matches we extended along diagonal
            let mut base_off = i32::MIN / 2;
            let mut prev_source = 'X'; // 'X' (mismatch), 'I', 'D'

            if curr_s >= x && !m_wf[curr_s - x].is_null() {
                let v = arena.get(&m_wf[curr_s - x], curr_k);
                if v >= 0 && v + 1 > base_off {
                    base_off = v + 1;
                    prev_source = 'X';
                }
            }
            let i_val = arena.get(&i_wf[curr_s], curr_k);
            if i_val > base_off {
                base_off = i_val;
                prev_source = 'I';
            }
            let d_val = arena.get(&d_wf[curr_s], curr_k);
            if d_val > base_off {
                base_off = d_val;
                prev_source = 'D';
            }

            let matches = (off - base_off).max(0) as u32;
            if matches > 0 {
                cigar.push(CigarOp::Match(matches));
            }

            match prev_source {
                'X' => {
                    cigar.push(CigarOp::Mismatch(1));
                    curr_s = curr_s.saturating_sub(x);
                }
                'I' => {
                    curr_state = 'I';
                }
                'D' => {
                    curr_state = 'D';
                }
                _ => break,
            }
        } else if curr_state == 'I' {
            cigar.push(CigarOp::Ins(1));
            // Came from M_(s - oe)[k + 1] or I_(s - e)[k + 1]
            let mut from_m = false;
            if curr_s >= oe && !m_wf[curr_s - oe].is_null() {
                let m_val = arena.get(&m_wf[curr_s - oe], curr_k + 1);
                let cur_i_val = arena.get(&i_wf[curr_s], curr_k);
                if m_val >= 0 && m_val + 1 == cur_i_val {
                    from_m = true;
                }
            }

            if from_m {
                curr_s = curr_s.saturating_sub(oe);
                curr_k += 1;
                curr_state = 'M';
            } else {
                curr_s = curr_s.saturating_sub(e);
                curr_k += 1;
            }
        } else if curr_state == 'D' {
            cigar.push(CigarOp::Del(1));
            let mut from_m = false;
            if curr_s >= oe && !m_wf[curr_s - oe].is_null() {
                let m_val = arena.get(&m_wf[curr_s - oe], curr_k - 1);
                let cur_d_val = arena.get(&d_wf[curr_s], curr_k);
                if m_val >= 0 && m_val == cur_d_val {
                    from_m = true;
                }
            }

            if from_m {
                curr_s = curr_s.saturating_sub(oe);
                curr_k -= 1;
                curr_state = 'M';
            } else {
                curr_s = curr_s.saturating_sub(e);
                curr_k -= 1;
            }
        }
    }

    if curr_s == 0 && curr_k == 0 && !m_wf.is_empty() {
        let init_m = arena.get(&m_wf[0], 0);
        if init_m > 0 && cigar.ops.is_empty() {
            cigar.push(CigarOp::Match(init_m as u32));
        }
    }

    cigar.reverse();
    AlignmentResult {
        score: final_score,
        cigar,
    }
}
