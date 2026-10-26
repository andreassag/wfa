use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),
    Mismatch(u32),
    Ins(u32),
    Del(u32),
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CigarBuffer {
    pub ops: Vec<CigarOp>,
}

impl CigarBuffer {
    pub fn new() -> Self {
        Self { ops: Vec::new() }
    }

    pub fn push(&mut self, op: CigarOp) {
        if let Some(last) = self.ops.last_mut() {
            match (last, op) {
                (CigarOp::Match(n1), CigarOp::Match(n2)) => { *n1 += n2; return; }
                (CigarOp::Mismatch(n1), CigarOp::Mismatch(n2)) => { *n1 += n2; return; }
                (CigarOp::Ins(n1), CigarOp::Ins(n2)) => { *n1 += n2; return; }
                (CigarOp::Del(n1), CigarOp::Del(n2)) => { *n1 += n2; return; }
                _ => {}
            }
        }
        self.ops.push(op);
    }

    pub fn reverse(&mut self) {
        self.ops.reverse();
    }

    pub fn to_sam_string(&self) -> String {
        let mut s = String::new();
        for op in &self.ops {
            match op {
                CigarOp::Match(len) => s.push_str(&format!("{}M", len)),
                CigarOp::Mismatch(len) => s.push_str(&format!("{}X", len)),
                CigarOp::Ins(len) => s.push_str(&format!("{}I", len)),
                CigarOp::Del(len) => s.push_str(&format!("{}D", len)),
            }
        }
        s
    }

    pub fn score(&self, x: i32, o: i32, e: i32) -> i32 {
        let mut total = 0;
        for op in &self.ops {
            match op {
                CigarOp::Match(_) => {}
                CigarOp::Mismatch(len) => total += (*len as i32) * x,
                CigarOp::Ins(len) => total += o + (*len as i32) * e,
                CigarOp::Del(len) => total += o + (*len as i32) * e,
            }
        }
        total
    }
}

impl fmt::Display for CigarBuffer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_sam_string())
    }
}
