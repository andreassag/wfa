#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Wavefront {
    pub lo: i32,
    pub hi: i32,
    pub slice_start: usize,
}

impl Wavefront {
    pub const NULL: Self = Self {
        lo: 0,
        hi: -1,
        slice_start: 0,
    };

    #[inline(always)]
    pub fn is_null(&self) -> bool {
        self.lo > self.hi
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        if self.is_null() {
            0
        } else {
            (self.hi - self.lo + 1) as usize
        }
    }
}

pub struct WavefrontArena {
    slab: Vec<i32>,
}

impl WavefrontArena {
    pub fn new() -> Self {
        Self {
            slab: Vec::with_capacity(1024 * 1024),
        }
    }

    pub fn with_capacity(cap: usize) -> Self {
        Self {
            slab: Vec::with_capacity(cap),
        }
    }

    pub fn allocate(&mut self, lo: i32, hi: i32, init_val: i32) -> Wavefront {
        if lo > hi {
            return Wavefront::NULL;
        }
        let len = (hi - lo + 1) as usize;
        let slice_start = self.slab.len();
        self.slab.resize(slice_start + len, init_val);
        Wavefront { lo, hi, slice_start }
    }

    #[inline(always)]
    pub fn get(&self, wf: &Wavefront, k: i32) -> i32 {
        if k < wf.lo || k > wf.hi {
            i32::MIN / 2
        } else {
            let idx = wf.slice_start + (k - wf.lo) as usize;
            self.slab[idx]
        }
    }

    #[inline(always)]
    pub fn set(&mut self, wf: &Wavefront, k: i32, val: i32) {
        if k >= wf.lo && k <= wf.hi {
            let idx = wf.slice_start + (k - wf.lo) as usize;
            self.slab[idx] = val;
        }
    }

    pub fn clear(&mut self) {
        self.slab.clear();
    }
}
