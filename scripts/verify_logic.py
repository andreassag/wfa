#!/usr/bin/env python3
"""
verify_logic.py - Standalone mathematical verification of WFA and Gotoh DP algorithms.
"""

def lcp(q, t, i, j):
    match = 0
    while i + match < len(q) and j + match < len(t) and q[i + match] == t[j + match]:
        match += 1
    return match

def gotoh_dp(q, t, x=4, o=6, e=2):
    n, m = len(q), len(t)
    inf = 10**9
    M = [[inf] * (m + 1) for _ in range(n + 1)]
    I = [[inf] * (m + 1) for _ in range(n + 1)]
    D = [[inf] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    for j in range(1, m + 1):
        D[0][j] = o + j * e
        M[0][j] = D[0][j]
    for i in range(1, n + 1):
        I[i][0] = o + i * e
        M[i][0] = I[i][0]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            I[i][j] = min(M[i - 1][j] + o + e, I[i - 1][j] + e)
            D[i][j] = min(M[i][j - 1] + o + e, D[i][j - 1] + e)
            cost = 0 if q[i - 1] == t[j - 1] else x
            diag = M[i - 1][j - 1] + cost
            M[i][j] = min(diag, I[i][j], D[i][j])

    return M[n][m]

def wfa_solve(q, t, x=4, o=6, e=2):
    n, m = len(q), len(t)
    target_k = m - n
    oe = o + e

    # Wavefront maps: score -> dict of {k: offset}
    M_wf = {}
    I_wf = {}
    D_wf = {}

    init_match = lcp(q, t, 0, 0)
    M_wf[0] = {0: init_match}
    I_wf[0] = {}
    D_wf[0] = {}

    if target_k == 0 and init_match >= n:
        return 0

    for s in range(1, 10000):
        cur_I = {}
        cur_D = {}
        cur_M = {}

        # Collect candidate diagonals
        diags = set()
        if s - oe in M_wf:
            for k in M_wf[s - oe]:
                diags.add(k - 1)
                diags.add(k + 1)
        if s - e in I_wf:
            for k in I_wf[s - e]:
                diags.add(k - 1)
        if s - e in D_wf:
            for k in D_wf[s - e]:
                diags.add(k + 1)
        if s - x in M_wf:
            for k in M_wf[s - x]:
                diags.add(k)

        for k in diags:
            i_val = -1
            if s - oe in M_wf and (k + 1) in M_wf[s - oe]:
                i_val = max(i_val, M_wf[s - oe][k + 1] + 1)
            if s - e in I_wf and (k + 1) in I_wf[s - e]:
                i_val = max(i_val, I_wf[s - e][k + 1] + 1)
            if i_val >= 0:
                cur_I[k] = i_val

            d_val = -1
            if s - oe in M_wf and (k - 1) in M_wf[s - oe]:
                d_val = max(d_val, M_wf[s - oe][k - 1])
            if s - e in D_wf and (k - 1) in D_wf[s - e]:
                d_val = max(d_val, D_wf[s - e][k - 1])
            if d_val >= 0:
                cur_D[k] = d_val

            m_val = -1
            if s - x in M_wf and k in M_wf[s - x]:
                m_val = max(m_val, M_wf[s - x][k] + 1)
            if k in cur_I:
                m_val = max(m_val, cur_I[k])
            if k in cur_D:
                m_val = max(m_val, cur_D[k])

            if m_val >= 0:
                ext = lcp(q, t, m_val, m_val + k)
                final_off = m_val + ext
                cur_M[k] = final_off
                if k == target_k and final_off >= n:
                    return s

        M_wf[s] = cur_M
        I_wf[s] = cur_I
        D_wf[s] = cur_D

    return -1

if __name__ == '__main__':
    # Read sample data
    def read_fa(p):
        with open(p) as f:
            return "".join(l.strip() for l in f if not l.startswith('>'))

    t = read_fa("data/sample_target.fa")
    q = read_fa("data/sample_query.fa")

    dp_score = gotoh_dp(q, t)
    wfa_score = wfa_solve(q, t)

    print(f"Gotoh Quadratic DP Score: {dp_score}")
    print(f"Wavefront Alignment Score: {wfa_score}")
    assert dp_score == wfa_score, "Score mismatch between DP and WFA!"
    print("Verification passed: 100% exact score equivalence!")
