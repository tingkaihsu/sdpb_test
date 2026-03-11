# Polynomial Matrix Programs (PMP)
Consider a **collection** of *symmetric* polynomial matrices $P_{j}^{n}(x)$:
$$
    M_{n}^{j}(x) = 
    \begin{pmatrix}
        P^n_{j, 11}(x) & \cdots & P^{n}_{j,1m_j}(x)\\
        \cdots & \cdots & \cdots\\
        P^{n}_{j, m_j1}(x) & \cdots & P^{n}_{j,m_jm_j}(x)
    \end{pmatrix}\,,
$$

# Input to SDPB

`PMP` must first be converted, using `pmp2sdp` program, into an internal `SDP` format that `SDPB` can quickly load.