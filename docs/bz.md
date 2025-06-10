# Momenta in the 1st Brillouin Zone

## Reciprocal Lattice Vectors

Define the matrix of lattice basis vectors and the matrix of reciprocal lattice basis vectors as

$$
\begin{equation}
    A_{ij} = a_{ij}, \quad B_{ij} = b_{ij}
\end{equation}
$$

where $a_{ij}$ ($b_{ij}$) is the $i$th component (under the Cartesian basis) of $a_j$ ($b_j$). By definition

$$
\begin{equation*}
    a_m \cdot b_n = a_{im} b_{in}
    = (A^\T B)_{mn} = 2\pi \delta_{mn}
\end{equation*}
$$

therefore

$$
\begin{equation}
    B = 2\pi (A^\T)^{-1}
\end{equation}
$$

## List of Momenta in 1st Brillouin Zone

A momenta in the 1st Brillouin zone ($d$-dimensional) can be expressed as (for the $i$th component)

$$
k_i = \sum_{j=1}^d \frac{n_j}{N_j} b_{ij}
\quad \text{or} \quad
k = B \frac{n}{N}
$$

- $\{b_j\}_{j=1}^d$ are reciprocal lattice vectors, and $b_{ij}$ is the $i$th component of $b_j$. 
- $\{n_j\}$ are integers for PBC, or half-integers for Anti-PBC.

We can find $n$ from $k$ as

$$
\frac{n}{N} = B^{-1} k
\ \Rightarrow \ n = N * (B^{-1} k)
$$

where $*$ means element-wise multiplication. 

## Change of Basis

For a $d$-dimensional vector space $V$, let $\{e_i\}_{i=1}^d$ be the standard Cartesian basis vectors, and $\{a_i\}_{i=1}^d$ be another set of basis vectors. A vector $v \in V$ can be expressed as (summation over paired indices is implied)

$$
\begin{equation}
    v = e_i v_i(e) = a_j v_j(a)
\end{equation}
$$

Let $a_{ij}$ be the $e_i$-component of $a_j$, i.e.

$$
\begin{equation}
    a_j = e_i a_{ij}
\end{equation}
$$

Let $A$ be the matrix of $a_{ij}$. Then

$$
\begin{equation}
    e_i v_i(e) = e_i a_{ij} v_j(a) 
    \ \Rightarrow \ 
    v(e) = A v(a)
\end{equation}
$$

Or equivalently, 

$$
\begin{equation}
    v(a) = A^{-1} v(e)
\end{equation}
$$

## Using Lattice Basis Vectors

In a $d$-dimensional Bravais lattice (can contain multiple sites in one unit cell), the lattice basis vectors and the reciprocal lattice basis vectors are related by

$$
\begin{equation}
    a_i \cdot b_j = 2\pi \delta_{ij}
\end{equation}
$$

The inner product $k \cdot R$ can be expressed in the components under the lattice basis as

$$
\begin{equation}
    k \cdot R
    = (k_i b_i) \cdot (R_j a_j)
    = 2 \pi k_i R_i
\end{equation}
$$

## From Discrete to Continuous<br>Brillouin Zone

For a periodic lattice with $N_i$ sites along the $i$-th dimension, the momenta in the 1st Brillouin zone are given by (in components under )

$$
\begin{equation}
    k = \left(
        \frac{m_1}{N_1}, \cdots,  
        \frac{m_d}{N_d}
    \right),
    \quad m_i \in \Z \text{ and } 
    -\frac{N_i}{2} \le m_i < \frac{N_i}{2}
\end{equation}
$$

As the number of unit cells $N \to \infty$, the momenta in the 1st Brillouin zone become continuous. The summation over them becomes an integration:

$$
\begin{equation}
    \frac{1}{N} \sum_k \to 
    \int_\mathrm{BZ} d^dk
\end{equation}
$$

where $k$ refers to components of $k$ under the reciprocal lattice basis; the region BZ is

$$
\begin{equation}
    \mathrm{BZ} = \left\{
        (k_1, ..., k_d) \,\middle|\,
        -\frac{1}{2} \le k_i < \frac{1}{2}
    \right\}
\end{equation}
$$
