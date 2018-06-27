## antiunitary symmetries for fermions

Following (*Ryu, 2010; Chiu, 2016*), by definition
$$
U_T^\dagger H U_T = H^*, \quad U_T^* U_T=\pm1
$$

$$
U_C^\dagger H U_C = -H^T, \quad U_C^* U_C=\pm1
$$

Combine these two, we have
$$
U_T^\dagger U_C^\dagger H U_C U_T =U_S^\dagger H U_S= -H^\dagger
$$

* In the eigenbasis of $U_S$, we have
  $$
  H=\left[\begin{matrix}0&h\\h^\dagger&0\end{matrix}\right]
  $$
  furthermore, if $h$ is real, we get (?)
  $$
  \left[\begin{matrix} I & 0 \\ 0 & -I \end{matrix}\right] H  \left[\begin{matrix} I & 0 \\ 0 & -I \end{matrix}\right] = -H^T
  $$
  which belongs to BDI class and is just a special case of sign problem free condition as shown in (*Lei Wang, 2015*).

* Ten-fold classes: $\mathcal{T}^2=0,\pm1; \mathcal{C}^2=0,\pm1; \mathcal{S}^2=0,1$

| $(\mathcal{T}^2, \mathcal{C}^2, \mathcal{S}^2)$ |     class      |               QMC examples               |
| :---------------------------------------------: | :------------: | :--------------------------------------: |
|                    (-1,0,0)                     |      AII       |     doped negative-U Hubbard model?      |
|                    (-1,1,1)                     |      DIII      |     doped negative-U Hubbard model?      |
|                    (-1,-1,1)                    |      CII       |                                          |
|                     (1,1,1)                     |      BDI       |  spinless-V model on bipartite lattice   |
|                    (0,-1,0)?                    |       C        | PH-symmetric Hubbard model (complex HS)? |
|                    (1,-1,1)?                    |       CI       |  PH-symmetric Hubbard model (real HS)?   |
|       (0,0,0), (0,0,1), (1,0,0), (0,1,0)        | A, AIII, AI, D |                                          |

(FAILED TRIES: too naive...)

## Majonara representation

It seems using Majorana representation, one can unify all known examples without sign problem: (See *Wei, 2016; Li, 2016; Wang, 2015; Huffman, 2014*) Majoranal class (Majorana reflection positivity) + Kramers class (Majorana Kramers positivity).

In the Majorana basis $[\gamma^1, \gamma^2]$, the noninteracting hamiltonian is
$$
H=\left[ \begin{matrix} A_1 & B \\ -B^T & A_2 \end{matrix} \right]
$$
where $A_i^T=-A_i$. Imposing a time reversal symmetry $\mathcal{T}^+=\tau_1\mathcal{K}$, $\tau_1 H \tau_1 = H^*$, we get
$$
A_2=A_1^*, \quad B=-B^\dagger
$$
which is the Majorana reflection positive condition EXCEPT the requirement of positive(or negative) semidefinite B. 

1. How to understand the missing condition? 
2. Where is the missing condision in Hong Yao's work? Is it a symmetry? Usually, a symmetry gives a degeneracy. But here, we may have no degeneracy.
3. Does some kind of chiral(sublattice) symmetry play any role in the positive condition?
