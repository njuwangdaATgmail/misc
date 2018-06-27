### harmonic oscillator, phonon, and path integral

Consider a harmonic oscillator
$$
H=\frac12kx^2+\frac1{2m}p^2=\Omega\left(b^\dagger b+\frac12\right)
$$
where $\Omega=\sqrt{k/m}$. By Hamiltonian canonical equation
$$
\dot{x}=\frac{\partial H}{\partial p}=\frac{p}{m}
$$

Its path integral form is
$$
U=\int \mathrm{D}x \exp\left[{\int (p\dot{x}-H)}\right]=\int \mathrm{D}x \exp\left[{\int \left(\frac12m\dot{x}^2-\frac12kx^2\right)}\right]
$$
which is just the Lagrangian form of the path integral. After Wick rotation for finite temperature field theory, we have
$$
Z=\int \mathrm{D}x \exp\left[-{\int \left(\frac12m\dot{x}^2+\frac12kx^2\right)}\right]
$$
from Gaussian integral or Green's function technique, we have
$$
\frac12k\langle x^2\rangle=\frac{T}{2}\sum_{i\omega_n}\frac{\Omega^2}{\omega_n^2+\Omega^2}
$$
The frequency summation can be completed as
$$
\frac12k\langle x^2 \rangle=\frac\Omega4T\sum_{i\omega_n}\left(\frac{1}{i\omega_n+\Omega}-\frac{1}{i\omega_n-\Omega}\right) \\=\frac\Omega4\left[-n_b(-\Omega)+n_b(\Omega)\right]=\frac\Omega2\left[n_b(\Omega)+\frac12\right]
$$

+ Using the path integral method, we have obtained the zero-point motion energy. However, where does the quantum effect come from? (Answer: $\beta\rightarrow\infty$)