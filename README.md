# Periodic Microstructure Model

Microstructures can repeat themselves in a structure in a periodic manner or in an aperiodic manner. Some analytical solutions exist in the literature for periodic microstructures. For instance, [Barbero et al.](#references) developed explicit formulas for a composite reinforced by isotropic, circular-cylindrical fibers, which are periodically arranged in a squared array (see Figure in the paper).

In the following text, we always assume that the fibers are aligned with the x-axis. These fibers are equally spaced with a distance of $2a_2$ in the direction of the y-axis and with a distance of $2a_3$ in the direction of the z-axis. For the case of square symmetry, we assume $2a_2 = 2a_3$, leading to 6 unique components for the stiffness tensor.

Let $λ_m$ and $μ_m$ be the Lamé parameters for the matrix, and $λ_f$ and $μ_f$ be the parameters for the fiber material. Let the volume fraction of the fiber be denoted by $V_f$. The coefficients are written as given by Barbero et al.:

```math
\begin{align}
    {C^*}_{1111} &= \lambda_m + 2 \mu_m - \frac{V_f}{D} 
                                                \left( 
                                                \frac{S^2_3}{\mu_m^2}
                                                - \frac{ 2 S_6 S_3}{\mu_m^2 g}
                                                - \frac{a S_3}{\mu_m c}
                                                + \frac{S^2_6 - S^2_7}{\mu_m^2 g^2}
                                                + \frac{ a S_6 + b S_7}{c  \mu_m  g}
                                                + \frac{ a^2 - b^2}{4 c^2}
                                                \right)
    \\
    {C^*}_{1122} &= \lambda_m + \frac{V_f b}{D} 
                                                \left( 
                                                \frac{S_3}{2 c \mu_m}
                                                -\frac{S_6 - S_7}{2 c \mu_m g}
                                                -\frac{a + b}{4 c^2}
                                                \right)
    \\
    {C^*}_{2233} &= \lambda_m + \frac{V_f }{D} 
                                                \left( 
                                                \frac{a S_7}{2 c \mu_m g }
                                                -\frac{b a + b^2 }{4 c^2}
                                                \right)
    \\
    {C^*}_{2222} &= \lambda_m + 2 \mu_m - \frac{V_f}{D} 
                                                \left( 
                                                - \frac{a S_3}{2\mu_m c}
                                                + \frac{ a S_6 }{2 c \mu_m  g}
                                                + \frac{ a^2 - b^2}{4 c^2}
                                                \right)
    \\
    {C^*}_{1212} &= \mu_m - V_f 
                                                \left( 
                                                - \frac{2 S_3}{\mu_m }
                                                + \frac{ 1 }{\mu_m -\mu_f}
                                                + \frac{ 2 S_7}{\mu_m \left(1-\nu_m\right)}
                                                \right)^{-1}
    \\
    {C^*}_{1313} &= \mu_m - V_f 
                                                \left( 
                                                - \frac{S_3}{\mu_m }
                                                + \frac{ 1 }{\mu_m -\mu_f}
                                                \right)^{-1}
\end{align}
```

where,

```math
\begin{equation*}
    D = \frac{a S_3^2}{2  \mu_m^2 c}
        - \frac{a S_6 S_3^2}{c  \mu_m^2 g}
        + \frac{a \left(S_6^2 S_7^2 \right)}{2 c \mu_m^2 g^2}
        + \frac{S_3 \left( b^2 - a^2\right)}{2  \mu_m c^2}
        + \frac{S_6 \left( a^2 - b^2\right) + S_7 \left( ab + b^2\right)}{2 c^2 \mu_m^2 g}
        + \frac{a^3 - 2 b^3 - 3 a b^2}{8 c^3}
\end{equation*}
```

and,
```math
\begin{align*}
    a &= \mu_f - \mu_m - 2 \mu_f \nu_m + 2 \mu_m \nu_f
    \\
    b &= - \mu_m  \nu_m + \mu_f  \nu_f + 2  \mu_m  \nu_m  \nu_f - 2  \mu_f  \nu_m  \nu_f
    \\
    c &= \left(\mu_m - \mu_f\right)  \left(    \mu_f - \mu_m 
            + \mu_f  \nu_f 
            - \mu_m  \nu_m 
            + 2  \mu_m  \nu_f 
            - 2  \mu_f  \nu_m 
            + 2  \mu_m  \nu_m  \nu_f  
            - 2  \mu_f  \nu_m  \nu_f
            \right) 
    \\
    g &= 2 - 2  \nu_m
\end{align*}
```

The constants $S_3, S_6$, and $S_7$ for the composite reinforced by long circular cylindrical fibers, periodically arranged in a square array, aligned with the x-axis with $a_2 = a_3$ are given in [Barbero et al.]:

```math
\begin{align}
S_3 &= 0.49247 - 0.47603 * V_f - 0.02748 * V_f^2
\\
S_6 &= 0.36844 - 0.14944 * V_f - 0.27152 * V_f^2
\\
S_7 &= 0.12346 - 0.32035 * V_f + 0.23517 * V_f^2
\end{align}
```

The tensor $C^*$ needs only 6 components because of the square symmetry. This tensor can be transformed into a more general result that would be applicable for transversely isotropic composite material.

Let $T$ be the coordinate transformation matrix, and rotation $θ$ about the x-axis yields:

```math
B(\theta) = T^T(\theta) C^* T(\theta)
```

Then the equivalent transversely isotropic tensor is obtained by averaging:

```math
\bar{B} = (1 / π) \int_0^\pi B (\theta) dθ
```

Then, using the relations between the engineering constants and the components of $\bar{B}$, the following engineering constants can be obtained:

```math
\begin{align*}
E_1 &= {C^*}_{1111} - \frac{2  {C^*}_{1122}}{ ({C^*}_{2222} + {C^*}_{2233})}
\\
\end{align*}
```

## References
- [Barbero et al.](#) (EJ Barbero and R Luciano. Micromechanical formulas for the relaxation tensor of linear viscoelastic
composites with transversely isotropic fibers. International Journal of Solids and Structures, 32(13):
1859–1872, 1995.)

