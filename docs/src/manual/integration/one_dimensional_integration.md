## One Dimensional Fractal Integration 

### Fractal Interpolation Coefficients
Let $N > 2$ be natural number. Let $x_0 < x_1 < \ldots < x_N$ be real numbers. Let $I = [x_0, x_N]$, and $I_i = [x_i, x_{i + 1}]$ for $i = 0, \ldots, N - 1$ be the subintervals of $I$. Consider homeomorphisms $L_i(x): I \mapsto I_i$ such that
```math
\begin{aligned}
|L_i(x_j) - L_i(x_k)| & \leq s_i |x_j - x_k| \\
L_i(x_0) = x_i, &\quad L_i(x_N) = x_{i + 1}  
\end{aligned}
```
for all $x_j, x_k \in I, \; 0 \leq s_i < 1$. Consider $N$ continuous transformations $F_i : I \times \mathbb{R} \mapsto \mathbb{R}$ such that 
```math
\begin{aligned}
|F_i(x_j, y_j) - F_i(x_j, y_k)| &\leq r_i |y_j - y_k| \\
F_i(x_0, y_0) = y_i &\quad F_i(x_N, y_N) = y_{i + 1} 
\end{aligned}
```
for $x_j \in I, \; y_j, y_k \in \mathbb{R}, \; 0 \leq r_i < 1$. Let $w_i : I \times \mathbb{R} \mapsto I_i \times \mathbb{R}$ such that
```math
w_i(x, y) = (L_i(x), F_i(x, y)), \; i = 1, \ldots, N
```
The followings hold for the IFS $\mathcal{I} = \{ I \times \mathbb{R}; w_i, \; i = 1, \ldots, N \}$:

* It has a unique attractor $G$ which is the graph of a continuous function $f:I\mapsto \mathbb{R}$ satisfying $f(x_i) = y_i$ for $i = 0, \ldots, N$.
* Let $(\mathcal{F}, d)$ be a complete metric space such that $\mathcal{F} = \{ \tilde{f} \in \mathcal{C}.[^1] (I): \tilde{f}(x_0) = y_0, \; \tilde{f}(x_N) = y_N \}$ with the metric $d(\tilde{f}_1, \tilde{f}_2) = max\{ |\tilde{f}_1(x) - \tilde{f}_2(x)| : x \in I \}$. Let $M: \mathcal{F} \mapsto \mathcal{F}$ be a mapping defined by $(M\tilde{f})(x) = F_i(L_i^{-1}(x), \tilde{f}(L_i^{-1}(x)), \; x \in I$ for $i = 0, \ldots, N$. Then, $M$ has unique fixed point $f$ such that $f = \lim_{n \mapsto \infty} M^n(\tilde{f})$ for any $\tilde{f} \in \mathcal{F}$ and $f$ satisfies the condition given above.

The function $f$ whose graph is the attractor of the IFS $\mathcal{I}$ and which is the fixed point of the mapping $M$ given in Theorem above is called a fractal interpolation function (FIF).

### One Dimensional Fractal Integration 
It is also possible to calculate the definite integrals of the fractal interpolation functions over the interpolation domain. The crucial result here is that the value of the integrals does not depend on the explicit formula of the interpolation function, but depends on the coefficients of the transformations that used to construct the interpolant. 

We first start with
```math 
\begin{aligned}
  I =  \int_{I} f(x)dx = \int_{x_0}^{x_N} f(x) dx
\end{aligned}
```
where ``f`` is the fractal interpolation function, i.e. the interpolant. Since ``f(x)`` is the fixed point of ``Tf``, we have ``f(x) = Tf(x)``. Then,
```math
\begin{aligned}
    I
    = \int_{x_0}^{x_N} Tf(x) dx
\end{aligned}
```
Since ``Tf(x) = F_n(L_n^{-1}(x), f(L_n^{-1}))``, where, in case one dimensional integration, 
```math
\begin{aligned}
    L_n(x) &= \alpha_1^n x + \beta_1^n \\
    F_n(x, y) &= \alpha_2^n x + \alpha_3^n y + \beta_2^n \\
\end{aligned}
```
we have,
```math
\begin{aligned}
    I
    = \int_{x_0}^{x_N} \left( \alpha_2^n L_n^{-1}(x) + \alpha_3^n f(L_n^{-1}(x)) + \beta_2^n \right) dx
\end{aligned}
```
We can divide integral ``I`` into multiple subintervasl ``I_n``,
```math
\begin{aligned}
    I  
    = \sum_{n=1}^N \int_{x_{n-1}}^{x_n} \left( \alpha_2^n L_n^{-1}(x) + \alpha_3^n f(L_n^{-1}(x)) + \beta_2^n \right) dx
    = \sum_{n=1}^N \int_{I_n} \left( \alpha_2^n L_n^{-1}(x) + \alpha_3^n f(L_n^{-1}(x)) + \beta_2^n \right) dx
\end{aligned}
```
To make the right hand side of the above equation look like the left hand side, we need to apply a suitable transformation. That is, we need to transform the integral ``\int_{I_n}`` over the subinterval ``I_n`` to the interpolation interval ``I``. 

!!! note

    Consider the transformation ``L`` shown in the figure below where $\Omega_{n} \subset \mathbb{R}^{2}$ and $\Omega \subset \mathbb{R}^{2}$. 
    ![1d_transformations](../../assests/1D_transformations.svg)

    We have the following equality for the definite integrals.
    ```math
    \begin{aligned}
        I
        = \int_{\Omega_n} h(x) dx
        = \int_{\Omega} h(L(\bar{x})) |J_L|d \bar{x} 
    \end{aligned}
    ```
    where $|J_L|$ is the determinant of the Jacobian of the transformation $L$.

Using the note given above, we can write 
```math 
\begin{aligned}
    I &= \sum_{n=1}^N \int_{I} \left( \alpha_2^n \bar{x} + \alpha_3^n f(\bar{x}) + \beta_2^n \right) d \bar{x} \\
        &= \sum_{n=1}^N \int_{I} \left( \alpha_2^n \bar{x} + \beta_2^n \right) d \bar{x}  +
        \sum_{n=1}^N \int_{I} \left( \alpha_3^n f(\bar{x}) \right) d \bar{x} \\
        &= \sum_{n=1}^N \int_{I} \left( \alpha_2^n \bar{x} + \beta_2^n \right) d \bar{x}  +
        \sum_{n=1}^N \left( \alpha_3^n \right) \int_{I} f(\bar{x}) d \bar{x} \\
        &= \sum_{n=1}^N \int_{I} \left( \alpha_2^n \bar{x} + \beta_2^n \right) d \bar{x}  +
        \sum_{n=1}^N \left( \alpha_3^n \right) I \\
\end{aligned}
```
from which we have, 
```math 
    I = \dfrac{K_1}{1-K_2}
   
```
where 
```math
    \begin{aligned}
    K_1 &= \sum_{n=1}^N \int_{I} \left( \alpha_2^n \bar{x} + \beta_2^n \right) d \bar{x} \\
    K_2 &= 1 - \sum_{n=1}^N \alpha_3^n 
    \end{aligned}
```

It worths pointing out that the value $I$ of the integration depends on just the coefficients of the transformations, not on the explicit expression of the interpolation function $f$.
