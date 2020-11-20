## Two Dimensional Fractal Integration 

Consider the data points $(x_i, y_i), \; i = 1, 2, \ldots, N$ triangulized over a triangular domain $\Omega \subset \mathbb{R}^n$. Associated with the data point $(xi, yi)$ we consider real values $z_i \in \mathbb{R}, \; i = 1, 2, \ldots, N$. In addition, we consider additional real data values $t_i \; i = 1, 2, \ldots, N$. Our aim is to find and interpolation function, i.e., interpolant $f: \Omega \mapsto \mathbb{R}^2$ such that $(z_i, t_i) = f(x_i, y_i), \; i = 1, 2, \ldots, N$.

Consider the iterated function system $\mathcal{I} = \{\Omega \times \mathbb{R}^n; w_n, n = 1, 2, \ldots, K\}$ where 
```math
w_n(x, y, z, t) = 
\begin{bmatrix} 
    \alpha_1^n & \alpha_2^n & 0 & 0 \\ 
    \alpha_3^n & \alpha_4^n & 0 & 0 \\ 
    \alpha_5^n & \alpha_6^n & \alpha_7^n & \alpha_8^n \\ 
    \alpha_9^n & \alpha_{10}^n & \alpha_{11}^n & \alpha_{12}^n 
\end{bmatrix}
\begin{bmatrix}
    x \\ y \\ z \\ t 
\end{bmatrix}
+ 
\begin{bmatrix} 
    \beta_1^n \\
    \beta_2^n \\
    \beta_3^n \\
    \beta_4^n
\end{bmatrix}
```
We can divide the transformation $w_n, \; n = 1, 2, \ldots, K$ as follows,
```math 
w_n(x, y, z, t) = (L_n(x, y), F_n(x, y, z, t))
```
where 
```math 
L_n(x, y) = 
\begin{bmatrix} 
    \alpha_1^n & \alpha_2^n \\ 
    \alpha_3^n & \alpha_4^n 
\end{bmatrix}
\begin{bmatrix}
    x \\ y 
\end{bmatrix}
+ 
\begin{bmatrix} 
    \beta_1^n \\
    \beta_2^n 
\end{bmatrix}
```
and 
```math 
F_n(x, y, z, t) = 
\begin{bmatrix} 
    \alpha_5^n & \alpha_6^n & \alpha_7^n \alpha_8^n \\
    \alpha_9^n & \alpha_{10}^n & \alpha_{11}^n \alpha_{12}^n \\
\end{bmatrix}
\begin{bmatrix}
    x \\ y \\ z \\ t 
\end{bmatrix}
+ 
\begin{bmatrix}
    \beta_n^3 \\ \beta_n^4
\end{bmatrix}
```
Given the free parameters $\alpha_7^n, \alpha_8^n, \alpha_{10}^n, \alpha_{12}^n, \; n = 1, \ldots, K$, the remaining coefficients of the transformations $w_n, \; n = 1, 2, \ldots, K$ can be found from the boundary conditions, 
```math 
L_n(\tilde{x}_j, \tilde{y}_j) = (x_k, y_k) \quad k = 1, 2, \ldots, K \quad j = 1, 2, 3
```
where the $\{(\tilde{x}_j, \tilde{y}_j), \; j = 1, 2, 3\}$ are the set of the vertex points of the interpolation $\Omega$ and $\{(x_j, y_j), \; j = 1, 2, 3\}$ are the set of vertex points of the triangle $\Omega_k \; k = 1, 2, \ldots, K$ where $K$ is the number of triangles in $\Omega$. Here the free variables $t_n, \; n = 1, 2, \ldots, N$ and the $\alpha_7^n, \alpha_8^n, \alpha_{10}^n, \alpha_{12}^n$ are called the \emph{hidden variables.}

Then the graph of the interpolant $f$ is the attractor of the IFS $\mathcal{I}$. Furthermore, the interpolant $f$ is the fixed of a mapping $Tf: \mathbb{F} \mapsto \mathbb{F}$ where $\mathbb{F}$ is the set of continuous functions $g: \Omega \mapsto \mathbb{R}^2$ such that 
```math 
Tf(x, y) = F_n(L_n^{-1}(x, u), f(L_n^{-1}(x, u)))
```

The definite integral of the interpolant $f$ over the interpolation domain $\Omega$ can be computed by using just the transformation coefficients. To this end, let us start, 
```math 
I = 
\int_\Omega f(x, y) dx dy = 
\begin{bmatrix}
    I_1 \\ 
    I_2 
\end{bmatrix}
```
Since $f$ is the fixed point of the transformation $Tf$, then $f(x,y) = Tf(x, y)$, then 
```math 
I = 
\int_\Omega Tf(x, y) dx dy = \sum_{n=1}^K \int_{\Omega^n} F_n(L_n^{-1}(x, y), f(L_n^{-1}(x, y))) dx dy
```
Using the transformation $(x, y) = L_n(\tilde{x}, \tilde{y})$, we can write 
```math 
\begin{aligned}
I &= 
\begin{bmatrix} 
    I_1 \\ 
    I_2 
\end{bmatrix}
\int_\Omega F_n(\tilde{x}, \tilde{y}, f(\tilde{x}, \tilde{y})) |J_{L_n}| d\tilde{x} d\tilde{y} \\ 
&= 
\begin{bmatrix}
\int_\Omega (\alpha_5^n \tilde{x} + \alpha_6^n \tilde{y} + \alpha_7^n f_1(\tilde{x}, \tilde{y}) + \alpha_8^n f_2(\tilde{x}, \tilde{y}) + \beta_3^n) d\tilde{x} d\tilde{y} \\
\int_\Omega (\alpha_9^n \tilde{x} + \alpha_{10}^n \tilde{y} + \alpha_{11}^n f_1(\tilde{x}, \tilde{y}) + \alpha_{12}^n f_2(\tilde{x}, \tilde{y}) + \beta_4^n) d\tilde{x} d\tilde{y} \\
\end{bmatrix} \\ 
&= 
\begin{bmatrix} 
\Lambda_1 + W_{11} I_1 + W_{12} I_2 \\ 
\Lambda_2 + W_{12} I_1 + W_{22} I_2  
\end{bmatrix}
\end{aligned}
```
where 
```math 
\begin{aligned}
\Lambda_{1} &= \sum_{n = 1}^K \int_\Omega (\alpha_5^n \tilde{x} + \alpha_6^n \tilde{y} + \beta_3^n) |J_{L_n}| d\tilde{x} d\tilde{y}  \\
\Lambda_{2} &= \sum_{n = 1}^K \int_\Omega (\alpha_9^n \tilde{x} + \alpha_{10}^n \tilde{y} + \beta_4^n) |J_{L_n}| d\tilde{x} d\tilde{y}  \\ 
W_{11} &= \sum_{n = 1}^K \alpha_7^n |J_{L_n}|  \\ 
W_{12} &= \sum_{n = 1}^K \alpha_8^n |J_{L_n}|  \\ 
W_{21} &= \sum_{n = 1}^K \alpha_{11}^n |J_{L_n}|  \\ 
W_{12} &= \sum_{n = 1}^K \alpha_{12}^n |J_{L_n}| 
\end{aligned}
```
Hence we can write, 
```math 
\begin{aligned}
\begin{bmatrix} 
1 - W_{11} & -W_{12} \\ 
-W_{12} & 1 - W_{22} 
\end{bmatrix} 
\begin{bmatrix}
    I_1 \\ 
    I_2
\end{bmatrix} 
&= 
\begin{bmatrix}
\Lambda_1 \\ 
\Lambda_2
\end{bmatrix} \\ 
\begin{bmatrix} 
\Gamma_{11} & \Gamma_{12} \\ 
\Gamma_{21} & \Gamma_{22} \\ 
\end{bmatrix}
\begin{bmatrix} 
I_1 \\ I_2
\end{bmatrix}
&= 
\begin{bmatrix} 
\Lambda_1 \\ 
\Lambda_2
\end{bmatrix}
\end{aligned}
```
By solving this system of equations,
```math 
\begin{aligned}
I_1 &= \dfrac{\Lambda_1 \Gamma_{22} - \Lambda_2 \Gamma_{12}}{\Gamma_{11} \Gamma_{22} - \Gamma_{12} \Gamma_{21}} \\ 
I_2 &= \dfrac{\Lambda_2 \Gamma_{11} - \Lambda_1 \Gamma_{21}}{\Gamma_{11} \Gamma_{22} - \Gamma_{12} \Gamma_{21}} 
\end{aligned}
```
We have, 
```math 
\begin{aligned} 
\Lambda_{1} &= \int_\Omega (\alpha_5^n \tilde{x} + \alpha_6^n \tilde{y} + \beta_3^n) |J_{L_n}| d\tilde{x} d\tilde{y} \\ 
\Lambda_{2} &= \int_\Omega (\alpha_9^n \tilde{x} + \alpha_{10}^n \tilde{y} + \beta_4^n) |J_{L_n}| d\tilde{x} d\tilde{y} 
\end{aligned}
```
These integrals are over the triangular integration domain $\Omega$. To ease the evaluation of these integral, we can transform the integral over a triangular domain $\Theta$ consisting of the vertex points $(0,0), (0, 1), (1, 0)$ using a transformation of the form 
```math 
T(\bar{x}, \bar{y}) = 
\begin{bmatrix} 
    \tilde{x} \\ \tilde{y} 
\end{bmatrix} = 
\begin{bmatrix}
a_{11} & a_{12} \\ 
a_{21} & a_{22} 
\end{bmatrix}
\begin{bmatrix} 
    \bar{x} \\ \bar{y} 
\end{bmatrix}
+ 
\begin{bmatrix} 
    b_1 \\ b_2 
\end{bmatrix}
```
Then we have, 
```math 
\begin{aligned} 
\Lambda_{11} 
    &= \sum_{n = 1}^K  \int_\Delta \left( \alpha_5^n (a_1 \bar{x} + a_2 \bar{y} + b_1) + \alpha_6^n (a_3 \bar{x} + a_4 \bar{y} + b_2) + \beta_3^n  \right) |J_{L_n}| |J_T| d\bar{x} d\bar{y}
\end{aligned}
```
Since, 
```math 
\begin{aligned} 
\int_\Delta (a_1 \bar{x} + a_2 \bar{y} + b_1) d\bar d\bar{x} \\
    &= \int_{0}^{1} \int_{0}^{-\bar{x} + 1}  (a_1 \bar{x} + a_2 \bar{y} + b_1) d\bar d\bar{x}  \\ 
    &= \int_{0}^{1} \left( a_1 \bar{x} (-\bar{x} + 1}) + a_2 \dfrac{(-\bar{x} + 1)^2}{2} + b_1(-\bar{x} + 1) \right)d\bar{x} \
    &= \dfrac{1}{6} (a_1 + a_2 + 3 b_1) 
\end{aligned}
```
we can write $\Lambda_{1}$
```math 
\begin{aligned} 
\Lambda_{1} = |J_T| \left( 
    \dfrac{a_1 + a_2 + 3 b_1}{6} \sum_{n = 1}^K \alpha_5^n |J_{L_n}| + 
    \dfrac{a_3 + a_4 + 3 b_2}{6} \sum_{n = 1}^K \alpha_6^n |J_{L_n}| + 
    \sum_{n = 1}^K \beta_3^n |J_{L_n}|  
    \right) 
\end{aligned}
```
and 
```math 
\begin{aligned} 
\Lambda_{2} = |J_T| \left( 
    \dfrac{a_1 + a_2 + 3 b_1}{6} \sum_{n = 1}^K \alpha_9^n |J_{L_n}| + 
    \dfrac{a_3 + a_4 + 3 b_2}{6} \sum_{n = 1}^K \alpha_{10}^n |J_{L_n}| + 
    \sum_{n = 1}^K \beta_4^n |J_{L_n}|  
    \right) 
\end{aligned}
```
Lastly, we have 
```math 
    |J_{L_n}| = |\alpha_1^n \alpha_4^n - \alpha_2^n \alpha_3^n|
```
Thus, the definite integral of the interpolant $f$ can be evaluated using the just the coefficients of the transformations. 