This page includes notes about one dimensional fractal interpolation and one dimensional hidden fractal interpolation.

## Preliminaries

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
Then following is the main theorem of fractal interpolation functions put forward by Barnsley.

The followings hold for the IFS $\mathcal{I} = \{ I \times \mathbb{R}; w_i, \; i = 1, \ldots, N \}$:

* It has a unique attractor $G$ which is the graph of a continuous function $f:I\mapsto \mathbb{R}$ satisfying $f(x_i) = y_i$ for $i = 0, \ldots, N$.
* Let $(\mathcal{F}, d)$ be a complete metric space such that $\mathcal{F} = \{ \tilde{f} \in \mathcal{C}.[^1] (I): \tilde{f}(x_0) = y_0, \; \tilde{f}(x_N) = y_N \}$ with the metric $d(\tilde{f}_1, \tilde{f}_2) = max\{ |\tilde{f}_1(x) - \tilde{f}_2(x)| : x \in I \}$. Let $M: \mathcal{F} \mapsto \mathcal{F}$ be a mapping defined by $(M\tilde{f})(x) = F_i(L_i^{-1}(x), \tilde{f}(L_i^{-1}(x)), \; x \in I$ for $i = 0, \ldots, N$. Then, $M$ has unique fixed point $f$ such that $f = \lim_{n \mapsto \infty} M^n(\tilde{f})$ for any $\tilde{f} \in \mathcal{F}$ and $f$ satisfies the condition given above.

The function $f$ whose graph is the attractor of the IFS $\mathcal{I}$ and which is the fixed point of the mapping $M$ given in Theorem above is called a fractal interpolation function (FIF).

## 1D Fractal Interpolation 

Consider the affine IFS $\mathcal{I} = \{ I \times \mathbb{R}; w_i, \; i = 1, \ldots, N \}$ with,
```math 
\begin{aligned}
L_i(x) &= a_i x + e_i \\
F_i(x, y) &= c_i x + d_i y + f_i
\end{aligned}
\quad i = 0, \ldots, N - 1
```
where $d_i, \; i = 1, \ldots, N$ are free parameters, also called as vertical scaling factors, such that $|d_i| < 1$. 

Given the vertical scaling factors $d_i, \; i = 1, \ldots, N$, the coefficients $a_i$ and $e_i$ of the affine transformation $L_i$, and $c_i$ and $f_i$ of the affine transformation $F_i$ are determined to satisfy the constraints, respectively.

From the boundary conditions, we have
```math
\begin{bmatrix} 
x_0  & 1 \\ 
x_N  & 1
\end{bmatrix}

\begin{bmatrix}
a_i \\
e_i
\end{bmatrix} = 

\begin{bmatrix}
x_i\\ 
x_{i + 1}
\end{bmatrix}

\quad i = 0, \ldots, N - 1
```
which implies 
```math
\begin{bmatrix}
a_i \\
e_i
\end{bmatrix} = 

\dfrac{1}{x_0 -  x_N}

\begin{bmatrix} 
1  & -1 \\ 
-x_N  & x_0 \\ 
\end{bmatrix}
\begin{bmatrix}
x_i \\ 
x_{i + 1}
\end{bmatrix}
```

Similarly, we have, 
```math 
\begin{bmatrix} 
x_0  & 1 \\ 
x_N  & 1 \\ 
\end{bmatrix}

\begin{bmatrix}
c_i \\
f_i
\end{bmatrix} = 

\begin{bmatrix}
y_i - d_i y_0 \\ 
y_{i + 1} - d_i y_N
\end{bmatrix}
\quad i = 0, \ldots, N - 1
```
which implies 
```math
\begin{bmatrix}
c_i \\
f_i
\end{bmatrix} = 

\dfrac{1}{x_0 -  x_N}

\begin{bmatrix} 
1  & -1 \\ 
-x_N  & x_0 \\ 
\end{bmatrix}
\begin{bmatrix}
y_i - d_i y_0 \\ 
y_{i + 1} - d_i y_N
\end{bmatrix}
```

[^1]: $\mathcal{C}(\bm{X})$ denotes the set of continuous functions defined over the set $\bm{X}$.

## 1D Hidden Fractal Interpolation 

Consider the following IFS, 
```math
w_i = \left(
\begin{bmatrix}
x \\ 
y \\
z
\end{bmatrix}
\right) = 

\begin{bmatrix}
a_i & 0 & 0 \\
c_i & d_i & h_i \\
k_i & l_i & m_i 
\end{bmatrix}

\begin{bmatrix}
x \\ 
y \\
z
\end{bmatrix}

+ 

\begin{bmatrix}
e_i \\ 
f_i \\
g_i
\end{bmatrix}
```
We assume that the following boundary conditions are satisfied. 
```math 
w_i = \left(
\begin{bmatrix}
x_0 \\ 
y_0 \\
z_0
\end{bmatrix}
\right) = 
\begin{bmatrix}
x_{i} \\ 
y_i \\
z_i
\end{bmatrix}
```
and 
```math 
w_i = \left(
\begin{bmatrix}
x_N \\ 
y_N \\
z_N
\end{bmatrix}
\right) = 

\begin{bmatrix}
x_{i + 1} \\ 
y_{i + 1} \\
z_{i + 1}
\end{bmatrix}
```
Given the free parameters  ``z_i, d_i, h_i, l_i, m_i``, which are called as hidden variables, the remaining parameter ``a_i, c_i, k_i, e_i, f_i, g_i`` can be calculated. 

Using the boundary conditions, we can write 
```math
\begin{bmatrix} 
x_0  & 1 \\ 
x_N  & 1
\end{bmatrix}

\begin{bmatrix}
a_i \\
e_i
\end{bmatrix} = 

\begin{bmatrix}
x_i\\ 
x_{i + 1}
\end{bmatrix}

\quad i = 0, \ldots, N - 1
```
which implies 
```math
\begin{bmatrix}
a_i \\
e_i
\end{bmatrix} = 

\dfrac{1}{x_0 -  x_N}

\begin{bmatrix} 
1  & -1 \\ 
-x_N  & x_0 \\ 
\end{bmatrix}
\begin{bmatrix}
x_i \\ 
x_{i + 1}
\end{bmatrix}
```

Similarly, we have 
```math
\begin{bmatrix} 
x_0  & 1 \\ 
x_N  & 1
\end{bmatrix}

\begin{bmatrix}
c_i \\
f_i
\end{bmatrix} = 

\begin{bmatrix}
y_i - d_i y_0 - h_i z_0 \\ 
y_{i + 1} - d_i y_N - h_i z_N
\end{bmatrix}

\quad i = 0, \ldots, N - 1
```
which implies 
```math
\begin{bmatrix}
c_i \\
f_i
\end{bmatrix} = 

\dfrac{1}{x_0 -  x_N}

\begin{bmatrix} 
1  & -1 \\ 
-x_N  & x_0 \\ 
\end{bmatrix}
\begin{bmatrix}
y_i - d_i y_0 - h_i z_0\\ 
y_{i + 1} - d_i y_N - h_i z_N
\end{bmatrix}
```

Again, 
```math
\begin{bmatrix} 
x_0  & 1 \\ 
x_N  & 1
\end{bmatrix}

\begin{bmatrix}
k_i \\
g_i
\end{bmatrix} = 

\begin{bmatrix}
z_i - d_i y_0 - h_i z_0 \\ 
z_{i + 1} - d_i y_N - h_i z_N
\end{bmatrix}

\quad i = 0, \ldots, N - 1
```
which implies 
```math
\begin{bmatrix}
k_i \\
g_i
\end{bmatrix} = 

\dfrac{1}{x_0 -  x_N}

\begin{bmatrix} 
1  & -1 \\ 
-x_N  & x_0 \\ 
\end{bmatrix}
\begin{bmatrix}
z_i - d_i y_0 - h_i z_0\\ 
z_{i + 1} - d_i y_N - h_i z_N
\end{bmatrix}
```