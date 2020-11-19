## Two Dimensional Fractal Interpolation

To generalize the concepts given in the previous section for FIFs to the FISs, consider that $\Omega = \{\tilde{P}_j = (\tilde{x}_j, \tilde{y}_j), \; j = 1, 2, 3\}$ is a triangular domain in the plane as shown in Fig. \ref{fig: triangulation}. Let $P = \{P_i = (x_i, y_i), \; i = 1, \ldots, N\}$ be given points in the plane scattered over $\Omega$ containing the points $\tilde{P}_j, \; j = 1, 2, 3$. A triangulation $T(\Omega)$ of $\Omega$ over the points $P$ is given with the set, \begin{equation} T(\Omega) = \{ \Omega^i = \{P_j^i\} : P_j^i \in P, \; j = 1, 2, 3, \;  i  = 1, \ldots, K \} \label{eq: triangulation} \end{equation} Note that $T(\Omega)$ consists of non-degenerate triangles such that $\Omega^i \cap \Omega^j = \emptyset, \; i \neq j$ and $\Omega = \cup_{i  = 1}^K \Omega^i$.

```@raw html
<center>
    <img src="../../assests/triangulation.svg" alt="model" width="30%"/>
</center>
```

Consider that each $P_i = (x_i, y_i) \in P$ is associated with $z_i \in \mathbb{R}$. An interpolation function--also called interpolant, corresponding to the data set $\{(x_i, y_i, z_i), \; i = 1, \ldots, N\}$  is a continuous function $f: \Omega \mapsto \mathbb{R}$ such that $f(x_i, y_i) = z_i$. The answer to the problem of constructing the interpolation function $f$ is in two folds:

1. Construction of an IFS $\mathcal{I} = \{\Omega \times \mathbb{R}; w_i; \; i = 1, \ldots, K \}$ whose attractor is the graph$G$ of the function $f$ satisfying the interpolation points, i.e.$f(x_i, y_i) = z_i, \; i = 1, \ldots, N$. 

2. Construction of a contraction mapping $M: \mathcal{F} \mapsto\mathcal{F}$ where $\mathcal{F} =\{ \tilde{f} \in \mathcal{C}(\Omega) : \tilde{f}(\tilde{x}_j, \tilde{y}_j) = \tilde{z}_j \}$ such that the fixed point $f$ of the mapping $M$ satisfies the interpolation points, i.e. $f(x_i, y_i) = z_i, \; i = 1, \ldots, N$

### Construction of the IFS $\mathcal{I}$
Addressing to the first step, we put forward the following IFS $I = \{ \Omega \times \mathbb{R}; w_i, \; i = 1, \ldots, K \}$ with affine transformations $w_i$,

```math
    \begin{aligned}
     w_i(x, y, z) &= (L_i(x, y), F_i(x, y, z))\\ 
    &=\begin{bmatrix}
    \alpha_1^i & \alpha_2^i & 0 \\
    \alpha_3^i & \alpha_4^i & 0 \\
    \alpha_5^i & \alpha_6^i & \alpha_7^i
    \end{bmatrix}
    \begin{bmatrix}
    x \\ y \\ z
    \end{bmatrix} + \begin{bmatrix}
    \beta_1^i \\ \beta_2^i \\ \beta_3^i
    \end{bmatrix} 
     \; i = 1, \ldots, K,

    \end{aligned}

```
where $L_i : \Omega \mapsto \Omega^i$

```math
    L_i(x, y) = \begin{bmatrix}
    \alpha_1^i & \alpha_2^i \\ \alpha_3^i & \alpha_4^i
    \end{bmatrix} \begin{bmatrix}
    x \\ y
    \end{bmatrix} + \begin{bmatrix}
    \beta_1^i \\ \beta_2^i
    \end{bmatrix}
```
are contraction mappings for the $z$ axis satisfying the boundary conditions,
```math
F_i(\tilde{x}_j, \tilde{y}_j, \tilde{z}_j) = z_j^i, \; j = 1, 2, 3.
```
where $\alpha_7^i$ are arbitrary contractivity factors satisfying $|\alpha_7^i| < 1, \; i = 1, \ldots, K$, also called as vertical scaling factors.

Given the vertical scaling factors $\alpha_7^i, \; i = 1, \ldots, K$, the coefficients $\alpha_k^i, \; k = 1 , \ldots, 6, \; i = 1, \ldots, K$ can be found using the boundary conditions in (\ref{eq: plane boundary conditions}) and (\ref{eq: z values boundary conditons}) which results in following system of equations,

```math
\begin{aligned}
\alpha_1^i \tilde{x}_j + \alpha_2^i \tilde{y}_j + \beta_1^i &= x_j^i \\
\alpha_3^i \tilde{x}_j + \alpha_4^i \tilde{y}_j + \beta_2^i &= y_j^i \\
\alpha_5^i \tilde{x}_j + \alpha_6^i \tilde{y}_j + \beta_3^i &= z_j^i - \alpha_7^i \tilde{z}_j 
\end{aligned}
```
for $i = 1, \ldots, K, \; j = 1, 2, 3$. This system can be rewritten in block diagonal matrix equation system,
```math
\begin{bmatrix}
\bm{P} & \bm{0} & \ldots & \bm{0} \\
\bm{0} & \bm{P} & \ldots & \bm{0} \\
\vdots & \vdots & \ddots & \vdots \\
\bm{0} & \bm{0} & \ldots & \bm{P} \\
\end{bmatrix}
\begin{bmatrix}
\bm{r}_1 \\ \bm{r}_2 \\ \bm{r}_3
\end{bmatrix} = \begin{bmatrix}
\bm{c}_1 \\ \bm{c}_2 \\ \bm{c}_3
\end{bmatrix}
```
where
```math
\begin{aligned}
\bm{P} &= 
\begin{bmatrix}
\tilde{x}_1 & \tilde{y}_1 & 1 \\ 
\tilde{x}_2 & \tilde{y}_2 & 1 \\ 
\tilde{x}_3 & \tilde{y}_3 & 1 \\ 
\end{bmatrix}\\
\bm{r}_1 &= [\bm{r}_1^1, \ldots, \bm{r}_1^i, \ldots, \bm{r}_1^K], \; \bm{r}_1^i = [\alpha_1^i, \alpha_2^i, \beta_1^i] \\
\bm{r}_2 &= [\bm{r}_2^1, \ldots, \bm{r}_2^i, \ldots, \bm{r}_2^K], \; \bm{r}_2^i = [\alpha_3^i, \alpha_4^i, \beta_2^i]  \\
\bm{r}_3 &= [\bm{r}_3^1, \ldots, \bm{r}_3^i, \ldots, \bm{r}_3^K], \; \bm{r}_3^i = [\alpha_5^i, \alpha_6^i, \beta_3^i] \\
\bm{c}_1 &= [\bm{c}_1^1, \ldots, \bm{c}_1^i, \ldots, \bm{c}_1^K], \; \bm{c}_1^i = [x_1^i, x_2^i, x_3^i] \\
\bm{c}_2 &= [\bm{c}_2^1, \ldots, \bm{c}_2^i, \ldots, \bm{c}_2^K], \; \bm{c}_2^i = [y_1^i, y_2^i, y_3^i] \\
\bm{c}_3 &= [\bm{c}_3^1, \ldots, \bm{c}_3^i, \ldots, \bm{c}_3^K], \; \bm{c}_3^i = [z_1^i, z_2^i, z_3^i] - \alpha_7^i [\tilde{z}_1^i, \tilde{z}_2^i, \tilde{z}_3^i]
\end{aligned}
```
We  have uncoupled system of equations,
```math
\bm{P} \bm{r}_j^i = \bm{c}_j^i, \; j = 1, 2, 3, \; i = 1, \ldots, N
```
Since the points $\{\tilde{P}_j, \; j = 1, 2, 3\}$ forms a non-degenerate triangular region $\Omega$, $\bm{P}^{-1}$ exists and can be solved for the coefficients of the IFS $\mathcal{I}$ which gives, 
```math
\bm{r}_j^i = \bm{P}^{-1} \bm{c}_j^i, \; j = 1, 2, 3, \; i = 1, \ldots, K
```

### Construction of the mapping $M$

Inspired by the reasoning given in Theorem 1, we propose the following conjecture for the second part of the problem.

!!! note 
    
    Consider the transformation ``L`` shown in the figure below where $\Omega_n$ and $\Omega \subset \mathbb{R}^{2}$ and the IFS $\mathcal{I}$ in (\ref{eq: IFS}), let $(\mathcal{F}, d)$ be a complete metric space such that

    
    
    

	
