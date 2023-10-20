#### Problem 1

$$
\begin{align*}
dX(t) &= \Delta_{1}(t)dS_{1}(t) + \Delta_{2}(t)dS_{2}(t) + \Gamma(t)dM(t) \\
&+ (S_{1}(t)d\Delta_{1}(t) + d\Delta_{1}(t)dS_{1}(t) + S_{2}(t)d\Delta_{2}(t) + d\Delta_{2}(t)dS_{2}(t) + M(t)d\Gamma(t) + d\Gamma(t)dM(t))\\
\end{align*}
$$

The self-financing condition
$$
\begin{align*}
S_{1}(t)d\Delta_{1}(t) + d\Delta_{1}(t)dS_{1}(t) + S_{2}(t)d\Delta_{2}(t) + d\Delta_{2}(t)dS_{2}(t) + M(t)d\Gamma(t) + d\Gamma(t)dM(t) = 0
\end{align*}
$$
Under this condition, we have
$$
\begin{align*}
dX(t) &= \Delta_{1}(t)dS_{1}(t) + \Delta_{2}(t)dS_{2}(t) + \Gamma(t)dM(t)
\end{align*}
$$
To derive (1.5), from
$$
\begin{align*}
dD(t)dX(t) &= - R(t)D(t)(\Delta_{1}(t)dS_{1}(t) + \Delta_{2}(t)dS_{2}(t) + \Gamma(t)dM(t))dt = 0\\
dS_{i}(t) &= \alpha_{i}S_{i}(t)dt + \sigma_{i}S_{i}(t)dB(t)\\
dM(t) &= R(t)M(t)dt
\end{align*}
$$
we have
$$
\begin{align*}
d(D(t)X(t)) &= D(t)dX(t) + X(t)dD(t) + dD(t)dX(t) \\
&= D(t)(\Delta_{1}(t)dS_{1}(t) + \Delta_{2}(t)dS_{2}(t) + \Gamma(t)dM(t))
- (\Delta_{1}(t)S_{1}(t) + \Delta_{2}(t)S_{2}(t) + \Gamma(t)M(t))R(t)D(t)dt\\
&= D(t)(\Delta_{1}(t)dS_{1}(t) + \Delta_{2}(t)dS_{2}(t)) - (\Delta_{1}(t)S_{1}(t) + \Delta_{2}(t)S_{2}(t))R(t)D(t)dt\\
&= D(t)(\Delta_{1}(t)(\alpha_{1}-R(t))S_{1}(t) + \Delta_{2}(t)(\alpha_{2}-R(t))S_{2}(t))dt + D(t)(\Delta_{1}(t)\sigma_{1}S_{1}(t) + \Delta_{2}(t)\sigma_{2}S_{2}(t))dB(t)
\end{align*}
$$
when $R(t) = r$ is a constant, we have
$$
\begin{align*}
d(D(t)X(t)) &= D(t)(\sigma_{1}\Delta_{1}(t)S_{1}(t) + \sigma_{2}\Delta_{2}(t)S_{2}(t))dB(t) + D(t)((\alpha_{1}-r)\Delta_{1}(t)S_{1}(t) + (\alpha_{2}-r)\Delta_{2}(t)S_{2}(t))dt
\end{align*}
$$


#### Problem 2

Given that the derivative security can be hedged, there exist a portfolio value process $X(t)$ with initial condition $X(0)$ satisfies
$$
\begin{align*}
X(T) = V(T)
\end{align*}
$$
Assume there are 2 risk-neutral probability measures $\mathbb {\tilde P}_{1}$ and $\mathbb {\tilde P}_{2}$, then the discounted portfolio value probability $D(t)X(t)$ is a martingale under both  $\mathbb {\tilde P}_{1}$ and $\mathbb {\tilde P}_{2}$. so we have
$$
\begin{align*}
V_{1}(t) &= \frac{1}{D(t)}\mathbb{\tilde E}_{1}[D(T)V(T)|\mathcal{F}(t)] = \frac{1}{D(t)}\mathbb{\tilde E}_{1}[D(T)X(T)|\mathcal{F}(t)] = X(t)\\
&=\frac{1}{D(t)}\mathbb{\tilde E}_{2}[D(T)X(T)|\mathcal{F}(t)] = \frac{1}{D(t)}\mathbb{\tilde E}_{2}[D(T)V(T)|\mathcal{F}(t)] = V_{2}(t)
\end{align*}
$$
So we can conclude that if $V(t)$ can be hedged, then its price under any risk-neutral probability measure is the same.