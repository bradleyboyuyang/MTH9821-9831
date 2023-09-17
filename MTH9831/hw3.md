### q2
the question suffice to show that
$$
\mathbb{E}(I(t) | \mathcal{F}(s)) = \mathbb{E}(I(t) - I(s) + I(s) | \mathcal{F}(s)) = \mathbb{E}(I(t) - I(s) | \mathcal{F}(s)) + I(s)= I(s),\quad \forall s < t.
$$
suppose $s \in [t_m, t_{m+1})$ and $t\in [t_n, t_{n+1})$, then
$$
I(t) - I(s) = \Delta_m (B(t_{m+1}) - B(s)) + \sum_{i=m+1}^{n}\Delta_{i}(B(t_{i+1}) - B(t_{i})) + \Delta_k(B(t)-B(t_{n})).
$$
for $i\ge m$, we have
$$
\mathbb{E}(\Delta_{i}(B(u) - B(v))|\mathcal{F}(s))=0, \text{ for } \forall \max(t_i,s)\le v\le u
$$
this is because $B(u) - B(v)$ is independent with $\Delta_i$ and $\mathcal{F}(s)$, so we have
$$
\begin{align*}
\mathbb{E}(\Delta_{i}(B(u)-B(v))|\mathcal{F}(s))&=\mathbb{E}(\mathbb{E}[\Delta_{i}(B(u)-B(v))|\mathcal{F}(s),\Delta_{i}]|\mathcal{F}(s))\\
&=\mathbb{E}(\Delta_{i}\mathbb{E}[B(u)-B(v)|\mathcal{F}(s),\Delta_{i}]|\mathcal{F}(s))=0
\end{align*}
$$
Therefore we can conclude that $I(t)$ is an $\mathcal{F}(t)$-martingale


### q6
#### (i)
Expectation:
$$
\begin{align*}
\mathbb E\int_0^tB^2(s)ds&=\int_0^t\mathbb E[B^2(s)]ds=\int_0^tsds=\frac{1}{2}t^2\\
\end{align*}
$$
Variance:
$$
\begin{align*}
\mathbb E(\int_0^tB^2(s)ds)^2&=\mathbb E(\int_0^t\int_0^tB^2(s)B^2(v)dsdv)\\
&=\int_0^t\int_0^t\mathbb E[B^2(s)B^2(v)]dsdv\\
&=\iint_{0<s<v<t}\mathbb E[B^2(s)B^2(v)]dsdv+\iint_{0<v<s<t}\mathbb E[B^2(s)B^2(v)]dsdv\\
&=2\int_0^t\int_0^v\mathbb E[B^2(s)B^2(v)]dsdv\\
&=2\int_0^t\int_0^v\mathbb E[B^2(s)(B(s)+(B(v)-B(s))^2)]dsdv\\
&=2\int_0^t\int_0^v(3s^2+s(v-s))dsdv\\
&=2\int_0^t\frac{7}{6}v^3dv=\frac{7}{12}t^4\\
Var(\int_0^tB^2(s)ds)&=\mathbb E(\int_0^tB^2(s)ds)^2-(\mathbb E\int_0^tB^2(s)ds)^2=\frac{1}{3}t^4
\end{align*}
$$
#### (ii)
$$
\mathbb E(X^2(t))=\mathbb E(e^{2B^2(t)})=\mathbb E(e^{2tY})
$$
where $Y\sim\chi^2_1$. so we obtain
$$
\mathbb E(X^2(t))=\mathbb E(e^{2tY})=(1-4t)^{-1/2}
$$
for the expectation to be finite, we require $t\in[0,\frac{1}{4})$