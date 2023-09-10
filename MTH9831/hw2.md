## q2

To start with
$$
\limsup_{t\to0^+} \frac{|B(t)|}{\sqrt{2t\ln\ln(1/t)}} = \limsup_{s\to+\infty} \frac{|B(s^{-1})|}{\sqrt{(2/s)\ln\ln(s)}} = \limsup_{s\to+\infty} \frac{|sB(s^{-1})|}{\sqrt{2s\ln\ln(s)}}
$$

According to part (iii) of Theorem 2.1, $B^*(t) = tB(t^{-1})$ is a Brownian motion. By using law of iterated logatithm
$$
\limsup_{s\to\infty}\frac{B(s)}{\sqrt{2s\ln\ln(s)}} = 1 \quad a.s.
$$

since we have
$$
|sB(s^{-1})| = max(B^*(s), -B^*(s))
$$

both of which are Brownian motion, so we can conclude that
$$
\begin{align*}
\limsup_{s\to+\infty} \frac{|sB(s^{-1})|}{\sqrt{2s\ln\ln(s)}} &= \limsup_{s\to+\infty} \frac{max(B^*(s), -B^*(s))}{\sqrt{2s\ln\ln(s)}} \\
&= max(\limsup_{s\to\infty}\frac{B^*(s)}{\sqrt{2s\ln\ln(s)}}, \limsup_{s\to\infty}\frac{-B^*(s)}{\sqrt{2s\ln\ln(s)}}) \\
&= 1 \quad a.s
\end{align*}
$$


## q4
### (i)
1) $\mathbb{E}[|Z(t)|] = mathbb{E}[Z(t)] = 1 < \infty$ 
2) $Z(t)$ is $\mathcal F$-measurable
3) for $s\le t$, we have
$$
\begin{align*}
\mathbb{E}[Z(t)|\mathcal F(s)] &= \mathbb{E}[exp\{\sigma X(t) - (\sigma\mu + \frac{1}{2}\sigma^2)t\}|\mathcal F(s)] \\
&= Z(s)\mathbb{E}[exp\{\sigma(X(t)-X(s)) - (\sigma\mu + \frac{1}{2}\sigma^2)(t-s)\}|\mathcal F(s)] \\
&= Z(s)\mathbb{E}[exp\{\sigma(X(t)-X(s)) - (\sigma\mu + \frac{1}{2}\sigma^2)(t-s)\}] \\
&= Z(s)exp\{-(\sigma\mu + \frac{1}{2}\sigma^2)(t-s)\}\mathbb{E}[exp\{\sigma(X(t)-X(s))\}] \\
&= Z(s)
\end{align*}
$$

Here we use the property: 
</br> (i) $X(t)-X(s)$ is independent with $\mathcal F(s)$
</br> (ii) $X(t)-X(s)$ ～ $N(\mu t, \sqrt{t-s})$

so we can conclude $Z(t)$ is a martingale

### (ii)
$Z(t\wedge\tau_{m})$ is a stopped martingale. so we have
$$
1 = Z(0) = Z(t\wedge\tau_{m}) = \mathbb{E}[exp\{\sigma X(t\wedge\tau_{m}) - (\sigma\mu + \frac{1}{2}\sigma^2)(t\wedge\tau_{m})\}], \quad t\ge0
$$

### (iii)
since we have 
$$
0\le exp\{\sigma X(t\wedge\tau_{m})\}\le e^{\sigma m}
$$
which is bounded. take limit on both sides of the equation in (ii), we get
$$
\lim_{t\to\infty}exp\{\sigma X(t\wedge\tau_{m}) - (\sigma\mu + \frac{1}{2}\sigma^2)(t\wedge\tau_{m})\} = \mathbb{I}_{(\tau < \infty)}exp\{\sigma m-(\sigma\mu + \frac{1}{2}\sigma^2)\tau_{m}\}
$$
so we have
$$
\mathbb{E}[exp\{\sigma m-(\sigma\mu + \frac{1}{2}\sigma^2)\tau_{m}\}\mathbb{I}_{(\tau < \infty)}] = 1 \\
\mathbb{E}[exp\{-(\sigma\mu + \frac{1}{2}\sigma^2)\tau_{m}\}\mathbb{I}_{(\tau < \infty)}] = e^{-\sigma m}
$$
take limit on both sides as $\sigma\to 0$, we get $\mathbb{E}[\tau_m < \infty] = 1$, equivalently
$$
\mathbb{P}[\tau_m < \infty] = 1
$$
so we have
$$
\mathbb{E}[exp\{-(\sigma\mu + \frac{1}{2}\sigma^2)\tau_{m}\}] = e^{-\sigma m}
$$
set $\alpha = (\sigma\mu + \frac{1}{2}\sigma^2)$, since $\sigma$ is positive, then we get $\sigma = -\mu + \sqrt{2\alpha + \mu^2}$. 
</br> so for all $\alpha > 0$, we have
$$
\mathbb{E}[e^{-\alpha\tau_{m}}] = e^{m\mu-m\sqrt{2\alpha + \mu^2}}
$$

### (iv)
take derivative with respect to $\alpha$, then we have
$$
\mathbb{E}[\tau_m e^{-\alpha\tau_{m}}] = \frac{m}{\sqrt{2\alpha+\mu^2}}e^{m\mu-m\sqrt{2\alpha + \mu^2}}
$$
let $\alpha\to0$, we get 
$$
\mathbb{E}[\tau_m] = \frac{m}{\mu}
$$

### (v)
for $\sigma > -2\mu$, we have $\sigma\mu + \frac{1}{2}\sigma^2 > 0$. so similar as before, we have
$$
\mathbb{E}[exp\{\sigma m-(\sigma\mu + \frac{1}{2}\sigma^2)\tau_{m}\}\mathbb{I}_{(\tau < \infty)}] = 1 \\
\mathbb{E}[exp\{-(\sigma\mu + \frac{1}{2}\sigma^2)\tau_{m}\}\mathbb{I}_{(\tau < \infty)}] = e^{-\sigma m}
$$
and then get


