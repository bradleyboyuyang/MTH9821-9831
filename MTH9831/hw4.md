### q11

Let $f(t,r(t)) = e^{\beta t} r(t) $, since $f_t=\beta e^{\beta t} $, $f_r = e^{\beta t}$, and $f_{rr}=0$ applying Ito's lemma, we have

$$
df(t,B(t)) = \frac{\partial f}{\partial t} dt + \frac{\partial f}{\partial r} dr + \frac{1}{2} \frac{\partial^2 f}{\partial r^2} drdr \\
d(e^{\beta t} r(t)) = \beta e^{\beta t}  r(t) dt + e^{\beta t}  dr + 0 \\
d(e^{\beta t} dr(t))  = \beta e^{\beta t}  r(t) dt + e^{\beta t}  (\alpha dt - \beta r(t)dt + \sigma d B(t)) \\
d(e^{\beta t} dr(t))  = e^{\beta t}  (\alpha dt + \sigma d B(t)) \\
e^{\beta t}r(t) = r(0) + \int_0^t e^{\beta s} (\alpha ds + \sigma d B(s)) \\
r(t) = e^{-\beta t} r(0) + \alpha \int_0^t e^{-\beta (t-s)} ds + \sigma \int_0^t e^{-\beta (t-s)} d B(s) \\
= e^{-\beta t} r(0) + \frac{\alpha}{\beta} (1-e^{-\beta t}) + \sigma \int_0^t e^{-\beta (t-s)} d B(s) \\
$$
Hence, $r(t)$ is a Gaussian process,  with mean $e^{-\beta t} r(0)  + \cfrac{\alpha}{\beta} (1-e^{-\beta t})$ and variance $\int_{0}^te^{2\beta(s-t)ds} = \frac{\sigma^2}{2\beta} (1-e^{-2\beta t})$, namely, $r(t) \sim \mathcal{N}(e^{-\beta t} r(0)  + \cfrac{\alpha}{\beta} (1-e^{-\beta t}), \frac{\sigma^2}{2\beta} (1-e^{-2\beta t}))$

### q12
(a) From question 11, we adopt exactly the same approach and can get:

$$
dr(t) = (\alpha-\beta r(t))dt + \sigma \sqrt{r(t)} d B(t) \\
r(t) = e^{-\beta t} r(0) + \frac{\alpha}{\beta} (1-e^{-\beta t}) + \sigma \int_0^t e^{-\beta (t-s)} \sqrt{r(s)} d B(s) \\
$$
Therefore the mean is given by
$$
E[r_t] = e^{-\beta t} r_0 + \frac{\alpha}{\beta} (1-e^{-\beta t})
$$
Then by Ito isometry, we have the variance
$$
Var(r_t) = \sigma^2 E[\int_0^t e^{-2\beta (t-s)} r(s) ds] = \sigma^2 \int_0^t e^{-2\beta (t-s)} E[r(s)] ds \\
$$

(b) Assume that $4\alpha=\sigma^2$, let $X(t)=\sqrt{r(t)}$, then $X(t)$ satisfies the following SDE:
$$
dX(t) = \frac{1}{2\sqrt{r(t)}}dr(t) - \frac{1}{4} \frac{1}{r(t)^{3/2}} dr(t) dr(t) \\
= \frac{1}{2\sqrt{r(t)}} (\alpha-\beta r(t))dt + \frac{1}{2} \sigma d B(t) - \frac{1}{4 r(t)^{3/2}} (\sigma^2 r(t)) dt \\ 
= \frac{1}{2\sqrt{r(t)}} (\alpha-\beta r(t) - \frac{1}{2} \sigma^2) dt + \frac{1}{2} \sigma d B(t) \\
= -\frac{1}{2X(t)}(\alpha+\beta X^2(t))dt + \frac{1}{2} \sigma d B(t) \\
$$

(c) Using part (b) to determine the distribution of r(t):