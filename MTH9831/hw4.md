### q6

First, we verify that $S(t)=Ae^{\sigma B(t)-\sigma^2t/2}$. 
</br>According to Ito's Formula, we have
$$
\begin{align*}
dS(t)=S_t(t)dt+S_x(t)dB(t)+\frac{1}{2}S_{xx}(t)dt=-\frac{\sigma^2}{2}S(t)dt+\sigma S(t)dB(t)+\frac{1}{2}\sigma^2S(t)dt=\sigma S(t)dB(t)\\
\end{align*}
$$
and we have $S(0)=A$. so $S(t)=Ae^{\sigma B(t)-\sigma^2t/2}$ solves the quation.

From $dS(t)=\sigma S(t)dB(t)$, we conclude $S(t)$ is a martingale. so we have
$$
\begin{align*}
\mathbb E[S(t)] = \mathbb E[S(0)] = A
\end{align*}
$$
Applying Ito's formula to $S^2(t)$ and take expectation, we obtain
$$
\begin{align*}
dS^2(t)&=2S(t)dS(t)+d[S,S](t)=2\sigma S^2(t)dB(t)+\sigma^2B^2(t)dt\\
S^2(t)&=S^2(0)+\int_0^t2\sigma S^2(r)dB(r)+\int_0^t\sigma^2B^2(r)dr\\
\mathbb E(S^2(t))&=\mathbb E(S^2(0))+\mathbb E(\int_0^t2\sigma S^2(r)dB(r)) +\mathbb E(\int_0^t\sigma^2B^2(r)dr)\\
&=A^2+\sigma^2\int_0^t\mathbb E(B^2(r))dr\\
&=A^2+\sigma^2\int_0^t\mathbb rdr=A^2+\frac{1}{2}\sigma^2t^2
\end{align*}
$$
therefore we obtain
$$
\begin{align*}
Var(S(t))=\mathbb E(S^2(t))-(\mathbb ES(t))^2=\frac{1}{2}\sigma^2t^2
\end{align*}
$$

### q7
based on using integration by parts, we have
$$
\begin{align*}
\int_0^tS(u)du&=tS(t)-\int_0^tudS(u)=tS(t)-\int_0^tu\sigma S(u)dB(u)=\int_0^t(t-u\sigma S(u))dB(u)
\end{align*}
$$
therefore we obtain
$$
\begin{align*}
\mathbb E\int_0^tS(u)du&=\mathbb E\int_0^t(t-u\sigma S(u))dB(u)=0\\
Var[\int_0^tS(u)du]&=\mathbb E\int_0^t(t-u\sigma S(u))^2du\\
&=\int_0^t(t^2-2u\sigma t\mathbb E(S(u))+u^2\sigma^2\mathbb E(S^2(u)))du\\
&=\int_0^tt^2du-2\sigma t\int_0^tAudu+\sigma^2\int_0^tu^3du\\
&=(1-A\sigma) t^3+\frac{1}{4}\sigma^2t^4

\end{align*}
$$