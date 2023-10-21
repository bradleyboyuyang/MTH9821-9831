#### Problem 3

##### (a) 

Log-normal SDE for dividend paying stock $S$
$$
\begin{align*}
dS_t = (\mu - q)S_tdt + \sigma S_tdW_t
\end{align*}
$$
Use Ito's Lemma to derive SDE for $F_{t}^{T}$
$$
\begin{align*}
dF_{t}^{T} &= d(S_{t}e^{-q(T-t)}) = e^{-q(T-t)}dS_{t} + S_{t}de^{-q(T-t)} + dS_{t}de^{-q(T-t)}\\
&= e^{-q(T-t)}((\mu - q)S_tdt + \sigma S_tdW_t) + qe^{-q(T-t)}S_{t}dt\\
&= e^{-q(T-t)}S_{t}(\mu dt + \sigma dW_t)\\
&= \mu F_{t}^{T}dt + \sigma F_{t}^{T}dW_{t}
\end{align*}
$$

##### (b)

From
$$
\begin{align*}
F_{t}^{T} &= S_{t}e^{-q(T-t)}\\
C(S_{t}, t) &= V(F_{t}^{T}, t)
\end{align*}
$$
We can obtain
$$
\begin{align*}
\frac{\partial V}{\partial t} &= \frac{\partial C}{\partial S}\frac{\partial S}{\partial t} + \frac{\partial C}{\partial t} = -qSC_S + C_t\\
\frac{\partial V}{\partial F} &= \frac{\partial C}{\partial S}\frac{\partial S}{\partial F} + \frac{\partial C}{\partial t}\frac{\partial t}{\partial F} = e^{q(T-t)}C_S\\
\frac{\partial^2 V}{\partial F^2} &= \frac{\partial V_F}{\partial S}\frac{\partial S}{\partial F} + \frac{\partial V_F}{\partial t}\frac{\partial t}{\partial F} = e^{2q(T-t)}C_{SS}
\end{align*}
$$
Replace $V$ and $F$ with $C$ and $S$, we obtain the BS PDE for a dividend paying asset
$$
\begin{align*}
-qSC_S + C_t +\frac{1}{2}\sigma^2S^2C_{SS} + rSC_{S} - rC &= 0\\
C_t + (r-q)SC_{S} + \frac{1}{2}\sigma^2S^2C_{SS} - rC &= 0
\end{align*}
$$
