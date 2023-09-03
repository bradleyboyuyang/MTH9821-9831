What is the distribution of

Find the distribution of 
$$
(X_{2} | X_{1} = X)
$$

Solution:
We know that we can decompose:

$$
X_{2} - \mu_{2} = \Lambda (X_{1} - \mu_{1}) + W
$$
We now calculate W: 
First we multiply by 
$$
X_{1} - \mu_{1}
$$
from the right to get.

$$
(X_{2} - \mu_{2} ) (X_{1} - \mu_{1})^{T} = (\Lambda (X_{1} - \mu_{1}) + W)(X_{1} - \mu_{1})^{T}
$$
Now we apply the expected value and get:
$$
C_{21} = \Lambda C_{11}
$$
as the expected Value of W is 0 and it is independent of X1 by construction.

Hence we derive:
$$
\Lambda = C_{21} C_{11}^{-1}
$$

We now want to find W.
We multiply by
$$
X_{2} - \mu_{2}
$$
on the right and 
$$
(\Lambda (X_{1} - \mu_{1}) + W)
$$
on the left to get 
$$
(X_{2} - \mu_{2} ) (X_{2} - \mu_{2})^{T} = (\Lambda (X_{1} - \mu_{1}) + W) (\Lambda (X_{1} - \mu_{1}) + W)^{T}
$$
$$
C_{22} = C_{21} C_{11}^{-1} (X_{1}-\mu_{1})(X_{1}-\mu_{1})^{T} (C_{21}C_{11}^{-1})^{T} + WW^{T}
$$
$$
C_{22} = C_{21} C_{11}^{-1}C_{12} + Var(W)
$$
Hence the distribution of 
$$ 
(X_{2}|X_{1} = X) = 
$$
$$
(\Lambda(X_{1} - \mu_{1}) + \mu_{2} + W | X_{1} =X)
$$
$$
\Lambda (X - \mu_{1}) + \mu_{2} + W
$$
$$ 
 C_{21} C_{11}^{-1} (X - \mu_{1}) + \mu_{2} + W
$$
where W is distributed like:
$$
\mathcal{N}( 0, C_{22} - C_{21}C_{11}^{-1}C_{12})
$$


