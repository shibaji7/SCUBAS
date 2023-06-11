<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->
#### Electrically Long Cable
---
Consider an Ocean-Earth section with physical length, $L$, and propagation constant, $\gamma$. The section has an adjustment distance $\frac{1}{\gamma}$. For an electrically-long transmission line, where $L>\frac{4}{\gamma}$ (i.e. length is greater than four times of adjustment distance) we have the following scenario,

$$
e^{\gamma L}>>1>>e^{-\gamma L}
$$

When the transmission line length is considerably shorter than the adjustment distance it is referred to as 'electrically-long' and the equivalent-$\pi$ components reduce to (Boteler et a, 2013):

$$
Y_E=\frac{2}{Z_0e^{\gamma L}}
$$

$$
\frac{Y'}{2}=\frac{1}{Z_0}
$$

$$
I_E=\frac{E}{Z}
$$

$$
V_i=-U_k=-\frac{E}{\gamma}
$$

$$
V(x)=Ve^{-\gamma (L-x)}-Ve^{-\gamma x}
$$

where $V_i=-V_k=V=\frac{E}{\gamma}$.

##### Example Codes
``` py

```

#### Electrically Short Cable
For an electrically-short section, where the physical length is less than the adjustment distance, i.e., $L<\frac{1}{\gamma}$. When the transmission line length is considerably shorter than the adjustment distance it is referred to as 'electrically-short' and the equivalent-$\pi$ components reduce to (Boteler et a, 2013):

$$
e^{\pm\gamma L}\approx 1\pm\gamma L
$$

$$
Y_E=\frac{1}{Z_0\gamma L}
$$

$$
\frac{Y'}{2}=0
$$

$$
I_E=\frac{E}{Z}
$$

$$
V_i=-V_k=-\frac{LE}{2}=V
$$

$$
V(x)=\frac{V}{2L}(2+\gamma L)(2x-L)
$$