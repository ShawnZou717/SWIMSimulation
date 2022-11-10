# CalculatingBackscatteringCoefficientsForEachGrid
This function uses a simplified Cox-Munk coefficient model to calculate BC (Backscattering Coefficient) for each grid of the splitted ocean. The model comes from **Eq(2.21)** in the Doctoral dissertation of Chu XiaoQing.

* 储小青. 海浪波谱仪海浪遥感方法及应用基础研究[D].中国科学院研究生院（海洋研究所）,2011.

# 0. Preliminary Knowledge
The simplied Cox-Munk model of BC is:
$$
\sigma^{0}(\theta,\phi)=
\frac{\left| R(0) \right|^{2}}{2\sigma_{u}\sigma_{c}cos^{4}\theta}\cdot
exp\left[ -\frac{tan^{2}\theta}{2\sigma_{\phi}^{2}} \right]\\
\frac{1}{\sigma^{2}_{\phi}}=\frac{cos^{2}\phi}{\sigma^{2}_{u}}+\frac{sin^{2}\phi}{\sigma_{c}^2}
$$
$\theta$ denotes the incident angle. $\phi$ denotes the relative scanning direction with respect to the peak wave direction. $\left| R(0) \right|^{2}$ represents the Fresnel reflection ratio, for $20^oC$ sea water, $\left| R(0) \right|^{2}=0.61$ for EM waves at Ku band. $\sigma_{u}^2$ denotes the variance of sea surface slope along wind direction. $\sigma_{c}^2$ denotes the variance of sea surface slope perpendicular to wind direction.

# 1. Determination of Parameters
## 1.1 Variance of Sea Surface Slope
Theoretically, the key to calculate the variance of sea slope is through the integration over its pdf:
$$
\sigma_{u}^{2}=\iint\xi_{x}^{2}\cdot p(\xi_{x},\xi_y)d\xi_{x}d\xi_{y}\\
\sigma_{c}^{2}=\iint\xi_{y}^{2}\cdot p(\xi_{x},\xi_y)d\xi_{x}d\xi_{y}
$$
Now, it is impossible to know the exact pdf of an actual sea surface. But let's review the definition of the DWS. Assuming the sea surface to be ergodic and a weak stationary process in space and time, the DWS (directional wave spectrum) is defined as the Fourier transform of autocorrelation of sea surface:
$$
F(\overrightarrow{k})=FT\left(\left< \xi(\overrightarrow{x_1})\xi(\overrightarrow{x_2}) \right>\right)
=FT\left(X\left( \overrightarrow{\Delta x} \right)\right)\\
\overrightarrow{\Delta x} = \overrightarrow{x_2}-\overrightarrow{x_1}
$$
Thus, the autocorrelation of sea surface denotes the inverse Fourier transform of DWS:
$$
\left< \xi(\overrightarrow{x_1})\xi(\overrightarrow{x_2}) \right>=\iint F(\overrightarrow{k}) e^{j\overrightarrow{k}\overrightarrow{\Delta x}} d \overrightarrow{k}
$$
Now, take two partial differentiations of the autocorrelation on x component of $\overrightarrow{x_1}$ and $\overrightarrow{x_2}$ respectively, assuming the switchable operation between partial differentiation and set average, we get:
$$
\left<\frac{
  \partial^{2} \xi(\overrightarrow{x_1})\xi(\overrightarrow{x_2}) 
}{\partial x_1\partial x_2}\right>=
\iint k_{x}^{2}F(\overrightarrow{k})e^{j\overrightarrow{k}\overrightarrow{\Delta x}}d \overrightarrow{k}
$$
Surprisingly, the autocorrelation of sea suraface slope can be easily obtained by the integration over slope spectrum. And by setting $\overrightarrow{\Delta x}=0$, the partial differential of the set average on x and y component equal to $\sigma_{u}^{2}$ and $\sigma_{c}^{2}$ respectively. i. e.
$$
\left<
\xi_{x}^{2}\right>=\sigma_{u}^{2}=\iint k_{x}^{2}F(\overrightarrow{k})d \overrightarrow{k}\\
\left<
\xi_{y}^{2}\right>=\sigma_{c}^{2}=\iint k_{y}^{2}F(\overrightarrow{k})d \overrightarrow{k}
$$
It's so neat that you may feel this urge to apply it immediately for BC calculation. However, **physics always set a boundary for mathematical equations**. In our case, the entire detection system works based on the **Kirchhoff approximation** who requires $kRcos^{3}\theta>>1$. $k$ denotes the working wavenumber of radar, $R$ denotes the sea surface curvature radius, $\theta$ denotes the incident angle. Limited by that, only a certain interval of slope spectrum whose wavenumber position satisfies the constraint condition contributes to the variance of sea surface slope. That certian interval of wavenumber is called cut-off wavenumber $k_d$. Thus with a given $k_d$, the variance can be figured by:
$$
\sigma_{i}^{2}=\iint_{|\overrightarrow{k}|<k_d} k_{i}^{2}F(\overrightarrow{k})d \overrightarrow{k}
$$

Luckly, some empirical cut-off wavenumber has been established as a function of wind speed at 10m altitude. You can check Chu Xiaoqing's doctoral dissertation page.29-30 for more details. Since our detection system works at Ku band, only the Ku band function is used:
$$
k_d=171.844-29.9715\cdot U_{10}+11.0462\cdot U_{10} ln(U_{10})-0.8727\cdot U_{10}^{2}+0.13932\cdot U_{10}^{2}ln(U_{10})
$$

$U_{10}$ represents the wind speed at 10m altitude. Other than integrating slope spectrum with cut-off wavenumber, there is a even simpler way to calculate the variance of sea surface: a straightforward empirical variance as a function of wind speed usually in a linear form:

$$
\sigma_{u}^{2}=0.00078545\times U_{10}+0.0092407\\
\sigma_{c}^{2}=0.00052799\times U_{10}+0.0097295
$$

 Since the simulation usually takes experiments for hundreds of time, the direct empirical equation is used to speed up the process **instead of integration method discussed above**. That is also why the simulation software requires at least one wind wave model during experiment. 

At the end, though the integration method is not coded in the software by default, feel free to add it into the source code if you tend to play with the statistics of sea surface slope and see what happens.

## 1.2 Incident Angle
The incident angle in the BC equation doesn't mean the incident angle of the radar but represents the incident angle with respect to normal vector of sea surface at each grid. Fig.1 illusrate this geometric relation.

<img src="C:\Users\ZJ\Pictures\GeometricRelation.png" alt="OceanOFSymmetricDWS.jpg"/>
<center><p><b>Fig.1 Definition of Incident Angle</b></p></center>

The incident angle $\theta$ can be easily acquired by vector product:
$$
cos\theta=-\frac{\overrightarrow{n}\cdot\overrightarrow{r}}{\left|\overrightarrow{n}\right|\cdot\left|\overrightarrow{r}\right|}
$$

Fig.2 illustrates how the local coordinate of each scanning is established. The origin of coordinate relies on the center point of 3dB beam projected on the horizontal plane. X axis extends along the scanning direction and y axis is perpendicular to x axis. Now, the radar's incident angle $\Xi$, flying height $H$, surface elevation $\xi(x,y)$ are all given. It is not hard to get the position vectoc of a given surface grid:

$$
\overrightarrow{r}=\left(x-H\cdot tan(\Xi), y, \xi-H \right)
$$

The normal vector is also easy to calculate. Let's recall that for a given curve function $F(x,y,z)=0$, the normal vector at a given position is $\left(\partial F/\partial x, \partial F/\partial y, \partial F/\partial z\right)$. If you want to learn more about this please refer to [Wolfram MathWorld content about normal vector](https://mathworld.wolfram.com/NormalVector.html). In our case, since $z=\xi(x,y)$, the normal vector of sea surface is:
$$\overrightarrow{n}=\left(\partial \xi/\partial x, \partial \xi/\partial y, -1\right)$$ Substituting the 2 vectors into the $cos\theta$'s equation, then the rest of work is ony for computer to do.

<img src="C:\Users\ZJ\Pictures\LocalCoordinate.png" alt="OceanOFSymmetricDWS.jpg"/>
<center><p><b>Fig.2 Establishment of Coordinate</b></p></center>

