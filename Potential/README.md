# 1. Pair potential
## 1.1 Lennard-Jones (12-6)
The Lenard-Jones (12-6) pair potential and its forces working on the atoms $i$ and $j$ are expressed as below:

$$
\begin{aligned}
    U_{LJ}&=4 \epsilon_{ij}
        \left[
            \left({\sigma_{ij} \over |\vec{r_{ij}}|}\right)^{12}-
            \left({\sigma_{ij} \over |\vec{r_{ij}}|}\right)^{6}
        \right]\\
    \vec{F}_i&=-{\partial U_{LJ} \over \partial \vec{r_{i}}}=
        48\epsilon_{ij}{\sigma_{ij}^{12} \over |\vec{r_{ij}}|^{13}}
        {\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}
        -24\epsilon_{ij}{\sigma_{ij}^{6} \over |\vec{r_{ij}}|^{7}}
        {\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}\\
    \vec{F}_j&=-{\partial U_{LJ} \over \partial \vec{r_{j}}}=
        48\epsilon_{ij}{\sigma_{ij}^{12} \over |\vec{r_{ij}}|^{13}}
        {\partial |\vec{r_{ij}}| \over \partial \vec{r_{j}}}
        -24\epsilon_{ij}{\sigma_{ij}^{6} \over |\vec{r_{ij}}|^{7}}
        {\partial |\vec{r_{ij}}| \over \partial \vec{r_{j}}}\\
\end{aligned}
$$

Where, $\epsilon_{ij}$ and $\sigma_{ij}$ are the Lennard-Jones potential parameters and $\vec{r_{ij}}$ is the distance between atoms $i$ and $j$.  And derivarivatives are

$$
{\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}=
{\partial |\vec{r_{ij}}| \over \partial \vec{r_{ij}}}{\partial \vec{r_{ij}} \over \partial \vec{r_{i}}}\\
{\partial |\vec{r_{ij}}| \over \partial \vec{r_{j}}}=
{\partial |\vec{r_{ij}}| \over \partial \vec{r_{ij}}}{\partial \vec{r_{ij}} \over \partial \vec{r_{j}}}
$$

Where, ${\partial |\vec{r_{ij}}| \over \partial \vec{r_{ij}}}=
{\vec{r_{ij}} \over |\vec{r_{ij}}|}$ (see appendix 1)
and $\vec{r_{ij}}=\vec{r_{j}}-\vec{r_{i}}$ make last derivatives to ${\partial \vec{r_{ij}} \over \partial \vec{r_{i}}}=-1$ and ${\partial \vec{r_{ij}} \over \partial \vec{r_{j}}}=1$, respectively.  In summary, the equation is

$$
\begin{aligned}
    \vec{F}_i&=-48\epsilon_{ij}{\sigma_{ij}^{12} \over |\vec{r_{ij}}|^{13}}
        {\vec{r_{ij}} \over |\vec{r_{ij}}|}
        +24\epsilon_{ij}{\sigma_{ij}^{6} \over |\vec{r_{ij}}|^{7}}
        {\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
    \vec{F}_j&=48\epsilon_{ij}{\sigma_{ij}^{12} \over |\vec{r_{ij}}|^{13}}
        {\vec{r_{ij}} \over |\vec{r_{ij}}|}
        -24\epsilon_{ij}{\sigma_{ij}^{6} \over |\vec{r_{ij}}|^{7}}
        {\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
\end{aligned}
$$

In the program, the parameters are defined as $A_{ij}=48\epsilon_{ij}\sigma_{ij}^{12}$ and $B_{ij}=24\epsilon_{ij}\sigma_{ij}^{6}$ and call the parameters in the pair interaction calculation loop to save the computation time.
## 1.2 Coulomb
The Coulombic pair potential and its forces are expressed as below:

$$
\begin{aligned}
    U_{Coul}&={q_iq_j \over 4 \pi \epsilon_0}{1 \over |\vec{r_{ij}}|}\\
    \vec{F}_i&=-{\partial U_{LJ} \over \partial \vec{r_{i}}}=
        {q_iq_j \over 4 \pi \epsilon_0}{1 \over |\vec{r_{ij}}|^2}
        {\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}\\
    \vec{F}_j&=-{\partial U_{LJ} \over \partial \vec{r_{j}}}=
        {q_iq_j \over 4 \pi \epsilon_0}{1 \over |\vec{r_{ij}}|^2}
        {\partial |\vec{r_{ij}}| \over \partial \vec{r_{j}}}\\
\end{aligned}
$$

Where, $\epsilon_0$ is the permittivity of vacuum, $q_i$ and $q_j$ are the partial charges of atoms $i$ and $j$.  Similar to Lennard-Jones potential, transforming the last derivatives yields:

$$
\begin{aligned}
    \vec{F}_i&=-{q_iq_j \over 4 \pi \epsilon_0}{1 \over |\vec{r_{ij}}|^2}
    {\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
    \vec{F}_j&={q_iq_j \over 4 \pi \epsilon_0}{1 \over |\vec{r_{ij}}|^2}
    {\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
\end{aligned}
$$

While the electrical interaction is able to directory calculate via above equation, the force working for long distance (proportional to $r^{-2}$ which is not converged immediately). Hence, in MD simulation the Ewaled method is generally applied for Coulomb interaction calculation. However, the Ewald method is not applicable if the total system charge is not neutral. Also, the Ewald method is not best choice for the infinity diluted particle calculation system, e.g., aerosol simulation.

## 1.3 Buckingham
The Buckingham pair potential is described below

$$
    U _{Buck}=A\exp \left(-B|\vec{r_{ij}|}\right)-{\frac {C}{|\vec{r_{ij}|}^{6}}}
$$

Where, A, B, and C are potential parameters.  The first and second terms of right hand show the repulsion and attractive potentials, respectively.  Its force is:

$$
\begin{aligned}
    \vec{F}_i&=-{\partial U_{Buck} \over \partial \vec{r_{i}}}=
        -AB\exp \left(-B\vec{r_{ij}}\right){\vec{r_{ij}} \over |\vec{r_{ij}}|}
        +6{\frac {C}{\vec{r_{ij}}^{7}}}{\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
    \vec{F}_j&=-{\partial U_{Buck} \over \partial \vec{r_{j}}}=
        AB\exp \left(-B\vec{r_{ij}}\right){\vec{r_{ij}} \over |\vec{r_{ij}}|}
        -6{\frac {C}{\vec{r_{ij}}^{7}}}{\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
\end{aligned}
$$


# 2. Many body interaction
## 2.1 EAM

$$
{U_{EAM}=F_{\alpha}\left(\sum _{j\neq i}\rho _{\beta }(r_{ij})\right)+{\frac {1}{2}}\sum _{j\neq i}\phi _{\alpha \beta }(r_{ij})}
$$

## 2.2 Tersoff
The Tersoff potential is described as

$$
\begin{aligned}
    U_{Tersoff}&={1 \over 2}\sum_i{\sum_{j\neq i}}
    f_C(r_{ij})
    [f_R(r_{ij})+b_{ij}f_A(r_{ij})]\\
    f_C(r_{ij})&=
    \begin{cases}
        1 & (r<R-D) \\
        {1 \over 2}-{1 \over 2}sin\left({\pi \over 2}{r-R \over D}\right) & (R-D<r<R+D) \\
        0 & (r>R+D)
    \end{cases} \\
    f_R(r_{ij})&=Aexp(-\lambda_1r_{ij})\\
    f_A(r_{ij})&=-Bexp(-\lambda_2r_{ij})\\
    b_{ij}&=(1+\beta^n\zeta_{ij}^{n})^{-{1 \over 2n}}\\
    \zeta_{ij}&=\sum_{k \neq j}f_C(r_{ik})g(\theta_{jik})exp\left[\lambda_3^3(r_{ij}-r_{ik})^3\right]\\
    g(\theta_{ijk})&=1+{c^2 \over d^2}-{c^2 \over d^2+\left[h-cos(\theta_{jik})\right]^2}\\
\end{aligned}
$$

$f_C(r_{ij})$ is the smoothing function. The first term of right hand of the potential equations $f_R(r_{ij})$ is the repulsion potential and second term $b_{ij}f_A(r_{ij})$ is the attractive potential working between atoms $i$ and $j$ but the attractive term is depended on the third atoms, $k$ locations. $R$, $D$, $A$, $B$, $\lambda_1$, $\lambda_2$, $\lambda_3$, $n$, $c$, $d$, and $h$ are the Tersoff potential parameters. In Tersoff potential, $r_{ij}=|\vec{r_{ij}}|$ and $r_{ik}=|\vec{r_{ik}}|$ are used.
## 2.2 Stillinger-Weber

# 3. Bond potential
## 3.1 Harmonic
The harmonic bond potential and forces working to two atoms are expressed as below:

$$
\begin{aligned}
U_{bond}&=k_{bond}(|\vec{r_{ij}|}-r_0)^2\\
\vec{F}_i&=-{\partial U_{bond} \over \partial \vec{r_{i}}}=-2k_{bond}(|\vec{r_{ij}}|-r_0){\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}\\
\vec{F}_j&=-{\partial U_{bond} \over \partial \vec{r_{j}}}=-2k_{bond}(|\vec{r_{ij}}|-r_0){\partial |\vec{r_{ij}}| \over \partial \vec{r_{j}}}\\
\end{aligned}
$$

Since ${\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}=-{\vec{r_{ij}} \over |\vec{r_{ij}}|}$ and ${\partial |\vec{r_{ij}}| \over \partial \vec{r_{j}}}={\vec{r_{ij}} \over |\vec{r_{ij}}|}$ (see Lennard-Jones potential),

$$
    \begin{aligned}
        \vec{F}_i&=2k_{bond}(\vec{r_{ji}}-r_0)
            {\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
        \vec{F}_j&=-2k_{bond}(\vec{r_{ji}}-r_0)
            {\vec{r_{ij}} \over |\vec{r_{ij}}|}\\
    \end{aligned}
$$

## 3.2 Morse

# 4. Angle potential
## 4.1 Harmonic
The harmonic angular potential is expressed as below:
$$
\begin{aligned}
    U_{angle}&=k_{angle}(\theta_{ijk}-\theta_0)^2\\
\end{aligned}
$$

The forces working to each atoms, $i$, $j$, and $k$ are

$$
\begin{aligned}
    \vec{F}_i&=-{\partial U_{angle} \over \partial \vec{r_{i}}}=-2k_{angle}(\theta_{ijk}-\theta_0){\partial \theta_{ijk} \over \partial \vec{r_{i}}}\\
    \vec{F}_k&=-{\partial U_{angle} \over \partial \vec{r}_{k}}=-2k_{angle}(\theta_{ijk}-\theta_0){\partial \theta_{ijk} \over \partial \vec{r_{k}}}\\
    \vec{F}_j&=-\vec{F}_i-\vec{F}_k\\
\end{aligned}
$$

$\theta_{ijk}$ is given as following equation.

$$
\begin{aligned}
    \theta_{ijk}&=cos^{-1}\left({\vec{r_{ji}} \cdot \vec{r_{jk}} \over |\vec{r_{ji}}||\vec{r_{jk}}|}\right)\\
\end{aligned}
$$

Its derivative of $\vec{r_{i}}$, ${\partial \theta_{ijk}\over \partial \vec{r_{i}}}$ is (there may be better way to derive this)

$$
\begin{aligned}
    {\partial \theta_{ijk} \over \partial \vec{r_{i}}}&=
    {\partial \theta_{ijk} \over \partial |r_{ji}|}
    {\partial |r_{ji}| \over \partial \vec{r_{i}}}\\
\end{aligned}
$$

Where, first component of right hand, ${\partial \theta_{ijk} \over \partial |r_{ji}|}$ is

$$
\begin{aligned}
    {\partial \theta_{ijk} \over \partial |r_{ji}|}
    &={-1 \over \sqrt{1-cos^2\theta_{ijk}}}
    \left({\partial \vec{r_{ji}} \over \partial |\vec{r_{ji}}|}
    \cdot {\vec{r_{jk}} \over |\vec{r_{ji}}||\vec{r_{jk}}|}
    +{-1 \over |\vec{r_{ji}}|^2}{\vec{r_{ji}} \cdot \vec{r_{jk}} \over |\vec{r_{jk}}|}\right)\\
    &={-1 \over sin\theta_{ijk}}
    \left({|\vec{r_{ji}}| \over \vec{r_{ji}}}
    \cdot {\vec{r_{jk}} \over |\vec{r_{ji}}||\vec{r_{jk}}|}
    +{-1 \over |\vec{r_{ji}}|}cos\theta_{ijk}\right)\\
\end{aligned}
$$

Second component is  ${\partial |\vec{r_{ij}}| \over \partial \vec{r_{i}}}=-{\vec{r_{ij}} \over |\vec{r_{ij}}|}$, hence, the differentiation of $\theta_{ijk}$ with respect to $\vec{r_{i}}$ become

$$
\begin{aligned}
    {\partial  \theta_{ijk} \over \partial \vec{r_{i}}}
    &={-1 \over sin\theta_{ijk}}
    \left({|\vec{r_{ji}}| \over \vec{r_{ji}}}
    \cdot {\vec{r_{jk}} \over |\vec{r_{ji}}||\vec{r_{jk}}|}
    -{1 \over |\vec{r_{ji}}|}cos\theta_{ijk}\right)
    {\vec{r_{ji}} \over |r_{ji}|}\\
    &={-1 \over |\vec{r_{ji}}|sin\theta_{ijk}}
    \left({\vec{r_{jk}} \over |\vec{r_{jk}}|}
    -{\vec{r_{ji}} \over |\vec{r_{ji}}|}cos\theta_{ijk}\right)\\
\end{aligned}
$$

The differentiation of $\theta_{ijk}$ with respect to $\vec{r_{k}}$, ${\partial \theta_{ijk}\over \partial \vec{r_{k}}}$ is similarly calculated

$$
\begin{aligned}
    {\partial  \theta_{ijk} \over \partial \vec{r_{k}}}
    &={-1 \over |\vec{r_{jk}}|sin\theta}
    \left({\vec{r_{ji}} \over |\vec{r_{ji}}|}
    -{\vec{r_{jk}} \over |\vec{r_{jk}}|}cos\theta_{ijk}\right)\\
\end{aligned}
$$

Therefore the forces are

$$
\begin{aligned}
    \vec{F}_i&=2k_{angle}\left(\theta_{ijk}-\theta_0\right)
    {1 \over |\vec{r_{ji}}|sin\theta}
    \left({\vec{r_{jk}} \over |\vec{r_{jk}}|}
    -{\vec{r_{ji}} \over |\vec{r_{ji}}|}cos\theta_{ijk}\right)\\
    \vec{F}_k&=2k_{angle}(\theta_{ijk}-\theta_0)
    {1 \over |\vec{r_{jk}}|sin\theta}
    \left({\vec{r_{ji}} \over |\vec{r_{ji}}|}
    -{\vec{r_{jk}} \over |\vec{r_{jk}}|}cos\theta_{ijk}\right)\\
    \vec{F}_j&=-\vec{F}_i-\vec{F}_k\\
\end{aligned}
$$

# 5. Dihedral potential

# appendix
## appendix 1: Derivation of ${\partial |\vec{r}| \over \partial \vec{r}}={\vec{r} \over |\vec{r}|}$
Since the definition, $\vec{r}=(x,y,z)$, the differentiation of $\vec{r}$ with respect to $|\vec{r}|$ is

$$
{\partial \vec{r} \over \partial |\vec{r}|}=
    \left({\partial x \over \partial |\vec{r}|}, 
    {\partial y \over \partial |\vec{r}|},
    {\partial z \over \partial |\vec{r}|}\right)\\
$$

According to $|\vec{r}|=\sqrt{x^2+y^2+z^2}$, $x=\sqrt{|\vec{r}|^2-y^2-z^2}$.  Hence, the differentiation of $x$ with respect to $|\vec{r}|$ is

$$
    {\partial x \over \partial |\vec{r}|}=
    {1 \over 2}(|\vec{r}|^2-y^2-z^2)^{-{1 \over 2}} \cdot 2|\vec{r}|=
    {|\vec{r}| \over x}
$$

${\partial y \over \partial |\vec{r}|}={|\vec{r}| \over y}$ and ${\partial z \over \partial |\vec{r}|}={|\vec{r}| \over z}$ are also able to be derived with same way and substituting it to first equation yields:

$$
    {\partial \vec{r} \over \partial |\vec{r}|}=
    \left({|\vec{r}| \over x}, 
    {|\vec{r}| \over y},
    {|\vec{r}| \over z}\right)=
    {|\vec{r}| \over \vec{r}}\\
$$

The inverses of this equation is desired form.
