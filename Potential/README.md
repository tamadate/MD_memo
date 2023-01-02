# 1. Pair potential
## 1.1 Lennard-Jones (12-6)
The Lenard-Jones (12-6) pair potential and its forces working on the atoms $i$ and $j$ are expressed as below:
$$
\begin{aligned}
    U_{LJ}&=4 \epsilon_{ij}
        \{({\sigma_{ij} \over \vec{r}_{ij}})^{12}-
        ({\sigma_{ij} \over \vec{r}_{ij}})^{6}\}\\
    \vec{F}_i&=-{\partial U_{LJ} \over \partial \vec{r}_i}=
        48\epsilon_{ij}{\sigma_{ij}^{12} \over \vec{r}_{ij}^{13}}
        {\partial \vec{r}_{ij} \over \partial \vec{r}_i}
        -24\epsilon_{ij}{\sigma_{ij}^{6} \over \vec{r}_{ij}^{7}}
        {\partial \vec{r}_{ij} \over \partial \vec{r}_i}\\
    \vec{F}_j&=-{\partial U_{LJ} \over \partial \vec{r}_j}=
        48\epsilon_{ij}{\sigma_{ij}^{12} \over \vec{r}_{ij}^{13}}
        {\partial \vec{r}_{ij} \over \partial \vec{r}_j}
        -24\epsilon_{ij}{\sigma_{ij}^{6} \over \vec{r}_{ij}^{7}}
        {\partial \vec{r}_{ij} \over \partial \vec{r}_j}\\
\end{aligned}
$$
Where, $\epsilon_{ij}$ and $\sigma_{ij}$ are the Lennard-Jones potential parameters and $\vec{r}_{ij}$ is the distance between atom $i$ and $j$. Since $\vec{r}_{ij}=\vec{r}_j-\vec{r}_i$, the derivatives can be calculated as ${\partial \vec{r}_{ij} \over \partial \vec{r}_i}=-1$ and ${\partial \vec{r}_{ij} \over \partial \vec{r}_j}=1$, hence,
$$
\begin{aligned}
    \vec{F}_i&=-48\epsilon_{ij}{\sigma_{ij}^{12} \over \vec{r}_{ij}^{13}}
        +24\epsilon_{ij}{\sigma_{ij}^{6} \over \vec{r}_{ij}^{7}}\\
    \vec{F}_j&=48\epsilon_{ij}{\sigma_{ij}^{12} \over \vec{r}_{ij}^{13}}
        -24\epsilon_{ij}{\sigma_{ij}^{6} \over \vec{r}_{ij}^{7}}\\
\end{aligned}
$$
In the program, the parameters are defined as $A_{ij}=48\epsilon_{ij}\sigma_{ij}^{12}$ and $B_{ij}=24\epsilon_{ij}\sigma_{ij}^{6}$ and call the parameters in the pair interaction calculation loop to save the computation time.
## 1.2 Coulomb
The Coulombic pair potential and its forces are expressed as below:
$$
\begin{aligned}
    U_{Coul}&={q_iq_j \over 4 \pi \epsilon_0}{1 \over \vec{r}_{ij}}\\
    \vec{F}_i&=-{\partial U_{LJ} \over \partial \vec{r}_i}=
        {q_iq_j \over 4 \pi \epsilon_0}{1 \over \vec{r}_{ij}^2}
        {\partial \vec{r}_{ij} \over \partial \vec{r}_i}\\
    \vec{F}_j&=-{\partial U_{LJ} \over \partial \vec{r}_j}=
        {q_iq_j \over 4 \pi \epsilon_0}{1 \over \vec{r}_{ij}^2}
        {\partial \vec{r}_{ij} \over \partial \vec{r}_j}\\
\end{aligned}
$$
Where, $\epsilon_0$ is the permittivity of vacuum, $q_i$ and $q_j$ are the partial charges of atoms $i$ and $j$, $\vec{r}_{ij}$ is the distance of atoms $i$ and $j$. Since $\vec{r}_{ij}=\vec{r}_j-\vec{r}_i$, the derivatives can be calculated as ${\partial \vec{r}_{ij} \over \partial \vec{r}_i}=-1$ and ${\partial \vec{r}_{ij} \over \partial \vec{r}_j}=1$, hence,
$$
\begin{aligned}
    \vec{F}_i&=-{q_iq_j \over 4 \pi \epsilon_0}{1 \over \vec{r}_{ij}^2}\\
    \vec{F}_j&={q_iq_j \over 4 \pi \epsilon_0}{1 \over \vec{r}_{ij}^2}\\
\end{aligned}
$$
While the electrical interaction is able to directory calculate via above equation, the force working for long distance (proportional to $r^{-2}$ which is not converged immediately). Hence, in MD simulation the Ewaled method is generally applied for Coulomb interaction calculation. However, the Ewald method is not applicable if the total system charge is not neutral. Also, the Ewald method is not best choice for the infinity diluted particle calculation system, e.g., aerosol simulation.

# 2. Many body interaction
## 2.1 EAM
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
        {1 \over 2}-{1 \over 2}sin({\pi \over 2}{r-R \over D}) & (R-D<r<R+D) \\
        0 & (r>R+D)
    \end{cases} \\
    f_R(r_{ij})&=Aexp(-\lambda_1r_{ij})\\
    f_A(r_{ij})&=-Bexp(-\lambda_2r_{ij})\\
    b_{ij}&=(1+\beta^n\zeta_{ij}^{n})^{-{1 \over 2n}}\\
    \zeta_{ij}&=\sum_{k \neq j}f_C(r_{ik})g(\theta_{jik})exp[\lambda_3^3(r_{ij}-r_{ik})^3]\\
    g(\theta_{ijk})&=1+{c^2 \over d^2}-{c^2 \over d^2+(h-cos(\theta_{jik})^2}\\
\end{aligned}
$$
$f_C(r_{ij})$ is the smoothing function. The first term of right hand of the potential equations $f_R(r_{ij})$ is the repulsion potential and second term $b_{ij}f_A(r_{ij})$ is the attractive potential working between atoms $i$ and $j$ but the attractive term is depended on the third atoms, $k$ locations. $R$, $D$, $A$, $B$, $\lambda_1$, $\lambda_2$, $\lambda_3$, $n$, $c$, $d$, and $h$ are the Tersoff potential parameters. In Tersoff potential, $r_{ij}=|\vec{r}_{ij}|$ and $r_{ik}=|\vec{r}_{ik}|$ is used.
## 2.2 Stillinger-Weber

# 3. Bond potential
## 3.1 Harmonic
The harmonic bond potential and forces working to two atoms are expressed as below:
$$
\begin{aligned}
U_{bond}&=k_{bond}(\vec{r}_{ji}-r_0)^2\\
\vec{F}_i&=-{\partial U_{bond} \over \partial \vec{r}_i}=-2k_{bond}(\vec{r}_{ji}-r_0){\partial \vec{r}_{ji} \over \partial \vec{r}_i}\\
\vec{F}_j&=-{\partial U_{bond} \over \partial \vec{r}_j}=-2k_{bond}(\vec{r}_{ji}-r_0){\partial \vec{r}_{ji} \over \partial \vec{r}_j}\\
\end{aligned}
$$
Since $\vec{r}_{ji}=\vec{r}_i-\vec{r}_j$, the derivatives can be calculated as ${\partial \vec{r}_{ji} \over \partial \vec{r}_i}=1$ and ${\partial \vec{r}_{ji} \over \partial \vec{r}_j}=-1$, hence,
$$
\begin{aligned}
\vec{F}_i&=-2k_{bond}(\vec{r}_{ji}-r_0)\\
\vec{F}_j&=2k_{bond}(\vec{r}_{ji}-r_0)
\end{aligned}
$$

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
    \vec{F}_i&=-{\partial U_{angle} \over \partial \vec{r}_i}=-2k_{angle}(\theta_{ijk}-\theta_0){\partial \theta_{ijk} \over \partial \vec{r}_i}\\
    \vec{F}_k&=-{\partial U_{angle} \over \partial \vec{r}_{k}}=-2k_{angle}(\theta_{ijk}-\theta_0){\partial \theta_{ijk} \over \partial \vec{r}_k}\\
    \vec{F}_j&=-\vec{F}_i-\vec{F}_k\\
\end{aligned}
$$

$\theta_{ijk}$ is given as following equation.
$$
\begin{aligned}
    \theta_{ijk}&=cos^{-1}({\vec{r}_{ji} \cdot \vec{r}_{jk} \over |\vec{r}_{ji}||\vec{r}_{jk}|})\\
\end{aligned}
$$
Its derivative of $\vec{r}_i$, ${\partial \theta_{ijk}\over \partial \vec{r}_i}$ is (there may be better way to derive this)
$$
\begin{aligned}
    {\partial \theta_{ijk} \over \partial \vec{r}_i}&=
    {\partial \theta_{ijk} \over \partial |r_{ji}|}
    {\partial |r_{ji}| \over \partial \vec{r}_{ji}}
    {\partial \vec{r}_{ji} \over \partial \vec{r}_i}\\

    {\partial \theta_{ijk} \over \partial |r_{ji}|}
    &={-1 \over \sqrt{1-cos^2\theta_{ijk}}}
    ({\partial \vec{r}_{ji} \over \partial |\vec{r}_{ji}|} 
    \cdot {\vec{r}_{jk} \over |\vec{r}_{ji}||\vec{r}_{jk}|}
    +{-1 \over |\vec{r}_{ji}|^2}{\vec{r}_{ji} \cdot \vec{r}_{jk} \over |\vec{r}_{jk}|})\\
    &={-1 \over sin\theta_{ijk}}
    ({|\vec{r}_{ji}| \over \vec{r}_{ji}} 
    \cdot {\vec{r}_{jk} \over |\vec{r}_{ji}||\vec{r}_{jk}|}
    +{-1 \over |\vec{r}_{ji}|}cos\theta_{ijk})\\
    {\partial |\vec{r}_{ji}| \over \partial \vec{r}_{ji}}
    &={\vec{r}_{ji} \over |r_{ji}|}\\
    {\partial \vec{r}_{ji} \over \partial \vec{r}_i}
    &={\partial (\vec{r}_i-\vec{r}_j) \over \partial \vec{r}_i}=1\\

    {\partial  \theta_{ijk} \over \partial \vec{r}_i}
    &={-1 \over sin\theta_{ijk}}
    ({|\vec{r}_{ji}| \over \vec{r}_{ji}} 
    \cdot {\vec{r}_{jk} \over |\vec{r}_{ji}||\vec{r}_{jk}|}
    -{1 \over |\vec{r}_{ji}|}cos\theta_{ijk})
    {\vec{r}_{ji} \over |r_{ji}|}\\
    &={-1 \over |\vec{r}_{ji}|sin\theta_{ijk}}
    ({\vec{r}_{jk} \over |\vec{r}_{jk}|}
    -{\vec{r}_{ji} \over |\vec{r}_{ji}|}cos\theta_{ijk})\\
\end{aligned}
$$

Similally, the derivative of $\vec{r}_k$, ${\partial \theta_{ijk}\over \partial \vec{r}_k}$ is calculatable
$$
\begin{aligned}
    {\partial  \theta_{ijk} \over \partial \vec{r}_k}
    &={-1 \over |\vec{r}_{jk}|sin\theta}
    ({\vec{r}_{ji} \over |\vec{r}_{ji}|}
    -{\vec{r}_{jk} \over |\vec{r}_{jk}|}cos\theta_{ijk})\\
\end{aligned}
$$

Therefore the forces are
$$
\begin{aligned}
    \vec{F}_i&=2k_{angle}(\theta_{ijk}-\theta_0)
    {1 \over |\vec{r}_{ji}|sin\theta}
    ({\vec{r}_{jk} \over |\vec{r}_{jk}|}
    -{\vec{r}_{ji} \over |\vec{r}_{ji}|}cos\theta_{ijk})\\
    \vec{F}_k&=2k_{angle}(\theta_{ijk}-\theta_0)
    {1 \over |\vec{r}_{jk}|sin\theta}
    ({\vec{r}_{ji} \over |\vec{r}_{ji}|}
    -{\vec{r}_{jk} \over |\vec{r}_{jk}|}cos\theta_{ijk})\\
    \vec{F}_j&=-\vec{F}_i-\vec{F}_k\\
\end{aligned}
$$

# 5. Dihedral potential