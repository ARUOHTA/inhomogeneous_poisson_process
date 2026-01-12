\title{
Estimation and Model Identification for Continuous Spatial Processes
}

\author{
By A. V. VECCHIA † \\ Colorado School of Mines, Golden, USA
}
[Received May 1987. Revised January 1988]

\section*{SUMMARY}

\begin{abstract}
Formal parameter estimation and model identification procedures for continuous domain spatial processes are introduced. The processes are assumed to be adequately described by a linear model with residuals that follow a second-order stationary Gaussian random field and data are assumed to consist of noisy observations of the process at arbitrary sampling locations. A general class of two-dimensional rational spectral density functions with elliptic contours is used to model the spatial covariance function. An iterative estimation procedure alleviates many of the computational difficulties of conventional maximum likelihood estimation for non-lattice data. The procedure is applied to several generated data sets and to an actual ground-water data set.
\end{abstract}

Keywords: ANISOTROPY; GAUSSIAN PROCESS; MEASUREMENT ERROR; NON-LATTICE DATA; SPATIAL CORRELATION; SPATIAL REGRESSION

\section*{1. INTRODUCTION}

Let \(\left\{Z(x, y),(x, y) \in R^{2}\right\}\) be a continuous domain spatial process governed by a linear regression model with spatially correlated errors,
\[
Z(x, y)=f^{\mathrm{T}}(x, y) \beta+\xi(x, y),
\]
where \(f^{\mathrm{T}}(x, y)\) is a \(1 \times r\) vector of known regressor functions, \(\beta\) is an \(r \times 1\) unknown parameter vector and \(\{\xi(x, y)\}\) is a real-valued, second-order stationary spatial stochastic process with zero mean and covariance function
\[
\Gamma(u, v)=\operatorname{cov}\{\xi(x, y), \xi(x+u, y+v)\} .
\]

This paper is concerned with the identification of a parametric model for \(\Gamma(u, v)\) and estimation of the parameters of the identified model along with the regression parameters \(\beta\). Data on which the estimation and model identification procedures are based consist of \(n\) observations from equation (1) that possibly include additive measurement error,
\[
z_{i}=Z\left(x_{i}, y_{i}\right)+\eta_{i}, \quad 1 \leqslant i \leqslant n,
\]
where the \(\eta_{i}\) are i.i.d. zero-mean random variables with common variance \(\sigma_{\eta}^{2}\) and ( \(x_{i}, y_{i}\) ) denotes the co-ordinates of the \(i\) th sampling location. Cases when sampling locations do not follow a systematic pattern are of primary concern for this study; more efficient methods that take advantage of patterns that occur in the covariance matrix of the observations (2) could be developed for systematic sampling patterns such as a rectangular lattice. In addition, lattice data would be amenable to a frequency analysis approach whereas the procedures presented here are based on a spatial domain analysis. Ripley (1981) gives a good review of procedures for analysing lattice data.

\footnotetext{
† Address for correspondence: Department of Mathematics, Colorado School of Mines, Golden, CO 80401, USA.
}

Statistical methods for the analysis of non-lattice data based on a nearest neighbour approach have been presented by Besag (1975). Discrete spatial autoregressive moving average schemes such as described in Cliff and Ord (1981) may also be extended to non-lattice data. However, such approaches do not directly model the continuous process from which the data are taken. In situations when the process being sampled represents a stochastic model of a physical system, it is often necessary to determine the local stochastic properties of the continuous process. For instance, in fluid mechanics \(\{Z(x, y)\}\) might represent a velocity potential field in which case \(\{-\partial Z / \partial x\}\) and \(\{-\partial Z / \partial y\}\) represent the velocity fields in the \(x\) and \(y\) directions. Knowledge of the stochastic properties of \(\{Z(x, y)\}\) over a continuum in space can be used to infer stochastic properties of the velocity fields, whereas the partial derivatives are not even defined for a discrete model.

Geostatistical methods originally developed by Matheron (1963, 1971) for application in the mining industry have gained wide acceptance as a general methodology for estimating spatial dependence and interpolating continuous spatial processes. The geostatistical interpolation technique called universal kriging (Journel and Huijbregts, 1978) assumes much the same model framework as equations (1) and (2), except that the errors \(\eta_{i}\) in equation (2) would then represent small-scale variability known as the nugget effect and \(\{\xi(x, y)\}\) in equation (1) would be assumed only to have stationary increments. The variance of the increments is modelled by using a variogram selected from a class of variograms developed largely from empirical considerations. In this study, \(\{\xi(x, y)\}\) is assumed to be second order stationary and a parametric class of models developed by Vecchia (1985) is used to model \(\Gamma(u, v)\). This parametric class is extremely flexible in modelling the local covariance structure of \(\{\xi(x, y)\}\) and the models are physically, rather than empirically, based. The models are reviewed in Section 2.

The primary tool in geostatistics for model selection and estimation is the empirical variogram (Journel and Huijbregts, 1978), which is basically a method of moments estimate of the true variogram. A theoretical variogram may be fitted to the empirical variogram by ad hoc methods or by a weighted least squares procedure such as described by Cressie (1985). Alternative approaches to estimation that do not involve the empirical variogram, such as maximum likelihood or minimum norm quadratic estimation, have been considered by Mardia and Marshall (1984), Marshall and Mardia (1985), Kitanidis (1983, 1987) and Stein (1987). In this paper, the maximum likelihood approach is used based on the assumption that \(\{\xi(x, y)\}\) is Gaussian. A new statistical methodology for obtaining the maximum likelihood estimates is presented in Section 3, namely a sequence of approximate likelihood functions \(L_{m}\) is defined, where \(L_{m}\) approaches the conventional likelihood function as \(m\) approaches \(n\) but is very easy to compute for small \(m\). An iterative estimation procedure is developed whereby estimates based on \(L_{1}\) are used as initial values for estimates based on \(L_{2}\), which are in turn used as initial values for estimates based on \(L_{3}\), and so on. Statistics computed at each step of the iterative procedure are used to determine when the iterative parameter estimates have converged. The estimation statistics also can be used to discriminate between specific models for \(\{\xi(x, y)\}\), thus allowing the best model from among a number of postulated models to be identified. In Section 4, the estimation and model identification procedures are applied to several generated data sets and to an actual ground-water data set.

Some computational details are presented throughout the paper to demonstrate
the feasibility of the new statistical methodology for analysing large data sets with irregularly spaced data. However, computational considerations are of secondary importance in this study; development of more efficient programs to perform the estimation will hopefully result from wider use of the new methodology.

\section*{2. PARAMETRIC CLASS OF MODELS FOR SPATIAL CORRELATION}

The models postulated here for the spatial covariance structure of \(\{\xi(x, y)\}\) are most easily stated in terms of the spectral density function
\[
S\left(k_{1}, k_{2}\right)=\int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \mathrm{e}^{\mathrm{i} u k_{1}+\mathrm{i} v k_{2}} \Gamma(u, v) \mathrm{d} u \mathrm{~d} v
\]
where it is assumed that \(\Gamma(u, v)\) is absolutely integrable so that equation (3) exists in the ordinary sense. The assumed form for the spectral density function is
\[
S\left(k_{1}, k_{2}\right)=S(\kappa)=\sigma^{2} \prod_{j=1}^{q}\left|\kappa^{2}+\theta_{j}\right|^{2 n_{j}} / \prod_{j=1}^{p}\left|\kappa^{2}+\phi_{j}\right|^{2 m_{j}}
\]
where
\[
\kappa^{2}=\left[\lambda^{-1}\left(k_{1} \cos \alpha-k_{2} \sin \alpha\right)\right]^{2}+\left[\lambda\left(k_{1} \sin \alpha+k_{2} \cos \alpha\right)\right]^{2},
\]
\(p\) is a positive integer, \(q\) is a non-negative integer, and \(n_{j}\) and \(m_{j}\) are positive integers satisfying \(\Sigma m_{j} \geqslant \Sigma n_{j}+1\). Equation (4) is the general form for a rational spectral density function of a two-dimensional process that can be transformed to a second-order isotropic, or directionally independent, process via rotation and scaling of the co-ordinate system. It is analogous to the class of rational spectral density functions for continuous parameter one-dimensional processes, which result in autoregressive moving average models (Priestley (1981), p. 283). Some of the \(\theta_{j} \mathrm{~s}\) and \(\phi_{j} \mathrm{~s}\) may be complex valued in the general case considered by Vecchia (1985). However, allowing complex-valued parameters presents computational problems with the estimation that render it impractical at this time, so for this study the parameter space of equations (4) and (5) is taken to be
\[
\theta_{j} \in R_{1}, \phi_{j} \in R_{1}^{+}, \sigma^{2} \in R_{1}^{+}, \lambda \in R_{1}^{+}, \alpha \in[0, \pi / 2)
\]
where \(R_{1}\) is the real line and \(R_{1}^{+}\)denotes the positive half-line. The \(\phi_{j}\) s need to be positive for \(S(\kappa)\) to be a valid spectral density function, and the restrictions on \(\lambda\) and \(\alpha\) result in no loss of generality in describing the shape of the elliptic form (5). Furthermore, there is no loss of generality in assuming a single scaling parameter \(\lambda\) in equation (5). One rotated axis is scaled by \(\lambda\) and the other by \(\lambda^{-1}\) so that the resulting transformation from \(\left(k_{1}, k_{2}\right)\) to the new co-ordinate system has unit Jacobian. This property is useful when making comparisons between isotropic models with \(\lambda=1\) and anisotropic models with \(\lambda \neq 1\). A similar method of scaling was used by Brewer and Mead (1986) to model spatial correlation.

The covariance function corresponding to equation (4) is given by
\[
\Gamma(u, v)=\Gamma(r)=\sigma^{2}(2 \pi)^{-1}(-1)^{M-1} \sum_{j=1}^{p}\left[\left(2 m_{j}-1\right)!\right]^{-1} \partial^{2 m_{j}-1}\left\{w_{j} K_{0}\left(r \sqrt{ } \phi_{j}\right)\right\} / \partial \phi_{j}^{2 m_{j}-1}
\]
where
\[
\begin{aligned}
& r^{2}=[\lambda(u \cos \alpha-v \sin \alpha)]^{2}+\left[\lambda^{-1}(u \sin \alpha+v \cos \alpha)\right]^{2}, \\
& M=\sum_{j=1}^{p} 2 m_{j}, \\
& w_{j}=\prod_{l=1}^{q}\left(\theta_{l}-\phi_{j}\right)^{2 n_{l}} / \prod_{\substack{l=1 \\
l \neq j}}^{p}\left(\phi_{j}-\phi_{l}\right)^{2 m_{l}},
\end{aligned}
\]
and \(K_{0}(\cdot)\) is the modified Bessel function of the second kind, order zero. To compute the process variance, \(\Gamma(0), K_{0}\left(r \sqrt{ } \phi_{j}\right)\) needs to be replaced by \(-\log \left(\sqrt{ } \phi_{j}\right)\) (Vecchia (1985), proposition 3). Methods for recursively evaluating the derivatives in equation (7) are given by Vecchia (1985) and need not be repeated here. The modified Bessel functions may be computed from any of a number of existing routines.

\section*{3. ESTIMATION}

\subsection*{3.1. Approximate Likelihood Functions}

Consider evaluation of the likelihood function of the observations (2) based on the assumptions that \(\{\xi(x, y)\}\) in equation (1) is a Gaussian process with spectral density (4) and that the measurement errors \(\eta_{i}\) are i.i.d. \(N\left(0, \sigma_{\eta}^{2}\right)\). Let the collection of observations in array form be given by
\[
z=F \beta+\xi+\eta
\]
where \(F\) is \(n \times r\) with \(i\) th row \(f_{i}^{\mathrm{T}}, z^{\mathrm{T}}=\left[z_{1}, \ldots, z_{n}\right], \xi^{\mathrm{T}}=\left[\xi_{1}, \ldots, \xi_{n}\right]\) and \(\eta^{\mathrm{T}}=\left[\eta_{1}, \ldots\right.\), \(\left.\eta_{n}\right]\). Any relationship between the order of the observations in the array (8) and their sampling locations is not important at this time; advantageous ways of ordering the observations are discussed further on. Under the Gaussian assumptions, the likelihood function of array (8) is
\[
L(z)=(2 \pi)^{-n / 2}\left(\sigma^{2} \gamma_{0}\right)^{-n / 2}\left|R+v^{2} I\right|^{-1 / 2} \exp \left[-\left(2 \sigma^{2} \gamma_{0}\right)^{-1} \varepsilon^{\mathrm{T}}\left(R+v^{2} I\right)^{-1} \varepsilon\right]
\]
where
\[
\begin{aligned}
\gamma_{0} & =\operatorname{var}\{\xi(x, y)\} / \sigma^{2} \\
R & =\operatorname{corr}(\xi) \\
v^{2} & =\sigma_{\eta}^{2} / \operatorname{var}\{\xi(x, y)\}
\end{aligned}
\]
and
\[
\varepsilon=z-F \beta .
\]

To develop the approximate likelihood functions, note that equation (9) is equivalent to
\[
L(z)=\prod_{i=1}^{n} p\left(z_{i} \mid z_{j}, 1 \leqslant j \leqslant i-1\right)
\]
where \(p(\cdot \mid \cdot)\) denotes the conditional normal probability density function. In general it might be expected that the random variables \(\left\{z_{j}, 1 \leqslant j \leqslant i-1\right\}\) would contain a
great deal of superfluous and/or redundant information in predicting \(z_{i}\) if \(i\) is large. This means that the \(i\) th conditional density in equation (10) may be nearly equivalent to \(p\left(z_{i} \mid z_{i m}\right)\) where \(z_{i m}\) is a vector consisting of a few, say at most \(m\), observations from among \(\left\{z_{j}, 1 \leqslant j \leqslant i-1\right\}\). Formally selecting which observations to include in \(z_{i m}\) according to some criterion such as maximising the multiple correlation coefficient between \(z_{i}\) and \(z_{i m}\) is not feasible because \(z_{i m}\) would then depend on both the form of the model (4) and on specific parameter values. This causes serious instabilities in iterative parameter estimation procedures and greatly increases computation times. The only logical non-model-dependent method for choosing \(z_{i m}\) is to select those observations whose sample locations are closest to \(z_{i}\) in some sense. Hence, define the approximate likelihood of order \(m\) to be
\[
L_{m}(z)=\prod_{i=1}^{n} p\left(z_{i} \mid z_{i m}\right)
\]
where \(z_{i m}\) is an array consisting of the \(\min (i-1, m)\) observations from among \(z_{1}, \ldots\), \(z_{i-1}\) that are closest to \(z_{i}\) in the sense of ordinary Euclidean distance, \(d_{i j}= \sqrt{ }\left[\left(x_{i}-x_{j}\right)^{2}+\left(y_{i}-y_{j}\right)^{2}\right]\). Ties among the \(d_{i j}, 1 \leqslant j \leqslant i-1\), may be resolved in any consistent manner. In cases when the anisotropy parameter \(\lambda\) in equation (5) is large, a better approximation to the full likelihood might be obtained by replacing \(z_{i m}\) in equation (11) with \(z_{i m}(\lambda, \alpha)\), where \(z_{i m}(\lambda, \alpha)\) consists of the \(\min (i-1, m)\) observations from among \(z_{1}, \ldots, z_{i-1}\) that lie within the smallest contour of the covariance function (7) centred at ( \(x_{i}, y_{i}\) ). However, unless \(\lambda\) and \(\alpha\) are fixed and known, the benefits of such an approach are outweighed by the added complexity of allowing \(z_{i m}\) to depend on \(\lambda\) and \(\alpha\).

Substituting in equation (11) formulae for the conditional normal density functions (Graybill, 1976) and expressing the result in terms of correlations yields
\[
L_{m}(z)=(2 \pi)^{-n / 2}\left(\sigma^{2} \gamma_{0}\right)^{-n / 2} \prod_{i=1}^{n} \omega_{i m}^{-1 / 2} \exp \left[-\left(2 \sigma^{2} \gamma_{0}\right)^{-1} \sum_{i=1}^{n} \omega_{i m}^{-1} e_{i m}^{2}\right]
\]
where
\[
\begin{aligned}
\omega_{i m} & =1+v^{2}-r_{i m}^{T}\left(R_{i m}+v^{2} I\right)^{-1} r_{i m}, \\
e_{i m} & =\varepsilon_{i}-r_{i m}^{T}\left(R_{i m}+v^{2} I\right)^{-1} \varepsilon_{i m}, \\
r_{i m}^{T} & =\operatorname{corr}\left(\xi_{i}, \xi_{i m}\right)
\end{aligned}
\]
and
\[
R_{i m}=\operatorname{corr}\left(\xi_{i m}, \xi_{i m}\right) .
\]

The \(\varepsilon_{i} \mathrm{~s}, \gamma_{0}\) and \(v\) are defined in equation (9) and the arrays \(\varepsilon_{i m}\) and \(\xi_{i m}\) are defined analogously to \(z_{i m}\).

Although equation (9) is invariant to permutations of the observations in \(z\), equation (12) is not, especially for small \(m\). Hence, it is clearly desirable to have a systematic ordering procedure so that equation (12) is a unique expression. From a purely computational standpoint it is advantageous to have small values of \(d_{i j}\) associated with small values of \(i-j, j<i\). For definiteness, we assume that the \(z_{i}\) are ordered with respect to increasing values of either \(x_{i}\) or \(y_{i}\). A plot of data locations may indicate which ordering would result in a greater positive correlation between \(d_{i j}\) and
\(i-j\). Otherwise, either ordering may be chosen. All the data sets analysed in this paper are ordered with respect to increasing \(y\) co-ordinates.

The approximate likelihood functions as expressed by equation (12) are in a suitable form for maximising with respect to the parameters, which is the subject of the next section, but first we consider some possible alternative approximations to the likelihood function. For instance, the conditioning set \(z_{i m}\) in equation (11) could be replaced with the set of \(m\) observations that are nearest to \(z_{i}\). In such an approach, which is along the lines of the pseudolikelihood technique of Besag (1975), the ordering of the observations in \(z\) would be irrelevant. However, by keeping track of order equation (12) has the advantage of approaching the true likelihood function as \(m\) increases. Yet another approach would be to let \(z_{i m}\) consist of all observations \(z_{j}\), \(1 \leqslant j \leqslant i-1\), that are within a fixed distance \(d_{m}\) of \(z_{i}\) and let the number of observations in \(z_{i m}\) vary. With lattice data, this approach would essentially lead to results equivalent to those presented later. However, with irregularly spaced data, the approximate likelihood functions based on the fixed distance approach tend to exhibit erratic fluctuations that are not conducive to the iterative estimation approach described here.

\subsection*{3.2. Iterative Maximum Likelihood Estimation}

Maximisation of equation (9) becomes impractical for moderate to large numbers of observations. Mardia and Marshall (1984) have developed a procedure for exact maximum likelihood estimation within the same model framework as equation (1). However, they mention that the procedure becomes impractical for data sets exceeding 150 points. There is a need, then, for a procedure that may be applied to large as well as to small data sets. Moreover, basing the estimation on the parametric class of covariance functions (7) causes the computational requirements for obtaining the full correlation matrix \(R\) in equation (9) to become excessive. The estimation procedure discussed in this section has been developed to alleviate these difficulties.

Consider maximisation of equation (12), for fixed \(m\), with respect to the unknown parameters. Such parameters may include \(\beta\), the spatial covariance parameters in equation (4), the scaling and rotation parameters in equation (5) and the measurement error parameter \(v\). If an isotropic model is assumed, \(\lambda=1\) and \(\alpha=0\) are fixed. Also, if there is no measurement error, \(v=0\) is fixed. Using differentiation to determine values \(\sigma_{m}^{2}\) and \(\beta_{m}\) that minimise \(-2 \log L_{m}\) with respect to \(\sigma^{2}\) and \(\beta\) yields
\[
\begin{gathered}
-2 \log L_{m}^{*}=n[1+\log (2 \pi)]+n \log \left(\sigma_{m}^{2} \gamma_{0}\right)+\sum_{i=1}^{n} \log \omega_{i m} \\
\sigma_{m}^{2}=\left(n \gamma_{0}\right)^{-1} \sum_{i=1}^{n} \omega_{i m}^{-1} \tilde{e}_{i m}^{2}
\end{gathered}
\]
and
\[
\beta_{m}=\left[\sum_{i=1}^{n} \omega_{i m}^{-1} g_{i m} g_{i m}^{\mathrm{T}}\right]^{-1}\left[\sum_{i=1}^{n} \omega_{i m}^{-1} g_{i m} h_{i m}\right]
\]
where
\[
g_{i m}=f_{i}-F_{i m}^{\mathrm{T}}\left(R_{i m}+v^{2} I\right)^{-1} r_{i m}
\]
and
\[
h_{i m}=z_{i}-r_{i m}^{\mathrm{T}}\left(R_{i m}+v^{2} I\right)^{-1} z_{i m} .
\]

In equation (14), \(\tilde{e}_{i m}\) denotes \(e_{i m}\) in equation (12) but with \(\beta_{m}\) replacing \(\beta\) and, in equation (16), \(F_{\text {im }}\) is \(\min (i-1, m) \times r\) with \(k\) th row \(f_{k_{i m}}^{\mathrm{T}}\) where \(k_{\text {im }}\) is the index corresponding to the \(k\) th element of \(z_{i m}\). In the iterative procedures described here, it is sometimes necessary to fix \(\beta\) at a value \(\widetilde{\beta}\) other than its conditional minimum value (15); for instance, at the ordinary least squares estimate
\[
\tilde{\beta}_{\mathrm{OLS}}=\left[\sum_{i=1}^{n} f_{i} f_{i}^{\mathrm{T}}\right]^{-1}\left[\sum_{i=1}^{n} f_{i} z_{i}\right] .
\]

In that case, equations (13) and (14) are still the basic estimation equations except that \(\tilde{e}_{i m}\) then denotes \(e_{i m}\) with \(\beta\) equal to its specified value \(\tilde{\beta}\).

Let the \(d\)-dimensional collection of unknown parameters in equation (13) be given by
\[
\underset{(d \times 1)}{\psi}=\left[\theta^{\mathrm{T}}, \phi^{\mathrm{T}}, \lambda, \alpha, v\right]^{\mathrm{T}}
\]
where \(\theta\) is \(q \times 1\) and \(\phi\) is \(p \times 1\), and it is understood that some of the parameters from among \(\lambda, \alpha\) and \(v\) may be omitted from \(\psi\) depending on the model specification. \(\gamma_{0}\) in equation (13) is simply a function of the elements of \(\psi\) obtained from equation (7) with \(r=0\) and \(\sigma^{2}=1\). Estimates obtained by numerically minimising equation (13) with respect to \(\psi\), denoted \(\hat{\psi}_{m}\), along with the resulting estimates \(\hat{\beta}_{m}\) and \(\hat{\sigma}_{m}^{2}\) obtained from equations (14) and (15) evaluated at \(\hat{\psi}_{m}\), are dubbed approximate maximum likelihood estimates of order \(m\). Parameter estimates obtained by minimising equation (13) with \(\beta\) fixed at \(\widetilde{\beta}\) are also denoted by \(\hat{\psi}_{m}\) and \(\hat{\sigma}_{m}^{2}\). To minimise equation (13), a Fortran optimisation program coded by R. B. Schnabel and based on quasi-Newton methods for unconstrained non-linear problems (Dennis and Schnabel, 1983) is used for all the applications in this paper. When \(\beta\) is included in the estimation, equation (13) is not minimised directly with \(\beta_{m}\) always given by equation (15). Rather, \(\beta_{m}\) is held fixed, say at \(\beta_{m 0}\), in equations (13) and (14) until a conditional minimum with respect to \(\psi\) is found, say at \(\psi_{m 1}\). Then, equation (15) is evaluated at \(\psi_{m 1}\) to obtain an updated value \(\beta_{m 1}\) of \(\beta_{m}\) which in turn is used in equations (13) and (14) to obtain a new estimate \(\psi_{m 2}\) of \(\psi\). This process may be repeated several times until the difference between \(\beta_{m k}\) and \(\beta_{m, k-1}\) becomes small. On the basis of application of the procedure to numerous data sets, \(k\) seldom needs to be larger than 2 or 3 .

The iterative estimation procedure starts with \(m=1\) and \(\beta\) fixed at the ordinary least squares estimate \(\widetilde{\beta}\) (equation (18)). Then \(\hat{\psi}_{1}\) is computed as described. A rough initial value for \(\psi\) may be determined by interactively evaluating \(-2 \log L_{1}^{*}\) at a few points. After \(\hat{\psi}_{1}\) is obtained, it may be used as the initial estimate for obtaining \(\hat{\psi}_{2}\), and so on, until the estimates are deemed to have converged. A statistic that is useful for determining when the estimates have converged is
\[
\Lambda_{m}=-2 \log \left[L_{m}^{*}\left(\hat{\psi}_{m}\right)\right]
\]
where the right-hand side is evaluated from equation (13) with \(\beta\) fixed at \(\widetilde{\beta}\). As \(m\) increases, \(\Lambda_{m}\) approaches \(-2 \log \hat{L}\) where \(\hat{L}\) is the maximised exact likelihood function (9). Our experience from the analysis of numerous simulated and actual data sets is that there is generally a small value of \(m\), say \(m^{\prime}\), such that fluctuations in \(\Lambda_{m}\) become negligible for \(m>m^{\prime}\). As indicated in the following applications, it is usually a simple matter to select \(m^{\prime}\) from inspection.
\(\beta\) is held fixed at \(\widetilde{\beta}\) throughout the process of selecting \(m^{\prime}\), after which \(\beta\) may be
included in the estimation, if desired, to obtain \(\widehat{\psi}_{m^{\prime}}\) and \(\widehat{\beta}_{m^{\prime}}\). Furthermore, anisotropy need be included in the model only after selecting \(m^{\prime}\). This results in considerable savings in computation time while usually having a negligible effect on the selection process. If \(\hat{\lambda}_{m}\), is considerably different from unity, it may be desirable to compare \(-2 \log \left[L_{m^{\prime}}^{*}\left(\hat{\psi}_{m^{\prime}}\right)\right]\) with \(-2 \log \left[L_{m}^{*}\left(\hat{\psi}_{m^{\prime}}\right)\right]\) for some values \(m>m^{\prime}\) to see whether a larger value of \(m^{\prime}\) is needed for the anisotropic model.

\section*{4. APPLICATIONS}

\subsection*{4.1. Simulated Data}

Six generated data sets will be analysed to demonstrate the iterative estimation procedure described in the previous sections and to demonstrate how the procedure may be used to identify the best model from among several postulated models. The particular configurations of models, sampling areas and sample sizes are summarised in Table 1. To simplify the discussion, all six data sets were generated with zero mean and no measurement error and the values \(\beta=0\) and \(v=0\) remain fixed throughout the analyses. An example with \(\beta\) and \(v\) included is presented in Section 4.2. Three different models were chosen to demonstrate the flexibility in the shape of the correlation functions arising from equations (4) and (5). The true correlation function for each of the models, as a function of scaled distance \(r\), is presented in Fig. 1 and the direction and degree of anisotropy is illustrated in Fig. 2. For each model, two different sampling areas are presented to give some idea of the effect that the density of sampling points in relation to the degree of correlation has on the estimation. Because emphasis is on unequally spaced observations, the sampling locations were generated independently from a uniform distribution over the respective sampling areas. The observations were obtained by generating i.i.d. standard normal variables and transforming them with the lower diagonal Cholesky factorisation of the appropriate covariance matrix.

Two different models will be fitted to each data set, one (model A) being the correct model and the other (model B) an imposter. The models are summarised in Table 2. A key statistic for model discrimination in time series analysis is Akaike's information criterion (Priestley (1981), p. 373). An analogous statistic is defined here for the sequential estimation procedure,
\[
A_{m}=\Lambda_{m}+2 d,
\]

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 1
Generated data sets}
\begin{tabular}{|l|l|l|l|l|l|l|}
\hline Data sets & \(S(\kappa)\) & \(\lambda^{2}\) & \(\alpha\) (deg) & \(n\) & Sampling areas † & Correlation function \\
\hline 1, 2 & \(100.0 /\left(\kappa^{2}+0.1\right)^{2}\) & 2.00 & 30 & 100 & \begin{tabular}{l}
\(30 \times 30\) \\
\(45 \times 45\)
\end{tabular} & Figs 1(a) and 2(a) \\
\hline 3, 4 & \(4.0\left(\kappa^{2}+0.0\right)^{2} /\left(\kappa^{2}+0.1\right)^{6}\) & 1.00 & - & 100 & \begin{tabular}{l}
\(\mathbf{4 0} \boldsymbol{\times} \mathbf{4 0}\) \\
\(60 \times 60\)
\end{tabular} & Figs 1(b) and 2(b) \\
\hline 5, 6 & \(1.0\left(\kappa^{2}-0.021\right)^{2} /\left(\kappa^{2}+0.015\right)^{4}\) & 2.64 & 60 & 100 & \begin{tabular}{l}
\(150 \times 150\) \\
\(225 \times 225\)
\end{tabular} & Figs 1(c) and 2(c) \\
\hline
\end{tabular}
\end{table}

\footnotetext{
† The first sampling area corresponds to the odd-numbered data set, the second to the even-numbered data set.
}

\section*{DATA SET 1}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2eacf06d-1d2d-46d3-ab54-4a081861b9d2-09.jpg?height=1322&width=1417&top_left_y=233&top_left_x=105}
\captionsetup{labelformat=empty}
\caption{Fig. 1. True and estimated scaled correlation functions for generated data sets 1,3 and 5 (the true correlation functions are on the left-hand side)}
\end{figure}
where \(\Lambda_{m}\) is defined by equation (19) and \(d\) is the number of elements in \(\hat{\psi}_{m}\). After the values of \(\boldsymbol{A}_{\boldsymbol{m}}\) begin to stabilise, the model with the smallest value is selected. Both \(\Lambda_{m}\) and \(d\) depend on the particular model being fitted to the data, although this dependence is not explicitly indicated in equation (20).

A simple yet effective discrimination procedure is first to assume an isotropic model and then to allow anisotropy after the best isotropic model has been selected. Table 3 summarises the results of the isotropic model discrimination analyses for the generated data sets. The odd-numbered data sets, which correspond to a higher density of sampling points than the even-numbered data sets, are characterised by larger fluctuations in \(A_{m}\) for small values of \(m\). This is to be expected, because neighbouring values of the process are more closely related with higher sampling density. The

\section*{DATA SET 1}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2eacf06d-1d2d-46d3-ab54-4a081861b9d2-10.jpg?height=1888&width=1417&top_left_y=215&top_left_x=107}
\captionsetup{labelformat=empty}
\caption{Fig. 2. Contour plots of true and estimated two-dimensional correlation functions for data sets 1, 3 and 5 (the true correlation functions are on the left-hand side): the contour levels are \(0.9,0.7,0.5,0.3\) and 0.1}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 2
Postulated models for generated data sets}
\begin{tabular}{|l|l|l|l|}
\hline \multirow[t]{2}{*}{Model} & \multirow[b]{2}{*}{1,2} & Data sets & \multirow[b]{2}{*}{5,6} \\
\hline & & 3, 4 & \\
\hline A & \(\sigma^{2} /\left(\kappa^{2}+\phi\right)^{2}\) & \(\sigma^{2}\left(\kappa^{2}+\theta\right)^{2} /\left(\kappa^{2}+\phi\right)^{6}\) & \(\sigma^{2}\left(\kappa^{2}+\theta\right)^{2} /\left(\kappa^{2}+\phi\right)^{4}\) \\
\hline B & \(\sigma^{2} /\left(\kappa^{2}+\phi\right)^{4}\) & \(\sigma^{2} /\left(\kappa^{2}+\phi\right)^{4}\) & \(\sigma^{2} /\left(\kappa^{2}+\phi\right)^{2}\) \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 3
Model discrimination statistics \(A_{m}\) for generated data sets}
\begin{tabular}{|l|l|l|l|l|l|l|l|}
\hline \multirow[t]{2}{*}{Data set} & \multicolumn{7}{|l|}{Model} \\
\hline & & & 4 & 6 & 8 & 10 & 100 \\
\hline \multirow[t]{2}{*}{1} & A & 501.70 & 495.06 & 493.59 & 492.47 & 492.83 & 493.11 \\
\hline & B & 510.99 & 500.62 & 500.92 & 500.75 & 500.80 & 500.73 \\
\hline \multirow[t]{2}{*}{2} & A & 565.02 & 565.06 & 565.12 & 565.33 & 564.87 & 563.91 \\
\hline & B & 573.09 & 573.99 & 574.10 & 574.78 & 574.70 & 574.82 \\
\hline \multirow[t]{2}{*}{3} & A & 290.02 & 234.05 & 233.05 & 237.46 & 234.87 & 233.94 \\
\hline & B & 288.40 & 236.41 & 235.62 & 237.29 & 235.62 & 236.59 \\
\hline \multirow[t]{2}{*}{4} & A & 378.53 & 369.49 & 364.83 & 363.51 & 364.66 & 365.25 \\
\hline & B & 376.46 & 368.44 & 366.43 & 365.31 & 365.58 & 366.20 \\
\hline \multirow[t]{2}{*}{5} & A & 215.11 & 213.86 & 207.53 & 202.01 & 205.49 & 204.92 \\
\hline & B & 222.98 & 225.97 & 222.17 & 221.74 & 221.82 & 221.51 \\
\hline \multirow[t]{2}{*}{6} & A & 243.71 & 245.44 & 245.56 & 245.37 & 245.94 & 246.91 \\
\hline & B & 251.35 & 253.70 & 253.95 & 253.87 & 253.89 & 254.89 \\
\hline
\end{tabular}
\end{table}
values of \(A_{m}\) remain relatively stable for \(m \geqslant 6\) in all the model runs except \(3 A\) and \(5 A\), which remain stable for \(m \geqslant 10\). The last column gives \(A_{100}\), which is the value obtained by maximising the full likelihood (9). The approximate likelihood functions are essentially equivalent to the full likelihood function for \(m \geqslant 10\) in all the runs and for much smaller values of \(m\) in many of the runs. There is no guarantee that \(A_{m}\) may not stabilise for several values of \(m\) and then jump. However, such behaviour is unlikely unless a plot of data locations indicates a great deal of clumping.

The discrimination statistics clearly indicate that model \(A\) is preferred in data sets 1, 2, 5 and 6 and they indicate a slight preference for model \(A\) in data sets 3 and 4. For each data set, selecting the wrong model could have serious implications when using the model for interpolation or simulation. First, consider data sets 1 and 2. It is straightforward to show that a process with spectral density function (4) is mean square differentiable up to order ( \(\Sigma 2 m_{j}-\Sigma 2 n_{j}-2\) ). Thus model \(A\) results in a process that has no mean-square derivatives, whereas model \(B\) results in a twice mean-square differentiable process. For data sets 3 and 4, both postulated models have the same smoothness properties and model \(B\) becomes identical with model \(A\) by letting \(\theta\) approach \(\phi\). Hence, discriminating between models \(A\) and \(B\) is a test for the significance of the moving average parameter \(\theta\). As evidenced by Fig. 1, \(\theta\) is important in allowing flexibility in the shape of the correlation function.

After selecting the best isotropic model, a more detailed analysis can be used to test for anisotropy. This could be done by starting over with \(m=1\) and going through
the entire iterative estimation procedure with \(\lambda\) and \(\alpha\) included in the estimation. However, an alternative procedure that is much simpler and has yielded equivalent results in all the test cases considered is to estimate the parameters of the anisotropic model with \(m\) fixed at a value \(m^{\prime}\) after which the isotropic parameter estimates are deemed to have converged. Estimates of \(\phi\) and \(\theta\) obtained from the isotropic analysis along with \(\hat{\lambda}=1.0\) can be used for initial values.

The anisotropic analyses were carried out with \(m^{\prime}=8\) for each of the generated data sets. The results presented in Table 4 show that the values of \(A_{8}\) for data sets 1, 2, 5 and 6 are significantly smaller than the corresponding values from Table 3, correctly indicating anisotropy in those data sets. Values of \(A_{8}\) for the anisotropic models fit to data sets 3 and 4 were larger than the corresponding values in Table 3, so \(\lambda\) and \(\alpha\) are not included in the final model for those data sets. Plots of the estimated correlation functions for data sets 1,3 and 5 are compared with the true functions in Figs 1 and 2.

To achieve some idea of the relative computation times involved in the estimation, a sampling of the central processor unit (CPU) times required to compute the approximate likelihood function (13) are compared in Table 5 with the CPU times required to evaluate the full likelihood function
\[
-2 \log L^{*}=n[1+\log (2 \pi)]+n \log \left(\tilde{\sigma}^{2} \gamma_{0}\right)+\log \left|R+v^{2} I\right|
\]
where
\[
\tilde{\sigma}^{2}=\left(n \gamma_{0}\right)^{-1} \tilde{\varepsilon}^{\mathbf{T}}\left(R+v^{2} I\right)^{-1} \tilde{\varepsilon} .
\]

Equations (21) and (22) were obtained from equation (9) in the same way that equations (13) and (14) were obtained from equation (12). All quadratic forms \(\boldsymbol{u}^{\mathbf{T}} \boldsymbol{M}^{-\mathbf{1}} v\) and determinants \(|M|\) required in the computations were computed through use of the Cholesky factorisation. There should be no problems in computing the approximate likelihood functions for very large data sets as long as \(m \leqslant 10\) is adequate, because the computation times and storage requirements involved in computing \(-2 \log \hat{L}_{m}\) increase in direct proportion to \(n\). Sorting costs involved in determining the indices of the nearest neighbours \(z_{i m}\) in equation (11), with \(m \leqslant 10\), are minimal because the sorting needs to be done only once and the indices saved along with the original data.

For comparisons, the likelihood functions for data set 1 were evaluated using a simple exponential covariance function, \(\Gamma(r)=\gamma_{0} \exp (-0.25 r)\), which is similar to the estimated covariance function given in the first row of Table 4. It can be seen from

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 4
Estimated models for generated data sets ( \(m=8\) )}
\begin{tabular}{|l|l|l|l|l|}
\hline Data set & \(\boldsymbol{A}_{\mathbf{8}}\) & \(\tilde{S}(\kappa)\) & \(\lambda^{2}\) & â (deg) \\
\hline 1 & 485.06 & \(108.35 /\left(\kappa^{2}+0.196\right)^{2}\) & 2.37 & 29 \\
\hline 2 & 546.54 & \(62.60 /\left(\kappa^{2}+0.043\right)^{2}\) & 2.98 & 37 \\
\hline 3 & 237.46 & \(2.71\left(\kappa^{2}-0.013\right)^{2} /\left(\kappa^{2}+0.069\right)^{6}\) & 1.00 & - \\
\hline 4 & 363.51 & \(4.52\left(\kappa^{2}-0.015\right)^{2} /\left(\kappa^{2}+0.085\right)^{6}\) & 1.00 & - \\
\hline 5 & 193.15 & \(0.72\left(\kappa^{2}-0.022\right)^{2} /\left(\kappa^{2}+0.013\right)^{4}\) & 2.58 & 71 \\
\hline 6 & 241.99 & \(0.86\left(\kappa^{2}-0.029\right)^{2} /\left(\kappa^{2}+0.023\right)^{4}\) & 5.73 & 66 \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 5
CPU times required to evaluate approximate and full likelihood functions †}
\begin{tabular}{|l|l|l|l|l|l|l|}
\hline Data set ‡ & \(m=2\) & \(m=4\) & Approximate (s) \(m=6\) & \(m=8\) & \(m=10\) & Full (s) \\
\hline 1 & 0.36 & 0.79 & 1.36 & 2.11 & 2.85 & 7.73 \\
\hline & 0.09 & 0.25 & 0.51 & 0.99 & 1.37 & 2.19 \\
\hline 3 & 0.52 & 1.20 & 2.03 & 3.02 & 4.10 & 12.01 \\
\hline 5 & 0.45 & 0.98 & 1.64 & 2.49 & 3.33 & 9.04 \\
\hline
\end{tabular}
\end{table}
† All times are on a Prime 9950 minicomputer and represent a single evaluation of the respective functions at the parameter values given in Table 4.
‡ The second line for data set 1 gives computation times assuming an exponential covariance function.

Table 5 that the disparity in computation times between the approximate and full likelihood functions is not as great for the exponential covariance. However, there is still a considerable advantage to the approximate likelihood functions. More efficient methods of evaluating equation (7), such as interpolation or approximation methods, could result in considerable savings. However, such methods are not considered in this study.

\subsection*{4.2. Actual Data}

A data set consisting of water levels in 93 observation wells from an aquifer in the Saratoga Valley in south-central Wyoming (Lenfest, 1986) is now analysed. The water levels were measured in October and November 1980 and can be assumed for all practical purposes to have been measured at a single point in time. A plot of the original data indicates a linear trend, suggesting the model
\[
Z(x, y)=\beta_{1}+\beta_{2} x+\beta_{3} y+\xi(x, y)
\]
where \(Z(x, y)\) is the true water level, in metres above sea level, and \((x, y)\) are expressed in kilometres. The ordinary least squares estimate of the trend surface is \(\widetilde{\beta}_{1}+\widetilde{\beta}_{2} x+ \tilde{\beta}_{3} y=2227.6-0.34 x-2.92 y\) and the residuals from this estimated trend surface are plotted in Fig. 3. There is a large range of about 160 m in the residuals and there apparently is a significant degree of non-randomness. With ground-water data of this type, there is generally a significant measurement error. Hence, the observations are assumed to be given by equation (2) and the measurement error parameter \(v>0\) is included in the estimation.

Proceeding as before, the first step is to identify an isotropic model for \(\{\xi(x, y)\}\). Table 6 summarises the model discrimination results assuming four different isotropic

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 6
Model discrimination statistics \(A_{m}\) for the Saratoga data}
\begin{tabular}{|l|l|l|l|l|l|l|}
\hline \multirow[t]{2}{*}{Model} & \multicolumn{2}{|c|}{\multirow{2}{*}{}} & & & & \\
\hline & \(I\) & 2 & 3 & 4 & 5 & 6 \\
\hline \(\sigma^{2} /\left(\kappa^{2}+\phi\right)^{2}\) & 694.11 & 677.21 & 665.53 & 658.55 & 654.68 & 655.21 \\
\hline \(\sigma^{2} /\left(\kappa^{2}+\phi\right)^{4}\) & 694.09 & 675.28 & 666.84 & 658.54 & 562.77 & 654.26 \\
\hline \(\sigma^{2}\left(\kappa^{2}+\theta\right)^{2} /\left(\kappa^{2}+\phi\right)^{4}\) & 694.86 & 676.01 & 667.49 & 660.84 & 653.19 & 654.02 \\
\hline \(\sigma^{2}\left(\kappa^{2}+\theta\right)^{2} /\left(\kappa^{2}+\phi\right)^{6}\) & 695.51 & 676.84 & 666.77 & 658.76 & 652.32 & 652.95 \\
\hline
\end{tabular}
\end{table}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2eacf06d-1d2d-46d3-ab54-4a081861b9d2-14.jpg?height=909&width=1404&top_left_y=196&top_left_x=107}
\captionsetup{labelformat=empty}
\caption{Fig. 3. Plot of residuals from the ordinary least squares trend surface for the Saratoga Valley data}
\end{figure}
models. The statistics for all four models stabilise by \(m=6\) after which the last model is selected as the best.

Next, the covariance parameters of the selected isotropic model were re-estimated with \(m=6\) and with \(\beta\) included, resulting in the last line of Table 7. As an illustration of the changes in parameter estimates as \(m\) increases, the approximate maximum likelihood estimates for smaller values of \(m\) are also included in Table 7. The penultimate column gives the estimated standard deviation of \(\xi(x, y)\), in metres, and the last column gives the estimated standard deviation of the measurement errors \(\eta_{i}\), also in metres. Inclusion of anisotropy in the estimation revealed that \(\lambda\) and \(\alpha\) did not significantly improve the model. A comparison of the estimated correlation functions for the isotropic and anisotropic models with \(m=6\), Fig. 4, indicates very little difference.

\begin{table}
\captionsetup{labelformat=empty}
\caption{TABLE 7
Parameter estimates for the final model for the Saratoga data}
\begin{tabular}{|l|l|l|l|l|l|l|l|l|}
\hline \(m\) & \(\boldsymbol{\beta}_{1}\) & \(\hat{\beta}_{2}\) & \(\ddot{\beta}_{3}\) & \(\ddot{\theta}\) & \(\hat{\phi}\) & \(\hat{\sigma}\) & \(\hat{\sigma} \hat{\gamma}_{o}^{1 / 2}\) & \(\hat{\nu} \hat{\sigma} \hat{\gamma}_{0}^{1 / 2}\) \\
\hline 1 & 2246.98 & -0.24 & -3.58 & -0.533 & 0.445 & 29.21 & 31.26 & 3.25 \\
\hline 2 & 2309.77 & -0.34 & -4.91 & -0.346 & 0.039 & 0.21 & 73.28 & 4.41 \\
\hline 3 & 2283.26 & -0.28 & -3.89 & -0.121 & 0.114 & 6.45 & 45.74 & 4.43 \\
\hline 4 & 2256.43 & -0.19 & -3.42 & -0.209 & 0.200 & 12.49 & 37.83 & 3.80 \\
\hline 5 & 2256.74 & -0.59 & -3.18 & -0.215 & 0.174 & 9.22 & 41.32 & 4.20 \\
\hline 6 & 2253.69 & -0.67 & -3.07 & -0.212 & 0.182 & 9.86 & 39.12 & 4.39 \\
\hline
\end{tabular}
\end{table}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2eacf06d-1d2d-46d3-ab54-4a081861b9d2-15.jpg?height=628&width=1403&top_left_y=170&top_left_x=113}
\captionsetup{labelformat=empty}
\caption{Fig. 4. Contour plots of estimated isotropic and anisotropic two-dimensional correlation functions for the Saratoga Valley data: the contour levels are \(0.9,0.7,0.5,0.3\) and 0.1}
\end{figure}

\section*{5. CONCLUDING REMARKS}

The estimation and model identification procedures of this paper are computationally feasible for non-lattice spatial data sets with virtually any sample size, provided that the data adhere reasonably closely to the assumptions of Section 1. The most restrictive assumption is that anisotropy in \(\{\xi(x, y)\}\) be completely specified via rotation and scaling of the co-ordinate system; investigation of other forms of anisotropy would require a richer class of models than equations (4) and (5). Another restrictive assumption is the choice of a specified linear model for the mean of \(Z(x, y)\) in equation (1). In situations when it is undesirable to specify such a rigid function, a resistant approach to removing the stationary trend in \(Z(x, y)\), such as the median polish approach described in Cressie (1986), may be applied. The residuals from the median polish may then be analysed using these methods with \(\beta=0\).

Although the procedures were developed largely on intuitive grounds, their performance on generated and actual data sets have consistently supported several important conclusions.
(a) For most reasonable sampling schemes, virtually all the information necessary for estimation and model identification is contained in the approximate likelihood functions \(L_{m}\) for \(1 \leqslant m \leqslant 10\).
(b) The model identification methods are robust to misspecification of anisotropy, thus allowing that the particular form of spectral density be selected assuming isotropy.
(c) The iterative estimation statistics are effective for identifying the smoothness properties of a continuous spatial process based on a sparse set of observations at point locations.
(d) The iterative estimation procedures are effective for identifying elliptical anisotropy in the data.

\section*{ACKNOWLEDGEMENTS}

This work was supported by the Water Resources Division of the US Geological Survey. The author is grateful to M. R. Karlinger, B. M. Troutman, R. L. Naff, the referees and the editor for their many helpful comments.

\section*{REFERENCES}

Besag, J. (1975) Statistical analysis of non-lattice data. Statistician, 24, 179-195.
Brewer, A. C. and Mead, R. (1986) Continuous second order models of spatial variation with application to the efficiency of field crop experiments. J. R. Statist. Soc. A, 149, 314-348.
Cliff, A. D. and Ord, J. K. (1981) Spatial Processes, Models and Applications. London: Pion.
Cressie, N. A. C. (1985) Fitting variogram models by weighted least squares. J. Int. Ass. Math. Geol., 17, 563-586. (1986) Kriging nonstationary data. J. Amer. Statist. Ass., 81, 625-634.

Dennis, J. E. and Schnabel, R. B. (1983) Numerical Methods for Unconstrained Optimization and Nonlinear Equations. Englewood Cliffs: Prentice-Hall.
Graybill, F. A. (1976) Theory and Application of the Linear Model, p. 106. North Scituate: Duxbury.
Journel, A. G. and Huijbregts, C. J. (1978) Mining Geostatistics. London: Academic Press.
Kitanidis, P. K. (1983) Statistical estimation of polynomial generalized covariance functions and hydrologic applications. Wat. Resour. Res., 19, 909-921.
(1987) Parametric estimation of covariances of regionalized variables. Wat. Resour. Res., 23, 557-567.

Lenfest, L. W., Jr (1986) Ground-water levels and use of water for irrigation in the Saratoga Valley, south-central Wyoming, 1980-81. Water-Resources Investigation Report 84-4040. Cheyenne: US Geological Survey.
Mardia, K. V. and Marshall, R. J. (1984) Maximum likelihood estimation of models for residual covariance in spatial regression. Biometrika, 71, 135-146.
Marshall, R. J. and Mardia, K. V. (1985) Minimum norm quadratic estimation of components of spatial covariance. J. Int. Ass. Math. Geol., 17, 517-525.

Matheron, G. (1963) Principles of geostatistics. Econom. Geol., 58, 1246-1266.
-(1971) The theory of regionalized variables and its applications. Cahiers du Centre de Morphologie Mathematique, No. 5. Fontainbleau: Centre de Morphologie Mathematique.
Priestley, M. B. (1981) Spectral Analysis and Time Series. London: Academic Press.
Ripley, B. D. (1981) Spatial Statistics, ch. 5. New York: Wiley.
Stein, M. L. (1987) Minimum norm quadratic estimation of spatial variograms. J. Amer. Statist. Ass., 82, 765-772.
Vecchia, A. V. (1985) A general class of models for stationary two-dimensional random processes. Biometrika, 72, 281-291.