\title{
Exact Bayesian inference in spatiotemporal Cox processes driven by multivariate Gaussian processes
}

\author{
Flávio B. Gonçalves \\ Universidade Federal de Minas Gerais, Brazil \\ and Dani Gamerman \\ Universidade Federal do Rio de Janeiro, Brazil
}
[Received March 2015. Revised April 2017]

\begin{abstract}
Summary. We present a novel inference methodology to perform Bayesian inference for spatiotemporal Cox processes where the intensity function depends on a multivariate Gaussian process. Dynamic Gaussian processes are introduced to enable evolution of the intensity function over discrete time. The novelty of the method lies on the fact that no discretization error is involved despite the non-tractability of the likelihood function and infinite dimensionality of the problem. The method is based on a Markov chain Monte Carlo algorithm that samples from the joint posterior distribution of the parameters and latent variables of the model. A particular choice of the dominating measure to obtain the likelihood function is shown to be crucial to devise a valid Markov chain Monte Carlo algorithm. The models are defined in a general and flexible way but they are amenable to direct sampling from the relevant distributions because of careful characterization of its components. The models also enable the inclusion of regression covariates and/or temporal components to explain the variability of the intensity function. These components may be subject to relevant interaction with space and/or time. Real and simulated examples illustrate the methodology, followed by concluding remarks.
\end{abstract}

Keywords: Augmented model; Dynamic Gaussian process; Intractable likelihood; Markov chain Monte Carlo sampling; Point pattern

\section*{1. Introduction}

A Cox process is an inhomogeneous Poisson process where the intensity function (IF) evolves stochastically. It is also referred to as a doubly stochastic process. Cox processes (Cox, 1955) have been extensively used in a variety of areas to model point process phenomena. Effects in Cox processes may present spatiotemporal variation to reflect the possibility of interaction between space-time and other model components. They can be traced back to log-Gaussian Cox processes (Møller et al., 1998), where a Gaussian process (GP) representation is used for the log-intensity (see also Diggle (2014), and references within).

The application of (GP-driven) Cox processes is closely related to two main problems: simulation and inference. These are difficult problems because of the infinite dimensionality of the process and the intractability of the likelihood function. Simulation is one of the main tools to tackle the inference problem, which primarily consists of estimating the unknown IF and

\footnotetext{
Address for correspondence: Flávio B. Gonçalves, Departamento de Estatística, Universidade Federal de Minas Gerais, Avda Antônio Carlos 6627, Belo Horizonte, Minas Gerais 31270-901, Brazil.
E-mail: fbgoncalves@est.ufmg.br
}
potential unknown parameters. However, prediction is often a concern, i.e. what should we expect in a future realization of the same phenomenon?

Solutions for the inference problem have required, until recently, the use of discrete approximations (see, for example, Møller et al. (1998), Brix and Diggle (2001) and Reis et al. (2013)). These represent a considerable source of error and, therefore, ought to be used with care (see Simpson et al. (2016)). Approximated (discretization-based) methods have as their main disadvantages
(a) producing biased results,
(b) it is difficult (sometimes impossible) to quantify the bias (error) involved,
(c) the level of discretization that is required to obtain good (sufficiently close to the exact) results is unknown and case specific and
(d) it is not always clear to which limit these approaches converge.

This motivates the development of exact methodologies, i.e. free from discretization errors, and helps to understand its advantages.

The importance of a given exact approach basically relies on its applicability in terms of modelling flexibility and computational cost. Exact solutions for inference on infinite dimensional processes with intractable likelihood can be found, for example, in Beskos et al. (2006), Sermaidis et al. (2013) and Goncalves et al. (2017). In this paper, the term exact refers to the fact that no discretization-based approximation is used. In particular, the methodology proposed here has Markov chain Monte Carlo (MCMC) error (Markov chain convergence plus Monte Carlo error) as the only source of inaccuracy, which is generally well understood and controlled.

One non-parametric exact approach to the analysis of spatial point patterns was proposed in Adams et al. (2009). They considered a univariate GP to describe the IF dynamics and an augmented model for the data and latent variables that simplifies the likelihood function. Theoretical concern about the algorithm of Adams et al. (2009) and empirical evidence to support it are provided in Section 5 of the on-line supplementary material. Another non-parametric exact approach was adopted in Kottas and Sansó (2007), where a particular factorization of the IF was proposed and Dirichlet processes priors were used. Their work was extended to the spatiotemporal context by Xiao et al. (2015).

The aim of this work is to propose an exact inference methodology for spatiotemporal Cox processes in which the IF dynamics are driven by a GP. The exactness feature stems from an augmented model approach as in Adams et al. (2009). However, we generalize their point pattern models by firstly considering spatiotemporal models and, secondly, by using multivariate (possibly dynamic) GPs to allow the inclusion of different model components (regression and temporal effects) in a flexible manner. Space and time may be considered continuous or discrete. In this paper we ultimately consider the general formulation of continuous space and discrete time. This is actually the most general formulation, since its continuous time version can be seen as a continuous space process where time is one of the dimensions.

Our methodology also introduces a particularly suited MCMC algorithm that enables direct simulation from the full conditional distributions of the GP and of other relevant latent variables. A particular choice of dominating measure to obtain the likelihood function is crucial to derive these sampling steps. Moreover, estimation of (possibly intractable) functionals of the IF and prediction based on the output of the MCMC algorithm are straightforward. We also provide formal proofs of the validity of the MCMC algorithm.

We believe that our methodology is the first valid exact alternative for inference in GPdriven Cox processes. It is robust under model complexity and specification, i.e. it is valid and has the same general formulation for any model where the prior measure on the main
component of the IF is a GP-univariate or multivariate, space or space-time, with or without covariates, etc. Furthermore, its complexity is mitigated by the good convergence properties of the MCMC algorithm proposed and computational time optimization strategies based on matrix algebra results. Finally, we present two real examples where we compare our methodology with a standard approximated methodology and demonstrates its applicability to reasonably sized problems.

The paper is organized as follows. Section 2 presents the class of models to be considered and the augmented model to derive the MCMC algorithm. Section 3 presents the general Bayesian approach and addresses some identifiability and implementation issues. Section 4 describes the MCMC algorithm for the spatial model and Section 5 presents its extension to the spatiotemporal case. Section 7 presents two real data examples to illustrate the methodology. Final remarks and possible directions for future work are presented in Section 7. The on-line supplementary material presents the aforementioned analysis of the algorithm of Adams et al. (2009) as well as a collection of simulated data analysis, proofs for the main results and other practical aspects of our methodology.

\section*{2. Model specification}

In this section we present the complete probabilistic model for spatiotemporal point processes with GP-driven intensities. We break the presentation in parts considering the different levels and generalizations of the model.

\subsection*{2.1. The general Cox process model}

We consider a Poisson process \(Y=\left\{Y_{t} ; t \in \mathcal{T}\right\}\) in \(S \times \mathcal{T}\), where \(S\) is some compact region in \(\mathbb{R}^{d}\) and \(\mathcal{T}\) is a finite set of \(\mathbb{N}\). It can be seen as a Poisson process in a region \(S\) that evolves in time. We assume that the Poisson process has an intensity function \(\lambda_{t}(s): S \times \mathcal{T} \rightarrow \mathbb{R}^{+}\).

We assume that the IF is a function of a (multivariate spatiotemporal) GP and covariates. A GP \(\beta\) is a stochastic process in some space such that the joint distribution of any finite collection of points in this space is Gaussian. This space may be defined so that we have spatial or spatiotemporal processes.

A detailed presentation of (dynamic) GPs is given in Section 2.3. For now, let \(\beta:=\left\{\beta_{0}, \beta_{1}, \ldots\right.\), \(\left.\beta_{q}\right\}\) be a collection of \(q+1\) independent GPs (a multivariate GP) in \(S \times \mathcal{T}\), for \(\mathcal{T}=\{0, \ldots, T\}\). We assume the following model for \(Y\) :
\[

\]
where \(\lambda_{\mathcal{T}}^{*}=\left(\lambda_{0}^{*}, \ldots, \lambda_{T}^{*}\right), \Phi\) is the distribution function of the standard Gaussian distribution, \(\mathrm{GP}_{\theta}\) is a GP indexed by (unknown) parameters \(\theta\) and \(W=\left\{W_{t}(s)\right\}\) is a set of covariates. We assume that \(f\) is linear in the co-ordinates of \(\beta\) and write \(W\) in a way such that \(f\left\{\beta_{t}(s), W_{t}(s)\right\}= W_{t}(s) \beta_{t}(s)\). Any distribution function of a continuous random variate or any general bounded function may be used instead of \(\Phi\). For example, a common choice is the logistic function, which was used by Adams et al. (2009) and is very similar to \(\Phi\) (the largest difference is 0.0095 ). The particular choice in expression (2) contributes to the construction of the MCMC algorithm
as discussed in Section 4. These links from the GP to the IF are bounded to \([0,1]\) in contrast with the unbounded exponential link of the log-Gaussian Cox processes that was proposed by Møller et al. (1998). We are interested in estimating not only the overall rate \(\lambda_{S, \mathcal{T}}\) but also the univariate GPs, separately, given their possibly meaningful interpretation.

Parameter \(\lambda_{t}^{*}\) ought to represent the supremum of the IF at time \(t\). One extreme possibility is to assume \(\lambda_{t}^{*}=\lambda^{*}\), which is a reasonable assumption in the case processes whose maximum intensity is time invariant. In this case, a common choice for the prior distribution of \(\lambda^{*}\) is \(\mathcal{G}\left(\alpha_{\lambda}, \beta_{\lambda}\right)\)-a gamma distribution. At the other extreme, unrelated parameters vary independently over time according to independent \(\mathcal{G}\left(\alpha_{\lambda_{t}}, \beta_{\lambda_{t}}\right)\) prior distributions. In between them, models allow for temporal dependence between (successive) \(\lambda_{t}^{*}\) s. One such formulation with attractive features is briefly described in Section 5.

The GPs may represent some relevant model features such as the effect of covariates and/or spatiotemporal components such as the spatially varying trend components or seasonality of the baseline intensity. They may also be space and/or time invariant. One common example is the model with \(q\) covariates:
\[
f\left\{\beta_{t}(s), W_{t}(s)\right\}=\beta_{0, t}(s)+\beta_{1, t}(s) W_{1, t}(s)+\ldots+\beta_{q, t}(s) W_{q, t}(s)
\]

This approach enables the use of extra prior information through covariates. The spatiotemporal variation of the effects is particularly relevant in applications where covariates present significant interaction with space and/or time. Examples are provided in Pinto et al. (2015).

The first advantage of the formulation in expressions (1)-(4) is that it allows exact simulation of data from the model, which is the key to develop exact inference methods. Exact simulation of the model is based on a key result from Poisson processes called Poisson thinning. This is a variant of rejection sampling for point processes proposed by Lewis and Shedler (1979) and is given in the following algorithm 1.

Step 1: simulate a Poisson process \(\left(s_{1}, \ldots, s_{K}\right)\) with constant intensity function \(\lambda_{t}^{*}\) on \(S\) -
(a) simulate \(K_{t} \sim \operatorname{Poisson}\left\{\lambda_{t}^{*} \mu(S)\right\}\), where \(\mu(S)\) is the volume of \(S\);
(b) distribute the \(K_{t}\) points uniformly on \(S\).

Step 2: simulate \(\beta_{t}\) and observe \(W_{t}\) at points \(\left\{s_{1}, \ldots, s_{K}\right\}\).
Step 3: keep each of the \(K_{t}\) points with probability \(\lambda_{t}\left(s_{k}\right) / \lambda_{t}^{*}\).
Step 4: output the points kept in the previous step.
To simulate the process at additional times \(t \in \mathcal{T}\), it is enough to perform algorithm 1 for each \(t\) and to simulate the GP conditionally on the points previously simulated. The idea of Poisson thinning has been applied in related contexts by Gonçalves and Roberts (2014).

\subsection*{2.2. The augmented model}

Performing exact inference for Cox processes is a challenging problem mainly because of the intractability of the likelihood function, which is given by
\[
L\left(\lambda_{S, \mathcal{T}}, y\right)=\exp \left\{-\sum_{t=0}^{T} \int_{S} \lambda_{t}(s) \mathrm{d} s\right\} \prod_{t=1}^{T} \prod_{n=1}^{N_{t}} \lambda_{t}\left(s_{t, n}\right)
\]
where \(s_{n, t}\) is the location of the \(n\)th event of \(Y_{t}\).
The crucial step to develop exact methods is to avoid dealing with likelihood (6). One possible solution is to define an augmented model for \(Y\) and some additional variable \(X\), such that the joint (pseudo)likelihood based on \((X, Y)\) is tractable. This poses the problem as a missing data
problem (where inference is based on the joint (pseudo)likelihood of the data and some latent variable) and enables us to use standard methods. The augmented model is constructed on the basis of the Poisson thinning that was presented in algorithm 1. The same strategy was employed in Adams et al. (2009) to circumvent the non-tractability problem.

Firstly, define \(X=\left\{X_{t} ; t \in \mathcal{T}\right\}\) where each \(X_{t}\) is a homogeneous Poisson process with intensity \(\lambda_{t}^{*}\) on \(S\) and the \(X_{t}\) s are mutually independent. Now let \(\left\{s_{t, k}\right\}_{k=1}^{K_{t}}\) be the locations of the \(K_{t}\) events of \(X_{t}\). We also define \(T\) vectors \(Z_{t}, t=1, \ldots, T\), with each co-ordinate taking values in \(\{0,1\}\) such that \(\left(Z_{t} \mid X, \beta_{K_{t}}, W_{K_{t}}\right)\) is a random vector \(\left(Z_{t, 1}, \ldots, Z_{t, K_{t}}\right)\), where the \(Z_{t, k} \mathrm{~s}\) are all independent with \(Z_{t, k} \sim \operatorname{Ber}\left[\Phi\left\{W_{t}\left(s_{t, k}\right) \beta_{t}\left(s_{t, k}\right)\right\}\right]\) and \(\left(\beta_{K_{t}}, W_{K_{t}}\right)\) is \((\beta, W)\) at the points from \(X_{t}\). Finally, define \(Y_{t}=h\left(Z_{t}, X_{t}\right)\) as the non-zero co-ordinates of the vector \(\left(Z_{t, 1} s_{t, 1}, \ldots, Z_{t, K_{t}} s_{t, K_{t}}\right)\), which leads to model (1).

Namely, the augmented model defines \(Y\) as the events remaining from performing the Poisson thinning to a Poisson process \(X\). It is important to note, however, that only \(Y\) is observed. We define \(\left\{s_{t, n}\right\}_{n=1}^{N_{t}}\) as the \(N_{t}\) events of \(Y_{t}\) and \(\left\{s_{t, m}\right\}_{m=1}^{K_{t}-N_{t}}\) as the \(M_{t}:=K_{t}-N_{t}\) thinned events. Most importantly, this approach leads to a tractable likelihood when the joint distribution of \(X\) and \(Y\) is considered, as is shown in Section 3.1.

The spatial model is a particular case where \(T=1\) which implies that \(X\) and \(Y\) are Poisson processes on \(S\) with IFs \(\lambda^{*}\) and \(\lambda(s)\) respectively. We observe \(\left\{s_{n}\right\}_{n=1}^{N}\) from \(Y\) and simplify the notation above accordingly. Note that the spatial model for unidimensional \(S\) is generally seen as the commonly used Cox process in (continuous) time.

\subsection*{2.3. Dynamic Gaussian processes}

GPs are a very flexible component to handle spatial variation, especially when smooth processes are expected. We say that \(\beta\) follows a stationary GP in \(S\) if \(\beta(s) \sim N\left(\mu, \sigma^{2}\right)\) and \(\operatorname{cov}\left\{\beta(s), \beta\left(s^{\prime}\right)\right\}= h\left(s, s^{\prime}\right)\), for \(s, s^{\prime} \in S\), constants \(\mu\) and \(\sigma^{2}\) and a differentiable (almost everywhere) function \(h\). Further simplification is obtained if isotropy can be assumed, leading to \(h\left(s, s^{\prime}\right)=\rho\left(\left|s-s^{\prime}\right|\right)\). In this case, the process is denoted by \(\beta \sim \operatorname{GP}\left(\mu, \sigma^{2}, \rho\right)\) and \(\rho\) is referred to as the correlation function.

Typical choices for \(h\) belong to the \(\gamma\)-exponential family of covariance functions:
\[
h\left(s, s^{\prime}\right)=\sigma^{2} \exp \left\{-1 /\left(2 \tau^{2}\right)\left|s-s^{\prime}\right|^{\gamma}\right\}, \quad 0<\gamma \leqslant 2 .
\]

The case where \(\gamma>1\) leads to almost surely differentiable paths (surfaces).
GPs can be extended in many directions. The most important here are extensions to handle multivariate GPs and extensions to cope with space and time. There are various ways to allow for multivariate responses. The main were reviewed in Gamerman et al. (2007) and include independent GPs, dependent processes with a common correlation function or linear mixtures of independent GPs, the co-regionalization models described by Wackernagel (2003). Similar ideas were also developed within the machine learning literature as multiple-output models (see, for example, Boyle and Frean (2004)). A review by Alvarez et al. (2012) discusses the relationship between co-regionalization and multiple-output models.

Extensions to cope with space and discrete time were considered by Gelfand et al. (2005) (see also Wikle and Cressie (1999)). A process \(\beta\) follows a dynamic GP in discrete time if it can be described by a difference equation
\[
\beta_{t^{\prime}}(\cdot)=G_{t^{\prime}, t} \beta_{t}(\cdot)+w_{t^{\prime}, t}(\cdot), \quad w_{t^{\prime}, t} \sim \mathrm{GP}
\]
where the multivariate GP disturbances \(w_{t^{\prime}, t}(\cdot)\) are zero mean and time independent; they are also taken as identically distributed in the equidistant case \(t^{\prime}=t+1\). The law of the process is
completed with a GP specification for \(\beta_{0}(\cdot)\). Similar processes were proposed in continuous time by Brix and Diggle (2001).

Several options are available for the temporal transition matrix \(G\), including the identity matrix. If additionally the disturbance processes \(w\) consist of independent GPs then the resulting process consists of independent univariate dynamic GPs. Gamerman (2010) presented some alternatives to model trend and seasonality of the IF.

It is also possible to consider the use of non-spatiotemporal covariates to explain the IF variation. This approach, however, requires adaptations to the original model and can be found in Pinto et al. (2015).

\section*{3. Inference for the spatial model}

We now focus on the inference problem of estimating the IF \(\lambda_{S, \mathcal{T}}\), parameter \(\lambda^{*}\) and potentially unknown parameters \(\theta\) from the GP, based on observations from the Poisson process \(Y\). We shall also discuss how to make prediction. To make the presentation of the methodology as clear as possible, we consider first the (purely) spatial process and then the generalization for the spatiotemporal case.

\subsection*{3.1. Posterior distribution}

Let \(\left\{s_{k}\right\}_{k=1}^{K}\) and \(\left\{s_{n}\right\}_{n=1}^{N}\) be the events from \(X\) and \(Y\) respectively, and \(\left\{s_{m}\right\}_{m=1}^{K-N}\) be the thinned events. Furthermore, \(\beta_{K}, \beta_{N}\) and \(\beta_{M}\) are \(\beta\) at \(\left\{s_{k}\right\}_{k=1}^{K},\left\{s_{n}\right\}_{n=1}^{N}\) and \(\left\{s_{m}\right\}_{m=1}^{K-N}\) respectively, and \(W_{K}, W_{N}\) and \(W_{M}\) are the respective subsets of \(W\). Note that \(\left\{s_{k}\right\}_{k=1}^{K}=\left(\left\{s_{n}\right\}_{n=1}^{N} \cup\left\{s_{m}\right\}_{m=1}^{K-N}\right)\) and \(\beta_{K}=\left(\beta_{N}, \beta_{M}\right)\). We express the vector of all random components of the model by \(\left(\left\{s_{n}\right\}_{n=1}\right.\), \(\left.\left\{s_{m}\right\}_{m=1}^{K-N}, \beta,\left\{s_{k}\right\}_{k=1}^{K}, K, \lambda^{*}, \theta, W\right)\).

Initially, assume that \(W\) is deterministic and consider only \(\beta_{K}\) instead of \(\beta\). Thus, define \(\psi:=\left(\left\{s_{n}\right\}_{n=1}^{N},\left\{s_{m}\right\}_{m=1}^{K-N}, \beta_{K},\left\{s_{k}\right\}_{k=1}^{K}, K, \lambda^{*}, \theta\right)\). This strategy simplifies the problem (making it finite dimensional) but still makes it possible to estimate the infinite dimensional remainder of \(\beta\)-to be discussed further ahead. Note that, because of redundancy issues, the components of the model could be specified in other ways.

We now specify the joint distribution of all the random components of the model- \((\psi \mid W, S)\). Note that the joint posterior density that we are aiming for is proportional to this. It is important to make an appropriate choice of a dominating measure with respect to which we write the density of \(\psi\). This choice depends on the specification of the components and is not unique. For example, one may choose to write the density of \(\left(K,\left\{s_{k}\right\}_{k=1}^{K} \mid \lambda^{*}\right)\) with respect to the measure of a unit rate Poisson process, resulting in
\[
\exp \left\{-\left(\lambda^{*}-1\right) \mu(S)\right\}\left(\lambda^{*}\right)^{K} \propto \exp \left\{-\lambda^{*} \mu(S)\right\}\left(\lambda^{*}\right)^{K}
\]
which is the usual Poisson process likelihood function. However, this is a valid likelihood for \(\lambda^{*}\) only, since the dominating measure is independent of this parameter. It is not a valid likelihood function for \(K\), which, in our case, is unknown and must be estimated. This gives good intuition about why this is not a good choice for the dominating measure. Although it is one possibility, it makes the derivation of the full conditional distributions (or the acceptance probability of potential Metropolis-Hastings steps) more difficult.

We choose to write the density of \((\psi \mid W, S)\) with respect to the dominating measure given by the product measure \(\mathbb{Q}:=\delta^{K} \otimes \mathbb{L}^{K} \otimes \mathbb{L}^{K} \otimes \delta \otimes \mathbb{L} \otimes \mathbb{L}^{d_{\theta}}\), where \(\delta\) is the counting measure, \(\mathbb{L}^{d}\) is the \(d\)-dimensional Lebesgue measure and \(d_{\theta}\) is the dimension of the parameter vector \(\theta\). This choice is related to the factorization that we choose. If we let \(\mathbb{P}\) be the probability measure of
our full model, the density \(\pi\) of \((\psi \mid W, S)\), which is defined as the Radon-Nikodym derivative of \(\mathbb{P}\) with respect to \(\mathbb{Q}\), is given by
\[
\begin{aligned}
\pi(\psi \mid W, S) & =\pi\left(\left\{s_{n}\right\},\left\{s_{m}\right\} \mid\left\{s_{k}\right\}, \beta_{K}, W\right) \pi\left(\beta_{K} \mid\left\{s_{k}\right\}, \theta\right) \pi\left(\left\{s_{k}\right\} \mid K, S\right) \pi\left(K \mid \lambda^{*}, S\right) \pi\left(\lambda^{*}\right) \pi(\theta) \\
& =\Phi_{N}\left(W_{N} \beta_{N} ; I_{N}\right) \Phi_{K-N}\left(-W_{M} \beta_{M} ; I_{M}\right) \pi_{\mathrm{GP}}\left(\beta_{K} \mid \theta\right) \exp \left\{-\lambda^{*} \mu(S)\right\}\left(\lambda^{*}\right)^{K} \frac{1}{K!} \pi\left(\lambda^{*}\right) \pi(\theta)
\end{aligned}
\]
where \(\left\{s_{n}\right\}=\left\{s_{n}\right\}_{n=1}^{N},\left\{s_{m}\right\}=\left\{s_{m}\right\}_{m=1}^{K-N}\) and \(\left\{s_{k}\right\}=\left\{s_{k}\right\}_{k=1}^{K} . \Phi_{k}\left(\cdot ; I_{k}\right)\) is the distribution function of the \(k\)-dimensional Gaussian distribution with mean vector 0 and covariance matrix \(I_{k}\left(k\right.\)-dimensional identity matrix) and \(\beta_{N}=\left(\beta_{0}\left(s_{1}\right) \ldots \beta_{0}\left(s_{N}\right) \ldots \beta_{q}\left(s_{1}\right) \ldots \beta_{q}\left(s_{N}\right)\right)^{\prime}\). Also, \(W_{N}=\) ( \(I_{N} W_{1} \ldots W_{q}\) ), where \(W_{i}\) is an \(N \times N\) diagonal matrix with the ( \(n, n\) )-entry \(W_{i}\left(s_{n}\right)\)-the \(i\) th covariate at location \(s_{n}\). Furthermore, \(\pi_{\mathrm{GP}}\left(\beta_{K} \mid \theta\right)\) is the density of the multivariate GP at locations \(\left\{s_{k}\right\}_{k=1}^{K}\) with respect to \(\mathbb{L}^{K}\). Finally, \(\pi\left(\lambda^{*}\right)\) and \(\pi(\theta)\) are the prior densities of \(\lambda^{*}\) and \(\theta\) respectively, with respect to \(\mathbb{L}\) and \(\mathbb{L}^{d_{\theta}}\).

\subsection*{3.2. Estimation of the intensity function}

The MCMC algorithm to be proposed in Section 4.2 outputs samples from the posterior distribution of the IF \(\lambda_{S}\) at the observed locations \(\left\{s_{n}\right\}_{n=1}^{N}\) and at another finite collection of locations \(\left\{s_{m}\right\}_{m=1}^{K-N}\) which varies between the iterations of the algorithm. Nevertheless, we need to have posterior estimates of \(\lambda_{S}\) over the whole space \(S\). It is quite simple, though possibly costly, to sample exactly from this posterior (at any finite collection of locations). It can be done by adding an extra step to the MCMC algorithm or by a sampling procedure after the MCMC runs. Both schemes may suffer from high computational cost but the former is considerably cheaper if well designed. Details are provided in Section 4.3. Finally, efficient solutions based on lower dimension or a nearest neighbour approximation may be employed when \(K\) is too large-see Section 4 of the on-line supplementary material.

\subsection*{3.3. Model identifiability and practical implementation}

The model proposed may suffer from identifiability problems concerning parameter \(\lambda^{*}\). The natural way to identify it is to have this parameter as the supremum of the IF which, under the Bayesian approach, should be achieved by an appropriate specification of the prior distribution. Any prior that identifies this parameter as equal to or greater than the supremum solves the identifiability problem. But the former makes the parameter interpretable and optimizes the computational cost (stochastically minimizes \(K\) ).

A reasonable choice for the prior of \(\lambda^{*}\) is a (truncated) gamma distribution for which the hyperparameter could be specified through an empirical analysis of the data set, more specifically, by obtaining an empirical estimate of the intensity in a small area with the highest concentration of points. This area should be reasonably chosen to give a good idea of the supremum of the IF. Note that the data are being used only to identify the model, which is different from using the data twice in a model which is already identified.

The prior distribution of \(\beta\) may also help in the identification of \(\lambda^{*}\). Note that, if \(\beta\) is estimated to be high (greater than 2) at any location, it implies that \(\lambda^{*}\) is (practically) identified as the supremum of the IF. In this sense, the prior on \(\beta\) may be specified in a way to favour such a scenario by, for example, fixing a positive mean parameter and/or a variance parameter that is coherent with the standard Gaussian distribution.

Generally speaking, identifiability is an important issue when estimating the IF of a non-
homogeneous Poisson process. It is well known that the reliability of the estimates relies on the amount of data that are available. In this sense, the higher is the actual IF the better. In a Bayesian framework, in particular, the prior on the IF plays an important role in the identification and estimation of this function. This is related to the fact that the data do not contain much information about the hyperparameters of GP priors. The simulated examples in section 1 of the on-line supplementary material explore this issue and provide some insight on how to proceed in general.

Another important issue is the computational cost from dealing with GPs. Despite their great flexibility in a variety of statistical modelling problems, GPs have a considerable practical limitation when it comes to computational cost. More specifically, simulating an \(n\)-dimensional GP has a cost which is typically of the order of \(n^{3}\). This means that in our case the cost would be \(O\left(K^{3}\right)\), without involving the procedures in Section 3.2.

Nevertheless, this issue is mitigated by several reasons. Firstly, our MCMC algorithm has very good properties in terms of speed of convergence and auto-correlation (see Fig. 3 of the on-line supplementary material) which, in turn, implies that not too many iterations are required to obtain good results. Secondly, the computational cost is feasible for reasonably sized problems because of
(a) matrix algebra strategies to avoid the computation of inverse matrices,
(b) the fast convergence of the embedded Gibbs sampler to sample from a truncated normal distribution and
(c) computational strategies to estimate the IF in a fine grid
(all three points are explained in detail in Section 4 of the supplementary material).
Alternative approaches based on discrete approximations are bound to suffer from similar dimensionality issues. Note, however, that the relevant order of magnitude there is defined by the number of grid points, which may be larger than the number of events. Also, non-trivial tuning of the Metropolis-Hasting algorithm (see, for example, Møller et al. (1998)) is crucial for devising an efficient MCMC algorithm to converge in feasible time. This is in sharp contrast with our algorithm that samples directly from the full conditional distribution of the GP and presents good convergence properties.

Furthermore, for the cases where the cost is still too high, some (approximating) strategies may be employed to reduce it. Most importantly, none of these defy the exactness of our methodology. Lower dimension approximations (see Banerjee et al. (2013) and Carlin et al. (2007)) may be used to deal with the GP prior and to speed the computation of covariance matrices, their inverses and their Cholesky decompositions-all required for the MCMC algorithm. This type of approximation is on the second level of the model and the first level (likelihood) is still dealt with in an exact set-up. Simpson et al. (2016) (section 3) argued that approximations in the first level have much more effect on the results. Furthermore, from a different perspective, which agrees with the argument from Simpson et al. (2016), the use of our proposed methodology combined with approximations for the GP may be seen as a fully exact set-up where the prior on \(\beta\) is given by the probability measure defined by the GP approximation. Therefore, as long as the approximation defines a Gaussian probability measure, this may be seen as an alternative model for which exact inference is carried out under our (exact) methodology. Another strategy that may reduce the computational cost is to use a carefully well-designed Metropolis-type step to sample the GP. Two possible solutions are Hamiltonian Monte Carlo sampling (see Adams et al. (2009)) and the Metropolis-adjusted Langevin algorithm (see Møller et al. (1998)). These will avoid the embedded Gibbs sampling that is required to sample directly from the full conditional distribution of the GP.

\section*{4. Computation for the spatial model}

In this section, we present the computational details to perform inference in the spatial model. The methodology proposed consists of an MCMC algorithm which has the exact joint posterior distribution of the unknown components of the model as its invariant distribution. More specifically, the algorithm is a Gibbs sampler. The derivation of the full conditional distributions is not straightforward for several reasons: intractability issues; the redundancy among some of the components; the hierarchical structure of the model, especially the fact that the observations are not (explicitly) on the first level, because of thinning. To sample directly from the full conditional distributions, it is essential to be able to simulate from a general class of multivariate skew normal distributions. We define such a class and propose an algorithm to sample from it.

\subsection*{4.1. A general class of multivariate skew normal distributions}

We consider a general class of skew normal distributions that was originally proposed in Arellano-Valle and Azzalini (2006) and present it here in a particularly useful way for the context of our work. Equally important, we also propose an algorithm to sample from this distribution.

For a \(d\)-dimensional vector \(\xi\), an \(m \times d\) matrix \(W\) and a \(d \times d\) matrix \(\Sigma\), we define
\[
\begin{gathered}
U=\binom{U_{0}}{U_{1}} \sim N_{m+d}\left(0, \Sigma^{*}\right) \\
\Sigma^{*}=\left(\begin{array}{cc}
\Gamma & \Delta^{\prime} \\
\Delta & \Sigma
\end{array}\right)
\end{gathered}
\]
where \(\Gamma=I_{m}+W \Sigma W^{\prime}\) and \(\Delta^{\prime}=W \Sigma\). Let \(a=\left(a_{1}, \ldots, a_{r}\right)>b=\left(b_{1}, \ldots, b_{r}\right)\) mean that \(a_{i}>b_{i}, \forall i\), and define \(\gamma=W \xi\). We say that ( \(U_{1}+\xi \mid U_{0}>-\gamma\) ) has an \(\operatorname{SN}(\xi, \Sigma, W)\) distribution whose density is given in the following proposition.

Proposition 1. The density of ( \(U_{1}+\xi \mid U_{0}>-\gamma\) ) is given by
\[
f(z)=\frac{1}{\Phi_{m}(\gamma ; \Gamma)} \phi_{d}(z-\xi ; \Sigma) \Phi_{m}\left(W z ; I_{m}\right)
\]
where \(\phi_{d}(\cdot ; \Omega)\) and \(\Phi_{d}(\cdot ; \Omega)\) are the density and distribution function respectively of the \(d\) dimensional Gaussian distribution with mean vector 0 and covariance matrix \(\Omega\).

For a proof of proposition 1, see the on-line supplementary material-section 2.
We propose the following algorithm to sample from the density in equation (12). Define \(U_{0}^{*}= A^{-1} U_{0}\), where \(A\) is the lower triangular matrix that is obtained from the Cholesky decomposition of \(\Gamma\), i.e. \(\Gamma=A A^{\prime}\). This implies that \(U_{0}^{*} \sim N_{m}\left(0, I_{m}\right)\) and \(U_{0}=A U_{0}^{*}\). We use the following results to construct our algorithm:
\[
f\left(U_{1}, U_{0} \mid U_{0}>-\gamma\right)=f\left(U_{1} \mid U_{0}, U_{0}>-\gamma\right) f\left(U_{0} \mid U_{0}>-\gamma\right) .
\]

Proposition 2. ( \(A U_{0}^{*} \mid U_{0}^{*} \in B\) ) has the same distribution as ( \(U_{0} \mid U_{0}>-\gamma\) ), where \(B=\left\{u_{0}^{*}:\right. \left.A u_{0}^{*}>-\gamma\right\}\).

For a proof of proposition 2, see the on-line supplementary material-section 2.
The decomposition in equation (13) suggests that simulation from density (12) may be performed by firstly simulating ( \(U_{0} \mid U_{0}>-\gamma\) ) and then using this value to simulate from ( \(U_{1} \mid U_{0}\) ). Moreover, the simulation of \(U_{0}\) is more efficient (as described in section 3.1 of the on-line supplementary material) if we first simulate \(U_{0}^{*}\) and then apply the appropriate transformation,
as suggested by proposition 2. The algorithm to simulate from density (12) is as follows (algorithm 2).

Step 1: simulate a value \(u^{*}\) from \(\left(U_{0}^{*} \mid U_{0}^{*} \in B\right)\).
Step 2: obtain \(u=A u^{*}\).
Step 3: simulate a value \(z^{*}\) from \(\left(U_{1} \mid U_{0}=u\right) \sim \mathcal{N}\left(\Delta \Gamma^{-1} u, \Sigma-\Delta \Gamma^{-1} \Delta^{\prime}\right)\).
Step 4: obtain \(z=z^{*}+\xi\).
Step 5: output \(z\).
The simulation of step 3 is trivial. Step 1 consists of the simulation of a truncated (by linear constraints) multivariate normal distribution and cannot be performed directly. The simulation from this distribution is described in the on-line supplementary material-section 3.1.

\subsection*{4.2. The Gibbs sampling algorithm}

Note that, given the data \(\left\{s_{k}\right\}\), the remaining unknown quantities are ( \(\left\{s_{m}\right\}, \beta_{K}, K, \lambda^{*}, \theta\) ). We block these quantities as \(\left(\left\{s_{m}\right\}, \beta_{M}, K\right), \beta_{K}, \lambda^{*}\) and \(\theta\). Note that \(\beta_{M}\) is sampled twice. That is mainly because updating \(B_{K}\) instead of only \(\beta_{N}\) significantly improves the mixing of the chain (by reducing the correlation between blocks) and because sampling the first block without \(\beta_{M}\) is virtually impossible. All full conditional densities are proportional to \(\pi(\psi \mid W, S)\). Thus, elimination of constant terms leads to
\[
\begin{gathered}
\pi\left(\left\{s_{m}\right\}, \beta_{M}, K \mid \cdot\right) \propto \Phi_{K-N}\left(-W_{M} \beta_{M} ; I_{K-N}\right) \pi_{\mathrm{GP}}\left(\beta_{M} \mid \beta_{N}, \theta\right)\left|\left(\lambda^{*}\right)^{K} \frac{1}{K!}\right| 1(K \geqslant N) \\
\pi\left(\beta_{K} \mid \cdot\right) \propto \Phi_{N}\left(W_{N} \beta_{N} ; I_{N}\right) \Phi_{K-N}\left(-W_{M} \beta_{M} ; I_{K-N}\right) \pi_{\mathrm{GP}}\left(\beta_{K} \mid \theta\right) \\
\pi\left(\lambda^{*} \mid \cdot\right) \propto \exp \left\{-\lambda^{*} \mu(S)\right\}\left(\lambda^{*}\right)^{K} \pi\left(\lambda^{*}\right) \\
\pi(\theta \mid \cdot) \propto \pi_{\mathrm{GP}}\left(\beta_{K} \mid \theta\right) \pi(\theta)
\end{gathered}
\]

The four densities above are written with respect to dominating measures: \(\mathbb{L}^{K-N} \otimes \mathbb{L}^{K-N} \otimes \delta, \mathbb{L}^{K}, \mathbb{L}\) and \(\mathbb{L}^{d_{\theta}}\), in accordance with the dominating measure that is used to write equation (10).

Define \(\pi_{0}\) as a Poisson \(\left\{\lambda^{*} \mu(S)\right\}\) distribution truncated to \(\{N, N+1, \ldots\}\). We propose the following rejection sampling algorithm to sample from distribution (14) (algorithm 3).

Step 1: simulate \(\dot{K}\) from \(\pi_{0}\).
Step 2: if \(\dot{K}=N\), make \(\left\{\dot{s}_{m}\right\}_{m=1}^{\dot{K}-N}=\dot{\beta}_{M}=\emptyset\) and go to step 8; otherwise go to step 3.
Step 3: make \(m=1\) and \(\beta_{1: m-1}=\emptyset\).
Step 4: make \(r_{m}=1\).
Step 5: simulate \(\dot{s}_{r_{m}} \sim \mathcal{U}(S)\) and \(\dot{\beta}_{r_{m}}\left(\dot{s}_{r_{m}}\right)\) from \(\pi_{\mathrm{GP}}\left\{\dot{\beta}_{r_{m}}\left(\dot{s}_{r_{m}}\right) \mid \beta_{N}, \dot{\beta}_{1: m-1}, \theta\right\}\).
Step 6: simulate \(Z_{r_{m}} \sim \operatorname{Ber}\left[\Phi\left\{-W\left(\dot{s}_{r_{m}}\right) \beta\left(\dot{s}_{r_{m}}\right)\right\}\right]\).
Step 7:
(a) if \(Z_{r_{m}}=1\) and \(m<K-N\), set \(\dot{s}_{m}=\dot{s}_{r_{m}}, \dot{\beta}\left(\dot{s}_{m}\right)=\dot{\beta}_{r_{m}}\left(\dot{s}_{r_{m}}\right), \dot{\beta}_{1: m-1}=\dot{\beta}_{1: m-1} \cup \dot{\beta}_{r_{m}}\left(\dot{s}_{r_{m}}\right)\) and \(m=m+1\) and go to step 4;
(b) if \(Z_{r_{m}}=1\) and \(m=K-N\), set \(\dot{s}_{m}=\dot{s}_{r_{m}}\) and \(\dot{\beta}\left(\dot{s}_{m}\right)=\dot{\beta}_{r_{m}}\left(\dot{s}_{r_{m}}\right)\) and go to step 8;
(c) if \(Z_{r_{m}}=0\), set \(r_{m}=r_{m}+1\) and go to step 5 .

Step 8: output \(\left(\dot{K},\left\{\dot{s}_{m}\right\}_{m=1}^{\dot{K}-N}, \dot{\beta}_{M}\right)\).
Lemma 1. The output of algorithm 3 is an exact draw from the full conditional distribution (14).

For a proof of lemma 1, see the on-line supplementary material-section 2.
Note that algorithm 3 takes advantage of the factorization of the global acceptance probability to perform the accept-reject procedure pointwise and to avoid a much higher cost. The straightforward version of this algorithm would propose and accept-reject the variables all at once, resulting in a possibly very small acceptance probability. Firstly, \(K\) is sampled from \(\pi_{0}\) then, for each of the \(K-N\) locations, a pair ( \(s, \beta(s)\) ) is proposed from a \(\mathcal{U}(S)\) and the (conditional) GP prior and accepted with probability \(\Phi\{-W(s) \beta(s)\}\). Metropolis-Hastings alternatives may sound like an attractive possibility because of the lower computational cost but the usual choices for the proposal distribution may lead to slower convergence. This option performed poorly even for simple examples in some simulated studies considering both dependent and independent proposals.

The choice of the Gaussian cumulative distribution function and the linearity in \(\beta\) in expression (2) is justified by the fact that it makes it possible to sample directly from the full conditional distribution in expression (15). This leads to an algorithm with a reasonable computational cost and good convergence properties. The algorithm is as follows (algorithm 4).

Step 1: obtain \(W_{K}\) from \(\left(W_{N}, W_{M}\right)\) such that
\[
\Phi_{K}\left(W_{K} \beta_{K} ; I_{K}\right)=\Phi_{N}\left(W_{N} \beta_{N} ; I_{N}\right) \Phi_{K-N}\left(-W_{M} \beta_{M} ; I_{K-N}\right) .
\]

Step 2: sample \(\beta_{K} \sim \mathrm{SN}\left(\mu_{K}, \Sigma_{K}, W_{K}\right)\) by using algorithm 2, where \(\mu_{K}\) and \(\Sigma_{K}\) are the mean vector and covariance matrix respectively of \(\pi_{\mathrm{GP}}\left(\beta_{K} \mid \theta\right)\).
Step 3: output \(\beta_{K}\).
Lemma 2. The output of algorithm 4 is an exact draw from the full conditional distribution in expression (15).

Proof. Simply note that distribution (15) is proportional to the density of an \(\operatorname{SN}\left(\mu_{K}, \Sigma_{K}, W_{K}\right)\) distribution.

The blocking and sampling schemes of our Gibbs sampling result in good convergence properties. In fact, the examples that are presented here suggest that convergence is attained after a few iterations. Algorithm 3 may suggest a high computational cost as every try of the rejection sampling algorithm requires the simulation of the GP at a location given the existent locations. However, the most expensive part of this simulation is the computation of the inverse covariance matrix of the existing points which can be made considerably fast by using strategies based on matrix algebra. The same ideas may be used throughout the algorithm to reduce the computational time considerably-all those strategies are described in detail in section 4 of the on-line supplementary material. Algorithm 4 also has a reasonable cost for the reasons that are explained in that section.

The next step of the Gibbs sampler draws \(\lambda^{*}\) from its full conditional distribution. This can be obtained by routine calculations: if a conjugate gamma prior \(\mathcal{G}\left(\alpha_{\lambda}, \beta_{\lambda}\right)\) is adopted for \(\lambda^{*}\), its full conditional is \(\mathcal{G}\left\{\alpha_{\lambda}+K, \beta_{\lambda}+\mu(S)\right\}\).

The fourth and last step from the Gibbs sampler draws \(\theta\) from its full conditional distribution. This task may be carried out ordinarily-using direct simulation when possible or via an appropriately tuned Metropolis-Hastings step. There is also the option of breaking \(\theta\) into smaller blocks if that is convenient for computational reasons. One attractive option is to use an adaptive Gaussian random-walk Metropolis-Hastings step where the covariance matrix of the proposal is based on the empirical covariance matrix of the previous steps, as proposed by Roberts and Rosenthal (2009).

\subsection*{4.3. Estimating functionals of the intensity function}

One of the purposes of fitting a Cox process to an observed point pattern is to estimate functionals of the IF. These functionals may include the intensity itself, the mean number of points at some subregion, etc. The estimation is performed by sampling such functionals to obtain Monte Carlo estimates. The sampling step is performed on the basis of the following lemma.

Lemma 3. Let \(S_{0}=\left\{\tilde{s}_{1}, \ldots, \tilde{s}_{G}\right\}\) be a finite set of locations in \(S\). Given ( \(\beta_{K}, \theta\) ), \(\beta_{S_{0}}\) is independent of \(\left(\left\{s_{n}\right\}_{n=1}^{N}, W\right)\) and its posterior distribution is given by
\[
\pi\left(\beta_{S_{0}} \mid\left\{s_{n}\right\}_{n=1}^{N}, W, S\right)=\int \pi\left(\beta_{S_{0}} \mid \beta_{K}, \theta\right) \pi\left(\beta_{K}, \theta \mid\left\{s_{n}\right\}_{n=1}^{N}, W, S\right) \mathrm{d} \beta_{K} \mathrm{~d} \theta
\]

For a proof of lemma 3, see the on-line supplementary material-section 2.
Lemma 3 implies that, to sample the GP at some arbitrary location from its posterior distribution, it is enough to sample from the GP prior conditionally on the GP sample from the MCMC algorithm at locations \(\left\{s_{k}\right\}_{k=1}^{K}\). The most efficient way to do this is by adding a sampling step to each iteration of the MCMC algorithm.

To obtain estimates of the IF in the finite subset \(S_{0}\) (a fine squared grid, say) of \(S\), we need posterior samples of \(\left(\lambda^{*}, \beta_{S_{0}}\right)\). That is achieved by sampling \(\left(\lambda^{*}, \beta_{K}, \theta\right)\) from \(\pi\left(\lambda^{*}, \beta_{K}, \theta \mid\left\{s_{n}\right\}_{n=1}^{N}\right.\), \(W, S\) ) and then \(\beta_{S_{0}}\) from \(\pi\left(\beta_{S_{0}} \mid \beta_{K}, \theta\right)\) at each step of the Markov chain after convergence is assumed to hold. This way, at iteration \(j\) of the chain, a draw from the posterior of \(\lambda_{S_{0}}\) is given by \(\left\{\lambda^{*(j)} \Phi\left\{W(\tilde{s}) \beta^{(j)}(\tilde{s})\right\}, \tilde{s} \in S_{0}\right\}\).

Another interesting functional to be estimated is the integrated intensity \(\Lambda(R)=\int_{R} \lambda(s) \mathrm{d} s\) for some region \(R \subseteq S\). This is the mean number of points in \(R\). Monte Carlo estimates of \(E[\Lambda(R) \mid y]\) may be obtained without any discretization error by introducing a random variate \(U \sim \mathcal{U}(R)\) and noting that
\[
E_{U}[\lambda(U)]=\frac{1}{\mu(R)} \int_{R} \lambda(s) \mathrm{d} s
\]
(see Beskos et al. (2006)), thus suggesting the estimator
\[
\hat{\Lambda}(R)=\mu(R) \frac{1}{J} \sum_{j=1}^{J} \lambda^{(j)}\left(U^{(j)}\right)
\]
which is a strongly consistent estimator of \(E[I \mid y]\) by the strong law of large numbers for Markov chains. The samples of \(\lambda\) come from the posterior distribution. The accuracy of the estimator may be improved defining a partition of \(R\) and using one uniform distribution to approximate the integral from each subregion of the partition.

\section*{5. Spatiotemporal model}

It is straightforward to generalize the MCMC algorithm from Section 4.2 to the spatiotemporal case. Remember that \(\left(X_{0}, \ldots, X_{T}\right)\) are conditionally mutually independent homogeneous Poisson processes on \(S\), given \(\lambda_{\mathcal{T}}^{*}\). The temporal dependence of the model is defined through \(\beta\) (and, possibly, \(\lambda_{\mathcal{T}}^{*}\) ).

We now write the density of ( \(\psi \mid W, S\) ) with respect to the dominating measure given by the product measure of the counting measure and the Lebesgue measure with corresponding dimensions and obtain
\[
\pi(\psi \mid W, S)=\prod_{t=0}^{T}\left\{\Phi_{N_{t}}\left(W_{N_{t}} \beta_{N_{t}} ; I_{N_{t}}\right) \Phi_{K_{t}-N_{t}}\left(-W_{M_{t}} \beta_{M_{t}} ; I_{K_{t}-N_{t}}\right)\right\} \pi_{\mathrm{GP}}\left(\beta_{K_{\mathcal{T}}} \mid \theta\right)
\]
\[
\times \prod_{t=0}^{T}\left[\frac{\exp \left\{-\lambda_{t}^{*} \mu(S)\right\}\left(\lambda_{t}^{*}\right)^{K_{t}}}{K_{t}!} \pi\left(\lambda_{t}^{*}\right)\right] \pi(\theta)
\]
where the new notation has a natural interpretation and \(\pi_{\mathrm{GP}}\) is the density of the dynamic GP in expression (8).

We have at least two options for the blocking scheme. The first samples ( \(K_{t},\left\{s_{t, m}\right\}, \beta_{M_{t}}\) ) and \(\beta_{K_{t}}\) separately, for each time. This algorithm may, however, lead to a chain with poor mixing properties if \(T\) is large because of the temporal dependence of \(\beta\) (see Carter and Kohn (1994), Fruhwirth-Schnatter (1994) and Gamerman (1998)). This problem is mitigated by a blocking scheme that makes \(\left\{\left(K_{t},\left\{s_{t, m}\right\}, \beta_{M_{t}}\right)\right\}_{t=0}^{T}\) and \(\left\{\beta_{K_{t}}\right\}_{t=0}^{T}\) one block each. This choice eliminates the mixing problem that was mentioned above and is particularly appealing in the dynamic GP context. The first block is sampled by applying algorithm 3 to each time \(t\) and considering the GP temporal dependence in step 5. The same idea extends algorithm 4 (and step 4 of this). Defining \(\beta_{K_{(t-1)}}:=\left(\beta_{K_{0}}, \ldots, \beta_{K_{t-1}}\right)\) and \(\beta_{K_{-1}}=\emptyset\), the one at a time simulation is possible because of the factorization
\[
\pi\left(\left\{\beta_{K_{t}}\right\}_{t=0}^{T} \mid \cdot\right) \propto \prod_{t=0}^{T}\left\{\Phi_{N_{t}}\left(W_{N_{t}} \beta_{N_{t}} ; I_{N_{t}}\right) \Phi_{K_{t}-N_{t}}\left(-W_{M_{t}} \beta_{M_{t}} ; I_{K_{t}-N_{t}}\right) \pi_{\mathrm{GP}}\left(\beta_{K_{t}} \mid \beta_{K_{(t-1)}}, \theta\right)\right\} .
\]

The full conditional distribution of \(\theta\) is carried out as before and particular blocking schemes may be motivated by the spatiotemporal structure. Finally, for a prior \(\mathcal{G}\left(\alpha_{\lambda_{t}}, \beta_{\lambda_{t}}\right)\), the full conditional of each \(\lambda_{t}^{*}\) is \(\mathcal{G}\left\{\alpha_{\lambda_{t}}+K_{t}, \beta_{\lambda_{t}}+\mu(S)\right\}\). In the case \(\lambda_{t}^{*}=\lambda^{*}, \forall t\), the full conditional of this parameter is \(\mathcal{G}\left\{\alpha_{\lambda_{t}}+\Sigma_{t=1}^{T} K_{t}, \beta_{\lambda_{t}}+T \mu(S)\right\}\).

Extensions of the spatiotemporal model above can be proposed by adding a temporal dependence structure to \(\lambda_{0: T}^{*}\)-this is particularly useful for prediction. One interesting possibility is the Markov structure that was proposed by Gamerman et al. (2013) (in a state space model context) where \(\lambda_{0}^{*} \sim \mathcal{G}\left(a_{0}, b_{0}\right), \lambda_{t}^{*} \mid K_{1: t-1}, \lambda_{t-1}^{*}=w^{-1} \lambda_{t-1}^{*} \varsigma_{t}\) and \(\varsigma_{t} \sim \operatorname{beta}\left\{w a_{t},(1-w) a_{t}\right\}\). The full conditional distribution of \(\lambda_{0: T}^{*}\) is available in section 3.2 of the on-line supplementary material.

\subsection*{5.1. Prediction}

Suppose that we want to predict \(Y\) at future times \(\mathcal{T}^{*}=(T+1, \ldots, T+J)\). The algorithm to sample from the predictive distribution of \(Y_{\mathcal{T}^{*}}\), proceeds iteratively in time from \(T+1\) onwards. Firstly, we sample \(\lambda_{t}^{*}\) (which depends on the structure that has been adopted); then we apply algorithm 1 with the GP being simulated from \(\pi_{\mathrm{GP}}\left(\beta_{K_{t}} \mid \beta_{K_{\mathcal{T}}}, \beta_{K_{T+1: t-1}}, \theta\right)\). This algorithm is supported by the following result:
\[
\begin{aligned}
\pi\left(y_{\mathcal{T}^{*}} \mid y_{\mathcal{T}}\right) \propto & \int \pi\left(y_{\mathcal{T}^{*}}, \lambda_{\mathcal{T}^{*}}^{*}, \beta_{\mathcal{T}^{*}}, \lambda_{T}^{*}, \beta_{\mathcal{T}}, \theta \mid y_{\mathcal{T}}\right) \mathrm{d} \lambda_{\mathcal{T}^{*}}^{*} \mathrm{~d} \beta_{\mathcal{T}^{*}} \mathrm{~d} \lambda_{T}^{*} \mathrm{~d} \beta_{\mathcal{T}} \mathrm{d} \theta \\
= & \int \prod_{t=T+1}^{T+J}\left[\pi\left(y_{t} \mid \lambda_{t}^{*}, \beta_{t}\right) \pi\left(\lambda_{t}^{*} \mid \lambda_{t-1}^{*}, y_{t-1}\right) \pi\left(\beta_{t} \mid \beta_{t-1}, \theta\right)\right] \\
& \times \pi\left(\lambda_{T}^{*}, \beta_{\mathcal{T}}, \theta \mid y_{\mathcal{T}}\right) \mathrm{d} \lambda_{\mathcal{T}^{*}}^{*} \mathrm{~d} \beta_{\mathcal{T}^{*}} \mathrm{~d} \lambda_{T}^{*} \mathrm{~d} \beta_{\mathcal{T}} \mathrm{d} \theta
\end{aligned}
\]
where \(y_{\mathcal{T}}\) are the observed data at times \(\mathcal{T}, \beta_{\mathcal{T}}\) is the GP at the locations of \(y_{\mathcal{T}}\) and \(\left(y_{\mathcal{T}^{*}}, \beta_{\mathcal{T}^{*}}\right)\) represent these components at a finite collection of locations at times \(\mathcal{T}^{*}\).

The predictive distribution may be explored in various ways, especially in a point process context, by choosing convenient functions of the observations to analyse. This issue is illustrated in a simulated example that is presented in section 1.2 of the on-line supplementary material. Note that the same algorithm provides prediction of the IF \(\lambda\).

\section*{6. Applications}

We apply the proposed methodology to real and synthetic data sets to investigate its performance. Three examples with synthetic data consider unidimensional and bidimensional spatial models and a spatiotemporal model. Detailed results and comments can be found in section 1 of the on-line supplementary material. Two real data sets are used in this section to fit a spatial and a spatiotemporal model. The spatial example is also analysed by using an existing approximated methodology. All the simulations concerning our methodology are coded in Ox (Doornik, 2007) and run in a \(3.50-\mathrm{GHz}\) Intel i7 processor with 16 Gbytes of random-access memory. Results obtained for synthetic data are considerably good, especially considering that data were not generated from our models for the IF. The estimates recovered well the model components, their functionals and prediction. This was achieved by fast converging well mixing chains, with low auto-correlation function, for all synthetic and real applications.

\subsection*{6.1. Example 1: Lansing Woods data}

We analyse a data set referring to the location of white oak trees in a \(924 \mathrm{ft} \times 924 \mathrm{ft}\) plot in Lansing Woods, Clinton County, Michigan, USA. These data are available in the R package spatstat (Baddeley et al., 2015) and have also been analysed by some other researchers (see, for example, Baddeley et al. (2015)). It consists of the location of 448 white oak trees in \(S\), a rescaled square of size 10 . We adopt a constant mean function \(\mu=0\) and the covariance function given in expression (7) with \(\gamma=\frac{3}{2}, \sigma^{2}=2\) and \(\tau^{2}=2\). We also adopt a gamma(76,6) prior truncated to be below 15 for \(\lambda^{*}\).

The estimated IF that was obatained by our methodology is presented in Fig. 1. It shows a smoothly varying pattern, while also respecting the rate of occurrence of events. Results are based on an MCMC run for 500 iterations and estimates of the IF are obtained for a burn-in of 100 iterations. Each iteration takes around 30 s, mostly consumed in the inversion and Choleski decomposition of highly dimensional (of order \(10^{3}\) ) covariance matrices and in the embedded

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/03a0ce15-057d-4aa7-9f78-458d1a29488c-14.jpg?height=772&width=1426&top_left_y=1352&top_left_x=124}
\captionsetup{labelformat=empty}
\caption{Fig. 1. Exact analysis: (a) map of the posterior mean IF and observed data for the white oak example; (b) auto-correlation function of \(\Lambda\left([0,4]^{2}\right)\)}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/03a0ce15-057d-4aa7-9f78-458d1a29488c-15.jpg?height=1517&width=1409&top_left_y=185&top_left_x=145}
\captionsetup{labelformat=empty}
\caption{Fig. 2. Discretized analysis-map of the posterior mean IF and observed data for the white oak example, for configurations with area sizes decreasing towards 0 (all areas of each configuration are of equal size; the hyperparameters are set at \(\left(\mu, \sigma^{2}, \phi\right)=(0,4,2.5)\) in Igcp): (a) 1600 areas; (b) 10000 areas; (c) 40000 areas}
\end{figure}

Gibbs sampling of \(\beta\) (three iterations are enough). Our MCMC algorithm enables appropriately precise Monte Carlo estimates with a small number of iterations. As an example, the Monte Carlo estimate of the posterior mean of \(\Lambda\left([0,4]^{2}\right)\) is 126.504, with a percentage Monte Carlo error of \(0.47 \%\) (the number of events in this region is 126 ).

A discretization-based approach for the log-Gaussian Cox process (Møller et al., 1998) is also presented for comparison. The discretized methodology was run in the R package lgcp (Taylor et al., 2013). A number of hyperparameters values were considered and we report the smoothest estimate that was obtained. Fig. 2 shows the estimate of the IF for this methodology.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/03a0ce15-057d-4aa7-9f78-458d1a29488c-16.jpg?height=500&width=523&top_left_y=191&top_left_x=569}
\captionsetup{labelformat=empty}
\caption{Fig. 3. Province of New Brunswick and the area considered in the analysis}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1. Number of fires per year for the New Brunswick data}
\begin{tabular}{|l|l|}
\hline Year & Number of fires \\
\hline 1987 & 216 \\
\hline 1989 & 120 \\
\hline 1990 & 102 \\
\hline 1991 & 211 \\
\hline 1992 & 155 \\
\hline 1993 & 123 \\
\hline 1994 & 136 \\
\hline 1995 & 169 \\
\hline 1996 & 122 \\
\hline 1997 & 94 \\
\hline 1998 & 86 \\
\hline 1999 & 224 \\
\hline 2000 & 140 \\
\hline 2001 & 194 \\
\hline 2002 & 127 \\
\hline 2003 & 94 \\
\hline
\end{tabular}
\end{table}

The results indicate convergence to an unknown limit, qualitatively similar (but not equal) to our results. This was to be expected as the underlying models are not the same. In contrast, our results are already the estimates of a continuously varying IF.

\subsection*{6.2 Example 2: New Brunswick fires}

The data set is also provided in the R package spatstat. It is provided by the Department of Natural Resources of the province of New Brunswick, Canada, and consists of all fires falling under their jurisdiction for the years 1987-2003 inclusively (with the year 1988 omitted). We consider fires notified in the rectangular area (rescaled to \(45.5 \times 54.5\) ) that is shown in Fig. 3. Table 1 gives the number of fires per year.

We fit the dynamic model with \(f\left\{\beta_{t}(s), W_{t}(s)\right\}=\beta_{0, t}(s)\) and \(\beta_{0, t}(s)=\beta_{0, t-1}(s)+w_{t}(s)\), where \(\beta_{0,0}\) and \(w_{t}\) are GPs with the covariance function given in expression (7) and hyperparameters \(\left(0,1.75^{2}, 10\right)\) and \(\left(0,0.5^{2}, 15\right)\) respectively. We also adopt a time-independent structure for the

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/03a0ce15-057d-4aa7-9f78-458d1a29488c-17.jpg?height=1221&width=1393&top_left_y=193&top_left_x=157}
\captionsetup{labelformat=empty}
\caption{Fig. 4. Maps of the IF posterior mean in the New Brunswick fires example: years are ordered from the top row to the bottom row, and then from left to right in each row in the sequence 1987, 1989, 1990, ..., 2003}
\end{figure}
\(\lambda_{t}^{*}\)-parameters with \(\lambda_{t}^{*} \sim \operatorname{gamma}(15,100) 1(<0.4), \forall t\). The estimated IF is shown in Fig. 4. The results clearly show not only spatial but also temporal smoothness of the IF. The panels associated with each year exhibit similarity of spatial patterns over consecutive years as a direct consequence of the temporal dependence of our dynamic GP. This effect is more vividly seen at the edges of the area of study, in a U-shaped region with higher values of the IF, especially at the beginning and the end of the time window that was considered but still bearing compromise with the observed events.

\section*{7. Final remarks}

This paper proposes a novel methodology to perform exact Bayesian inference in spatiotemporal Cox processes in which the IF dynamics are described by a multivariate GP. We showed how the usual components of spatiotemporal point patterns such as trend, seasonality and covariates can be incorporated, with flexibility of their effects warranted by the GP prior.

The methodology is exact in the sense that no discrete approximation of the process is used and Monte Carlo error is the only source of inaccuracy. Inference is performed via MCMC
sampling; more specifically, a Gibbs sampler whose particular choice of blocking and sampling scheme leads to fast convergence. The validity of the methodology is established through the proofs of the main results. Finally, simulated and real data studies illustrate the methodology and provide empirical evidence of its efficiency.

This work may give rise to new problems and possibilities that may be considered in future work. An immediate extension of our models involves consideration of marks to the Poisson events. These marks may be described with a variety of components, whose effects are allowed to vary smoothly, in line with the models that are used for the IF. The introduction of non-spatiotemporal covariates is another extension that may be a very useful contribution to various areas of application. Finally, computational developments are still required to deal with very large data sets, which is a general problem when working with GPs. Computation with GPs is an area of very active and promising research that is likely to bring computational gains to our methodology against other existing methodologies in the near future. As a result, we can envision the development of software implementing our methods in the near future.

\section*{Acknowledgements}

The authors thank the referees and the Associate Editor for many useful comments that led to a much more improved version of the paper throughout. They also thank Gareth Roberts and Krzysztof Łatusziński for insightful discussions about MCMC sampling, Piotr Zwiernik for insightful discussions about matrices, Ryan Adams for providing the MATLAB code to run the algorithm of Adams et al. (2009), Jony A. Pinto, Jr, and Jesus E. Gamboa for their help in the data analyses with the discretized versions of the models and Daniel Simpson for drawing our attention to spatstat. The first author thanks the Fundação de Amparo à Pesquisa do Estado de Minas Gerais for financial support and the second author thanks Conselho Nacional de Desenvolvimento Científico e Tecnológico Brazil and the Fundação de Amparo à Pesquisa do Estado do Rio de Janeiro for financial support.

\section*{References}

Adams, R. P., Murray, I. and Mackay, D. J. C. (2009) Tractable nonparametric Bayesian inference in Poisson processes with Gaussian process intensities. In Proc. 26th Int. Conf. Machine Learning (eds L. Bottou and M. Littman). Montreal: Omnipress.

Alvarez, M. A., Rosasco, L. and Lawrence, N. D. (2012) Kernels for vector-valued functions: a review. Foundns Trends Mach. Learn., 4, 195-266.
Arellano-Valle, R. B. and Azzalini, A. (2006) On the unification of families of skew-normal distributions. Scand. J. Statist., 33, 561-574.

Baddeley, A., Rubak, E. and Turner, R. (2015) Spatial Point Patterns: Methodology and Applications with R. Boca Raton: Chapman and Hall-CRC.
Banerjee, A., Dunson, D. B. and Tokdar, S. T. (2013) Efficient Gaussian process regression for large datasets. Biometrika, 100, 75-89.
Beskos, A., Papaspiliopoulos, O., Roberts, G. O. and Fearnhead, P. (2006) Exact and computationally efficient likelihood-based inference for discretely observed diffusion processes (with discussion). J. R. Statist. Soc. B, 68, 333-382.
Boyle, P. and Frean, M. (2004) Dependent Gaussian processes. In Advances in Neural Information Processing Systems, vol. 17 (eds L. K. Saul, Y. Weiss and L. Bottou). Cambridge: MIT Press.
Brix, A. and Diggle, P. J. (2001) Spatiotemporal prediction for log-Gaussian Cox processes. J. R. Statist. Soc. B, 63, 823-841.
Carlin, B. P., Banerjee, S. and Finley, A. O. (2007) spBayes: an R package for univariate and multivariate hierarchical point-referenced spatial models. J. Statist. Softwr, 19, 1-24.
Carter, C. K. and Kohn, R. (1994) On Gibbs sampling for state space models. Biometrika, 81, 541-553.
Cox, D. R. (1955) Some statistical methods connected with series of events (with discussion). J. R. Statist. Soc. B, 17, 129-164.

Diggle, P. J. (2014) Statistical Analysis of Spatial and Spatio-temporal Point Patterns, 3rd edn. Boca Raton: Chapman and Hall.
Doornik, J. A. (2007) Object-oriented Matrix Programming using Ox, 3rd edn. London: Timberlake Consultants.
Fruhwirth-Schnatter, S. (1994) Data augmentation and dynamic linear models. J. Time Ser. Anal., 15, 183-202.
Gamerman, D. (1998) Markov chain Monte Carlo for dynamic generalised linear models. Biometrika, 85, 215-227.
Gamerman, D. (2010) Dynamic spatial models including spatial time series. In Handbook of Spatial Statistics (eds A. E. Gelfand, P. J. Diggle, M. Fuentes and P. Guttorp), pp. 437-448. Boca Raton: Chapman and HallCRC.
Gamerman, D., Salazar, E. and Reis, E. A. (2007) Dynamic Gaussian process priors with applications to the analysis of space-time data (with discussion). In Bayesian Statistics 8 (eds J. M. Bernardo, M. J. Bayarri, J. O. Berger, A. P. Dawid, D. Heckerman, A. F. M. Smith and M. West), pp. 149-174. Oxford: Oxford University Press.
Gamerman, D., Santos, T. R. and Franco, G. C. (2013) A non-Gaussian family of state-space models with exact marginal likelihood. J. Time Ser. Anal., 34, 625-645.
Gelfand, A. E., Banerjee, S. and Gamerman, D. (2005) Spatial process modelling for univariate and multivariate dynamic spatial data. Environmetrics, 16, 465-479.
Gonçalves, F. B. and Roberts, G. O. (2014) Exact simulation problems for jump-diffusions. Methodol. Comput. Appl. Probab., 16, 907-930.
Gonçalves, F. B., Roberts, G. O. and Łatuszyński, K. G. (2017) Exact Monte Carlo likelihood-based inference for jump-diffusion processes. To be published.
Kottas, A. and Sansó, B. (2007) Bayesian mixture modeling for spatial Poisson process intensities, with applications to extreme value analysis. J. Statist. Planng Inf., 137, 3151-3163.
Lewis, P. A. W. and Shedler, G. S. (1979) Simulation of a nonhomogeneous Poisson process by thinning. Navl Res. Logist. Q., 26, 403-413.
Møller, J., Syversveen, A. R. and Waagepetersen, R. P. (1998) Log Gaussian Cox processes. Scand. J. Statist., 25, 451-482.
Pinto, Jr, J. A., Gamerman, D., Paez, M. S. and Alves, R. H. F. (2015) Point pattern analysis with spatially varying covariate effects, applied to the study of cerebrovascular deaths. Statist. Med., 34, 1214-1226.
Reis, E. A., Gamerman, D., Paez, M. S. and Martins, T. G. (2013) Bayesian dynamic models for space-time point processes. Computnl Statist. Data Anal., 60, 146-156.
Roberts, G. O. and Rosenthal, J. S. (2009) Examples of adaptive MCMC. J. Computnl Graph. Statist., 18, 349-367.
Sermaidis, G., Papaspiliopoulos, O., Roberts, G. O., Beskos, A. and Fearnhead, P. (2013) Markov chain Monte Carlo for exact inference for diffusions. Scand. J. Statist., 40, 294-321.
Simpson, D., Illian, J. B., Lindgren, F., Sørbye, S. H. and Rue, H. (2016) Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103, 49-70.
Taylor, B. M., Davis, T. M., Rowlingson, B. S. and Diggle, P. J. (2013) lgcp: an R package for inference with spatial and spatio-temporal log-Gaussian Cox processes. J. Statist. Softwr., 52, 1-40.
Wackernagel, H. (2003) Multivariate Geostatistics. Berlin: Springer.
Wikle, C. and Cressie, N. (1999) A dimension-reduced approach to space-time Kalman filtering. Biometrika, 86, 815-829.
Xiao, S., Kottas, A. and Sansó, B. (2015) Modeling for seasonal marked point processes: an analysis of evolving hurricane occurrences. Ann. Appl. Statist., 9, 353-382.

\section*{Supporting information}

Additional 'supporting information' may be found in the on-line version of this article:
'Supplementary material'.