\title{
Dependent Bayesian nonparametric modeling of compositional data using random Bernstein polynomials*
}

\author{
Claudia Wehrhahn \\ Department of Statistics, University of California Santa Cruz, Santa Cruz, U.S.A. \\ e-mail: cwehrhah@ucsc.edu \\ Andrés F. Barrientos \\ Department of Statistics, Florida State University, Florida, U.S.A. \\ e-mail: abarrientos@fsu.edu \\ and
}

\author{
Alejandro Jara \\ MiDaS-Center for the Discovery of Structures in Complex Data and Department of Statistics, Pontificia Universidad Católica de Chile, Santiago, Chile \\ e-mail: atjara@uc.cl
}

\begin{abstract}
We discuss Bayesian nonparametric procedures for the regression analysis of compositional responses, that is, data supported on a multivariate simplex. The procedures are based on a modified class of multivariate Bernstein polynomials and on the use of dependent stick-breaking processes. A general model and two simplified versions of the general model are discussed. Appealing theoretical properties such as continuity, association structure, support, and consistency of the posterior distribution are established. Additionally, we exploit the use of spike-and-slab priors for choosing the version of the model that best adapts to the complexity of the underlying true data-generating distribution. The performance of the proposed model is illustrated in a simulation study and in an application to solid waste data from Colombia.
\end{abstract}

Keywords and phrases: Fully nonparametric regression, density regression, Dirichlet process, dependent Dirichlet processes.

Received August 2021.

\section*{1. Introduction}

The statistical modeling of compositional data plays a key role in many scientific areas, economics, political sciences, and engineering, among others, with

\footnotetext{
*C. Wehrhahn's research was supported by the "Programa de Becas de Postgrado de Chile, CONICYT", NSF-DMS 1738053 and ATD-DMS 1441433. A. Jara's work was supported by the Agencia Nacional de Investigación y Desarrollo (ANID) through the Fondo Nacional de Desarrollo Científico y Tecnológico (FONDECYT) grant No 1220907 and through grant NCN17_059 from Millennium Science Initiative Program, Millennium Nucleus Center for the Discovery of Structures in Complex Data (MIDAS).
}
applications ranging from human microbiome analyses to topic modeling in large text corpora, where the focus is on the modeling of vectors containing information on the relative frequencies in which the different components occur. The need for appropriate methodologies for the analysis of these type of data can be originated by the nature of the scientific questions or limitations associated with the measurement methods. For instance, recent advances in biological high-throughput sequencing technologies can only provide relative abundance information because they can only capture a limited total number of transcripts or sequences in a sample or do not control for the total number of microbes entering the measurement process. Thus, the resulting sequencing count data carry only relative abundance information about the different transcripts or taxa in a given sample. On the other hand, in our motivating problem, we are interested in describing and understanding the solid waste composition generated in a residential area of a city in Colombia.

From a mathematical point of view, compositional data can be defined as multivariate data supported on the \(m\)-dimensional simplex, \(\Delta_{m}\), given by
\[
\Delta_{m}=\left\{\left(y_{1}, \ldots, y_{m}\right) \in[0,1]^{m}: \sum_{i=1}^{m} y_{i} \leq 1\right\} .
\]

Since Aitchison (1982), several parametric regression models for compositional responses have been proposed. Common approaches transform the compositional responses from \(\Delta_{m}\) to \(\mathbb{R}^{m}\), and use the well known and familiar battery of statistical models for normally distributed responses (see, e.g., Aitchison, 1982; Atchison \& Shen, 1980; Shimizu et al., 2021; Wang et al., 2010). Other proposals use the Dirichlet distribution to model the compositional responses and link the Dirichlet parameters to covariates (see, e.g., Gueorguieva et al., 2008; Hijazi, 2003; Hijazi \& Jernigan, 2009; Van der Merwe, 2019). These models can be easily extended to allow for non-parametric functional forms in the relationship between the model parameters and the predictors (see, e.g., Di Marzio et al., 2015; Tsagris et al., 2020). However, they rely on particular parametric distributional forms which limits the type of inferences that can be obtained. Modeling approaches where the complete distribution of the compositional responses can flexibly vary as a function of the predictors are scarce in the literature. We aim to fill this gap by proposing a class of Bayesian nonparametric (BNP) predictordependent mixture models that enjoys appealing theoretical properties and is easy to use.

Most BNP approaches for collections of predictor-dependent probability distributions employ mixtures of densities from parametric families (see, e.g., Müller et al., 2015, and references therein). Mixture models are convenient for density estimation because they induce a prior distribution on densities by placing a prior distribution on the mixing measure. Dependent Dirichlet processes (MacEachern, 1999, 2000; Quintana et al., 2022) are often used as priors for the mixing distributions. Other extensions and alternative constructions for dealing with predictor-dependent probability distributions include the orderedcategory probit regression model (Karabatsos \& Walker, 2012), the dependent beta process (Trippa et al., 2011), the dependent tail-free processes (Jara \&

Hanson, 2011), the dependent neutral to the right processes and correlated twoparameter Poisson-Dirichlet processes (Epifani \& Lijoi, 2010; Leisen \& Lijoi, 2011), and the general class of dependent normalized completely random measures (Lijoi et al., 2014). Due to their flexibility and ease in computation, these models are routinely implemented in a wide variety of applications (see, e.g., Müller et al., 2015, and references therein).

BNP approaches for collections of predictor-dependent probability distributions have mainly focused on responses defined on the real line. Although those approaches can be applied to compositional responses, by transforming the responses from \(\Delta_{m}\) to \(\mathbb{R}^{m}, \mathrm{i}\) ) the resulting density in the simplex could not be well defined at the edges or ii) the resulting density in the simplex could be equal to zero at the edges. This can cause important problems if zeros are observed in the data because either the likelihood could not be defined, if i) holds true, or the likelihood would always be equal to zero, if ii) holds true. Also, other problem associated with the use of transformations is that it is not very clear that the resulting density is flexible at the edges of the simplex (please see Appendix A for more details on this).

We propose modeling compositional responses using a particular class of mixtures of Dirichlet probability density functions that naturally emerges from the theoretical properties and extensions of Bernstein polynomials (BP). Motivated by their uniform approximation properties, frequentist and Bayesian methods based on univariate BP have been proposed for the estimation of probability distributions supported on bounded intervals, unit hyper-cubes, and simplex spaces (see, e.g. Petrone, 1999a,b; Petrone \& Wasserman, 2002; Tenbusch, 1994; Ouimet, 2021; Babu \& Chaubey, 2006; Zheng et al., 2010). For example, Babu \& Chaubey (2006) studied a general multivariate version of the bivariate estimator proposed by Tenbusch (1994), while Zheng et al. (2010) constructed a multivariate Bernstein polynomial (MBP) prior for the spectral density of a random field. Key for our approach, Tenbusch (1994) considered multivariate extensions of Bernstein polynomials defined on \(\Delta_{2}\) to propose and study a density estimator. Tenbusch's approach is easy to extend to the \(m\)-dimensional case. The approach is based on the class of MBP associated with \(G\), a cumulative distribution function (CDF) on \(\Delta_{m}\), given by
\[
\widetilde{B}_{k, G}(\boldsymbol{y})=\sum_{\boldsymbol{j} \in \mathscr{J}_{k, m}} G\left(\frac{j_{1}}{k}, \ldots, \frac{j_{m}}{k}\right) \operatorname{Mult}(\mathbf{j} \mid k, \boldsymbol{y}), \boldsymbol{y} \in \Delta_{m}
\]
where \(k \in \mathbb{N}\) is the degree of the MBP, \(\boldsymbol{j}=\left(j_{1}, \ldots, j_{m}\right)\),
\[
\mathscr{J}_{k, m}=\left\{\left(j_{1}, \ldots, j_{m}\right) \in\{0, \ldots, k\}^{m}: \sum_{l=1}^{m} j_{l} \leq k\right\},
\]
and Mult ( \(\cdot \mid k, \boldsymbol{y}\) ) stands for the probability mass function of a multinomial distribution with parameters \((k, \boldsymbol{y})\).

Tenbusch's estimator arises by replacing \(G\) in (1) by the empirical CDF of observed data and it is not difficult to show that, if \(G\) is a CDF on \(\Delta_{m}\), then
\(\widetilde{B}_{k, G}(\cdot)\) is not a CDF on \(\Delta_{m}\) for a finite \(k\). In this case, \(\widetilde{B}_{k, G}(\cdot)\) can be expressed as a linear combination of CDFs of probability measures defined on \(\Delta_{m}\), where the coefficients are nonnegative but do not add up to 1 . Tenbusch's estimator is defined as the derivative of \(\widetilde{B}_{k, G}(\cdot)\) and, although it is consistent and optimal at the interior points of the simplex, it is not a valid density function for finite \(k\) and finite sample size. To avoid this problem, Barrientos et al. (2015) proposed a modified class of MBP by changing the set \(\mathscr{J}_{k, m}\). The class retains the well known approximation properties of the original version. Furthermore, when \(G\) is a CDF on \(\Delta_{m}\), the modified MBP is a genuine CDF with density function defined by a mixture of Dirichlet densities.

We propose a class of fully nonparametric regression models for compositional responses, by extending the class of MBP priors of Barrientos et al. (2015). The extension relies on predictor-dependent stick-breaking processes (see, Barrientos et al., 2017, for a similar extension for responses defined on the unit-interval). An important property of the considered model class is that the densities are well-defined in scenarios with compositional data containing zero values. When the response vector contains zero values, either models based on the Dirichlet distribution without restrictions on the parameter space or approaches based on the log-ratio transformation are not properly defined and cannot be employed unless a zero value imputation is applied first (please see Appendix A). We study theoretical properties of the proposed model class, such as continuity, association structure, support, and consistency of the posterior distribution. These properties are non-trivial extensions of the results obtained by Barrientos et al. (2017) for the unit-interval and their proofs are provided in the Appendix. The use of the dependent stick-breaking process raises the question of where to introduce the predictor dependency: on weights, atoms, or both. Each selection leading to a different version of the model. Rather than fitting all versions of the model, as done by Barrientos et al. (2017), we use spike-and-slab mixtures (George \& McCulloch, 1993) to define a prior that automatically chooses the version of the model that best accommodates to the complexity of the underlying data-generating mechanism. We evaluate the performance of the proposed approach using simulated data. The proposed approach is also applied for the analysis of solid waste data from Colombia.

The rest of the paper is organized as follows. The modified class of MBP and its main properties are summarized in Section 2. The proposed model class and its theoretical properties are discussed in Sections 3 and 4, respectively. Section 5 describes the main computational aspects. Section 6 illustrates the performance of the model using simulated data and in an application to solid waste in Colombia. A final discussion concludes the article.

\section*{2. Random multivariate Bernstein polynomials}

Based on Tenbusch's MBP, Barrientos et al. (2015) defined a modified class of MBP on the \(m\)-dimensional simplex and proposed a BNP density estimation model for compositional data. The modified class increases the domain
of function \(G\) and the size of the set \(\mathscr{J}_{k, m}\) from the original class of MBP on the \(m\)-dimensional simplex provided in Equation (1). For a given function \(G: \mathbb{R}^{m} \longrightarrow \mathbb{R}\), the associated modified class of MBP of degree \(k \in \mathbb{N}\) on \(\Delta_{m}\) is given by
\[
B(\boldsymbol{y} \mid k, G)=\sum_{\boldsymbol{j} \in \mathscr{H}_{k, m}} G\left(\frac{j_{1}}{k}, \ldots, \frac{j_{m}}{k}\right) \operatorname{Mult}(\mathbf{j} \mid k+m-1, \boldsymbol{y}), \boldsymbol{y} \in \Delta_{m}
\]
where \(\mathscr{H}_{k, m}=\left\{\left(j_{1}, \ldots, j_{m}\right) \in\{0, \ldots, k\}^{m}: \sum_{l=1}^{m} j_{l} \leq k+m-1\right\}\).
As shown by Barrientos et al. (2015), this class of MBP retains the appealing approximation properties of univariate BP and the standard class of MBP given in Equation (1). Specifically, if \(G\) is a real-valued function defined on \(\mathbb{R}^{m}\) and \(\left.G\right|_{\Delta_{m}}\) is its restriction on \(\Delta_{m}\), then \(B(\cdot \mid k, G)\) converges pointwise to \(\left.G\right|_{\Delta_{m}}\), as \(k\) goes to infinity, and the relation holds uniformly on \(\Delta_{m}\) if \(\left.G\right|_{\Delta_{m}}\) is a continuous function.

It is also possible to show that if \(G\) is the CDF of a probability measure defined on \(\Delta_{m}\), then \(B(\cdot \mid k, G)\) is also the restriction of the CDF of a probability measure defined on \(\Delta_{m}\). Furthermore, if \(G\) is the CDF of a probability measure defined on \(\Delta_{m}^{0}=\left\{\boldsymbol{y} \in \Delta_{m}: y_{j}>0, j=1, \ldots, m\right\}\), then \(B(\cdot \mid k, G)\) is the restriction of the CDF of a probability measure with density function given by the following mixture of Dirichlet distributions,
\[
b(\boldsymbol{y} \mid k, G)=\sum_{\boldsymbol{j} \in \mathscr{H}_{k, m}^{0}} G\left(\left(\frac{j_{1}-1}{k}, \frac{j_{1}}{k}\right] \times \ldots \times\left(\frac{j_{m}-1}{k}, \frac{j_{m}}{k}\right]\right) \operatorname{dir}(\boldsymbol{y} \mid \alpha(k, \boldsymbol{j})),
\]
where \(\mathscr{H}_{k, m}^{0}=\left\{\left(j_{1}, \ldots, j_{m}\right) \in\{1, \ldots, k\}^{m}: \sum_{l=1}^{m} j_{l} \leq k+m-1\right\}, \alpha(k, \boldsymbol{j})= \left(\boldsymbol{j}, k+m-\|\mathbf{j}\|_{1}\right),\|\cdot\|_{1}\) denotes the \(l_{1}\)-norm, and \(\operatorname{dir}\left(\cdot \mid\left(\alpha_{1}, \ldots, \alpha_{m+1}\right)\right)\) denotes the density function of an \(m\)-dimensional Dirichlet distribution with parameters \(\left(\alpha_{1}, \ldots, \alpha_{m+1}\right)\).

By considering the density function given by Equation (2), a random function \(G\), and a random degree \(k\), Barrientos et al. (2015) defined a BNP prior for densities defined on \(\Delta_{m}\). The model corresponds to a DP mixture model of specific Dirichlet densities given by
\[
\begin{aligned}
b(\boldsymbol{y} \mid k, G) & =\int_{\Delta_{m}} \operatorname{dir}(\boldsymbol{y} \mid \alpha(k,\lceil k \boldsymbol{\theta}\rceil)) G(d \boldsymbol{\theta}) \\
G \mid M, G_{0} & \sim D P\left(M, G_{0}\right) \\
k \mid \lambda & \sim p(\cdot \mid \lambda)
\end{aligned}
\]
where \(D P\left(\alpha, G_{0}\right)\) denotes a Dirichlet process with concentration parameter \(M>\) 0 and base distribution \(G_{0}\) on \(\Delta_{m}^{0}, p(\cdot \mid \lambda)\) is the probability mass function of a distribution on \(\mathbb{N}\) parameterized by \(\lambda\), and \(\lceil\cdot\rceil\) denotes the ceiling function.

\section*{3. The model}

Suppose that we observe regression data \(\left\{\left(\boldsymbol{y}_{i}, \boldsymbol{x}_{i}\right): i=1, \ldots, n\right\}\), where \(\boldsymbol{y}_{i}\) is a continuous \(\Delta_{m}\)-valued outcome vector and \(\boldsymbol{x}_{i} \in \mathscr{X} \subseteq \mathbb{R}^{p}\) is a \(p\)-dimensional
vector of exogenous predictors. We define the regression model for compositional responses by introducing predictor-dependency in the mixture model given in (3), which allows the complete shape of the conditional densities to flexibly vary with values of \(\boldsymbol{x}\). To this end, we replace the mixing measure \(G\) by a predictordependent mixing measure \(G_{\boldsymbol{x}}\). Under this approach, the random conditional densities are given by
\[
f_{\boldsymbol{x}}\left(\boldsymbol{y} \mid k, G_{\boldsymbol{x}}\right)=\int_{\Delta_{m}} \operatorname{dir}(\boldsymbol{y} \mid \alpha(k,\lceil k \boldsymbol{\theta}\rceil)) G_{\boldsymbol{x}}(d \boldsymbol{\theta})
\]
where the set of mixing distributions \(\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) follows a dependent stickbreaking process, with elements of the form \(G_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}) \delta_{\boldsymbol{\theta}_{j}(\boldsymbol{x})}(\cdot)\), where \(w_{j}(\boldsymbol{x})=V_{j}(\boldsymbol{x}) \prod_{l<j}\left[1-V_{l}(\boldsymbol{x})\right]\), and where \(V_{j}(\boldsymbol{x})\) and \(\boldsymbol{\theta}_{j}(\boldsymbol{x})\) are transformations of underlying stochastic processes.

\subsection*{3.1. The formal definition}

Let \(\mathscr{V}=\left\{v_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) and \(\mathscr{H}=\left\{h_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) be two sets of known bijective continuous functions, such that for every \(\boldsymbol{x} \in \mathscr{X}, v_{\boldsymbol{x}}: \mathbb{R} \longrightarrow[0,1]\) and \(h_{\boldsymbol{x}}: \mathbb{R}^{m} \longrightarrow \Delta_{m}^{0}\), are such that for every \(a \in \mathbb{R}\) and \(\boldsymbol{b} \in \mathbb{R}^{m}, v_{\boldsymbol{x}}(a)\) and \(h_{\boldsymbol{x}}(\boldsymbol{b})\) are continuous functions of \(\boldsymbol{x}\). Let \(\mathscr{P}\left(\Delta_{m}\right)\) be the set of all probability measures defined on \(\Delta_{m}\).

Definition 1. Let \(\mathscr{V}\) and \(\mathscr{H}\) be two sets of functions as before. Let \(F= \left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) be a \(\mathscr{P}\left(\Delta_{m}\right)\)-valued stochastic process such that:
(i) \(\eta_{j}=\left\{\eta_{j}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}, j \geq 1\), are independent and identically distributed real-valued stochastic processes with law indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{\eta}\).
(ii) \(\boldsymbol{z}_{j}=\left\{\boldsymbol{z}_{j}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}, j \geq 1\), are independent and identically distributed real-valued stochastic processes with law indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{\boldsymbol{z}}\).
(iii) \(k \in \mathbb{N}\) is a discrete random variable with distribution indexed by a finitedimensional parameter \(\lambda\).
(iv) For every \(\boldsymbol{x} \in \mathscr{X}\), the density function of \(F_{\boldsymbol{x}}\), w.r.t. Lebesgue measure, is given by the following dependent mixture of Dirichlet densities,
\[
f_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}) \operatorname{dir}\left(\cdot \mid \alpha\left(k,\left\lceil k \boldsymbol{\theta}_{j}(\boldsymbol{x})\right\rceil\right)\right),
\]
where \(\boldsymbol{\theta}_{j}(\boldsymbol{x})=h_{\boldsymbol{x}}\left(\boldsymbol{z}_{j}(\boldsymbol{x})\right),\left\lceil k \boldsymbol{\theta}_{j}(\boldsymbol{x})\right\rceil=\left(\left\lceil k \theta_{j 1}(\boldsymbol{x})\right\rceil, \ldots,\left\lceil k \theta_{j m}(\boldsymbol{x})\right\rceil\right)\), and \(w_{j}(\boldsymbol{x})=V_{j}(\boldsymbol{x}) \prod_{l<j}\left[1-V_{l}(\boldsymbol{x})\right]\), with \(V_{j}(\boldsymbol{x})=v_{\boldsymbol{x}}\left(\eta_{j}(\boldsymbol{x})\right)\).

The process \(F=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) will be referred to as dependent MBP process with parameters \(\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \mathscr{H}\right)\), and denoted by \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \mathscr{H}\right)\) and DMBPP for short.

In the search of parsimonious models, it is of interest to study two special cases of the general construction given by Definition 1. The case involving dependent stick-breaking processes with common weights and predictor-dependent support points is referred to as 'single-weights' DMBPP, while the case involving dependent stick-breaking processes with common support points and predictor-dependent weights is referred to as 'single-atoms' DMBPP. In what follows, we briefly discuss the definition of these special cases. Their formal definitions, needed to provide the proofs of the theoretical properties discussed in the following sections, are provided in Appendix B.

In the definition of the 'single-weights' DMBPP, the real-valued stochastic processes of condition (i) in Definition \(1, \eta_{j}=\left\{\eta_{j}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\), are replaced by \([0,1]\)-valued independent and identically distributed random variables, \(v_{j}\), with common distribution indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{v}\). In this special case, the density function of \(F_{\boldsymbol{x}}\) is given by
\[
f_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{\infty} w_{j} \operatorname{dir}\left(\cdot \mid \alpha\left(k,\left\lceil k \boldsymbol{\theta}_{j}(\boldsymbol{x})\right\rceil\right)\right),
\]
where \(\boldsymbol{\theta}_{j}(\boldsymbol{x})\) and \(\left\lceil k \boldsymbol{\theta}_{j}(\boldsymbol{x})\right\rceil\) are defined as in Definition 1 and \(w_{j}=v_{j} \prod_{l<j}\left[1-v_{l}\right]\). The process \(F=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) will be referred to as single-weight dependent MBP process with parameters \(\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\), and denoted by \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}\right.\), \(\left.\boldsymbol{\Psi}_{z}, \mathscr{H}\right)\) and \(w\) DMBPP for short.

In the definition of the 'single-atoms' DMBPP, the real-valued stochastic processes of condition (ii) in Definition \(1, \boldsymbol{z}_{j}=\left\{\boldsymbol{z}_{j}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\), are replaced by independent and identically distributed \(\Delta_{m}^{0}\)-valued random vectors, \(\boldsymbol{\theta}_{j}\), with common distribution indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{\boldsymbol{\theta}}\). In this case, the density function of \(F_{\boldsymbol{x}}\) is given by
\[
f_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}) \operatorname{dir}\left(\cdot \mid \alpha\left(k,\left\lceil k \boldsymbol{\theta}_{j}\right\rceil\right)\right),
\]
where \(w_{j}(\boldsymbol{x})\) are defined as in Definition 1 and \(\left\lceil k \boldsymbol{\theta}_{j}\right\rceil=\left(\left\lceil k \theta_{j 1}\right\rceil, \ldots,\left\lceil k \theta_{j m}\right\rceil\right)\). The process \(F=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) will be referred to as single-atoms dependent MBP process with parameters \(\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\), and denoted by \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}\right.\), \(\mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\) ) and \(\theta\) DMBPP for short.

Notice that the DMBPP, including its special cases, is well defined if the mapping induced by (iv) in Definition 1 is measurable, which is discussed in detail in Section 3.2. Notice also that expressions (5), (6), and (7) are indeed a density w.r.t. Lebesgue measure since, for every \(\boldsymbol{x} \in \mathscr{X}\),
\[
\sum_{j=1}^{\infty} \log \left[1-\mathrm{E}\left(v_{\boldsymbol{x}}\left\{\eta_{j}(\boldsymbol{x})\right\}\right)\right]=-\infty, \quad \text { and } \quad \sum_{j=1}^{\infty} \log \left[1-\mathrm{E}\left(v_{j}\right)\right]=-\infty
\]
which are sufficient and necessary conditions for the corresponding weights to add up to one with probability one. It is important to emphasize that DMBPP
generates dependent mixture of Dirichlet densities with constant support points and covariate-dependent weights,
\[
f_{\boldsymbol{x}}(\cdot)=\sum_{\boldsymbol{j} \in \mathcal{H}_{k, m}^{0}} W_{k, \boldsymbol{j}, \boldsymbol{x}} \times \operatorname{dir}(\cdot \mid \alpha(k, \boldsymbol{j})),
\]
where
\[
W_{k, \boldsymbol{j}, \boldsymbol{x}}=\left\{\begin{array}{l}
\sum_{l=1}^{\infty} w_{l}(\boldsymbol{x}) \delta_{\boldsymbol{\theta}_{l}(\boldsymbol{x})}\left(\left(\frac{j_{1}-1}{k}, \frac{j_{1}}{k}\right] \times \ldots \times\left(\frac{j_{m}-1}{k}, \frac{j_{m}}{k}\right]\right), \\
\sum_{l=1}^{\infty} w_{l} \delta_{\boldsymbol{\theta}_{l}(\boldsymbol{x})}\left(\left[\frac{j_{1}-1}{k}, \frac{j_{1}}{k}\right] \times \ldots \times\left(\frac{j_{m}-1}{k}, \frac{j_{m}}{k}\right]\right), \\
\sum_{l=1}^{\infty} w_{l}(\boldsymbol{x}) \delta_{\boldsymbol{\theta}_{l}}\left(\left(\frac{j_{1}-1}{k}, \frac{j_{1}}{k}\right] \times \ldots \times\left(\frac{j_{m}-1}{k}, \frac{j_{m}}{k}\right]\right),
\end{array}\right.
\]
for the DMBPP, \(w \mathrm{DMBPP}\), and \(\theta \mathrm{DMBPP}\), respectively.

\subsection*{3.2. The measurability of the processes}

In this section we show that the corresponding mappings defining the trajectories of DMBPP, \(w \mathrm{DMBPP}\), and \(\theta \mathrm{DMBPP}\) are measurable under the Borel \(\sigma\)-field generated by the weak product topology, \(L_{\infty}\) product topology, and \(L_{\infty}\) topology, which correspond to generalizations of standard topologies for spaces of single probability measures. The topologies considered here are formally defined in Appendix C.

Let \(\mathscr{D}\left(\Delta_{m}\right) \subset \mathscr{P}\left(\Delta_{m}\right)\) be the space of all probability measures defined on \(\Delta_{m}\) that are absolutely continuous w.r.t. Lebesgue measure and with continuous density function and consider the spaces \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}=\prod_{\boldsymbol{x} \in \mathscr{X}} \mathscr{P}\left(\Delta_{m}\right)\) and \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}=\prod_{\boldsymbol{x} \in \mathscr{X}} \mathscr{D}\left(\Delta_{m}\right)\). Theorem 1, which proof is provided in Appendix D.1, summarizes the measurability results for the different versions of the proposed model.

Theorem 1. Let \(\mathscr{B}_{1}, \mathscr{B}_{2}\), and \(\mathscr{B}_{3}\) be the Borel \(\sigma\)-field generated by the weak product topology, \(L_{\infty}\) product topology, and \(L_{\infty}\) topology, respectively. If \(F\) is a DMBPP, \(w \mathrm{DMBPP}\) or \(\theta \mathrm{DMBPP}\), defined on the appropriate measurable space \((\Omega, \mathscr{A})\), then the following mappings are measurable:
- \(F:(\Omega, \mathcal{A}) \longrightarrow\left(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}, \mathscr{B}_{1}\right)\).
- \(F:(\Omega, \mathcal{A}) \longrightarrow\left(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}, \mathscr{B}_{2}\right)\).
- \(F:(\Omega, \mathcal{A}) \longrightarrow\left(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}, \mathscr{B}_{3}\right)\).

\section*{4. The main properties}

We establish basic properties of the proposed class of models in this section. They include the characterization of the topological support, the continuity and association structure of the models, and the asymptotic behavior of the posterior distribution. Detailed proofs of Theorems 2-12 are provided in Appendix D. 2 - D.12, respectively.

\subsection*{4.1. The support of the processes}

Full support is a "necessary" property for a Bayesian model to be considered "nonparametric". In a fully nonparametric regression model setting, full support implies that the prior probability model assigns positive mass to any neighborhood of every collection of probability measures \(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\). Therefore, the definition of support strongly depends on the choice of a "distance" defining the basic neighborhoods. We provide sufficient conditions for \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) and \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) to be the support of DMBPPs under the weak product topology and the \(L_{\infty}\) product topology, respectively.

Theorem 2. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\), or a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). If \(F\) is defined such that:
(i) for every \(\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{L}\right) \in \mathscr{X}^{L}, L \geq 1\), the joint distribution of \(\left(\eta_{j}\left(\boldsymbol{x}_{1}\right), \ldots\right.\), \(\left.\eta_{j}\left(\boldsymbol{x}_{L}\right)\right), j \geq 1\), has full support on \(\mathbb{R}^{L}\),
(ii) for every \(\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{L}\right) \in \mathscr{X}^{L}, L \geq 1\), the joint distribution of \(\left(\boldsymbol{z}_{j}\left(\boldsymbol{x}_{1}\right), \ldots\right.\), \(\left.\boldsymbol{z}_{j}\left(\boldsymbol{x}_{L}\right)\right), j \geq 1\), has full support on \(\mathbb{R}^{m \times L}\),
(iii) \(k\) has full support on \(\mathbb{N}\),
(iv) \(v_{j}, j \geq 1\), has full support on \([0,1]\),
(v) \(\boldsymbol{\theta}_{j}, j \geq 1\), has full support on \(\Delta_{m}^{0}\),
then \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) and \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the weak product topology and the \(L_{\infty}\) product topology, respectively.

If stronger assumptions on the parameter space are imposed, a stronger support property can be obtained. Specifically, consider the sub-space \(\tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}} \subset \mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\), where
\(\tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}=\left\{\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}:(\boldsymbol{y}, \boldsymbol{x}) \longmapsto q_{\boldsymbol{x}}(\boldsymbol{y})\right.\) is continuous \(\}\),
and \(q_{\boldsymbol{x}}\) denotes the density function of \(Q_{\boldsymbol{x}}\) w.r.t. Lebesgue measure. The following theorem provides sufficient conditions for \(\tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) to be in the support of DMBPPs under the \(L_{\infty}\) topology.

Theorem 3. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) or a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). Assume that \(\boldsymbol{x} \in \mathscr{X}\) contains only continuous components and that \(\mathscr{X}\) is compact. If \(F\) is defined such that:
(i) for every \(B \in \mathscr{B}\left(\Delta_{m}\right)\), every \(\Delta_{m}^{0}\)-valued continuous mapping \(\boldsymbol{x} \mapsto f(\boldsymbol{x})\), and every \(j \geq 1\),
\[
\operatorname{Pr}\left\{\sup _{\boldsymbol{x} \in \mathscr{X}}\left|h_{\boldsymbol{x}}\left(\boldsymbol{z}_{j}(\boldsymbol{x})\right)-f_{\boldsymbol{x}}\right| \in B\right\}>0
\]
(ii) for every \(\epsilon>0\), every \([0,1]\)-valued continuous mapping \(\boldsymbol{x} \mapsto f(\boldsymbol{x})\), and every \(j \geq 1\),
\[
\operatorname{Pr}\left\{\sup _{\boldsymbol{x} \in \mathscr{X}}\left|v_{\boldsymbol{x}}\left(\eta_{j}(\boldsymbol{x})\right)-f_{\boldsymbol{x}}\right|<\epsilon\right\}>0
\]
(iii) \(k\) has full support on \(\mathbb{N}\),
(iv) \(v_{j}, j \geq 1\), has full support on \([0,1]\),
(v) \(\boldsymbol{\theta}_{j}, j \geq 1\), has full support on \(\Delta_{m}^{0}\),
then \(\tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) is contained in the support of \(F\) under the \(L_{\infty}\) topology.
An important consequence of the previous theorem is that the proposed processes can assign positive mass to arbitrarily small neighborhoods of any collection of probability measures \(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\), based on the supremum over the predictor space of Kullback-Leibler (KL) divergences between the predictor-dependent probability measures.

Theorem 4. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) or a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). Under the same assumptions of Theorem 3, it follows that
\[
\operatorname{Pr}\left\{\sup _{\boldsymbol{x} \in \mathscr{X}} \int_{\Delta_{m}} q_{\boldsymbol{x}}(\boldsymbol{y}) \log \left(\frac{q_{\boldsymbol{x}}(\boldsymbol{y})}{f_{\boldsymbol{x}}(\boldsymbol{y})}\right) d \boldsymbol{y}<\epsilon\right\}>0
\]
for every \(\epsilon>0\), and every \(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) with density functions \(\left\{q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\).

\subsection*{4.2. The continuity and association structure of the processes}

The characteristics of the stochastic processes used in the definitions of aDMBPP determine important properties of the resulting model. Regardless of the specific choice of the stochastic processes used in its definition, the use of almost surely (a.s.) continuous stochastic processes ensures that DMBPP and \(w\) DMBPP have a.s. a limit.

Theorem 5. Let \(F\) be \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\) or \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\), defined such that \(\mathscr{V}\) and \(\mathscr{H}\) are sets of equicontinuous functions of \(\boldsymbol{x}\), and for every \(i \geq 1\), the stochastic processes \(\eta_{i}\) and \(\boldsymbol{z}_{i}\) have a.s. continuous trajectories. Then, for every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), with \(\boldsymbol{x}_{l} \in \mathscr{X}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l}=\boldsymbol{x}_{0}, F_{\boldsymbol{x}}\) has a.s. a limit with the total variation norm.

An interesting property of the \(\theta\) DMBPP compared to the other version, and the general model, is that the use of a.s. continuous stochastic processes in the weights guarantees a.s. continuity of the 'single-atoms' DMBPP.

Theorem 6. Let \(F\) be a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\), defined such that \(\mathscr{V}\) is a set of equicontinuous functions, and such that for every \(j \geq 1\), the stochastic process \(\eta_{j}\) is a.s. continuous. Then, for every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), with \(\boldsymbol{x}_{l} \in \mathscr{X}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l}=\boldsymbol{x}_{0}\),
\[
\lim _{l \rightarrow \infty} \sup _{B \in \mathscr{B}\left(\Delta_{m}\right)}\left|F_{\boldsymbol{x}_{l}}(B)-F_{\boldsymbol{x}_{0}}(B)\right|=0, \text { a.s.. }
\]

That is, \(F_{\boldsymbol{x}_{l}}\) converges a.s. in total variation norm to \(F_{\boldsymbol{x}_{0}}\), when \(\boldsymbol{x}_{l} \longrightarrow \boldsymbol{x}_{0}\).

The dependence structure of DMBPPs is completely determined by the association structure of the stochastic processes used in their definition. For instance, under mild conditions on the stochastic processes defining the DMBPPs, the correlation between the corresponding random measures approaches to one as the predictor values get closer.

Theorem 7. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) or a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\), defined such that \(\mathscr{V}\) and \(\mathscr{H}\) are sets of equicontinuous functions, and such that for every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), with \(\boldsymbol{x}_{l} \in \mathscr{X}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l}=\boldsymbol{x}_{0}\), we have \(\eta_{j}\left(\boldsymbol{x}_{l}\right) \xrightarrow{\mathscr{L}} \eta_{j}\left(\boldsymbol{x}_{0}\right)\) and \(\boldsymbol{z}_{j}\left(\boldsymbol{x}_{l}\right) \xrightarrow{\mathscr{L}} \boldsymbol{z}_{j}\left(\boldsymbol{x}_{0}\right)\), as \(l \rightarrow \infty\), \(j \geq 1\). Then, for every \(\boldsymbol{y} \in \tilde{\Delta}_{m}\),
\[
\lim _{l \rightarrow \infty} \rho\left[F_{\boldsymbol{x}_{l}}\left(B_{\boldsymbol{y}}\right), F_{\boldsymbol{x}_{0}}\left(B_{\boldsymbol{y}}\right)\right]=1,
\]
where \(\rho(A, B)\) denotes the Pearson correlation between \(A\) and \(B, B_{\boldsymbol{y}}=\left[0, y_{1}\right] \times \ldots \times\left[0, y_{m}\right]\).

If the stochastic processes defining the DMBPP and \(w \mathrm{DMBPP}\) are such that the pairwise finite-dimensional distributions converge to the product of the corresponding marginal distributions as the Euclidean distance between the predictors grows larger, then under mild conditions the correlation between the corresponding random measures can approach zero. The following theorem shows that under the assumptions previously discussed, the marginal covariance between the random measures is equal to the covariance between the conditional expectations of the random measures, given the degree of the MBP.
Theorem 8. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\) or a \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}\right.\), \(\boldsymbol{\Psi}_{z}, \mathscr{H}\) ), defined such that \(\mathscr{V}\) and \(\mathscr{H}\) are sets of equicontinuous functions and there exists a constant \(\gamma>0\) such that if \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) and \(\left\|\boldsymbol{x}_{1}-\boldsymbol{x}_{2}\right\|> \gamma\), then \(\operatorname{Cov}\left[\mathbb{I}_{\left\{\eta_{j}\left(\boldsymbol{x}_{1}\right) \in A_{1}\right\}}, \mathbb{I}_{\left\{\eta_{j}\left(\boldsymbol{x}_{2}\right) \in A_{2}\right\}}\right]=0\), for every \(A_{1}, A_{2} \in \mathscr{B}(\mathbb{R})\), and \(\operatorname{Cov}\left[\mathbb{I}_{\left\{\boldsymbol{z}_{j}\left(\boldsymbol{x}_{1}\right) \in A_{3}\right\}}, \mathbb{I}_{\left\{\boldsymbol{z}_{j}\left(\boldsymbol{x}_{2}\right) \in A_{4}\right\}}\right]=0\), for every \(A_{3}, A_{4} \in \mathscr{B}\left(\mathbb{R}^{m}\right), j \geq 1\). Assume also that for every \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) and for every sequence \(\left\{\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)\right\}_{l=1}^{\infty}\), with \(\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right) \in \mathscr{X}^{2}\) and such that \(\lim _{l \rightarrow \infty}\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)=\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right)\), we have that
\[
\left(\eta_{j}\left(\boldsymbol{x}_{1 l}\right), \eta_{j}\left(\boldsymbol{x}_{2 l}\right)\right) \xrightarrow{\mathscr{L}}\left(\eta_{j}\left(\boldsymbol{x}_{1}\right), \eta_{j}\left(\boldsymbol{x}_{2}\right)\right),
\]
and
\[
\left(\boldsymbol{z}_{j}\left(\boldsymbol{x}_{1 l}\right), \boldsymbol{z}_{j}\left(\boldsymbol{x}_{2 l}\right)\right) \xrightarrow{\mathscr{L}}\left(\boldsymbol{z}_{j}\left(\boldsymbol{x}_{1}\right), \boldsymbol{z}_{j}\left(\boldsymbol{x}_{2}\right)\right),
\]
\(j \geq 1\), as \(l \rightarrow \infty\). Then, for every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\lim _{l \rightarrow \infty} \operatorname{Cov}\left[F_{\boldsymbol{x}_{1 l}}\left(B_{\boldsymbol{y}}\right), F_{\boldsymbol{x}_{2 l}}\left(B_{\boldsymbol{y}}\right)\right]=\operatorname{Cov}\left[E\left\{F_{\boldsymbol{x}_{1}}\left(B_{\boldsymbol{y}}\right) \mid k\right\}, E\left\{F_{\boldsymbol{x}_{2}}\left(B_{\boldsymbol{y}}\right) \mid k\right\}\right],
\]
with
\[
E\left\{F_{\boldsymbol{x}}\left(B_{\boldsymbol{y}}\right) \mid k\right\}=\sum_{\mathbf{j} \in \mathcal{H}_{k, m}} G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}, k}\right) \operatorname{Mult}(\mathbf{j} \mid k+m-1, \boldsymbol{y}),
\]
where \(B_{\boldsymbol{y}}=\left[0, y_{1}\right] \times \ldots \times\left[0, y_{m}\right], A_{\mathbf{j}, k}=\left[0, j_{1} / k\right] \times \ldots \times\left[0, j_{m} / k\right]\) and \(G_{0, \boldsymbol{x}}\) is the marginal probability measure of \(\boldsymbol{\theta}_{j}(\boldsymbol{x})\) defined on \(\Delta_{m}^{0}\).

Notice that the assumption \(\operatorname{Cov}\left[\mathbb{I}_{\left\{\eta_{j}\left(\boldsymbol{x}_{1}\right) \in A_{1}\right\}}, \mathbb{I}_{\left\{\eta_{j}\left(\boldsymbol{x}_{2}\right) \in A_{2}\right\}}\right]=0\), for every \(A_{1}\), \(A_{2} \in \mathscr{B}(\mathbb{R})\) is equivalent to assuming that \(\eta_{j}\left(\boldsymbol{x}_{1}\right)\) and \(\eta_{j}\left(\boldsymbol{x}_{2}\right)\) are independent. This also applies for the process \(\boldsymbol{z}_{j}\). Notice also that an example of a process meeting the conditions of Theorem 8 is the Gaussian process with spherical covariance function (see Banerjee et al., 2003, Chapter 2). From Theorem 8 it is easy to see that if DMBPP or \(w \mathrm{DMBPP}\) are specified such that the marginal distribution of \(k\) is degenerate, then the correlation between the corresponding random measures goes to zero, since \(\lim _{l \rightarrow \infty} \operatorname{Cov}\left[F_{\boldsymbol{x}_{1 l}}\left(B_{\boldsymbol{y}}\right), F_{\boldsymbol{x}_{2 l}}\left(B_{\boldsymbol{y}}\right)\right]=0\). For \(\theta \mathrm{DMBPP}\) the correlation between the associated random measures when the predictor values are far apart reaches a different limit. In such case, it is difficult to establish conditions on the prior specification ensuring that the limit is zero.
Theorem 9. Let \(F\) be a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\). Assume that \(\mathscr{V}\) is a set of equicontinuous functions and that there exists a constant \(\gamma>0\), such that if \(\boldsymbol{x}_{1}, \boldsymbol{x}_{2} \in \mathscr{X}\) and \(\left\|\boldsymbol{x}_{1}-\boldsymbol{x}_{2}\right\|>\gamma\), then \(\operatorname{Cov}\left[\mathbb{I}_{\left\{\eta_{j}\left(\boldsymbol{x}_{1}\right) \in A_{1}\right\}}, \mathbb{I}_{\left\{\eta_{j}\left(\boldsymbol{x}_{2}\right) \in A_{2}\right\}}\right]=0\), for every \(A_{1}, A_{2} \in \mathscr{B}(\mathbb{R}), j \geq 1\). Assume also that for every \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) and for every sequence \(\left\{\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)\right\}_{l=1}^{\infty}\), with \(\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right) \in \mathscr{X}^{2}\), such that \(\lim _{l \rightarrow \infty}\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)= \left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right)\), we have \(\left(\eta_{j}\left(\boldsymbol{x}_{1 l}\right), \eta_{j}\left(\boldsymbol{x}_{2 l}\right)\right) \xrightarrow{\mathscr{L}}\left(\eta_{j}\left(\boldsymbol{x}_{1}\right), \eta_{j}\left(\boldsymbol{x}_{2}\right)\right), j \geq 1\), as \(l \rightarrow \infty\). Then, for every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\begin{aligned}
\lim _{l \rightarrow \infty} \operatorname{Cov}\left[F_{\boldsymbol{x}_{1 l}}\left(B_{\boldsymbol{y}}\right)\right. & \left., F_{\boldsymbol{x}_{2 l}}\left(B_{\boldsymbol{y}}\right)\right]=\sum_{k_{1}=1}^{\infty} \operatorname{Pr}\left\{k=k_{1}\right\} \sum_{\mathbf{j}_{1}, \mathbf{j}_{2} \in \mathcal{H}_{k_{1}, m}} \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid k_{1}+m-1, \boldsymbol{y}\right) \\
& \times \sum_{j=1}^{\infty} E\left[w_{j}\left(\boldsymbol{x}_{1}\right)\right] E\left[w_{j}\left(\boldsymbol{x}_{2}\right)\right] \operatorname{Cov}\left[\mathbb{I}_{\left\{\boldsymbol{\theta}_{j} \in A_{\mathbf{j}_{1}, k_{1}}\right\}}, \mathbb{I}_{\left\{\boldsymbol{\theta}_{j} \in A_{\mathbf{j}_{2}, k_{1}}\right\}}\right] \\
& +\operatorname{Cov}\left[E\left\{F_{\boldsymbol{x}_{1}}\left(B_{\boldsymbol{y}}\right) \mid k\right\}, E\left\{F_{\boldsymbol{x}_{2}}\left(B_{\boldsymbol{y}}\right) \mid k\right\}\right]
\end{aligned}
\]
with \(E\left\{F_{\boldsymbol{x}}\left(B_{\boldsymbol{y}}\right) \mid k\right\}, B_{\boldsymbol{y}}, A_{\mathbf{j}, k}\), and \(G_{0, \boldsymbol{x}}\) as defined in Theorem 8 and \(\bar{M}\left(\mathbf{j}, \mathbf{j}_{1} \mid\right. k+m-1, \boldsymbol{y})=\operatorname{Mult}(\mathbf{j} \mid k+m-1, \boldsymbol{y}) \times \operatorname{Mult}\left(\mathbf{j}_{1} \mid k+m-1, \boldsymbol{y}\right)\).

Finally, although the trajectories of the DMBPP and \(w \mathrm{DMBPP}\) have a.s. a limit only, the autocorrelation function of all versions of the model are continuous under mild conditions on the elements defining the processes.
Theorem 10. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}\right.\), \(\left.\boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) or a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\), defined such that \(\mathscr{V}\) and \(\mathscr{H}\) are sets of equicontinuous functions. Assume that for every \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) and for every sequence \(\left\{\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)\right\}_{l=1}^{\infty}\), with \(\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right) \in \mathscr{X}^{2}\), such that \(\lim _{l \rightarrow \infty}\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)= \left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right)\), we have that
\[
\left(\eta_{j}\left(\boldsymbol{x}_{1 l}\right), \eta_{j}\left(\boldsymbol{x}_{2 l}\right)\right) \xrightarrow{\mathscr{L}}\left(\eta_{j}\left(\boldsymbol{x}_{1}\right), \eta_{j}\left(\boldsymbol{x}_{2}\right)\right),
\]
and
\[
\left(\boldsymbol{z}_{j}\left(\boldsymbol{x}_{1 l}\right), \boldsymbol{z}_{j}\left(\boldsymbol{x}_{2 l}\right)\right) \xrightarrow{\mathscr{L}}\left(\boldsymbol{z}_{j}\left(\boldsymbol{x}_{1}\right), \boldsymbol{z}_{j}\left(\boldsymbol{x}_{2}\right)\right),
\]
as \(l \rightarrow \infty, j \geq 1\). Then, for every \(\boldsymbol{y} \in \Delta_{m}^{0}\),
\[
\lim _{l \rightarrow \infty} \rho\left[F_{\boldsymbol{x}_{1 l}}\left(B_{\boldsymbol{y}}\right), F_{\boldsymbol{x}_{2 l}}\left(B_{\boldsymbol{y}}\right)\right]=\rho\left[F_{\boldsymbol{x}_{1}}\left(B_{\boldsymbol{y}}\right), F_{\boldsymbol{x}_{2}}\left(B_{\boldsymbol{y}}\right)\right],
\]
where \(B_{\boldsymbol{y}}=\left[0, y_{1}\right] \times \ldots \times\left[0, y_{m}\right]\).

\subsection*{4.3. The asymptotic behavior of the posterior distribution}

We study the asymptotic behavior of the posterior distribution of the proposed model class in this section. Here we assume that we observe a random sample \(\left(\boldsymbol{y}_{i}, \boldsymbol{x}_{i}\right), i=1, \ldots, n\). Recall that, as is common in regression settings, we assume that the predictor vector \(\boldsymbol{x}_{i}\) contains only exogenous predictors. Notice that the exogeneity assumption allows us to focus on the conditional density estimation problem, regardless of the data generating mechanism of the predictors, that is, if they are randomly generated or fixed by design (see, e.g. Barndorff-Nielsen, 1973, 1978; Florens et al., 1990). Let \(Q\) be the true probability measure generating the predictors, with density w.r.t. a corresponding \(\sigma\)-additive measure denoted by \(q\). By the exogeneity assumption, the true probability model for the response variable and predictors takes the form \(h_{0}(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) q_{0}(\boldsymbol{y} \mid \boldsymbol{x})\), where both \(q\) and \(\left\{q_{0}(\cdot \mid \boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\) are in free variation, with \(q_{0}(\boldsymbol{y} \mid \boldsymbol{x})\) denoting a conditional density defined on \(\Delta_{m}\), and \(\boldsymbol{x} \in \mathscr{X}\).

Theorem 11. Let \(F\) be a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}\right.\), \(\left.\boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) or a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). If the assumptions of Theorem 3 are satisfied, then the posterior distribution associated with the random joint distribution induced by the corresponding DMBPP model, \(h(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) f_{\boldsymbol{x}}(\boldsymbol{y})\), where \(q\) is the density generating the predictors, is weakly consistent at any joint distribution of the form \(h_{0}(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) q_{0}(\boldsymbol{y} \mid \boldsymbol{x})\), where \(\left\{q_{0}(\cdot \mid \boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\).

Although Theorem 11 assumes that \(\boldsymbol{x}\) contains only continuous predictors, a similar result can be obtained when \(\boldsymbol{x}\) contains only predictors with finite support (e.g., categorical, ordinal and discrete predictors) or a combination of continuous predictors and predictors with finite support.

The following theorem states a stronger posterior consistency result when a specific probit stick-breaking process is assumed in the definition of the \(\theta\) DMBPP.

Theorem 12. Let \(F\) be a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\). If \(\mathscr{X}=[0,1]^{p}\) and the \(\theta \mathrm{DMBPP}\) is defined such that
(i) for every \(j \in \mathbb{N}, \eta_{j}\) is a Gaussian process with zero mean function and covariance kernel given by \(c_{j}\left(\boldsymbol{x}, \boldsymbol{x}^{\prime}\right)=\tau^{2} \exp \left\{-A_{j}\left\|\boldsymbol{x}-\boldsymbol{x}^{\prime}\right\|^{2}\right\}\), where \(\left(\boldsymbol{x}, \boldsymbol{x}^{\prime}\right) \in \mathscr{X}^{2}\) and \(A_{j}\) is a random variable, such that for some positive constants \(\kappa\) and \(\kappa_{0}\), and some sequence \(r_{n} \uparrow \infty\), such that \(r_{n}^{p} n^{\kappa}(\log n)^{p+1}=o(n)\),
\[
\operatorname{Pr}\left\{A_{j}>\delta_{n}\right\} \leq \exp \left\{-n^{-\kappa_{0}} j^{\left(\kappa_{0}+2\right) / \kappa} \log j\right\}
\]
and
\[
\operatorname{Pr}\left\{A_{n}>r_{n}\right\} \leq \exp \{-n\}
\]
where \(\delta_{n}=O\left((\log n)^{2} / n^{5 / 2}\right)\),
(ii) for every \(v_{\boldsymbol{x}} \in \mathscr{V}, v_{\boldsymbol{x}} \equiv \Phi\), where \(\Phi\) denotes the CDF of a standardnormal distribution.
(iii) \(G_{0}\) has full support on \(\Delta_{m}^{0}\), where \(G_{0}\) is the distribution of \(\boldsymbol{\theta}_{j}, j \geq 1\).
(iv) \(k\) has full support on \(\mathbb{N}\),
(v) there exists a sequence \(k_{n} \in \mathbb{N}\) such that \(\log \left(\frac{k_{n}\left(k_{n}+m\right)!}{k_{n}!(m+1)!}\right) \preceq O(n)\) and \(\operatorname{Pr}\left\{k>k_{n}\right\} \preceq O(\exp \{-n\})\), where \(\preceq\) stands for inequality up to a constant.

Then, the posterior distribution associated with the random joint distribution induced by the θDMBPP model, \(h(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) f_{\boldsymbol{x}}(\boldsymbol{y})\), where \(q\) is the density generating the predictors, is \(L_{1}\)-consistent at any joint distribution of the form \(h_{0}(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) q_{0}(\boldsymbol{y} \mid \boldsymbol{x})\), where \(\left\{q_{0}(\cdot \mid \boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\).

For an example of how to construct the sequence of random variables \(A_{j}\), see Remark 5.12 in Pati et al. (2013).

\section*{5. Computational aspects}

As can be noted from the definitions of the proposed models, predictors can be included in different manners. In what follows, we consider special definitions by exploiting the relation between Gaussian processes and Bayesian linear regression models. We also make use of spike-and-slab prior distributions on the regression coefficients that allow us for an automatic selection of the version of the model that best accommodates to the complexity of the underlying datagenerating mechanism. Note that such a prior avoids the need to fit each version of the model, as done by Barrientos et al. (2017).

We specify the predictor dependent weights and atoms of the dependent stick-breaking process in the DMBPP by means of transformations of a linear predictor. To define the weights of the DMBPP we consider \(v_{\boldsymbol{x}}(a)=e^{a} /(1+ \left.e^{a}\right), a \in \mathbb{R}\), and the stochastic process \(\eta_{j}(\boldsymbol{x})=\beta_{0 j}^{\eta}+\boldsymbol{x}^{t} \boldsymbol{\beta}_{j}^{\eta}\), where \(\beta_{0 j}^{\eta} \in \mathbb{R}\) and \(\boldsymbol{\beta}_{j}^{\eta} \in \mathbb{R}^{p}\) are independent and identically distributed for \(j \geq 1\), and \(\boldsymbol{x}= \left(x_{1}, \ldots, x_{p}\right) \in \mathscr{X}^{p}\) denotes the vector of covariates. Similarly, to define the atoms of the dependent stick-breaking process in the DMBPP we consider the transformation
\[
h_{\boldsymbol{x}}(\mathbf{b})=\left(e^{b_{1}}, \ldots, e^{b_{m}}\right) /\left(1+\sum_{l=1}^{m} e^{b_{l}}\right), \quad \mathbf{b} \in \mathbb{R}^{m}
\]
and the stochastic process \(\boldsymbol{z}_{j}(\boldsymbol{x})=\left(\boldsymbol{z}_{j 1}(\boldsymbol{x}), \ldots, \boldsymbol{z}_{j m}(\boldsymbol{x})\right)\), where \(\boldsymbol{z}_{j l}(\boldsymbol{x})=\beta_{0 j l}^{\boldsymbol{z}}+ \boldsymbol{x}^{t} \boldsymbol{\beta}_{j l}^{z}\) and \(\beta_{0 j l}^{\boldsymbol{z}} \in \mathbb{R}\) and \(\boldsymbol{\beta}_{j l}^{\boldsymbol{z}} \in \mathbb{R}^{p}\) are independent and identically distributed for \(j \geq 1, l=1, \ldots, m\).

In order to choose the version of the DMBPP model that best adapts to the data and following George \& McCulloch (1993), we consider a two-components mixture of normal distributions with different variances as a prior distribution on the coefficients of the linear predictor associated with the covariates, that is, on \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j l}^{\boldsymbol{z}}\). For the intercepts of the linear predictors we assume \(\beta_{0 j}^{\eta} \stackrel{i i d}{\sim} N\left(0, \sigma_{\eta}^{2}\right)\) and \(\beta_{0 j l}^{\boldsymbol{z}} \stackrel{i i d}{\sim} N\left(0, \sigma_{\boldsymbol{z}}^{2}\right)\). For \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j l}^{\boldsymbol{z}}\) we introduce latent binary variables \(\gamma^{\eta}\) and \(\gamma^{z}\) and assume
\[
\boldsymbol{\beta}_{j}^{\eta} \mid \gamma^{\eta} \stackrel{i i d}{\sim} N_{p}\left(\mathbf{0}, \boldsymbol{\Sigma}_{1}^{\eta}\right)^{1-\gamma^{\eta}} \times N_{p}\left(\mathbf{0}, \boldsymbol{\Sigma}_{2}^{\eta}\right)^{\gamma^{\eta}}, \quad j \geq 1
\]
\[
\boldsymbol{\beta}_{j l}^{z} \mid \gamma^{z} \stackrel{i i d}{\sim} N_{p}\left(\mathbf{0}, \boldsymbol{\Sigma}_{1}^{z}\right)^{1-\gamma^{z}} \times N_{p}\left(\mathbf{0}, \boldsymbol{\Sigma}_{2}^{z}\right)^{\gamma^{z}}, \quad j \geq 1, \quad l=1, \ldots, m,
\]
where \(N_{p}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) denotes the \(p\)-dimensional multivariate normal distribution with mean vector \(\boldsymbol{\mu} \in \mathbb{R}^{p}\) and \(p \times p\) positive definite covariance matrix \(\boldsymbol{\Sigma}\). The covariance matrices \(\boldsymbol{\Sigma}_{1}^{\eta}\) and \(\boldsymbol{\Sigma}_{1}^{\boldsymbol{z}}\) define the "spike" component of the prior and are set such that define normal distributions that are highly concentrated around zero, while \(\boldsymbol{\Sigma}_{2}^{\eta}\) and \(\boldsymbol{\Sigma}_{2}^{z}\) define the "slab" component of the prior and are set such that the resulting normal distributions are less concentrated around zero. Therefore, binary parameters \(\gamma^{\eta}\) and \(\gamma^{\boldsymbol{z}}\), which are common for every \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j}^{z}\), control the predictor dependency structure of the model.

When the vector of binary variables \(\left(\gamma^{\eta}, \gamma^{z}\right)\) is equal to \((1,1),(0,1),(1,0)\), or \((0,0)\), then the chosen model is fully dependent, single-weight, single-atom, or predictor independent, respectively. To complete the prior for \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j l}^{\boldsymbol{z}}\), we consider
\[
\left(\gamma^{\eta}, \gamma^{z}\right) \sim \pi_{1} \delta_{(1,1)}+\pi_{2} \delta_{(0,1)}+\pi_{3} \delta_{(1,0)}+\pi_{4} \delta_{(0,0)}
\]
where \(\pi_{i} \geq 0\) and \(\pi_{1}+\pi_{2}+\pi_{3}+\pi_{4}=1\). Finally, to complete the prior specification for the DMBPP model we assume \(k \mid \lambda \sim \operatorname{Poisson}(\lambda) \mathbb{I}_{\{k \geq 1\}}\).

Posterior sampling of the DMBPP model can be based on any conditional algorithm designed for BNP mixture models. The specific implementation employed in the simulation study and in the application was based on a finite representation of the dependent stick-breaking process to a level \(N\) (Ishwaran \& James, 2001). We use Gibbs sampling algorithms to generate samples from the posterior distribution. To sample the non conjugate full conditional distributions of the coefficients in the linear predictors we use the slice sampler algorithm (Neal, 2003). We use a Metropolis-Hastings step (Tierney, 1994) to update the degree of the polynomial. The binary parameters are sampled from their conjugate categorical posterior distribution. More details are provided in Appendix E. The code employed here to fit the proposed model is available in GitHub.

\section*{6. Illustrations}

In this section we illustrate the performance of the model in a simulation study and in an application to solid waste recycling in the city of Santiago de Cali, Colombia. In the simulation study, we show the ability of the model to estimate the true conditional densities as well as its capacity to choose the version of the DMBPP model (fully-dependent, single-atoms, single-weights, or independent) that best accommodates to the complexity of the true data-generating mechanism. In the application, we compare the performance of our proposed model with the performance of a parametric Dirichlet regression model on a transformed version of the data.

\subsection*{6.1. Simulation study}

We consider four simulation scenarios representing varying degrees of complexity and shapes as the predictor varies, based on mixtures of predictor-dependent Dirichlet densities that are not particular cases of the implemented model. In all cases, the predictor is univariate and uniformly distributed on the \((0,1)\) interval. For Scenario I, both the weights and the parameters of the Dirichlet distributions depend on the predictor. Under this scenario, for small values of \(x\) the conditional density has one mode which splits into two and later merges into one again as the value of the predictor increases. For Scenario II, only the parameters of the Dirichlet densities depend on the predictor. Under this scenario, for small values of \(x\) the conditional density has three well separated modes, one at each corner of the simplex, which merge into two and later into only one irregularly shaped as the value of \(x\) increases. For Scenario III only weights depend on the predictor. Under this scenario, for small values of \(x\) the conditional density has only one mode which splits into two and later merges into one mode centered roughly in the middle of the simplex as the value of \(x\) increases. Finally, Scenario IV is given by a predictor-independent Dirichlet density. The specification of true conditional densities for Scenarios I - IV is given in Table 1.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1
Simulation Study: true conditional density functions considered in the simulation study.
Here \(w_{1}(x)=\frac{x}{4-3 x}, \boldsymbol{\theta}_{1}(x)=(25-20 x, 5+25 x, 3), \boldsymbol{\theta}_{2}(x)=(5,5+15 x, 30-17 x)\), \(\boldsymbol{\theta}_{3}(x)=(5+9 x, 30+9 x, 3+9 x)\), and \(x \in(0,1)\).}
\begin{tabular}{cl}
\hline \hline Scenario & \(f_{0}(\boldsymbol{y} \mid x)\) \\
\hline I & \(w_{1}(x) \operatorname{dir}\left(\boldsymbol{y} \mid \boldsymbol{\theta}_{1}(x)\right)+\left(1-w_{1}(x)\right) \operatorname{dir}\left(\boldsymbol{y} \mid \boldsymbol{\theta}_{2}(x)\right)\) \\
II & \(0.6 \operatorname{dir}\left(\boldsymbol{y} \mid \boldsymbol{\theta}_{1}(x)\right)+0.2 \operatorname{dir}\left(\boldsymbol{y} \mid \boldsymbol{\theta}_{2}(x)\right)+0.2 \operatorname{dir}\left(\boldsymbol{y} \mid \boldsymbol{\theta}_{3}(x)\right)\) \\
III & \(w_{1}(x) \operatorname{dir}(\boldsymbol{y} \mid(10,12,12))+\left(1-w_{1}(x)\right) \operatorname{dir}(\boldsymbol{y} \mid(24,6,6))\) \\
IV & \(\operatorname{dir}(\boldsymbol{y} \mid(35,25,40))\) \\
\hline \hline
\end{tabular}
\end{table}

We consider three sample sizes, \(n=250, n=500\), and \(n=1,000\) for each scenario. A Monte Carlo sample size of 100 was considered for each scenario and sample size. Following Zellner (1983), we consider \(\boldsymbol{\Sigma}_{l}^{\eta}=\tau_{l}^{\eta}\left(\mathbb{X}^{t} \mathbb{X}\right)^{-1}\) and \(\boldsymbol{\Sigma}_{l}^{z}= \tau_{l}^{z}\left(\mathbb{X}^{t} \mathbb{X}\right)^{-1}\), for \(l=1,2\), where \(\mathbb{X}\) denotes the design matrix without including the intercept, \(\tau_{1}^{\eta}\) and \(\tau_{1}^{z}\) are small positive values, while \(\tau_{2}^{\eta}\) and \(\tau_{2}^{z}\) are large positive values. In Appendix F we provide the justification for the particular choices. Finally, we set \(\sigma_{\eta}^{2}=\sigma_{\boldsymbol{z}}^{2}=100\). For the prior of the binary latent variables, \(\left(\gamma^{\eta}, \gamma^{z}\right)\), we set \(\pi_{1}=1 / t^{2}, \pi_{2}=\pi_{3}=(t-1) / 2 t^{2}\), and \(\pi_{4}=(t-1) / t\), with \(t>1\). A priori, larger values of \(t\) favor more parsimonious models. Two prior specifications were employed by setting \(t=2\) (Prior I) and \(t=10\) (Prior II). Under Prior I, the prior probability of the covariate independent model is \(\pi_{4}=0.50\), followed by the prior probability of the fully covariate dependent model, which is \(\pi_{1}=0.25\). Prior II strongly favors parsimonious models. Under this specification, the prior probability of the covariate independent model is \(\pi_{4}=0.90\), while the prior probability of its fully covariate dependent counterpart
is only \(\pi_{1}=0.01\). Finally, to complete the prior specification we consider \(\lambda=25\).
A single Markov chain was generated for each simulated data set. For \(n=250\) and \(n=500\) a chain of length 110,000 was generated and the posterior inference was based on a reduced chain of 10,000 samples obtained after a burn-in period of 10,000 and keeping 1 every 10 samples. A similar specification was considered for \(n=1,000\), but considered a burn-in period of 50,000 samples in such cases. To assess the performance of the proposed model in estimating the true data generating mechanism, we compute an estimate of the integrated- \(L_{1}\) and \(L_{\infty}\) distances, denoted by \(\widehat{I L}_{1}\) and \(\widehat{L}_{\infty}\), respectively. Specifically, we compute
\[
\begin{aligned}
& \widehat{I L}_{1}=\frac{1}{L} \frac{1}{M} \sum_{j=1}^{L} \sum_{i=1}^{M}\left|\hat{f}\left(\boldsymbol{y}_{i} \mid \boldsymbol{x}_{j}\right)-f_{0}\left(\boldsymbol{y}_{i} \mid \boldsymbol{x}_{j}\right)\right| \\
& \widehat{L}_{\infty}=\max _{i} \max _{j}\left|\hat{f}\left(\boldsymbol{y}_{i} \mid \boldsymbol{x}_{j}\right)-f_{0}\left(\boldsymbol{y}_{i} \mid \boldsymbol{x}_{j}\right)\right|
\end{aligned}
\]
where \(\widehat{f}(\cdot \mid x)\) denotes the posterior mean of the conditional density, \(f_{0}(\cdot \mid x)\) denotes the true conditional density, and \(\left\{\boldsymbol{y}_{i}\right\}_{i=1}^{M}\) and \(\left\{\boldsymbol{x}_{j}\right\}_{j=1}^{L}\) define an equally spaced grid of \(\Delta_{m}\) and \(\mathscr{X}\), respectively.

To assess the model's ability to choose the version that best accommodates to the complexity of the underlying true data-generating distribution, we select the combination of \(\left(\gamma^{\eta}, \gamma^{z}\right)\) that concentrates the highest posterior probability and compare it to the true predictor dependency structure of the simulation scenario. Recall that \(\left(\gamma^{\eta}, \gamma^{z}\right)\) control which part of the model depends on the predictor and that each of the simulation scenarios depend on the predictor in different ways. Scenario I involves predictors in weights and in Dirichlet densities, Scenario II only in Dirichlet densities, Scenario III only in weights, and Scenario IV does not depend on predictors at all.

Table 2 shows the mean, across replicates, of the integrated- \(L_{1}\) distance between the true and the posterior mean for each simulation scenario, sample size, and spike-and-slab prior. As expected, the integrated- \(L_{1}\) distance decreases as the sample size increases for each simulation scenario under both spike-andslab priors. For small samples sizes ( \(n=250,500\) ), the smallest integrated- \(L_{1}\) distances are observed for Scenario III, the single-atoms true model, while for \(n=1000\), the smallest integrated- \(L_{1}\) distance is observed for Scenario IV, the predictor independent true model. The largest integrated- \(L_{1}\) distances for small sample sizes are observed for Scenario IV, while for \(n=1000\) the largest distance is observed for Scenario II, the single-weights true model. The model seems to be robust regarding the choice of the spike-and-slab prior. Similar results are obtained when the \(L_{\infty}\) distance is considered, which are shown in Appendix G.

Table 3 shows the proportion of times, across Monte Carlo replicates, that the selected predictor dependency structure of the fit model agrees with the one of the true model. The proportion increases as the sample size increases for each simulation scenario and spike-and-slab prior specification. Remarkably, for Scenarios I and IV, the DMBPP model is able to choose the version of the model that is in agreement with the predictor dependency structure of the true model

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2
Mean, across Monte Carlo replicates, of the integrated \(L_{1}\) distance between the truth and the posterior mean of the conditional densities for each simulation scenario, spike-and-slab prior (Prior I and Prior II), and sample size ( \(n\) ).}
\begin{tabular}{cccccccc}
\hline \hline & \multicolumn{3}{c}{ Prior I } & & \multicolumn{3}{c}{ Prior II } \\
\cline { 2 - 4 } \cline { 6 - 8 } Scenario & \(n=250\) & \(n=500\) & \(n=1,000\) & & \(n=250\) & \(n=500\) & \(n=1,000\) \\
\hline I & 0.413 & 0.345 & 0.326 & & 0.414 & 0.349 & 0.322 \\
II & 0.479 & 0.426 & 0.411 & & 0.482 & 0.434 & 0.413 \\
III & 0.411 & 0.301 & 0.271 & & 0.410 & 0.301 & 0.267 \\
IV & 0.599 & 0.380 & 0.234 & & 0.599 & 0.380 & 0.230 \\
\hline \hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3
Proportion of times, across Monte Carlo replicates, in which the true predictor dependency structure is selected for each simulation scenario, spike-and-slab prior (Prior I and Prior II), and sample size ( \(n\) ).}
\begin{tabular}{cccccccc}
\hline \hline & \multicolumn{3}{c}{ Prior I } & & \multicolumn{3}{c}{ Prior II } \\
\cline { 2 - 4 } \cline { 6 - 8 } Scenario & \(n=250\) & \(n=500\) & \(n=1000\) & & \(n=250\) & \(n=500\) & \(n=1000\) \\
\hline I & 1.000 & 1.000 & 1.000 & & 1.000 & 1.000 & 1.000 \\
II & 0.440 & 0.670 & 0.870 & & 0.990 & 1.000 & 0.980 \\
III & 0.960 & 0.960 & 0.980 & & 0.960 & 0.970 & 0.980 \\
IV & 1.000 & 1.000 & 1.000 & & 1.000 & 1.000 & 1.000 \\
\hline \hline
\end{tabular}
\end{table}
for every replicated data set and sample size. For Scenario III, the proportion of times that the chosen version of the model agrees with the true model increases from 0.96 to 0.98 as the sample size increases from 250 to 1000 . The true model for which it is most difficult to choose the version of the model that agrees with the predictor dependency structure of the true model, is Scenario II (singleweights model) and under Prior I. Interestingly, it seems that the ability of the model to choose the version of the model that best fits the data is not completely related to the capacity of the model to estimate the conditional densities. For example, for sample sizes 250 and 500 , the smallest integrated \(L_{1}\) mean distances are observed for Scenario III, while the binary latent variable estimates agree with the predictor dependency structure of the true model the most for Scenarios I and IV. Again, the results are robust regarding the model selection prior distribution.

Figures 1 to 4 display the contour plot of the mean, across simulation, of the posterior mean of the conditional density for selected values of the predictor, each sample size and simulation scenario, under Prior I for \(\left(\gamma^{\eta}, \gamma^{z}\right)\). The results are consistent with the previous discussion. For the different values of the predictor, the figures show how close the estimates are to the true model and how they improve as the sample size increases.

\subsection*{6.2. Application to solid waste in Colombia}

In this section, we analyze data about solid waste generated in a residential area in the city of Santiago de Cali, Colombia. The data set was collected to

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-19.jpg?height=1649&width=1095&top_left_y=418&top_left_x=603}
\captionsetup{labelformat=empty}
\caption{Fig 1. Simulation study - Scenario I: contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for ( \(\gamma^{\eta}, \gamma^{z}\) ). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-20.jpg?height=1643&width=1094&top_left_y=419&top_left_x=416}
\captionsetup{labelformat=empty}
\caption{Fig 2. Simulation study - Scenario II: contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for ( \(\gamma^{\eta}, \gamma^{z}\) ). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-21.jpg?height=1641&width=1090&top_left_y=418&top_left_x=603}
\captionsetup{labelformat=empty}
\caption{Fig 3. Simulation study - Scenario III: contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for ( \(\gamma^{\eta}, \gamma^{z}\) ). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-22.jpg?height=1640&width=1101&top_left_y=421&top_left_x=416}
\captionsetup{labelformat=empty}
\caption{Fig 4. Simulation study - Scenario IV: contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for ( \(\gamma^{\eta}, \gamma^{z}\) ). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}
estimate the per capita daily production and characterization of solid waste in the city. The data set contains information about 261 block sides, for which solid waste was separated in different kinds of materials, including food, hygienic, and others. Finally, the proportion of these materials were registered for each block side. Additionally, the socioeconomic level of the residents in the area was recorded including the categories "low-low", "low", "medium-low", "medium", "medium-high", and "high". See Klinger et al. (2009) for more details regarding this data set.

In this analysis, the proportion of food, hygienic solid, and other type of waste were considered as the response vector on the 2 -dimensional simplex. The socioeconomic level was considered as a categorical covariate with a dummy variable representation leading to \(p=6\) predictors. Here, we consider a model specification similar to the one detailed in Section 5. Here we set \(\sigma_{\eta}^{2}=\sigma_{\boldsymbol{z}}^{2}=100\) and \(\lambda=25\). We refer the reader to Appendix H for a description regarding the selection of \(\tau_{l}^{\eta}\) and \(\tau_{l}^{z}\).

A single Markov chain with 300,000 samples was generated. Posterior inference was based on a reduced chain with 10,000 samples obtained after a 100,000 burn-in period and keeping 1 every 20 samples. To assess the performance of the proposed model, we also fit a parametric Dirichlet regression (PDR) model to the data. Due to presence of zero-coordinate vectors in the data, we transform the response vectors, \(\boldsymbol{y}\), using the transformation proposed by Smithson \& Verkuilen (2006), given by \(\boldsymbol{y}^{*}=[\boldsymbol{y}(n-1)+1 /(m+1)] / n\), where \(n\) is the size of the sample and \(m\) is the dimension of the simplex. Thus the parametric model was applied to the transformed responses \(\boldsymbol{y}^{*}\), such that
\[
\boldsymbol{y}_{i}^{*} \mid \boldsymbol{x}_{i}, \boldsymbol{\gamma} \sim \operatorname{dir}\left(\boldsymbol{\gamma}\left(\boldsymbol{x}_{i}\right)\right),
\]
where \(\boldsymbol{\gamma}\left(\boldsymbol{x}_{i}\right)=\left(\gamma_{1}\left(\boldsymbol{x}_{i}\right), \ldots, \gamma_{m}\left(\boldsymbol{x}_{i}\right)\right)\), with \(\log \left(\gamma_{l}(\boldsymbol{x})\right)=\boldsymbol{x}^{t} \boldsymbol{\beta}_{l}, l=1, \ldots, m\). We complete the model specification by assuming \(\boldsymbol{\beta}_{l} \sim N_{p+1}(\boldsymbol{m}, \boldsymbol{\Sigma})\), with \(\boldsymbol{m}=\mathbf{0}\), \(\boldsymbol{\Sigma}=\sigma^{2} \boldsymbol{I}_{p+1}, \sigma^{2}=100\), and \(\boldsymbol{I}_{p}\) being a \(p \times p\) identity matrix. The models were compared by means of their posterior predictive abilities, quantified by the log pseudo marginal likelihood (LPML) and the widely applicable information criterion (WAIC) (Watanabe \& Opper, 2010). The LPML, developed by Geisser \(\& \operatorname{Eddy}(1979)\), is given by \(\sum_{i=1}^{n} \log p_{M}\left(\boldsymbol{y}_{i} \mid \boldsymbol{Y}_{-i}\right)\), where \(p_{M}\left(\boldsymbol{y}_{i} \mid \boldsymbol{Y}_{-i}\right)\) is the posterior predictive distribution for observation \(\boldsymbol{y}_{i}\), based on the data \(\boldsymbol{Y}_{-i}\), under model \(M\), with \(\boldsymbol{Y}_{-i}\) denoting the observed data matrix after removing the \(i\) th observation. The \(p_{M}\left(\boldsymbol{y}_{i} \mid \boldsymbol{Y}_{-i}\right)\) is also known as the conditional predictive ordinate of observation \(i\) under model \(M\) and the method of Gelfand \& Dey (1994) was used in its computation. The WAIC is given by
\[
\mathrm{WAIC}=-\frac{1}{n} \sum_{i=1}^{n} \log E_{\text {post }}\left[p_{M}\left(\boldsymbol{y}_{i} \mid \boldsymbol{\theta}\right)\right]+\frac{1}{n} \sum_{i=1}^{n} \operatorname{Var}_{\text {post }}\left[\log p_{M}\left(\boldsymbol{y}_{i} \mid \boldsymbol{\theta}\right)\right]
\]
where \(p_{M}\left(\boldsymbol{y}_{i} \mid \boldsymbol{\theta}\right)\) is the density function for observation \(\boldsymbol{y}_{i}\), given parameter \(\boldsymbol{\theta}\), under model \(M\), and \(E_{\text {post }}\) and \(\operatorname{Var}_{\text {post }}\) denote the posterior mean and posterior variance, respectively. The second term on the right hand side of expression (12)

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-24.jpg?height=1063&width=1137&top_left_y=386&top_left_x=397}
\captionsetup{labelformat=empty}
\caption{Fig 5. Solid waste data: Posterior mean of the conditional density. The results for the for the DMBPP model under spike-and-slab Prior I are presented in the first and third columns. The results for the PDR model are presented in the second and fourth columns. Panels (a) and (b), (c) and (d), (e) and (f), (i) and (j), and (k) and (l) present the results for socioeconomic level low-low, low, medium-low, medium, medium-high, and high, respectively. The \(x\)-axis and \(y\)-axis denote the proportion of food and hygienic waste, respectively.}
\end{figure}
is a penalty for overfitting. In what follows, we compute WAIC as described by Gelman et al. (2013), page 173, and report \(-n W A I C\). Models with greater values of LPML and \(-n W A I C\) are to be preferred. For the parametric model, we consider the transformed data \(\boldsymbol{y}_{i}=\boldsymbol{y}_{i}^{*}\) in the computation of the LPML and WAIC criteria.

Figure 5 displays the conditional density estimates for the DMBPP model and the PDR model, for each socioeconomic level. The results for the DMBPP model re shown under spike-and-slab Prior I. They suggest that different socioeconomic levels show different recycling behaviors and the advantages of the nonparametric model are evident. The DMBPP model shows to be more flexible in estimating the conditional densities than the parametric model for varying values of the predictor, specially when the socioeconomic levels are "medium" and "high". We highlight that the parametric fit was only possible due to a pre-transformation of the data.

The LPML for the DMBPP and PDR model was 778 and 649, respectively. The \(-n W A I C\) the DMBPP and PDR model was 778 and 649 was 778 and 650, respectively. These goodness-of-fit criteria support and agree in that the DMBPP model provides a better fit for this data set than the PDR model. For DMBPP model, the posterior probability of \(\left(\gamma^{\eta}, \gamma^{z}\right)=(0,0)\) was approximately equal to zero, which is a formal test for the hypothesis that the densities for solid waste are the same across the socioeconomic levels. A probability close to zero is interpreted as little evidence in favor of this hypothesis. Additional results for the DMBPP model under spike-and-slap Prior II can be found in Appendix I, which are robust regarding the prior specification.

\section*{7. Discussion}

We have proposed a novel and general class of probability models for sets of predictor-dependent probability distributions supported on simplex spaces. The proposal corresponds to an extension of dependent univariate Bernstein polynomial processes proposed by Barrientos et al. (2017) and is based on the modified class of MBP proposed by Barrientos et al. (2015).

The proposed model class has appealing theoretical properties such as full support, well behaved correlation function, and consistent posterior distribution. The incorporation of a spike-and-slab prior for the predictor dependent stochastic processes involved in the model adapts well to the complexity of the underlying data-generating distribution. The approach also allows the user to formally test whether all predictors are simultaneously related to the compositional response. The study of the theoretical properties of the model selection component of our approach is the subject of ongoing research.

\section*{Appendix A: Dealing with compositional data and zero-valued entries}

In this appendix, we explain how the proposed approach is able to handle compositional observations with zero-valued entries. To simplify the discussion, we limit ourselves to the scenario where no predictors are present, but the arguments extend naturally to the predictor-dependent case. In some circumstances, such as in the Colombian solid waste application, some compositional observations may have entries equal to zero. As a result, analysts must determine whether such zeros are systematic and, if so, whether they should be handled using a zero-inflated modeling approach. Analysts might use standard modeling techniques if there is no indication that the statistical model must be degenerated to account for zero-valued entries (as we assume throughout this paper). Analysts should be aware that using standard modeling techniques can cause potential problems. In what follows, we illustrate some of those potential problems in a specific, yet common, scenario.

Suppose we will model compositional data with zero-valued entries using either a single Dirichlet distribution or a standard mixture of Dirichlet densities.

Assume that, out of a sample of \(n\) observations \(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\), only the first entry of \(\boldsymbol{y}_{1}\) is exactly equal to zero, i.e., \(y_{1,1}=0\). For the remaining entries of \(\boldsymbol{y}_{1}\) and all entries of \(\boldsymbol{y}_{2}, \ldots, \boldsymbol{y}_{n}\), assume they take values on ( 0,1 ). If we use the Dirichlet distribution to model \(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\) directly, we need first to notice that the likelihood is equal to
\[
y_{1,1}^{\gamma_{1}-1} \times \mathcal{L}_{-y_{1,1}}\left(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) ;\left(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\right)\right)
\]
where
\[
\begin{aligned}
& \mathcal{L}_{-y_{1,1}}\left(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) ;\left(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\right)\right) \\
& \quad=\left\{\frac{\Gamma\left(\gamma_{1}+\ldots+\gamma_{m+1}\right)}{\Gamma\left(\gamma_{1}\right) \times \ldots \times \Gamma\left(\gamma_{m+1}\right)} y_{1,2}^{\gamma_{2}-1} \times \ldots \times y_{1, m}^{\gamma_{m}-1}\left(1-\sum_{j=1}^{m} y_{i, j}\right)^{\gamma_{m+1}-1}\right\} \\
& \quad \times \prod_{i=2}^{n} \operatorname{dir}\left(\boldsymbol{y}_{i} \mid\left(\gamma_{1}, \ldots, \gamma_{m+1}\right)\right)
\end{aligned}
\]
and
\[
\begin{aligned}
\operatorname{dir}\left(\boldsymbol{y}_{i} \mid\right. & \left.\left(\gamma_{1}, \ldots, \gamma_{m+1}\right)\right)= \\
& \frac{\Gamma\left(\gamma_{1}+\ldots+\gamma_{m+1}\right)}{\Gamma\left(\gamma_{1}\right) \times \ldots \times \Gamma\left(\gamma_{m+1}\right)} y_{i, 1}^{\gamma_{1}-1} \times \ldots \times y_{i, m}^{\gamma_{m}-1}\left(1-\sum_{j=1}^{m} y_{i, j}\right)^{\gamma_{m+1}-1}
\end{aligned}
\]

Since we are assuming that \(y_{1,1}=0\), we have to consider the following three cases:
- If \(\gamma_{1}<1\), then \(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) \mapsto y_{1,1}^{\gamma_{1}-1} \mathcal{L}_{-y_{1,1}}\left(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) ;\left(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\right)\right) \equiv \infty\).
- If \(\gamma_{1}>1\), then \(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) \mapsto y_{1,1}^{\gamma_{1}-1} \mathcal{L}_{-y_{1,1}}\left(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) ;\left(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\right)\right) \equiv 0\).
- If \(\gamma_{1}=1\), then \(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) \mapsto y_{1,1}^{\gamma_{1}-1} \mathcal{L}_{-y_{1,1}}\left(\left(\gamma_{1}, \ldots, \gamma_{m+1}\right) ;\left(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\right)\right)\) is not constant and will take values on ( \(0, \infty\) ). For this case, we adopt the convention that \(0^{0}=1\).
Notice that, only in the case \(\gamma_{1}=1\), we will be able to make inferences about
\[
\left(\gamma_{2}, \ldots, \gamma_{m+1}\right)
\]

Since the constraint \(\gamma_{1}=1\) is entirely driven by the fact that the first entry of \(\boldsymbol{y}_{1}\) is equal to zero, adopting this modeling approach is inadequate.

We might also consider using a mixture model for it is a more general and flexible approach. For example, consider a standard DP mixture of Dirichlet densities of the form
\[
\sum_{j=1}^{\infty} w_{j} \operatorname{dir}\left(\cdot \mid \gamma_{j}=\left(\gamma_{1, j}, \ldots, \gamma_{m+1, j}\right)\right)
\]

Using this mixture model will lead to issues similar to those described above if the prior distribution for the atoms, \(\gamma_{j}\), is absolutely continuous with respect
to Lebesgue. We could overcome this issue if the prior distribution for \(\gamma_{j}\) is such that \(P\left\{\gamma_{1, j}=1\right\}>0\), which will allow making inferences while remaining flexible. For a more general scenario where we assume that zeros might occur in multiple entries and observations, we will be able to make inferences as long as the assumption \(P\left\{\gamma_{l, j}=1: j \in \mathcal{J}\right\}>0\) is met for each \(\mathcal{J} \subseteq\{1, \ldots, m\}\). This particular assumption is naturally met (under mild conditions) when using the class of mixture of Dirichlet densities derived from the MBP described in Section 2, that is \(\sum_{j=1}^{\infty} w_{j} \operatorname{dir}\left(\cdot \mid \alpha\left(k,\left\lceil k \boldsymbol{\theta}_{j}\right\rceil\right)\right)\). The assumption is met when the prior for \(\boldsymbol{\theta}_{j}\) has full support on \(\Delta_{m}^{0}\) or, equivalently, the prior for \(\boldsymbol{\theta}_{j}\) has positive density (with respect to Lebesgue) on \(\Delta_{m}^{0}\).

\section*{Appendix B: Formal definition of special cases of the general model}

Definition 2. Let \(\mathscr{V}\) and \(\mathscr{H}\) be two sets of functions as defined before. Let \(F=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) be a \(\mathscr{P}\left(\Delta_{m}\right)\)-valued stochastic process such that:
(i) \(v_{1}, v_{2}, \ldots\), are independent \([0,1]\)-valued random variables with common distribution indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{v}\).
(ii) \(\boldsymbol{z}_{j}=\left\{\boldsymbol{z}_{j}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}, j \geq 1\), are independent and identically distributed real-valued stochastic processes with law indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{z}\).
(iii) \(k \in \mathbb{N}\) is a discrete random variable with distribution indexed by a finitedimensional parameter \(\lambda\).
(iv) For every \(\boldsymbol{x} \in \mathscr{X}\), the density function of \(F_{\boldsymbol{x}}\), w.r.t. Lebesgue measure, is given by the following dependent mixture of Dirichlet densities,
\[
f_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{\infty} w_{j} \operatorname{dir}\left(\cdot \mid \alpha\left(k,\left\lceil k \boldsymbol{\theta}_{j}(\boldsymbol{x})\right\rceil\right)\right),
\]
where \(\alpha(k, \boldsymbol{j})=\left(\boldsymbol{j}, k+m-\|\mathbf{j}\|_{1}\right), \boldsymbol{\theta}_{j}(\boldsymbol{x})\) and \(\left\lceil k \boldsymbol{\theta}_{j}(\boldsymbol{x})\right\rceil\) are defined as in Definition 1, and \(w_{j}=v_{j} \prod_{l<j}\left[1-v_{l}\right]\).
The process \(F=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) will be referred to as single-weight dependent MBP process with parameters \(\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\), and denoted by \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}\right.\), \(\boldsymbol{\Psi}_{z}, \mathscr{H}\) ) and \(w \mathrm{DMBPP}\) for short.

Definition 3. Let \(\mathscr{V}\) and \(\mathscr{H}\) be two sets of functions as defined before. Let \(F=\{F(\boldsymbol{x}, \cdot): \boldsymbol{x} \in \mathscr{X}\}\) be a \(\mathscr{P}\left(\Delta_{m}\right)\)-valued stochastic process such that:
(i) \(\eta_{j}=\left\{\eta_{j}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}, j \geq 1\), are independent and identically distributed real-valued stochastic processes with law indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{\eta}\).
(ii) \(\boldsymbol{\theta}_{1}, \boldsymbol{\theta}_{2}, \ldots\), are independent \(\Delta_{m}^{0}\)-valued random vectors with common distribution indexed by a finite-dimensional parameter \(\boldsymbol{\Psi}_{\boldsymbol{\theta}}\).
(iii) \(k \in \mathbb{N}\) is a discrete random variable with distribution indexed by a finitedimensional parameter \(\lambda\).
(iv) For every \(\boldsymbol{x} \in \mathscr{X}\), the density function of \(F_{\boldsymbol{x}}\), w.r.t. Lebesgue measure, is given by the following dependent mixture of Dirichlet densities,
\[
f_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}) \operatorname{dir}\left(\cdot \mid \alpha\left(k,\left\lceil k \boldsymbol{\theta}_{j}\right\rceil\right)\right),
\]
where \(w_{j}(\boldsymbol{x})\) are defined as in Definition 1 and
\[
\left\lceil k \boldsymbol{\theta}_{j}\right\rceil=\left(\left\lceil k \theta_{j 1}\right\rceil, \ldots,\left\lceil k \theta_{j m}\right\rceil\right)
\]

The process \(F=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) will be referred to as single-atoms dependent \(M B P\) process with parameters \(\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\), and denoted by \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}\right.\), \(\mathscr{V}, \mathbf{\Psi}_{\boldsymbol{\theta}}\) ) and \(\theta\) DMBPP for short.

\section*{Appendix C: Topological bases and sub-bases}

A sub-base for the weak product topology for the space
\[
\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}=\prod_{\boldsymbol{x} \in \mathscr{X}} \mathscr{P}\left(\Delta_{m}\right)
\]
is given by sets of the form \(B_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\right)=\prod_{\boldsymbol{x} \in \mathscr{X}} \Delta_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(Q_{\boldsymbol{x}}\right)\), where
\[
\Delta_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(Q_{\boldsymbol{x}}\right)= \begin{cases}\mathscr{P}\left(\Delta_{m}\right), & \text { if } \boldsymbol{x} \neq \boldsymbol{x}_{0} \\ \left\{M_{\boldsymbol{x}} \in \mathscr{P}\left(\Delta_{m}\right):\left|\int_{\Delta_{m}} f d M_{\boldsymbol{x}}-\int_{\Delta_{m}} f d Q_{\boldsymbol{x}}\right|<\epsilon\right\}, & \text { if } \boldsymbol{x}=\boldsymbol{x}_{0}\end{cases}
\]
for every \(f: \Delta_{m} \longrightarrow \mathbb{R}\) bounded continuous function, \(\epsilon>0, \boldsymbol{x}_{0} \in \mathscr{X}\) and \(Q_{\boldsymbol{x}} \in \mathscr{P}\left(\Delta_{m}\right)\).

Let \(\mathscr{D}\left(\Delta_{m}\right) \subset \mathscr{P}\left(\Delta_{m}\right)\) be the space of all probability measures defined on \(\Delta_{m}\) that are absolutely continuous w.r.t. Lebesgue measure and with continuous density function. A sub-base for the \(L_{\infty}\) product topology for the space
\[
\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}=\prod_{\boldsymbol{x} \in \mathscr{X}} \mathscr{D}\left(\Delta_{m}\right)
\]
is given by sets of the form \(B_{\epsilon, \boldsymbol{x}_{0}}^{L_{\infty}}\left(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\right)=\prod_{\boldsymbol{x} \in \mathscr{X}} \Delta_{\epsilon, \boldsymbol{x}_{0}}^{L_{\infty}}\left(Q_{\boldsymbol{x}}\right)\), where
\[
\Delta_{\epsilon, \boldsymbol{x}_{0}}^{L_{\infty}}\left(Q_{\boldsymbol{x}}\right)= \begin{cases}\mathscr{D}\left(\Delta_{m}\right), & \text { if } \boldsymbol{x} \neq \boldsymbol{x}_{0} \\ \left\{M_{\boldsymbol{x}} \in \mathscr{D}\left(\Delta_{m}\right): \sup _{\boldsymbol{y} \in \Delta_{m}}\left|m_{\boldsymbol{x}}(\boldsymbol{y})-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\epsilon\right\}, & \text { if } \boldsymbol{x}=\boldsymbol{x}_{0}\end{cases}
\]
for every \(\epsilon>0, \boldsymbol{x}_{0} \in \mathscr{X}\) and \(Q_{\boldsymbol{x}} \in \mathscr{D}\left(\Delta_{m}\right)\), where \(m_{\boldsymbol{x}}\) and \(q_{\boldsymbol{x}}\) denote the density function of \(M_{\boldsymbol{x}}\) and \(Q_{\boldsymbol{x}}\), respectively.

Now, assume that the predictor vector \(\boldsymbol{x}\) contains only continuous predictors and that the predictor space \(\mathscr{X}\) is compact. A base for the \(L_{\infty}\) topology for the space \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}=\prod_{\boldsymbol{x} \in \mathscr{X}} \mathscr{D}\left(\Delta_{m}\right)\), is given by sets of the form
\[
\begin{aligned}
B_{\epsilon}^{L_{\infty}}\left(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\right)= & \\
& \left\{\left\{M_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}: \sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|m_{\boldsymbol{x}}(\boldsymbol{y})-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\epsilon\right\},
\end{aligned}
\]
for every \(\epsilon>0\) and \(Q_{\boldsymbol{x}} \in \mathscr{D}\left(\Delta_{m}\right)\).

\section*{Appendix D: Proof of theoretical results}

\section*{D.1. Proof of Theorem 1}

The proof of this theorem follows the same reasoning as the proof of Theorem 1 in Barrientos et al. (2015). First, we re-define the three versions of the model by means of a mapping \(S\) and the measurability of the process is then proved by showing that the mapping \(S\) is continuous. This is stated in Lemma 1 below, which is an extension of Lemma B.1.1 in the supplementary material of Barrientos et al. (2015).

The proofs of parts (i), (ii), and (iii) in Lemma 1 follow the same reasoning as the corresponding proofs for Lemma B.1.1. The main difference in the following proofs comes from the number of weights that the multivariate Bernstein polynomial on the \(m\)-dimensional symplex has. For completeness, the proof of Theorem 1 is provided below.

Let \(\mathbf{T}, \mathbf{T}^{\theta}\) and \(\mathbf{T}^{w}\) be dependent stick-breaking processes of the form:
- \(\mathbf{T}=\left\{T_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\), where \(T_{\boldsymbol{x}}(\omega, \cdot)=\sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}, \omega) \delta_{\boldsymbol{\theta}_{j}(\boldsymbol{x}, \omega)}(\cdot)\), where \(w_{j}(\boldsymbol{x}, \omega)\) and \(\boldsymbol{\theta}_{j}(\boldsymbol{x}, \omega)\) are defined as in Definition 1.
- \(\mathbf{T}^{\theta}=\left\{T_{\boldsymbol{x}}^{\theta}: \boldsymbol{x} \in \mathscr{X}\right\}\), where \(T_{\boldsymbol{x}}^{\theta}(\omega, \cdot)=\sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}, \omega) \delta_{\boldsymbol{\theta}_{j}(\omega)}(\cdot)\), where \(w_{j}(\boldsymbol{x}, \omega)\) and \(\boldsymbol{\theta}_{j}(\omega)\) are defined as in Definition 2.
- \(\mathbf{T}^{w}=\left\{T_{\boldsymbol{x}}^{w}: \boldsymbol{x} \in \mathscr{X}\right\}\), where \(T_{\boldsymbol{x}}^{w}(\omega, \cdot)=\sum_{j=1}^{\infty} w_{j}(\omega) \delta_{\boldsymbol{\theta}_{j}(\boldsymbol{x}, \omega)}(\cdot)\), where \(w_{j}(\omega)\) and \(\boldsymbol{\theta}_{j}(\boldsymbol{x}, \omega)\) are defined as in Definition 3.

Let \(S\) be a mapping defined on \(\mathbb{N} \times \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\) of the form
\[
S\left(k_{0}, \mathcal{Q}\right):=\left\{H\left(k_{0}, Q_{\boldsymbol{x}}\right): \boldsymbol{x} \in \mathscr{X}\right\},
\]
where \(k_{0} \in \mathbb{N}, \mathcal{Q}=\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\) and \(H\left(k_{0}, Q_{\boldsymbol{x}}\right)\) is the probability measure associated to the Bernstein polynomial of degree \(k_{0}\) of the measure \(Q_{\boldsymbol{x}}\). \(F\) can be expressed as \(S(k, \mathbf{T}), S\left(k, \mathbf{T}^{\theta}\right)\) or \(S\left(k, \mathbf{T}^{w}\right)\), when \(F\) corresponds to \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \mathscr{H}\right), \theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) and \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{\boldsymbol{z}}\right.\), \(\mathscr{H}\) ), respectively. Since \(\mathbf{T}, \mathbf{T}^{\theta}\) and \(\mathbf{T}^{w}\) are well-defined stochastic processes, to prove the measurability of \(F\), it suffices to prove the measurability of \(S\) which is proven by showing that mapping \(S\) is continuous. For this, it is necessary to consider some topologies in the space where the mapping is valued and defined. This topologies and spaces are described below.

Let \(\mathscr{T}_{1}\) be the weak product topology for the space \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) and let \(\mathscr{T}_{2}\) and \(\mathscr{T}_{3}\) be the \(L_{\infty}\) product topology and \(L_{\infty}\) topology for the space \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\), respectively. A sub-base for the weak product topology, \(\mathscr{T}_{4}\), for the space \(\mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}=\prod_{\boldsymbol{x} \in \mathscr{X}} \mathscr{P}\left(\tilde{\Delta}_{m}\right)\) is given by sets of the form \(\tilde{B}_{f, \epsilon, \boldsymbol{x}_{0}}^{W}(\mathcal{Q})= \prod_{\boldsymbol{x} \in \mathscr{X}} \tilde{\Delta}_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(Q_{\boldsymbol{x}}\right)\), where \(\tilde{\Delta}_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(Q_{\boldsymbol{x}}\right)=\Delta_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(Q_{\boldsymbol{x}}\right) \bigcap \mathscr{P}\left(\tilde{\Delta}_{m}\right)\), with \(\mathcal{Q} \in \mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}, f: \Delta_{m} \longrightarrow \mathbb{R}\) a bounded continuous function, \(\epsilon>0\) and \(\boldsymbol{x}_{0} \in \mathscr{X}\). A sub-base for the product topology, \(\mathscr{L}_{1}\), for the space \(\mathbb{N} \times \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\) is given by sets of the form \(B_{f, \epsilon, \boldsymbol{x}_{0}}^{D \times W}(\mathcal{Q})=\prod_{\boldsymbol{x} \in \mathscr{X}}\left[\left\{k_{0}\right\} \times \tilde{\Delta}_{f, \epsilon, \boldsymbol{x}_{0}}^{W}\left(Q_{\boldsymbol{x}}\right)\right]\). Finally, a sub-base
for the product topology, \(\mathscr{L}_{2}\), for the space \(\mathbb{N} \times \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\) is given by sets of the form \(B_{\epsilon, N}^{D \times L_{\infty}}\left(k_{0}, \mathcal{Q}\right)=\left\{k_{0}\right\} \times \tilde{\Delta}_{\epsilon, N}^{L_{\infty}}(\mathcal{Q})\), where \(\tilde{\Delta}_{\epsilon, N}^{L_{\infty}}(\mathcal{Q})\) is given by
\(\left\{\left\{M_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}: \max _{\mathbf{j} \in \mathcal{H}_{N, m}^{0}} \sup _{\boldsymbol{x} \in \mathscr{X}}\left|M_{\boldsymbol{x}}\left(A_{\mathbf{j}, N}\right)-Q_{\boldsymbol{x}}\left(A_{\mathbf{j}, N}\right)\right|<\epsilon\right\}\),
where \(k_{0} \in \mathbb{N}, N \in \mathbb{N}, \epsilon>0, A_{\mathbf{j}, N}=\left(\frac{j_{1}-1}{N}, \frac{j_{1}}{N}\right] \times \ldots, \times\left(\frac{j_{m}-1}{N}, \frac{j_{m}}{N}\right]\) and \(\mathcal{Q} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\).

The following lemma states that mapping \(S\) defined by expression (15), is continuous under \(\mathscr{T}_{1}, \mathscr{T}_{2}\) and \(\mathscr{T}_{3}\) in the space where \(S\) is valued, thus ensuring that \(F\) is measurable under \(\mathscr{B}_{1}, \mathscr{B}_{2}\) and \(\mathscr{B}_{3}\), respectively.

Lemma 1. Let \(S\) be a mapping defined as in Equation (15), then
(i) \(S:\left(\mathbb{N} \times \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}, \mathscr{L}_{1}\right) \longrightarrow\left(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}, \mathscr{B}_{1}\right)\),
(ii) \(S:\left(\mathbb{N} \times \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}, \mathscr{L}_{1}\right) \longrightarrow\left(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}, \mathscr{B}_{2}\right)\),
(iii) \(S:\left(\mathbb{N} \times \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}, \mathscr{L}_{2}\right) \longrightarrow\left(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}, \mathscr{B}_{3}\right)\),
are continuous.
The proof of each part of Lemma 1 is given below:
(i) Let \(\mathcal{Q} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}, k_{0} \in \mathbb{N}\) and
\[
V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon\right)=\bigcap_{i=1}^{L} \bigcap_{j=1}^{K_{i}} B_{f_{i j}, \epsilon, \boldsymbol{x}_{i}}^{W}\left(S\left(k_{0}, \mathcal{Q}\right)\right)
\]
where \(L, K_{i}, i \in\{1, \ldots, L\}\), are positive integers, \(f_{i j}, j=1, \ldots, K_{i}, i= 1, \ldots, L\), are bounded continuous functions, \(\epsilon>0\) and \(\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{L}\right) \in \mathscr{X}^{L}\). The proof is based on finding and open set \(U \in \mathscr{L}_{1}\) such that \(\left(k_{0}, \mathcal{Q}\right) \in U\) and \(S(U) \subseteq V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon\right)\).
Notice that for every \(\mathcal{M}=\left\{M_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\),
\[
\begin{aligned}
& \left|\int_{\Delta_{m}} f_{i j} d H\left(k_{0}, M_{\boldsymbol{x}_{i}}\right)-\int_{\Delta_{m}} f_{i j} d H\left(k_{0}, Q_{\boldsymbol{x}_{i}}\right)\right| \\
& \quad \leq \int_{\Delta_{m}}\left|f_{i j}(\boldsymbol{y})\right| \sum_{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}}\left|M_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)\right| \operatorname{dir}\left(\boldsymbol{y} \mid \alpha\left(k_{0}, \mathbf{j}\right)\right) \\
& \quad \leq \frac{M_{0}\left(k_{0}+m-1\right)!}{m!\left(k_{0}-1\right)!} N_{k_{0}}(\mathcal{M}, \mathcal{Q})
\end{aligned}
\]
where \(\alpha(k, \boldsymbol{j})=\left(\boldsymbol{j}, k+m-\|\mathbf{j}\|_{1}\right)\) and
\[
N_{k_{0}}(\mathcal{M}, \mathcal{Q})=\max _{i \in\{1, \ldots, L\}} \max _{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}}\left|M_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)\right|
\]
\(M_{0}=\max _{i \in\{1, \ldots, L\}} \max _{j \in\left\{1, \ldots, K_{i}\right\}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|f_{i j}(\boldsymbol{y})\right|,\|\cdot\|_{1}\) denotes the \(l_{1}\)-norm, and \(A_{\mathbf{j}, k_{0}}=\left(\frac{j_{1}-1}{k_{0}}, \frac{j_{1}}{k_{0}}\right] \times \ldots \times\left(\frac{j_{m}-1}{k_{0}}, \frac{j_{m}}{k_{0}}\right]\). From Lemma 1 in

Barrientos et al. (2012), there exists \(\mathcal{Q}^{\prime}=\left\{Q_{\boldsymbol{x}}^{\prime}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\) such that for every \(\boldsymbol{x} \in \mathscr{X}, Q_{\boldsymbol{x}}^{\prime}\) is absolutely continuous w.r.t Lebesgue measure and such that,
\[
N_{k_{0}}\left(\mathcal{Q}^{\prime}, \mathcal{Q}\right) \leq \frac{m!\left(k_{0}-1\right)!}{2 M_{0}\left(k_{0}+m-1\right)!} \epsilon
\]

Since \(Q_{\boldsymbol{x}_{i}}^{\prime}, i=1, \ldots, L\), is an absolutely continuous measure, w.r.t. Lebesgue measure, then \(A_{\mathbf{j}, k_{0}}, \mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}\), are sets of \(Q_{\boldsymbol{x}_{i}}^{\prime}\) continuity, i.e., the boundaries of \(A_{\mathbf{j}, k_{0}}\) have null \(Q_{\boldsymbol{x}_{i}}^{\prime}\) measure, for every \(\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}\) and every \(i=1, \ldots, L\). Thus, the set
\[
\begin{aligned}
U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right) & =\bigcap_{i=1}^{L}\left\{M_{\boldsymbol{x}_{i}} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right): \max _{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}}\left|M_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{\mathbf{j}, k_{0}}\right)\right| \leq \tilde{\epsilon}\right\} \\
& =\left\{\mathcal{M} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}: N_{k_{0}}\left(\mathcal{M}, \mathcal{Q}^{\prime}\right) \leq \tilde{\epsilon}\right\}
\end{aligned}
\]
belongs to \(\mathscr{T}_{4}\). Notice that if \(\tilde{\epsilon}=\frac{m!\left(k_{0}-1\right)!}{2 M_{0}\left(k_{0}+m-1\right)!} \epsilon\), then
\[
\left|\int_{\Delta_{m}} f_{i j} d H\left(k_{0}, M_{\boldsymbol{x}_{i}}\right)-\int_{\Delta_{m}} f_{i j} d H\left(k_{0}, Q_{\boldsymbol{x}_{i}}\right)\right|<\epsilon
\]
where \(H\left(k_{0}, Q_{\boldsymbol{x}}\right)\) is the probability measure associated to the multivariate Bernstein polynomial of measure \(Q_{\boldsymbol{x}}\) of degree \(k_{0}\). Therefore, if \(U=\left\{k_{0}\right\} \times U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right)\), then \(U \in \mathscr{L}_{1},\left(k_{0}, \mathcal{Q}\right) \in U\) and \(S(U) \subseteq V\left(S\left(k_{0}, \mathcal{Q}\right), \epsilon\right)\), which completes the proof of (i) in Lemma 1.
(ii) Let \(\mathcal{Q} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}, k_{0} \in \mathbb{N}\) and \(V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon\right)=\bigcap_{i=1}^{L} B_{\epsilon, \boldsymbol{x}_{i}}^{L_{\infty}}\left(S\left(k_{0}, \mathcal{Q}\right)\right)\), where \(L\) is a positive integer, \(\epsilon>0\) and \(\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{L}\right) \in \mathscr{X}^{L}\). The proof is based on finding and open set \(U \in \mathscr{L}_{1}\) such that \(\left(k_{0}, \mathcal{Q}\right) \in U\) and \(S(U) \subseteq V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon\right)\).
Notice that for every \(\mathcal{M}=\left\{M_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\),
\[
\begin{aligned}
\sup _{\boldsymbol{y} \in \Delta_{m}} \mid b(\boldsymbol{y} \mid & \left.k_{0}, M_{\boldsymbol{x}_{i}}\right)-b\left(\boldsymbol{y} \mid k_{0}, Q_{\boldsymbol{x}_{i}}\right) \mid \\
& \leq \sup _{\boldsymbol{y} \in \Delta_{m}} \sum_{\mathbf{j} \in \mathcal{H}_{k_{0}}^{0}, m}\left|M_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)\right| \operatorname{dir}\left(\boldsymbol{y} \mid \alpha\left(k_{0}, \mathbf{j}\right)\right) \\
& \leq \frac{M_{0}\left(k_{0}+m-1\right)!}{m!\left(k_{0}-1\right)!} N_{k_{0}}(\mathcal{M}, \mathcal{Q})
\end{aligned}
\]
where \(M_{0}=\max _{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}} \sup _{\boldsymbol{y} \in \Delta_{m}} \operatorname{dir}\left(\boldsymbol{y} \mid \mathbf{j}, k_{0}+m-\|\mathbf{j}\|_{1}\right), b\left(\boldsymbol{y} \mid k_{0}, M_{\boldsymbol{x}_{i}}\right)\) stands for the density function of the multivariate Bernstein polynomial of function \(M_{\boldsymbol{x}_{i}}\) of degree \(k_{0}\), and \(N_{k_{0}}(\mathcal{M}, \mathcal{Q})\) and \(A_{\mathbf{j}, k_{0}}\) are defined as in part (i) of the proof. By the same arguments from part (i), it follows that if \(U=\left\{k_{0}\right\} \times U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right)\), where \(U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right)\) is defined as in (17), with \(\tilde{\epsilon}=\frac{m!\left(k_{0}-1\right)!}{2 M_{0}\left(k_{0}+m-1\right)!} \epsilon\), then \(U \in \mathscr{L}_{1},\left(k_{0}, \mathcal{Q}\right) \in U\) and \(S(U) \subseteq V\left(S\left(k_{0}, \mathcal{Q}\right), \epsilon\right)\), which completes the proof of (ii) in Lemma 1.
(iii) Let \(\mathcal{Q} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}, k_{0} \in \mathbb{N}\) and \(V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon\right)=B_{\epsilon}^{L_{\infty}}\left(S\left(k_{0}, \mathcal{Q}\right)\right)\). The proof is based on finding and open set \(U \in \mathscr{L}_{2}\) such that \(\left(k_{0}, \mathcal{Q}\right) \in U\) and \(S(U) \subseteq V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon\right)\).
Notice that for every \(\mathcal{M}=\left\{M_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}\),
\[
\begin{aligned}
& \sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|b\left(\boldsymbol{y} \mid k_{0}, M_{x}\right)-b\left(\boldsymbol{y} \mid k_{0}, Q_{x}\right)\right| \\
& \quad \leq \frac{M_{0}\left(k_{0}+m-1\right)!}{m!\left(k_{0}-1\right)!} \sup _{\boldsymbol{x} \in \mathscr{X}} \max _{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}}\left|M_{\boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right)\right|
\end{aligned}
\]
where \(M_{0}=\max _{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}} \sup _{\boldsymbol{y} \in \Delta_{m}}, \operatorname{dir}\left(\boldsymbol{y} \mid \mathbf{j}, k_{0}+m-\|\mathbf{j}\|_{1}\right)\), and \(A_{\mathbf{j}, k_{0}}\) are defined as in the proof of (i). Then, if \(U=\left\{k_{0}\right\} \times \tilde{\Delta}_{\tilde{\epsilon}, k_{0}}^{L_{\infty}}(\mathcal{Q})\), where \(\tilde{\Delta}_{\tilde{\epsilon}, k_{0}}^{L_{\infty}}(\mathcal{Q})\) is defined as in (16), with \(\tilde{\epsilon}=\frac{m!\left(k_{0}-1\right)!}{M_{0}\left(k_{0}+m-1\right)!} \epsilon\), then \(U \in \mathscr{L}_{2}\), \(\left(k_{0}, \mathcal{Q}\right) \in U\) and \(S(U) \subseteq V\left(S\left(k_{0}, \mathcal{Q}\right), \epsilon\right)\), which completes the proof of (iii) in Lemma 1.

\section*{D.2. Proof of Theorem 2}

Proving that \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) and \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) are the support of \(F\) under the weak product and \(L_{\infty}\) product topology, are direct extensions of the proofs of Theorems 2 and 3 in Barrientos et al. (2017), respectively. For completeness, the proof of this theorem is provided in what follows. First, we prove that \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the weak product topology. Then we prove that \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the \(L_{\infty}\) product topology. In each case all three versions of \(F\) are considered.

To prove that \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the weak product topology, it suffices to prove that any open set of the weak product topology has positive \(P \circ F^{-1}\)-measure. Let \(\mathcal{Q} \in \mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) and \(V(\mathcal{Q} ; \epsilon)=\bigcap_{i=1}^{L} \bigcap_{j=1}^{K_{i}} B_{f_{i j}, \epsilon, \boldsymbol{x}_{i}}^{W}(\mathcal{Q})\), where \(L, K_{i}, i=1, \ldots, L\), are positive integers, \(f_{i j}, j=1, \ldots, K_{i}, i=1, \ldots, L\), are bounded continuous functions, \(\epsilon>0\) and \(\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{L}\right) \in \mathscr{X}^{L}\). From Lemma 1 in Barrientos et al. (2012), there exists \(\mathcal{Q}^{\prime}=\left\{Q_{\boldsymbol{x}}^{\prime}: \boldsymbol{x} \in \mathscr{X}\right\} \in \mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\), such that for every \(\boldsymbol{x} \in \mathscr{X}, Q_{\boldsymbol{x}}^{\prime}\) is absolutely continuous w.r.t Lebesgue measure and such that \(Q_{\boldsymbol{x}}^{\prime}=Q_{\boldsymbol{x}}\) if \(\boldsymbol{x} \neq \boldsymbol{x}_{i}\) and
\[
\left|\int_{\Delta_{m}} f_{i j} Q_{\boldsymbol{x}_{i}}-\int_{\Delta_{m}} f_{i j} d Q_{\boldsymbol{x}_{i}}^{\prime}\right|<\frac{\epsilon}{2}
\]
if \(\boldsymbol{x}=\boldsymbol{x}_{i}, i=1, \ldots, L\). Then, \(V\left(\mathcal{Q}^{\prime} ; \epsilon / 2\right) \subset V(\mathcal{Q} ; \epsilon)\). Since for every \(\boldsymbol{x} \in \mathscr{X}\), \(H\left(k, Q_{\boldsymbol{x}}^{\prime}\right)\) converges weakly to \(Q_{\boldsymbol{x}}^{\prime}\) as \(k \rightarrow \infty\), for every \(\epsilon>0\), there exists large enough \(k_{0} \in \mathbb{N}\) such that
\[
\left|\int_{\Delta_{m}} f_{i j} d H\left(k_{0}, Q_{\boldsymbol{x}_{i}}^{\prime}\right)-\int_{\Delta_{m}} f_{i j} d Q_{\boldsymbol{x}_{i}}^{\prime}\right|<\frac{\epsilon}{4}
\]
then \(V\left(S\left(k_{0}, \mathcal{Q}^{\prime}\right) ; \epsilon / 4\right) \subset V\left(\mathcal{Q}^{\prime} ; \epsilon / 2\right)\). By Lemma 1 part ( \(i\) ), there exists \(U= \left\{k_{0}\right\} \times U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right) \in \mathscr{L}_{1}\), with \(\tilde{\epsilon}=\frac{m!\left(k_{0}-1\right)!}{4 M_{0}\left(k_{0}+m-1\right)!} \epsilon, k_{0} \in \mathbb{N}\) and \(U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right) \in \mathscr{T}_{4}\),
such that \(S(U) \subset V\left(S\left(k_{0}, \mathcal{Q}^{\prime}\right) ; \epsilon / 4\right)\). Thus, to prove this theorem, it suffices to prove that \(P \circ F^{-1}(V(\mathcal{Q} ; \epsilon)) \geq P\{\omega \in \Omega:(k(\omega), \overline{\mathbf{T}}) \in U\}>0\), where \(U=\left\{k_{0}\right\} \times U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right)\), with \(U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right)\) defined as in (17) and \(\overline{\mathbf{T}}\) is either \(\mathbf{T}, \mathbf{T}^{\theta}\) or \(\mathbf{T}^{w}\). Before considering each case, note that there are \(N=\frac{\left(k_{0}+m-1\right)!}{m!\left(k_{0}-1\right)!}\) disjoint sets in \(\mathcal{H}_{k_{0}, m}^{0}\). Each of these sets is denoted by \(A_{[l], N}, l=1, \ldots, N\).

When \(\overline{\mathbf{T}}\) is \(\mathbf{T}\), the assumption that the stochastic processes \(\eta_{j}\) and \(\boldsymbol{z}_{j}\) are well defined and have full support implies that
\[
\begin{aligned}
& P\{\omega \in \Omega:( k(\omega), \mathbf{T}) \in U\} \geq P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& \times P\left\{\omega \in \Omega: \max _{i \in\{1, \ldots, L\} \mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}}\left|M_{\boldsymbol{x}_{i}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{\mathbf{j}, k_{0}}\right)\right| \leq \tilde{\epsilon}\right\}, \\
& \geq P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& \times \prod_{l=1}^{N} P\left\{\omega \in \Omega:\left(\boldsymbol{\theta}_{l}\left(\boldsymbol{x}_{1}, \omega\right), \ldots, \boldsymbol{\theta}_{l}\left(\boldsymbol{x}_{L}, \omega\right)\right) \in A_{[l], N}^{L}\right\} \\
& \times P\left\{\omega \in \Omega:\left(V_{l}\left(\boldsymbol{x}_{1}, \omega\right), \ldots, V_{l}\left(\boldsymbol{x}_{L}, \omega\right)\right) \in B_{l}^{L}, l \in\{1, \ldots, N\}\right\} \\
& \times \prod_{l=N+1}^{\infty} P\left\{\omega \in \Omega:\left(\boldsymbol{\theta}_{l}\left(\boldsymbol{x}_{1}, \omega\right), \ldots, \boldsymbol{\theta}_{l}\left(\boldsymbol{x}_{L}, \omega\right)\right) \in \tilde{\Delta}_{m}^{L}\right\} \\
& \times \prod_{l=N+1}^{\infty} P\left\{\omega \in \Omega:\left(V_{l}\left(\boldsymbol{x}_{1}, \omega\right), \ldots, V_{l}\left(\boldsymbol{x}_{L}, \omega\right)\right) \in[0,1]^{L}\right\}, \\
&>0,
\end{aligned}
\]
where
\(B_{1}^{L}=\bigotimes_{i=1}^{L}\left\{Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{[1], N}\right)-\frac{\tilde{\epsilon}}{4(N-1)} ; Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{[1], N}\right)+\frac{\tilde{\epsilon}}{4(N-1)}\right\}\),
\(B_{l}^{L}=\bigotimes_{i=1}^{L}\left\{\frac{Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{[l], N}\right)-\frac{\tilde{\epsilon}}{4(N-1)}}{\prod_{l_{1}<l}\left[1-V_{l_{1}}\left(\boldsymbol{x}_{i}, \omega\right)\right]} ; \frac{Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{[l], N}\right)+\frac{\tilde{\epsilon}}{4(N-1)}}{\prod_{l_{1}<l}\left[1-V_{l_{1}}\left(\boldsymbol{x}_{i}, \omega\right)\right]}\right\}, l=2, \ldots, N-1\),
\(B_{N}^{L}=\bigotimes_{i=1}^{L}\left\{\frac{Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{[N], N}\right)-\frac{\tilde{\epsilon}}{3}}{\prod_{l_{1}<N}\left[1-V_{l_{1}}\left(\boldsymbol{x}_{i}, \omega\right)\right]} ; \frac{Q_{\boldsymbol{x}_{i}}^{\prime}\left(A_{[N], N}\right)-\frac{\tilde{\epsilon}}{4}}{\prod_{l_{1}<N}\left[1-V_{l_{1}}\left(\boldsymbol{x}_{i}, \omega\right)\right]}\right\}\),
\(A_{[l], N}^{L}=\bigotimes_{i=1}^{L} A_{[l], N}, l=1, \ldots, N, \tilde{\Delta}_{m}^{L}=\bigotimes_{i=1}^{L} \tilde{\Delta}_{d}\) and \([0,1]^{L}=\bigotimes_{i=1}^{L}[0,1]\).
This completes the proof that \(F\) considered as \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \mathscr{H}\right)\) has weak product support.

When \(\overline{\mathbf{T}}\) is \(\mathbf{T}^{\theta}\), the assumption that the stochastic processes \(\eta_{j}\) and the random vectors \(\boldsymbol{\theta}_{j}\) are well defined and have full support imply that
\[
\begin{aligned}
P\{\omega \in \Omega:(k(\omega), & \left.\left.\mathbf{T}^{\theta}\right) \in U\right\} \geq P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& \quad \times \prod_{l=1}^{N} P\left\{\omega \in \Omega:\left(\boldsymbol{\theta}_{l}(\omega), \ldots, \boldsymbol{\theta}_{l}(\omega)\right) \in A_{[l], N}^{L}\right\}
\end{aligned}
\]
\[
\begin{aligned}
\times P\{\omega \in \Omega & \left.:\left(V_{l}\left(\boldsymbol{x}_{1}, \omega\right), \ldots, V_{l}\left(\boldsymbol{x}_{L}, \omega\right)\right) \in B_{l}^{L}, l \in\{1, \ldots, N\}\right\} \\
& >0
\end{aligned}
\]
where \(B_{1}^{L}, B_{l}^{L}, l=2, \ldots, N-1, B_{N}^{L}\) and \(A_{[l], N}^{L}, l=1, \ldots, N\), are defined as above. This completes the proof that \(F\) considered as \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) has weak product support.

Finally, when \(\overline{\mathbf{T}}\) is \(\mathbf{T}^{w}\). Since \(\Delta_{m}\) is a separable space and \(\tilde{\Delta}_{m}\) is dense in \(\Delta_{m}\), then the space of measures whose support points are finite subsets of \(\tilde{\Delta}_{m}\) is dense in \(\mathscr{P}\left(\Delta_{m}\right)\) (Parthasarathy, 1967). Then, for each \(\boldsymbol{x} \in \mathscr{X}\), there exists a probability measure \(\tilde{Q}_{\boldsymbol{x}}(\cdot)=\sum_{j=1}^{R} \tilde{w}_{j} \delta_{\tilde{\boldsymbol{\theta}}_{j}(\boldsymbol{x})}(\cdot)\), defined on \(\tilde{\Delta}_{m}\), where \(R\) is an integer, \(\tilde{w}_{j} \in[0,1], j=1, \ldots, R, \sum_{j=1}^{R} \tilde{w}_{j}=1\), and \(\tilde{\boldsymbol{\theta}}_{j}(\boldsymbol{x}) \in \tilde{\Delta}_{m}\) are continuous functions of \(\boldsymbol{x}, j=1, \ldots, R\), such that, for every \(\boldsymbol{x} \in \mathscr{X}, \mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}\),
\[
\left|\tilde{Q}_{\boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}}^{\prime}\left(A_{\mathbf{j}, k_{0}}\right)\right|<\frac{\tilde{\epsilon}}{2} .
\]

Then \(U^{\prime}(\tilde{Q} ; \tilde{\epsilon} / 2) \subset U^{\prime}\left(Q^{\prime} ; \tilde{\epsilon}\right)\), where \(U^{\prime}(\mathcal{Q} ; \epsilon)\) is defined as in (17). Thus, it suffices to prove that
\[
P\left\{\omega \in \Omega:\left(k(\omega), \mathbf{T}^{w}\right) \in\left\{k_{0}\right\} \times U^{\prime}(\tilde{Q} ; \tilde{\epsilon} / 2)\right\}>0
\]

Consider \(\left\{\tilde{A}_{[\tilde{l}], M}\right\}_{\tilde{l}=1}^{M}\), a finer partition of \(\mathcal{H}_{k_{0}, m}^{0}\) than \(\left\{A_{[l], N}\right\}_{l=1}^{N}\), such that \(A_{[1], N}=\bigcup_{\tilde{l}=1}^{n_{1}} \tilde{A}_{[\tilde{l}], M}\) and \(A_{[l], N}=\bigcup_{\tilde{l}=n_{l-1}+1}^{n_{l}} \tilde{A}_{\tilde{l}, M}, l=1, \ldots, N\), where \(\sum_{l=1}^{N} n_{l}= M\). Then, the assumption that the stochastic processes \(\boldsymbol{z}_{j}\) and the random variables \(v_{j}\) are well defined and have full support imply that
\(P\left\{\omega \in \Omega:\left(k(\omega), \mathbf{T}^{w}\right) \in U\right\} \geq P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\}\)
\(\times P\left\{\omega \in \Omega:\left(\left\lceil k_{0} \boldsymbol{\theta}_{1}\left(\boldsymbol{x}_{i}, \omega\right)\right\rceil-\left\lceil k_{0} \tilde{\boldsymbol{\theta}}_{j}\left(\boldsymbol{x}_{i}\right)\right\rceil\right)=\mathbf{0}, m=1, \ldots, L, j=1, \ldots, n_{1}\right\}\)
\(\times \prod_{l=2}^{N} P\left\{\omega \in \Omega:\left(\left\lceil k_{0} \boldsymbol{\theta}_{l}\left(\boldsymbol{x}_{i}, \omega\right)\right\rceil-\left\lceil k_{0} \tilde{\boldsymbol{\theta}}_{j}\left(\boldsymbol{x}_{i}\right)\right\rceil\right)=\mathbf{0}, i=1, \ldots, L, j=n_{l-1}+1, \ldots, n_{l}\right\}\)
\(\times P\left\{\omega \in \Omega: v_{l}(\omega) \in B_{l}^{L}, l=1, \ldots, N\right\}\),
> 0,
where
\[
\begin{aligned}
B_{1}^{L} & =\left\{\sum_{j=1}^{n_{1}} \tilde{w}_{j}-\frac{\tilde{\epsilon}}{8(N-1)} ; \sum_{j=1}^{n_{1}} \tilde{w}_{j}+\frac{\tilde{\epsilon}}{8(N-1)}\right\}, \\
B_{l}^{L} & =\left\{\frac{\sum_{j=n_{l-1}+1}^{n_{l}} \tilde{w}_{j}-\frac{\tilde{\epsilon}}{8(N-1)}}{\prod_{l_{1}<l}\left[1-V_{l_{1}}(\omega)\right]} ; \frac{\sum_{j=n_{l-1}+1}^{n_{l}} \tilde{w}_{j}+\frac{\tilde{\epsilon}}{8(N-1)}}{\prod_{l_{1}<l}\left[1-V_{l_{1}}(\omega)\right]}\right\}, l=2, \ldots, N-1,
\end{aligned}
\]
\[
B_{N}^{L}=\left\{\frac{1-\sum_{j=1}^{n_{N-1}} \tilde{w}_{j}-\frac{\tilde{\epsilon}}{6}}{\prod_{l_{1}<N}\left[1-V_{l_{1}}(\omega)\right]} ; \frac{1-\sum_{j=1}^{n_{N-1}} \tilde{w}_{j}-\frac{\tilde{\epsilon}}{8}}{\prod_{l_{1}<N}\left[1-V_{l_{1}}(\omega)\right]}\right\},
\]
which completes the proof that \(F\) considered as \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{H}\right)\) has weak product support. Thus the proof that \(\mathscr{P}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the weak product topology is completed.

Now we will prove that \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the \(L_{\infty}\) product toplogy. Note that it suffices to prove that any open set of the \(L_{\infty}\) product topology has positive \(P \circ F^{-1}\)-measure. Let \(\mathcal{Q} \in \mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) and \(V(\mathcal{Q} ; \epsilon)= \bigcap_{i=1}^{L} \Delta_{\epsilon, \boldsymbol{x}_{i}}^{L_{\infty}}(\mathcal{Q})\), where \(L\) is a positive integer, \(\epsilon>0\), and \(\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{L}\right) \in \mathscr{X}^{L}\). Recall that for every \(\boldsymbol{x} \in \mathscr{X}, Q_{\boldsymbol{x}} \in \mathscr{D}\left(\Delta_{m}\right)\) is an absolutely continuous measures, w.r.t. Lebesgue measure, with continuous density, \(q_{\boldsymbol{x}}\). By Theorem 1 in Barrientos et al. (2015), for every \(\epsilon>0\), there exists large enough \(k_{0} \in \mathbb{N}\), such that for every \(\boldsymbol{x} \in \mathscr{X}\),
\[
\sup _{\boldsymbol{y} \in \Delta_{m}}\left|b\left(\boldsymbol{y} \mid k_{0}, Q_{\boldsymbol{x}}\right)-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\frac{\epsilon}{2}
\]
where \(b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)\) stands for the density function of the multivariate Bernstein polynomial of degree \(k\) of function \(Q_{\boldsymbol{x}}\). Then \(V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon / 2\right) \subset V(\mathcal{Q} ; \epsilon)\). By Lemma 1 part ( \(i i\) ), there exists \(U=\left\{k_{0}\right\} \times U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right) \in \mathscr{L}_{1}\), with \(\tilde{\epsilon}= \frac{m!\left(k_{0}-1\right)!}{2 M_{0}\left(k_{0}+m-1\right)!} \epsilon, k_{0} \in \mathbb{N}\) and \(U^{\prime}\left(\mathcal{Q}^{\prime} ; \tilde{\epsilon}\right) \in \mathscr{T}_{4}\), such that \(S(U) \subset V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon / 2\right)\). In analogy with the weak product support proof, it suffices to prove that \(P \circ F^{-1}(V(\mathcal{Q} ; \epsilon)) \geq P\left\{\omega \in \Omega:(k(\omega), \overline{\mathbf{T}}) \in\left\{k_{0}\right\} \times U^{\prime}(\mathcal{Q} ; \tilde{\epsilon})\right\}>0\), where \(\overline{\mathbf{T}}\) is either \(\mathbf{T}, \mathbf{T}^{\theta}\) or \(\mathbf{T}^{w}\). By the same arguments used to prove the weak product support of \(F\), it follows that \(P \circ F^{-1}(V(\mathcal{Q} ; \epsilon))>0\), when \(F\) is considered as DMBPP, \(\theta \mathrm{DMBPP}\), or \(w \mathrm{DMBPP}\). This completes the proof that \(\mathscr{D}\left(\Delta_{m}\right)^{\mathscr{X}}\) is the support of \(F\) under the \(L_{\infty}\) product topology, and thus completes the proof of the theorem.

\section*{D.3. Proof of Theorem 3}

The proof of Theorem 3 follows the same reasoning of the proof of Theorem 4 in Barrientos et al. (2017). First, we state and prove Lemma 2 below, which is used in the proof of this theorem and is an extension of Lemma B.4.1 in the supplementary material of Barrientos et al. (2017), and then we prove that \(\tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) is contained in the support of \(F\) under the \(L_{\infty}\) topology.
Lemma 2. Let \(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) be an absolutely continuous measure, w.r.t. Lebesgue measure, such that the mapping \((\boldsymbol{x}, \boldsymbol{y}) \mapsto q_{\boldsymbol{x}}(\boldsymbol{y})\) is continuous, and consider \(\mathscr{X}\) a compact space on \(\mathbb{R}^{p}\). Denote \(b_{k, Q_{\boldsymbol{x}}}(\boldsymbol{y})\), the density function, w.r.t. Lebesgue measure, of the multivariate Bernstein polynomial of degree \(k\) of function \(Q_{\boldsymbol{x}}\). Then for every \(\epsilon>0\), there exists \(k_{0} \in \mathbb{N}\) such that for every
\[
\sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\epsilon .
\]

\section*{Proof of Lemma 2}

Without loss of generality, consider \(\mathscr{X}=[0,1]^{p}\), and a uniform marginal distribution for \(\boldsymbol{X}\) on \(\mathscr{X}\). Then, \(q_{\boldsymbol{x}}(\boldsymbol{y})\) denotes a joint density function on \(\Delta_{m} \times \mathscr{X}\). Note that \(b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)\) can be written as
\[
b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)=\sum_{\boldsymbol{j} \in \mathcal{H}_{k, m}^{0}}\left[\int_{A_{\mathbf{j}, k}} q_{\boldsymbol{x}}(\boldsymbol{y}) d \boldsymbol{y}\right] \times \operatorname{dir}(\boldsymbol{y} \mid \alpha(k, \boldsymbol{j}))
\]

Now, consider \(r \in \mathbb{N}, \boldsymbol{l}=\left(l_{1}, \ldots, l_{r}\right)\), where \(l_{s} \in \mathbb{N}, s=1, \ldots, r\), are positive integers, and define
\[
\begin{aligned}
p_{k, \boldsymbol{l}, Q_{\boldsymbol{x}}}(\boldsymbol{y})=\sum_{\boldsymbol{j} \in \mathcal{H}_{k, m}^{0}} \sum_{i_{1}=1}^{l_{1}} & \ldots \sum_{i_{r}=1}^{l_{r}}\left[\int_{B_{i_{1}}} \ldots \int_{B_{i_{r}}} \int_{A_{\mathbf{j}, k}} q_{\boldsymbol{x}}(\boldsymbol{y}) d \boldsymbol{y} d x_{r} \ldots d x_{1}\right] \\
\times & \prod_{s=1}^{r} \beta\left(x_{s} \mid a_{s}, b_{s}\right) \operatorname{dir}(\boldsymbol{y} \mid \alpha(k, \boldsymbol{j}))
\end{aligned}
\]
where \(B_{i_{s}}=\left(\frac{i_{s}-1}{l_{s}}, \frac{i_{s}}{l_{s}}\right], a_{s}=i_{s}, b_{s}=l_{s}-i_{s}+1, s=1, \ldots, r\), and \(\beta(\cdot \mid a, b)\) stands for a beta density with parameters \(a\) and \(b\). Since \((\boldsymbol{x}, \boldsymbol{y}) \mapsto q_{\boldsymbol{x}}(\boldsymbol{y})\) is a continuous mapping, it is easy to show that \(p_{k, \boldsymbol{l}, Q_{\boldsymbol{x}}}(\boldsymbol{y})\) can uniformly approximate any continuous density function defined on \(\Delta_{m} \times \mathscr{X}\). Thus, every \(k>k_{0}, l_{s}>l_{s, 0}\), \(s=1, \ldots, r\), it follows that
\[
\sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|p_{k, \boldsymbol{l}, Q_{\boldsymbol{x}}}(\boldsymbol{y})-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\epsilon / 2
\]

Now, noting that
\[
\sum_{i_{1}=1}^{l_{1}} \ldots \sum_{i_{r}=1}^{l_{r}}\left[\int_{B_{i_{1}}} \ldots \int_{B_{i_{r}}} \int_{A_{\mathbf{j}, k}} q_{\boldsymbol{x}}(\boldsymbol{y}) d \boldsymbol{y} d x_{r} \ldots d x_{1}\right] \times \prod_{s=1}^{r} \beta\left(x_{s} \mid a_{s}, b_{s}\right)
\]
is the density function of the multivariate Bernstein polynomial of degree \(l_{1}, \ldots\), \(l_{r}\), of the mapping
\[
\boldsymbol{x} \mapsto \int_{A_{\mathbf{j}, k}} q_{\boldsymbol{x}}(\boldsymbol{y}) d \boldsymbol{y}
\]
defined on \(\mathscr{X}\), it follows that (18) converges uniformly to \(\int_{A_{\mathbf{j}, k}} q_{\boldsymbol{x}}(\boldsymbol{y}) d \boldsymbol{y}\), as \(\left(l_{1}, \ldots, l_{r}\right) \rightarrow \infty\), component-wise. Therefore, for every \(l_{s}>l_{s, 1}, s=1, \ldots, r\),
\[
\sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)-p_{k, \boldsymbol{l}, Q_{\boldsymbol{x}}}(\boldsymbol{y})\right|<\sum_{\boldsymbol{j} \in \mathcal{H}_{k, m}^{0}} \frac{\tilde{\epsilon}}{2} \operatorname{dir}(\boldsymbol{y} \mid \alpha(k, \boldsymbol{j}))<\frac{\epsilon}{2},
\]
where \(\tilde{\epsilon}=\frac{m!(k-1)!}{M_{0}(k+m-1)!} \epsilon\), with \(M_{0}=\max _{\mathbf{j} \in \mathcal{H}_{k, m}^{0}} \sup _{\boldsymbol{y} \in \Delta_{m}} \operatorname{dir}(\boldsymbol{y} \mid \alpha(k, \boldsymbol{j}))\). Finally, for \(k>k_{0}, l_{s}>\max \left\{l_{s, 0}, l_{s, 1}\right\}, s=1, \ldots, r\), and an application of the triangle inequality, it follows that
\[
\sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\epsilon,
\]
which completes to proof of the lemma.
Now, note that to prove that \(\tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) is contained in the support of \(F\) under the \(L_{\infty}\) topology, it suffices to prove that any open set of the \(L_{\infty}\) topology has positive \(P \circ F^{-1}\)-measure. Let \(\mathcal{Q} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) and \(V(\mathcal{Q} ; \epsilon)=B_{\epsilon}^{L_{\infty}}(\mathcal{Q})\), \(\epsilon>0\). Recall that \(\mathscr{X}\) is compact, and \(Q_{\boldsymbol{x}} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)\) is an absolutely continuous measures, w.r.t. Lebesgue measure, with continuous density, \(q_{\boldsymbol{x}}\), sucth that \((\boldsymbol{x}, \boldsymbol{y}) \mapsto q_{\boldsymbol{x}}(\boldsymbol{y})\) is continuous. From Lemma 2, there exists large enough \(k_{0}\), such that,
\[
\sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|b\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\frac{\epsilon}{2},
\]
where \(\operatorname{bp}\left(\boldsymbol{y} \mid k, Q_{\boldsymbol{x}}\right)\) stands for the density function of the multivariate Bernstein polynomial of function \(Q_{\boldsymbol{x}}\) of degree \(k\). Then, \(V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon / 2\right) \subset V(\mathcal{Q} ; \epsilon)\). By Lemma 1 part (iii), there exists \(U=\left\{k_{0}\right\} \times \tilde{\Delta}_{\tilde{\epsilon}, k_{0}}^{L_{\infty}}(\mathcal{Q}) \in \mathscr{L}_{2}\), with \(\tilde{\epsilon}= \frac{m!\left(k_{0}-1\right)!}{2 M_{0}\left(k_{0}+m-1\right)!} \epsilon, k_{0} \in \mathbb{N}\), such that \(S(U) \subset V\left(S\left(k_{0}, \mathcal{Q}\right) ; \epsilon / 2\right)\). Thus, to prove this theorem, it suffices to prove that \(P \circ F^{-1}(V(\mathcal{Q} ; \epsilon)) \geq P\{\omega \in \Omega:(k(\omega), \overline{\mathbf{T}}) \in \left.\left\{k_{0}\right\} \times U^{\star}(\mathcal{Q} ; \tilde{\epsilon})\right\}>0\), where
\[
U^{\star}(\mathcal{Q} ; \tilde{\epsilon})=\left\{\mathcal{M} \in \mathscr{P}\left(\tilde{\Delta}_{m}\right)^{\mathscr{X}}: \sup _{\boldsymbol{x} \in \mathscr{X}} \max _{\mathbf{j} \in \mathcal{H}_{k_{0}, m}^{0}}\left|M_{\boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right)-Q_{\boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right)\right| \leq \tilde{\epsilon}\right\},
\]
and \(\overline{\mathbf{T}}\) is either \(\mathbf{T}, \mathbf{T}^{\theta}\) or \(\mathbf{T}^{w}\).
First we assume that \(\overline{\mathbf{T}}\) is \(\mathbf{T}\) and following a similar reasoning as in the proof of Theorem 2. Since the stochastic processes \(\eta_{j}\) and \(\boldsymbol{z}_{j}\) are well defined and have full support, \(A_{[l], N} \in \mathscr{B}\left(\Delta_{m}\right)\) and the mappings
\[
\begin{aligned}
\boldsymbol{x} & \mapsto Q_{\boldsymbol{x}}\left(A_{[l], N}\right), \\
\boldsymbol{x} & \mapsto \frac{Q_{\boldsymbol{x}}\left(A_{[l], N}\right) / 2}{\prod_{l_{1}<l}\left[1-V_{l_{1}}(\boldsymbol{x}, \omega)\right]},
\end{aligned}
\]
are continuous, it follows that,
\[
\begin{aligned}
P \circ & F^{-1}(V(\mathcal{Q} ; \epsilon)) \\
\geq & P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& \times \prod_{l=1}^{N} P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}}\left|\boldsymbol{\theta}_{l}(\boldsymbol{x}, \omega)\right| \in A_{[l], N}\right\}
\end{aligned}
\]
\[
\begin{aligned}
& \times P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}}\left|V_{1}(\boldsymbol{x}, \omega)-Q_{\boldsymbol{x}}\left(A_{[l], N}\right) / 2\right|<\frac{\tilde{\epsilon}}{2 N}\right\} \\
& \times P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}}\left|V_{l}(\boldsymbol{x}, \omega)-\frac{Q_{\boldsymbol{x}}\left(A_{[l], N}\right) / 2}{\prod_{l_{1}<l}\left[1-V_{l_{1}}(\boldsymbol{x}, \omega)\right]}\right|<\frac{\tilde{\epsilon}}{2 N}, l=2, \ldots, N\right\}, \\
&>0
\end{aligned}
\]
which completes the proof that when \(F\) is a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \mathscr{H}\right)\) it has full \(L_{\infty}\) support.

Now assume that \(\overline{\mathbf{T}}\) is \(\mathbf{T}^{\theta}\). Given the above proof, it is straightforward to prove that \(F\) considered as \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) has \(L_{\infty}\) support.

Finally, assume that \(\overline{\mathbf{T}}\) is \(\mathbf{T}^{w}\). Consider the partition \(\left\{A_{[l], N}\right\}_{l=1}^{N}\) of \(\Delta_{m}\) and for each \(l=1, \ldots, N\), consider \(\left\{\mathscr{X}_{l, j}\right\}_{j=1}^{N_{l}}\), a partition of space \(\mathscr{X}, N_{l} \in \mathbb{N}\), \(N_{l}>N\). Since \(Q_{\boldsymbol{x}} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)\) are such that \((\boldsymbol{y}, \boldsymbol{x}) \mapsto q_{\boldsymbol{x}}(\boldsymbol{y})\) are continuous, then \((\boldsymbol{y}, \boldsymbol{x}) \mapsto Q_{\boldsymbol{x}}(\boldsymbol{y})\) are continuous and can be approximated by functions of the form,
\[
\bar{Q}_{\boldsymbol{x}}(\boldsymbol{y})=\sum_{l=1}^{N} \sum_{j=1}^{N_{l}} a_{l, j} \mathbb{I}(\boldsymbol{x}, \boldsymbol{y})_{\left\{\mathscr{X}_{l, j} \times A_{[l], N}\right\}},
\]
where \(\left\{a_{l, j}\right\}_{j=1}^{N_{l}}, l=1, \ldots, N\) are positive constants, \(\mathbb{I}_{A}\) denotes the indicator function of set \(A, \boldsymbol{x} \in \mathscr{X}\) and \(\boldsymbol{y} \in \Delta_{m}\). Now, for each \(l=1, \ldots, N\), consider the mapping \(\left(a_{l, 1}, \ldots, a_{l, N_{l}}\right) \mapsto \tilde{w}_{l, j}=a_{l, j} / \sum_{j=1}^{N_{l}} a_{l, j}\) and the continuous mappings \(\boldsymbol{x} \mapsto \tilde{\boldsymbol{\theta}}_{l, j}(\boldsymbol{x})\), where \(\tilde{w}_{l, j} \in[0,1], \sum_{j=1}^{N_{l}} \tilde{w}_{l, j}=1, \tilde{\boldsymbol{\theta}}_{l, j}(\boldsymbol{x}) \in \tilde{\Delta}_{m}\) and \(\tilde{\boldsymbol{\theta}}\left(\mathscr{X}_{l, 1}, \ldots, \mathscr{X}_{l, N_{l}}\right)=\left\{\tilde{A}_{[l, j], N_{l}}\right\}_{j=1}^{N_{l}}\) is a finer partition of \(\mathcal{H}_{k_{0}, m}^{0}\) than \(\left\{A_{[l], N}\right\}_{l=1}^{N}\), such that \(A_{[l], N}=\bigcup_{j=1}^{n_{l}} \tilde{A}_{[l, j], N_{l}}, n_{l}<N_{l}\). Thus, for each \(l= 1, \ldots, N, \bar{Q}_{\boldsymbol{x}}\left(A_{[l], N}\right)\) can be written as a measure of the form
\[
\tilde{Q}_{\boldsymbol{x}}\left(A_{[l], N}\right)=\sum_{j=1}^{N_{l}} \tilde{w}_{l, j} \mathbb{I}\left\{\tilde{\boldsymbol{\theta}}_{l, j}(\boldsymbol{x})\right\}_{\left\{\tilde{A}_{[l, j], N_{l}}\right\}}
\]
such that, for every \(l=1, \ldots, N\),
\[
\sup _{\boldsymbol{x} \in \mathscr{X}}\left|\tilde{Q}_{\boldsymbol{x}}\left(A_{[l], N}\right)-Q_{\boldsymbol{x}}\left(A_{[l], N}\right)\right|<\frac{\tilde{\epsilon}}{2}
\]

Then \(U^{\star}(\tilde{\mathcal{Q}} ; \tilde{\epsilon} / 2) \subset U^{\star}(\mathcal{Q} ; \tilde{\epsilon})\), where \(U^{\star}(\mathcal{Q} ; \epsilon)\) is defined as (19). Thus, in analogy with the previous proofs, it suffices to prove that
\[
P\left\{\omega \in \Omega:\left(k(\omega), \mathbf{T}^{w}\right) \in\left\{k_{0}\right\} \times U^{\star}(\tilde{\mathcal{Q}} ; \tilde{\epsilon} / 2)\right\}>0
\]

Following a similar reasoning as in the proof of Theorem 2 when \(\overline{\mathbf{T}}\) was considered as \(\mathbf{T}^{w}\) and by the assumption that the stochastic processes \(\eta_{j}\) and \(\boldsymbol{z}_{j}\) are well defined and have full support, \(A_{[l], N} \in \mathscr{B}\left(\Delta_{m}\right)\), and the mappings
\[
\boldsymbol{x} \mapsto k_{0} \tilde{\boldsymbol{\theta}}_{l, j}(\boldsymbol{x}), \quad j=1, \ldots, n_{l}, \quad l=1, \ldots, N,
\]
are continuous, it follows that,
\[
\begin{aligned}
P \circ & F^{-1}(V(\mathcal{Q} ; \epsilon)) \\
\geq & P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& \times \prod_{l=1}^{N} P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}}\left(\left\lceil k_{0} \boldsymbol{\theta}_{l}(\boldsymbol{x}, \omega)\right\rceil-\left\lceil k_{0} \tilde{\boldsymbol{\theta}}_{l, j}(\boldsymbol{x}, \omega)\right\rceil\right)=\mathbf{0}, j=1, \ldots, n_{l}\right\} \\
& \times P\left\{\omega \in \Omega:\left|v_{1}(\omega)-\sum_{j=1}^{n_{1}} \tilde{w}_{1, j}(\omega) / 2\right|<\frac{\tilde{\epsilon}}{2 N}\right\} \\
& \times P\left\{\omega \in \Omega:\left|v_{l}(\omega)-\frac{\sum_{j=1}^{n_{l}} \tilde{w}_{l, j}(\omega) / 2}{\prod_{l_{2}<l}\left[1-v_{l_{2}}(\omega)\right]}\right|<\frac{\tilde{\epsilon}}{2 N}, l=2, \ldots, N\right\} \\
> & 0
\end{aligned}
\]
which completes the proof that when \(F\) is a \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{H}\right)\) it has full \(L_{\infty}\) support.

\section*{D.4. Proof of Theorem 4}

The proof of this Theorem is an extension of the proof of Corollary 1 in Barrientos et al. (2017). For completeness, the proof is given below.

Let \(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \in \tilde{\mathscr{D}}\left(\Delta_{m}\right)^{\mathscr{X}}\) with continuous density function \(\left\{q_{\boldsymbol{x}}: \boldsymbol{x} \in\right. \mathscr{X}\}\). Here we will prove that, for every \(\delta>0\), any Kullback-Leibler neighborhood of \(\left\{Q_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) has positive \(P \circ F^{-1}\)-measure. This is,
\[
P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}} K L\left(q_{\boldsymbol{x}}, f(\boldsymbol{x}, \omega)\right)<\delta\right\}>0
\]
where \(K L(q, p)=\int_{\Delta_{m}} q(\boldsymbol{y}) \log \left(\frac{q(\boldsymbol{y})}{p(\boldsymbol{y})}\right) d \boldsymbol{y}\). Since \(\mathscr{X}\) and \(\Delta_{m}\) are compact sets and \((\boldsymbol{x}, \boldsymbol{y}) \mapsto q_{\boldsymbol{x}}(\boldsymbol{y})\) is a continuous mapping, it follows that \(\inf _{\boldsymbol{x} \in \mathscr{X}} \inf _{\boldsymbol{y} \in \Delta_{m}} q_{\boldsymbol{x}}(\boldsymbol{y})\) exists and is bounded.

First, suppose that \(\inf _{\boldsymbol{x} \in \mathscr{X}} \inf _{\boldsymbol{y} \in \Delta_{m}} q_{\boldsymbol{x}}(\boldsymbol{y})>0\). If for every \(\epsilon>0\), \(\sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|f(\boldsymbol{x}, \omega)(\boldsymbol{y})-q_{\boldsymbol{x}}(\boldsymbol{y})\right|<\epsilon\), then \(\inf _{\boldsymbol{x} \in \mathscr{X}} \inf _{\boldsymbol{y} \in \Delta_{m}} f(\boldsymbol{x}, \omega)(\boldsymbol{y})>0\) and for every \(\epsilon^{\prime}>0\), there exists \(\epsilon>0\) such that for every \(\boldsymbol{x} \in \mathscr{X}\) and every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\log \left(\frac{q_{\boldsymbol{x}}(\boldsymbol{y})}{f(\boldsymbol{x}, \omega)(\boldsymbol{y})}\right)<\epsilon^{\prime} .
\]

This in turn implies that \(\sup _{\boldsymbol{x} \in \mathscr{X}} K L\left(q_{\boldsymbol{x}}, f(\boldsymbol{x}, \omega)\right)<\epsilon^{\prime}\). From Theorem 3, it follows that
\[
\begin{array}{r}
P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}} K L\left(q_{\boldsymbol{x}}, f(\boldsymbol{x}, \omega)\right)<\epsilon^{\prime}\right\}>P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}} \sup _{\boldsymbol{y} \in \Delta_{m}}\left|f(\boldsymbol{x}, \omega)-q_{\boldsymbol{x}}\right|<\epsilon\right\} \\
>0 .
\end{array}
\]

Now, suppose that \(\inf _{\boldsymbol{x} \in \mathscr{X}} \inf _{\boldsymbol{y} \in \Delta_{d}} q_{\boldsymbol{x}}(\boldsymbol{y}) \approx 0\). Here we use a similar reasoning as in the proof of Theorem 2 of Petrone \& Wasserman (2002). Consider \(a>0\) and
\[
q_{\boldsymbol{x}}^{1}(\boldsymbol{y})=\frac{q_{\boldsymbol{x}}(\boldsymbol{y}) \vee a}{\int_{\Delta_{m}} q_{\boldsymbol{x}}(\boldsymbol{y}) \vee a d \boldsymbol{y}},
\]
where \(a \vee b\) stands for the maximum between \(a\) and \(b\). Clearly \(q_{\boldsymbol{x}}^{1}(\boldsymbol{y})\) is a density function such that \(q_{\boldsymbol{x}}(\boldsymbol{y}) \leq C q_{\boldsymbol{x}}^{1}(\boldsymbol{y})\), with \(C=\int_{\Delta_{m}} q_{\boldsymbol{x}}(\boldsymbol{y}) \vee a d \boldsymbol{y}\), and is greater than zero. Hence \(\sup _{\boldsymbol{x} \in \mathscr{X}} K L\left(q_{\boldsymbol{x}}^{1}, f(\boldsymbol{x}, \omega)\right)<\epsilon^{\prime}\). Considering \(a\) and \(\epsilon^{\prime}\) sufficiently small, it follows that there exists \(\tilde{\epsilon}>0\), sucht that
\[
K L\left(q_{\boldsymbol{x}}, f(\boldsymbol{x}, \omega)\right) \leq(C+1) \log (C)+C\left\{K L\left(q_{\boldsymbol{x}}^{1}, f(\boldsymbol{x}, \omega)\right)+\sqrt{K L\left(q_{\boldsymbol{x}}^{1}, f(\boldsymbol{x}, \omega)\right)}\right\}<\tilde{\epsilon}
\]

Thus, from the first part of this proof, it follows that
\[
\begin{array}{r}
P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}} K L\left(q_{\boldsymbol{x}}, f(\boldsymbol{x}, \omega)\right)<\tilde{\epsilon}\right\} \geq P\left\{\omega \in \Omega: \sup _{\boldsymbol{x} \in \mathscr{X}} K L\left(q_{\boldsymbol{x}}^{1}, f(\boldsymbol{x}, \omega)\right)<\epsilon^{\prime}\right\} \\
>0
\end{array}
\]
which completes the proof of the theorem.

\section*{D.5. Proof of Theorem 5}

The following Lemma is used in the proofs of continuity and association structure of the processes. This Lemma states that equicontinuous families of functions preserve a.s. continuity and convergence in distribution of stochastic processes.

Lemma 3. Let \(\mathscr{F}=\left\{f_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) be a set of known bijective continuous functions such that for every \(\boldsymbol{x} \in \mathscr{X}, f_{\boldsymbol{x}}: \mathbb{R}^{m} \longrightarrow \mathbb{R}^{n}\) is such that for every \(\boldsymbol{a} \in \mathbb{R}^{m}, f_{\boldsymbol{x}}(\boldsymbol{a})\) is a continuous functions of \(\boldsymbol{x}\). In addition asume that \(\mathscr{F}\) is an equicontinuous family of functions of \(\boldsymbol{a}\) or \(\left\{\boldsymbol{x} \mapsto f_{\boldsymbol{x}}(\boldsymbol{a}): \boldsymbol{a} \in \mathbb{R}^{m}, f_{x} \in \mathscr{F}\right\}\) is an equicontinuous family of functions of \(\boldsymbol{x}\). Let \(g_{i}: \mathscr{X} \times \Omega \longrightarrow \mathbb{R}^{m}\), \(i \geq 1\), be stochastic processes defined on an appropiate probability space ( \(\Omega, \mathscr{A}, P\) ).
(i) If for every \(i \in \mathbb{N}\), the stochastic process \(g_{i}\) is \(P\)-a.s. continuous, then \(\boldsymbol{x} \mapsto f_{\boldsymbol{x}}\left\{g_{i}(\boldsymbol{x}, \cdot)\right\}, i \in \mathbb{N}\) is \(P\)-a.s. continuous.
(ii) Consider \(\left\{\boldsymbol{x}_{j}\right\}_{j=1}^{\infty} \subset \mathscr{X}\), such that \(\lim _{j \rightarrow \infty} \boldsymbol{x}_{j}=\boldsymbol{x}_{0} \in \mathscr{X}\). If \(g_{i}\left(\boldsymbol{x}_{j}, \cdot\right) \xrightarrow{\mathscr{L}} g_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\), as \(j \longrightarrow \infty\), then \(f_{\boldsymbol{x}_{j}}\left\{g_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \xrightarrow{\mathscr{L}} f_{\boldsymbol{x}_{0}}\left\{g_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\}\), as \(j \longrightarrow \infty\).

\section*{Proof of Lemma 3}
(i) First consider, for every \(\boldsymbol{x} \in \mathscr{X}, f_{\boldsymbol{x}}\) an equicontinuous of function of \(\boldsymbol{a}\). Consider \(\boldsymbol{x}_{0} \in \mathscr{X}\). Since \(f_{\boldsymbol{x}}\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right)\) is a continuous function of \(\boldsymbol{x}\), there exists \(\delta_{1}>0\) such that for all \(\boldsymbol{x} \in B\left(\boldsymbol{x}_{0}, \delta_{1}\right), \mid f_{\boldsymbol{x}}\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right)-\)
\(f_{\boldsymbol{x}_{0}}\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right) \mid<\epsilon / 2\). By assumption, \(g_{i}\), being a \(P\)-a.s. continuous stochastic process, implies that, for almost every \(\omega \in \Omega\), and for every \(\epsilon_{2}>0\), there exists \(\delta_{2}>0\) such that for all \(\boldsymbol{x} \in B\left(\boldsymbol{x}_{0}, \delta_{2}\right),\left|g_{i}(\boldsymbol{x}, \omega)-g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right|<\epsilon_{2}\). Hence, by the equicontinuity of \(f_{\boldsymbol{x}}\), for almost every \(\omega \in \Omega\), and every \(g_{i}(\boldsymbol{x}, \omega) \in B\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right), \epsilon_{2}\right),\left|f_{\boldsymbol{x}}\left(g_{i}(\boldsymbol{x}, \omega)\right)-f_{\boldsymbol{x}}\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right)\right|<\epsilon / 2\). Finally, considering \(\delta=\min \left\{\delta_{1}, \delta_{2}\right\}\) which does not depend on \(f_{\boldsymbol{x}}\), by the triangle inequality, it follows that for every \(\omega \in \Omega\) and every \(\boldsymbol{x} \in B\left(\boldsymbol{x}_{0}, \delta\right)\), \(\left|f_{\boldsymbol{x}}\left\{g_{i}(\boldsymbol{x}, \omega)\right\}-f_{\boldsymbol{x}_{0}}\left\{g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right\}\right|<\epsilon\), which completes this part of the proof. Now consider \(\left\{\boldsymbol{x} \mapsto f_{\boldsymbol{x}}(\boldsymbol{a}): \boldsymbol{a} \in \mathbb{R}^{m}, f_{x} \in \mathscr{F}\right\}\) an equicontinuous family of functions of \(\boldsymbol{x}\). The proof is similar to the previous. By the equicontinuity consideration, there exists \(\delta_{1}>0\) such that for every \(\boldsymbol{x} \in B\left(\boldsymbol{x}_{0}, \delta_{1}\right)\), \(\left|f_{\boldsymbol{x}}\left(g_{i}(\boldsymbol{x}, \omega)\right)-f_{\boldsymbol{x}_{0}}\left(g_{i}(\boldsymbol{x}, \omega)\right)\right|<\epsilon / 2\). Since \(g_{i}\) is a \(P\)-a.s. continuous stochastic process, for almost every \(\omega \in \Omega\), and for every \(\epsilon_{2}>0\), there exists \(\delta_{2}>0\) such that for every \(\boldsymbol{x} \in B\left(\boldsymbol{x}_{0}, \delta_{2}\right),\left|g_{i}(\boldsymbol{x}, \omega)-g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right|<\epsilon_{2}\). Due to continuity of \(f_{\boldsymbol{x}}\) as a function for \(\boldsymbol{a}\), it follows that for almost every \(\omega \in \Omega\), and every \(g_{i}(\boldsymbol{x}, \omega) \in B\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right), \epsilon_{2}\right),\left|f_{\boldsymbol{x}_{0}}\left(g_{i}(\boldsymbol{x}, \omega)\right)-f_{\boldsymbol{x}_{0}}\left(g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right)\right|<\epsilon / 2\). Finally, considering \(\delta=\min \left\{\delta_{1}, \delta_{2}\right\}\) which does not depend on \(f_{\boldsymbol{x}}\), by the triangle inequality, it follows that for every \(\omega \in \Omega\) and every \(\boldsymbol{x} \in B\left(\boldsymbol{x}_{0}, \delta\right)\), \(\left|f_{\boldsymbol{x}}\left\{g_{i}(\boldsymbol{x}, \omega)\right\}-f_{\boldsymbol{x}_{0}}\left\{g_{i}\left(\boldsymbol{x}_{0}, \omega\right)\right\}\right|<\epsilon\), which completes the proof of the first part of the lemma.
(ii) Consider \(\mathscr{F}\) an equicontinuous family of functions of \(\boldsymbol{a}\) or \(\left\{\boldsymbol{x} \mapsto f_{\boldsymbol{x}}(\boldsymbol{a})\right.\) : \(\left.\boldsymbol{a} \in \mathbb{R}^{m}, f_{x} \in \mathscr{F}\right\}\) an equicontinuous family of functions of \(\boldsymbol{x}\). If \(g_{i}\left(\boldsymbol{x}_{j}, \cdot\right) \xrightarrow{\mathscr{L}} g_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\), as \(j \longrightarrow \infty\), then by baby Skorohod's theorem (Resnick, 2019), there exist random variables \(\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\}_{j \geq 0}\) defined on the Lebesgue probability space \(([0,1], \mathscr{B}([0,1]), \lambda)\), where \(\lambda\) is the Lebesgue measure, such that for each fixed \(j \geq 0, g_{i}\left(\boldsymbol{x}_{j}, \cdot\right) \stackrel{d}{=} \tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\), and \(\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right) \longrightarrow \tilde{g}_{i}\left(\boldsymbol{x}_{0}, \cdot\right) \lambda\)-a.s. as \(j \longrightarrow \infty\). Since \(f_{\boldsymbol{x}}(\boldsymbol{a})\) is a continuous function of \(\boldsymbol{a}\), it follows that for \(\boldsymbol{x} \in \mathscr{X}, f_{\boldsymbol{x}}\left\{g_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \stackrel{d}{=} f_{\boldsymbol{x}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\}\). In particular, \(f_{\boldsymbol{x}_{j}}\left\{g_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \stackrel{d}{=} f_{\boldsymbol{x}_{j}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\}\) and \(f_{\boldsymbol{x}_{0}}\left\{g_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \stackrel{d}{=} f_{\boldsymbol{x}_{0}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\}\). Since \(\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right) \longrightarrow \tilde{g}_{i}\left(\boldsymbol{x}_{0}, \cdot\right) \lambda\)-a.s. as \(j \longrightarrow \infty\) then \(\tilde{g}_{i}\) is \(\lambda\)-a.s continuous. Therefore, by Lemma 3 part ( \(i\) ), \(f_{\boldsymbol{x}_{j}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \longrightarrow f_{\boldsymbol{x}_{0}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\} \lambda\)-a.s. as \(j \longrightarrow \infty\) which implies that \(f_{\boldsymbol{x}_{j}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \xrightarrow{\mathscr{L}} f_{\boldsymbol{x}_{0}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\}\). Thus
\[
f_{\boldsymbol{x}_{j}}\left\{g_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \stackrel{d}{=} f_{\boldsymbol{x}_{j}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{j}, \cdot\right)\right\} \stackrel{\mathscr{L}}{\longrightarrow} f_{\boldsymbol{x}_{0}}\left\{\tilde{g}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\} \stackrel{d}{=} f_{\boldsymbol{x}_{0}}\left\{g_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\},
\]
as \(j \longrightarrow \infty\), which completes proof of the lemma.
Now we provide the proof of the theorem. Firstly, assume that \(F\) is a DMBPP \(\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\). Since the elements of \(\mathscr{V}\) and \(\mathscr{H}\) are equicontinuous functions of \(\boldsymbol{x}\), and for every \(i \geq 1\), the stochastic processes \(\eta_{i}\) and \(\boldsymbol{z}_{i}\) are \(P\)-a.s. continuous, by Lemma 3 and continuous mapping theorem, it follows that \(\boldsymbol{x} \mapsto v_{\boldsymbol{x}}\left(\eta_{i}(\boldsymbol{x}, \cdot)\right)\), \(\boldsymbol{x} \mapsto w_{i}(\boldsymbol{x}, \cdot)\), and \(\boldsymbol{x} \mapsto \boldsymbol{\theta}_{i}(\boldsymbol{x}, \cdot)\) are \(P\)-a.s. continuous mappings. Now, the ceiling function being continuous from the left and having a limit from the right implies that, for \(i \geq 1\), and almost every \(\omega \in \Omega,\left\lceil k(\omega) \boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{l}, \omega\right)\right\rceil\) has a limit, as \(l \longrightarrow \infty\). Note that there exists \(M>0\) such that, for every \(\boldsymbol{y} \in \Delta_{m}, i \geq 1, \boldsymbol{x} \in \mathscr{X}\) and
\(\omega \in \Omega, \mathrm{d}\left(\boldsymbol{y} \mid \alpha\left(k(\omega),\left\lceil k(\omega) \boldsymbol{\theta}_{i}(\boldsymbol{x}, \omega)\right\rceil\right)\right) \leq M\), where \(\alpha(k, \mathbf{j})=\left(\mathbf{j}, k+m-\sum_{l=1}^{m} j_{l}\right)\), and that for every \(\boldsymbol{x} \in \mathscr{X}\) and \(\omega \in \Omega, \sum_{i=1}^{\infty} w_{i}(\boldsymbol{x}, \omega)=1\). Then by dominated convergence theorem for series, the density, w.r.t. Lebesgue measure, of \(F(\boldsymbol{x}, \cdot)\), \(f\left(\boldsymbol{x}_{l}, \omega\right)\), has a.s. a limit, say \(\tilde{f}\left(\boldsymbol{x}_{0}, \cdot\right)\). This is, for every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\operatorname{Pr}\left\{\omega \in \Omega: \lim _{l \rightarrow \infty} f\left(\boldsymbol{x}_{l}, \omega\right)(\boldsymbol{y})=\tilde{f}\left(\boldsymbol{x}_{0}, \omega\right)(\boldsymbol{y}),\right\}=1
\]

Let \(\tilde{F}(\boldsymbol{x}, \omega)\) be a probability measure with density function \(\tilde{f}(\boldsymbol{x}, \omega)\). A direct application of Scheffe's theorem implies that \(F\left(\boldsymbol{x}_{l}, \cdot\right)\) converges in total variation norm to \(\tilde{F}\left(\boldsymbol{x}_{0}, \cdot\right)\) as \(l \longrightarrow \infty\), a.s., this is,
\[
P\left\{\omega \in \Omega: \lim _{l \rightarrow \infty} \sup _{B \in \mathscr{B}\left(\Delta_{m}\right)}\left|F\left(\boldsymbol{x}_{l}, \omega\right)(B)-\tilde{F}\left(\boldsymbol{x}_{0}, \omega\right)(B)\right|=0,\right\}=1
\]
which completes the proof of the theorem when \(F\) is a DMBPP.
Now, assume that \(F\) is \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). This proof follows the same arguments as in the previous part, but the arguments related to the weights of the process are not needed. Thus, there exists a probability measure \(\tilde{F}(\boldsymbol{x}, \omega)\), such that
\[
P\left\{\omega \in \Omega: \lim _{l \rightarrow \infty} \sup _{B \in \mathscr{B}\left(\Delta_{m}\right)}\left|F\left(\boldsymbol{x}_{l}, \omega\right)(B)-\tilde{F}\left(\boldsymbol{x}_{0}, \omega\right)(B)\right|=0,\right\}=1
\]
which completes the proof of the theorem.

\section*{D.6. Proof of Theorem 6}

The proof of this theorem follows the reasoning of the proof of Theorem 6 in Barrientos et al. (2017). For completeness the proof is stated below. Here, we prove that for every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), with \(\boldsymbol{x}_{l} \in \mathscr{X}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l}=\boldsymbol{x}_{0}\),
\[
P\left\{\omega \in \Omega: \lim _{l \rightarrow \infty} \sup _{B \in \mathscr{B}\left(\Delta_{m}\right)}\left|F\left(\boldsymbol{x}_{l}, \omega\right)(B)-F\left(\boldsymbol{x}_{0}, \omega\right)(B)\right|=0,\right\}=1
\]

By assumption, for every \(i \geq 1\), the stochastic processes \(\eta_{i}\) are a.s. continuous, i.e., for every \(i \geq 1, \boldsymbol{x} \longmapsto \eta_{i}(\boldsymbol{x}, \cdot)\) is an a.s. continuous function. By Lemma 3, the equicontinuity assumption of \(\mathscr{V}\) as a function of \(\boldsymbol{x}\), and continuous mapping theorem, it follows that for every \(i \geq 1, \boldsymbol{x} \longmapsto w_{i}(\boldsymbol{x}, \cdot)\) is an a.s. continuous function. Therefore for every \(i \geq 1\) and every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l}=\boldsymbol{x}_{0}\), we have that \(\lim _{l \rightarrow \infty} w_{i}\left(\boldsymbol{x}_{l}, \cdot\right)=w_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\), a.s. Noting that there exists \(M>0\) such that, for every \(\boldsymbol{y} \in \Delta_{m}, i \geq 1\), and \(\omega \in \Omega, \mathrm{d}\left(\boldsymbol{y} \mid \alpha\left(k(\omega),\left\lceil k(\omega) \boldsymbol{\theta}_{i}(\omega)\right\rceil\right)\right) \leq M\), and that for every \(\boldsymbol{x} \in \mathscr{X}\) and \(\omega \in \Omega, \sum_{i=1}^{\infty} w_{i}(\boldsymbol{x}, \omega)=1\), dominated convergence theorem for series implies that the density, w.r.t. Lebesgue measure, of \(F(\boldsymbol{x}, \cdot)\) is a.s. continuos, i.e., for every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\operatorname{Pr}\left\{\omega \in \Omega: \lim _{l \rightarrow \infty} f\left(\boldsymbol{x}_{l}, \omega\right)(\boldsymbol{y})=f\left(\boldsymbol{x}_{0}, \omega\right)(\boldsymbol{y}),\right\}=1
\]

Finally, a direct application of Scheffe's theorem implies that \(F\left(\boldsymbol{x}_{j}, \cdot\right)\) converges in total variation norm to \(F\left(\boldsymbol{x}_{0}, \cdot\right)\) as \(j \longrightarrow \infty\), a.s., which completes the proof of the theorem.

\section*{D.7. Proof of Theorem 7}

Here we prove that for every \(\boldsymbol{y} \in \tilde{\Delta}_{m}\), every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), with \(\boldsymbol{x}_{l} \in \mathscr{X}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l}=\boldsymbol{x}_{0}\),
\[
\frac{E\left\{F\left(\boldsymbol{x}_{l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\}-E\left\{F\left(\boldsymbol{x}_{l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\} E\left\{F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\}}{\sqrt{\operatorname{Var}\left\{F\left(\boldsymbol{x}_{l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\} \operatorname{Var}\left\{F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\}}} \longrightarrow 1,
\]
as \(l \longrightarrow \infty\), where \(B_{\boldsymbol{y}}=\left[0, y_{1}\right] \times \ldots \times\left[0, y_{m}\right]\), and expectations are obtained by the law of total expectation conditioning on the degree of the polynomial. We state the complete proof for the general definition of \(F\). The proof for the simplified versions of \(F\) are straightforward. In order to reduce the notation, \(k(\omega)\) is denoted by \(k\) when necessary.

First, assume that \(F\) is a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\). Notice that for every \(\boldsymbol{y} \in \Delta_{m}\) and every \(\boldsymbol{x} \in \mathscr{X}\),
\(E\left\{F(\boldsymbol{x}, \cdot)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\}=\sum_{\mathbf{j} \in \mathcal{H}_{k_{0}, m}} E\left\{\left.F^{*}(\boldsymbol{x}, \cdot)\left(\frac{\mathbf{j}}{k}\right) \right\rvert\, k=k_{0}\right\} \operatorname{Mult}\left(\mathbf{j} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right)\),
where \(\mathcal{H}_{k, m}=\left\{\left(j_{1}, \ldots, j_{m}\right) \in\{0, \ldots, k\}^{m}: \sum_{l=1}^{m} j_{l} \leq k+m-1\right\}, \beta\left(k_{0}, \boldsymbol{y}\right)= \left(k_{0}+m-1, \boldsymbol{y}\right),(\mathbf{j} / k)=\left(j_{1} / k, \ldots, j_{m} / k\right)\), Mult \((\cdot \mid k, \boldsymbol{y})\) denotes the probability mass function of a multinomial distribution with parameters \((k, \boldsymbol{y})\), and
\[
F^{*}(\boldsymbol{x}, \cdot)\left(\frac{\mathbf{j}}{k}\right)=\sum_{i=1}^{\infty} w_{i}(\boldsymbol{x}, \cdot) \mathbb{I}\left\{\boldsymbol{\theta}_{i}(\boldsymbol{x}, \cdot)\right\}_{\left\{A_{\mathbf{j}, k}\right\}},
\]
where \(A_{\mathbf{j}, k}=\left[0, j_{1} / k\right] \times \ldots \times\left[0, j_{m} / k\right]\). Since the stochastic processes \(\left\{\eta_{i}\right\}_{i \geq 1}\) and \(\left\{\boldsymbol{z}_{i}\right\}_{i \geq 1}\) are independent and identically distributed, it follows that,
\[
\begin{aligned}
E\left\{\left.F^{*}(\boldsymbol{x}, \cdot)\left(\frac{\mathbf{j}}{k}\right) \right\rvert\, k=k_{0}\right\} & =\sum_{i=1}^{\infty} E\left\{w_{i}(\boldsymbol{x}, \cdot)\right\} E\left\{\mathbb{I}\left\{\boldsymbol{\theta}_{1}(\boldsymbol{x}, \cdot)\right\}_{\left\{A_{\mathbf{j}, k_{0}}\right\}}\right\}, \\
& =G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right),
\end{aligned}
\]
where \(G_{0, \boldsymbol{x}}(A)=G_{0}(\boldsymbol{x}, \cdot)(A)\) denotes the distribution function of \(\boldsymbol{\theta}_{1}(\boldsymbol{x}, \cdot)\) defined on \(\tilde{\Delta}_{m}\). Thus
\[
E\left\{F(\boldsymbol{x}, \cdot)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\}=\sum_{\mathbf{j} \in \mathcal{H}_{k_{0}, m}} G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right) \operatorname{Mult}\left(\mathbf{j} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right) .
\]

Noting that for every \(\boldsymbol{x}, \boldsymbol{x}_{0} \in \mathscr{X}\) and every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\begin{aligned}
& E\left\{F(\boldsymbol{x}, \cdot)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} \\
& \quad=\sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k, m}, \mathbf{j}_{2} \in \mathcal{H}_{k, m}}} E\left\{\left.F^{*}\left(\boldsymbol{x}, \boldsymbol{x}_{0}, \cdot\right)\left(\frac{\mathbf{j}_{1}}{k}, \frac{\mathbf{j}_{2}}{k}\right) \right\rvert\, k=k_{0}\right\} \times \overline{\mathrm{M}}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right),
\end{aligned}
\]
where \(\bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta(k, \boldsymbol{y})\right)=\operatorname{Mult}\left(\mathbf{j}_{1} \mid \beta(k, \boldsymbol{y})\right) \times \operatorname{Mult}\left(\mathbf{j}_{2} \mid \beta(k, \boldsymbol{y})\right)\), and
\[
\begin{aligned}
F^{*}\left(\boldsymbol{x}, \boldsymbol{x}_{0}, \cdot\right) & \left(\frac{\mathbf{j}_{1}}{k}, \frac{\mathbf{j}_{2}}{k}\right)=\sum_{i=1}^{\infty} w_{i}(\boldsymbol{x}, \cdot) w_{i}\left(\boldsymbol{x}_{0}, \cdot\right) \mathbb{I}\left\{\boldsymbol{\theta}_{i}(\boldsymbol{x}, \cdot)\right\}_{\left\{A_{\mathbf{j}_{1}}, k\right\}} \mathbb{I}\left\{\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\}_{\left\{A_{\mathbf{j}_{2}, k}\right\}}, \\
+ & \sum_{\substack{i, i_{1}=1, i \neq i_{1}}}^{\infty} w_{i}(\boldsymbol{x}, \cdot) w_{i_{1}}\left(\boldsymbol{x}_{0}, \cdot\right) \mathbb{I}\left\{\boldsymbol{\theta}_{i}(\boldsymbol{x}, \cdot)\right\}_{\left\{A_{\mathbf{j}_{1}, k}\right\}} \mathbb{I}\left\{\boldsymbol{\theta}_{i_{1}}\left(\boldsymbol{x}_{0}, \cdot\right)\right\}_{\left\{A_{\mathbf{j}_{2}, k}\right\}} \cdot
\end{aligned}
\]

Applying a similar reasoning as before, it follows that, for every \(\boldsymbol{x}, \boldsymbol{x}_{0} \in \mathscr{X}\) and every \(\boldsymbol{y} \in \Delta_{m}\),
\[
\begin{aligned}
E\left\{F(\boldsymbol{x}, \cdot)\left(B_{\boldsymbol{y}}\right) F\right. & \left.\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} \\
= & \sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}, m} \\
, \mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}}\left\{\sum_{i=1}^{\infty} E\left\{w_{i}(\boldsymbol{x}, \cdot) w_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\} G_{0, \boldsymbol{x}, \boldsymbol{x}_{0}}\left(A_{\mathbf{j}_{1}, k_{0}} \times A_{\mathbf{j}_{2}, k_{0}}\right),\right. \\
& \left.\quad+\sum_{\substack{i, i_{1}=1, i \neq i_{1}}}^{\infty} E\left\{w_{i}(\boldsymbol{x}, \cdot) w_{i_{1}}\left(\boldsymbol{x}_{0}, \cdot\right)\right\} G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}_{1}, k_{0}}\right) G_{0, \boldsymbol{x}_{0}}\left(A_{\mathbf{j}_{2}, k_{0}}\right)\right\}, \\
& \quad \times \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right),
\end{aligned}
\]
where \(G_{0, \boldsymbol{x}, \boldsymbol{x}_{0}}(A)=G_{0}\left(\left(\boldsymbol{x}, \boldsymbol{x}_{0}\right), \cdot\right)(A)\) denotes the joint distribution function of \(\left(\boldsymbol{\theta}_{i}(\boldsymbol{x}, \cdot), \boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right)\) defined on \(\tilde{\Delta}_{m}^{2}\). In particular, for \(\boldsymbol{x}=\boldsymbol{x}_{0}\),
\[
\begin{aligned}
E\left\{F(\boldsymbol{x}, \cdot)\left(B_{\boldsymbol{y}}\right)^{2} \mid\right. & \left.k=k_{0}\right\}=\sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}, m} \\
\mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}}\left\{\sum_{i=1}^{\infty} E\left\{w_{i}(\boldsymbol{x}, \cdot)^{2}\right\} G_{0, \boldsymbol{x}}\left(A_{\min \left\{\mathbf{j}_{1}, \mathbf{j}_{2}\right\}, \mathrm{k}_{0}}\right)\right. \\
& \left.+\sum_{\substack{i, i_{1}=1, i \neq i_{1}}}^{\infty} E\left\{w_{i}(\boldsymbol{x}, \cdot) w_{i_{1}}(\boldsymbol{x}, \cdot)\right\} G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}_{1}, k_{0}}\right) G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}_{2}, k_{0}}\right)\right\} \\
& \times \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right.
\end{aligned}
\]
where \(A_{\text {min }\left\{\mathbf{j}_{1}, \mathbf{j}_{2}\right\}, k}=\left[0, \min \left\{j_{11}, j_{21}\right\} / k\right] \times \ldots \times\left[0, \min \left\{j_{1 m}, j_{2 m}\right\} / k\right]\). By assumption, for every \(i \geq 1\), and every \(\left\{\boldsymbol{x}_{l}\right\}_{l=0}^{\infty}\), with \(\boldsymbol{x}_{l} \in \mathscr{X}\), such that \(\lim _{l \rightarrow \infty} \boldsymbol{x}_{l} =\boldsymbol{x}_{0}\), the processes \(\eta_{i}\left(\boldsymbol{x}_{l}, \cdot\right)\) and \(\boldsymbol{z}_{i}\left(\boldsymbol{x}_{l}, \cdot\right)\) converge in distribution to \(\eta_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\) and \(\boldsymbol{z}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\), respectively, as \(l \rightarrow \infty\). Since \(\mathscr{V}\) and \(\mathscr{H}\) are sets of equicontinuous functions of \(\boldsymbol{x}\), by Lemma 3, and continuous mapping theorem, it follows that \(w_{i}\left(\boldsymbol{x}_{l}, \cdot\right)\) converges in distribution to \(w_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\) and \(\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{l}, \cdot\right)\) converges in distribution to \(\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\), as \(l \rightarrow \infty\). Thus, for every \(a \in \tilde{\Delta}_{m}, \lim _{l \rightarrow \infty} G_{0, \boldsymbol{x}_{l}}(a)= G_{0, \boldsymbol{x}_{0}}(a)\). Noting that \(w_{i}(\boldsymbol{x}, \cdot)\) are bounded variables, Portmanteau's theorem implies that the mappings \(\boldsymbol{x} \mapsto E\left\{w_{i}(\boldsymbol{x}, \cdot)\right\}, \boldsymbol{x} \mapsto E\left\{w_{i}(\boldsymbol{x}, \cdot)^{2}\right\}\) and \(\boldsymbol{x} \mapsto E\left\{w_{i}(\boldsymbol{x}, \cdot) w_{i}\left(\boldsymbol{x}_{0}, \cdot\right)\right\}\), are continuous. Now, considering \(\boldsymbol{y} \in \tilde{\Delta}_{m}\), the above expressions and few applications of dominated convergence theorem for series, it
follows that,
\[
\begin{aligned}
\lim _{j \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{j}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\} & =\sum_{k_{0}=1}^{\infty} \lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right\} P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& =E\left\{F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\} \\
\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)^{2}\right\} & =\sum_{k_{0}=1}^{\infty} \lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)^{2} \mid k_{0}\right\} P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\} \\
& =E\left\{F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)^{2}\right\}
\end{aligned}
\]
and
\[
\begin{aligned}
\lim _{j \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{j}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\}= & \sum_{k_{0}=1}^{\infty} \lim _{j \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{j}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right\}, \\
& \times P\left\{\omega \in \Omega: k(\omega)=k_{0}\right\}, \\
= & E\left\{F\left(\boldsymbol{x}_{0}, \cdot\right)\left(B_{\boldsymbol{y}}\right)^{2}\right\} .
\end{aligned}
\]

Thus the proof is completed when \(F\) is a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{\boldsymbol{z}}, \mathscr{V}, \mathscr{H}\right)\).

\section*{D.8. Proof of Theorem 8}

The proof of this theorem is a straightforward extension of the proof of Theorem 8 in Barrientos et al. (2017). For completeness we state the proof in what follows. Here we use the law of total covariance conditioning on the degree of the polynomial.

Assume that \(F\) is a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\). By assumption, for every \(i \geq\) 1 , and every \(\left\{\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)\right\}_{l=1}^{\infty}\) with \(\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right) \in \mathscr{X}^{2}\) and \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\), such that \(\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right) \longrightarrow\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right)\), as \(l \rightarrow \infty\), the joint processes \(\left(\eta_{i}\left(\boldsymbol{x}_{1 l}, \cdot\right), \eta_{i}\left(\boldsymbol{x}_{2 l}, \cdot\right)\right)\) and \(\left(\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1 l}, \cdot\right), \boldsymbol{z}_{i}\left(\boldsymbol{x}_{2 l}, \cdot\right)\right)\) converge in distribution to the processes \(\left(\eta_{i}\left(\boldsymbol{x}_{1}, \cdot\right), \eta_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\), and \(\left(\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1}, \cdot\right), \boldsymbol{z}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\), as \(l \rightarrow \infty\), respectively. Since \(\mathscr{V}\) and \(\mathscr{H}\) are sets of equicontinuous functions of \(\boldsymbol{x}\), by Lemma 3 and continuous mapping theorem, it follows that for every \(i \geq 1,\left(w_{i}\left(\boldsymbol{x}_{1 l}, \cdot\right), w_{i}\left(\boldsymbol{x}_{2 l}, \cdot\right)\right)\) and \(\left(\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{1 l}, \cdot\right), \boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{2 l}, \cdot\right)\right)\) converge in distribution to \(\left(w_{i}\left(\boldsymbol{x}_{1}, \cdot\right), w_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\) and \(\left(\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{1}, \cdot\right), \boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\), as \(l \rightarrow \infty\), respectively. Thus, for every \(\boldsymbol{a}_{1} \in \tilde{\Delta}_{m}\) and \(\boldsymbol{a}_{2} \in \tilde{\Delta}_{m}, \lim _{l \rightarrow \infty} G_{0, \boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}}\left(\boldsymbol{a}_{1}, \boldsymbol{a}_{2}\right)= G_{0, \boldsymbol{x}_{1}, \boldsymbol{x}_{2}}\left(\boldsymbol{a}_{1}, \boldsymbol{a}_{2}\right)\), where \(G_{0, \boldsymbol{x}_{1}, \boldsymbol{x}_{2}}\) denotes the joint distribution function of \(\left(\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right.\), \(\left.\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\). Noting that for every \(\boldsymbol{x}, w_{i}(\boldsymbol{x}, \cdot)\) are bounded variables, Portmanteau's theorem implies that mappings \(\boldsymbol{x} \mapsto E\left\{w_{i}(\boldsymbol{x}, \cdot)\right\}\) and \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \mapsto E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right) w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}, i, i_{1} \in \mathbb{N}\), are continuous. In addition, for \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) such that \(\left\|\boldsymbol{x}_{1}-\boldsymbol{x}_{2}\right\|>\gamma\), the assumption
\[
\operatorname{Cov}\left[\mathbb{I}_{\left\{A_{1}\right\}}\left\{\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\}, \mathbb{I}_{\left\{A_{2}\right\}}\left\{\boldsymbol{z}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}\right]=0
\]
and
\[
\operatorname{Cov}\left[\mathbb{I}\left\{\eta_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\}_{\left\{A_{1}\right\}}, \mathbb{I}\left\{\eta_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}_{\left\{A_{2}\right\}}\right]=0
\]
imply that
\[
\begin{aligned}
& E\left\{\mathbb{I}\left\{\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\}_{\left\{A_{1}\right\}} \mathbb{I}\left\{\boldsymbol{z}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}_{\left\{A_{2}\right\}}\right\}=E\left\{\mathbb{I}\left\{\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\}_{\left\{A_{1}\right\}}\right\} E\left\{\mathbb{I}\left\{\boldsymbol{z}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}_{\left\{A_{2}\right\}}\right\} \\
& E\left\{\mathbb{I}\left\{\eta_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\}_{\left\{A_{1}\right\}} \mathbb{I}\left\{\eta_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}_{\left\{A_{2}\right\}}\right\}=E\left\{\mathbb{I}\left\{\eta_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\}_{\left\{A_{1}\right\}}\right\} E\left\{\mathbb{I}\left\{\eta_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}_{\left\{A_{2}\right\}}\right\}
\end{aligned}
\]

Therefore, for every \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) such that \(\left\|\boldsymbol{x}_{1}-\boldsymbol{x}_{2}\right\|>\gamma\), it follows that \(G_{0, \boldsymbol{x}_{1}, \boldsymbol{x}_{2}}\left(\boldsymbol{a}_{1}, \boldsymbol{a}_{2}\right)=G_{0, \boldsymbol{x}_{1}}\left(\boldsymbol{a}_{1}\right) G_{0, \boldsymbol{x}_{2}}\left(\boldsymbol{a}_{2}\right)\), and that for \(i, i_{1} \in \mathbb{N}\),
\[
E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right) w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\}=E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\} E\left\{w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\} .
\]

Now, considering the expressions from the proof of Theorem 7, for every \(\boldsymbol{y} \in \Delta_{m}\), \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) such that \(\left\|\boldsymbol{x}_{1}-\boldsymbol{x}_{2}\right\|>\gamma\), and an application of dominated convergence theorem, it follows that
\[
\begin{aligned}
\lim _{l \rightarrow \infty} E & \left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} \\
= & \sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}, m} \\
\mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}} \sum_{\substack{i=1 \\
i_{1}=1}}^{\infty} E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\} E\left\{w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\} G_{0, \boldsymbol{x}_{1}}\left(A_{\mathbf{j}_{1}, k_{0}}\right) G_{0, \boldsymbol{x}_{2}}\left(A_{\mathbf{j}_{2}, k_{0}}\right), \\
& \times \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right), \\
= & \lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} E\left\{F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} .
\end{aligned}
\]

Thus,
\[
\lim _{l \rightarrow \infty} \operatorname{Cov}\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\}=0
\]

Finally, by dominated convergence theorem for series, it follows that
\[
\begin{aligned}
\lim _{l \rightarrow \infty} \operatorname{Cov} & {\left[F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right] } \\
& =E\left\{\lim _{l \rightarrow \infty} \operatorname{Cov}\left[F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right]\right\}, \\
& +\operatorname{Cov}\left[\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right\}, \lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k(\cdot)\right\}\right], \\
& =\operatorname{Cov}\left[E\left\{F\left(\boldsymbol{x}_{1}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right\}, E\left\{F\left(\boldsymbol{x}_{2}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k(\cdot)\right\}\right],
\end{aligned}
\]
where for every \(\boldsymbol{x} \in \mathscr{X}\),
\[
E\left\{F(\boldsymbol{x}, \cdot)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\}=\sum_{\mathbf{j} \in \mathcal{H}_{k_{0}, m}} G_{0, \boldsymbol{x}}\left(A_{\mathbf{j}, k_{0}}\right) \operatorname{Mult}\left(\mathbf{j} \mid k_{0}+m-1, \boldsymbol{y}\right)
\]
which completes this part of the proof.
Assume now that \(F\) is a \(w \mathrm{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). By the same arguments as when \(F\) is the general model, for every \(\boldsymbol{y} \in \Delta_{m}\) and \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) and an application of dominated convergence theorem, it follows that
\[
\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\}
\]
\[
=\sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}, m} \\ \mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}} \sum_{\substack{i=1, i_{1}=1}}^{\infty} E\left\{w_{i}(\cdot) w_{i_{1}}(\cdot)\right\} G_{0, \boldsymbol{x}_{1}}\left(A_{\mathbf{j}_{1}, k_{0}}\right) G_{0, \boldsymbol{x}_{2}}\left(A_{\mathbf{j}_{2}, k_{0}}\right) \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right),
\]
and
\[
\begin{aligned}
\lim _{l \rightarrow \infty} & E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} E\left\{F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} \\
= & \sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}}, m \\
\mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}} \sum_{i=1,}^{\infty} E\left\{w_{i}(\cdot)\right\} E\left\{w_{i_{1}}(\cdot)\right\} G_{0, \boldsymbol{x}_{1}}\left(A_{\mathbf{j}_{1}, k_{0}}\right) G_{0, \boldsymbol{x}_{2}}\left(A_{\mathbf{j}_{2}, k_{0}}\right) \\
& \times \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right)
\end{aligned}
\]

Since \(\operatorname{Cov}\left[\sum_{i=1}^{\infty} w_{i}(\omega), \sum_{i_{1}=1}^{\infty} w_{i_{1}}(\omega)\right]=0\), it follows that
\[
\lim _{l \rightarrow \infty} \operatorname{Cov}\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right),\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\}=0
\]

Finally, the proof is completed using the same arguments as in the first part.

\section*{D.9. Proof of Theorem 9}

The proof of this theorem is a straightforward extension of the proof of Theorem 9 in Barrientos et al. (2017). For completeness we state the proof in what follows. We use the law of total covariance conditioning on the degree of the polynomial. Assume that \(F\) is a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\). By the same arguments as in the proof of the first part of Theorem 8, for every \(\boldsymbol{y} \in \Delta_{m}\) and \(\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\) such that \(\left\|\boldsymbol{x}_{1}-\boldsymbol{x}_{2}\right\|>\gamma\), and few applications of dominated convergence theorem, it follows that
\[
\begin{aligned}
\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\right. & \left.\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} \\
= & \sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}}, m \\
\mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}} \sum_{\substack{i=1 \\
i_{1}=1}}^{\infty} E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\} E\left\{w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\} \\
& \quad \times E\left\{\mathbb{I}\left\{\boldsymbol{\theta}_{i}(\cdot)\right\}_{\left\{A_{\mathbf{j}_{1}, k_{0}}\right\}} \mathbb{I}\left\{\boldsymbol{\theta}_{i 1}(\cdot)\right\}_{\left\{A_{\mathbf{j}_{2}, k_{0}}\right\}}\right\} \times \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right)
\end{aligned}
\]
and
\[
\begin{aligned}
\lim _{l \rightarrow \infty} & E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} E\left\{F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k=k_{0}\right\} \\
= & \sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}}, m \\
\mathbf{j}_{2} \in \mathcal{H}_{k_{0}}, m}} \sum_{i=1}^{\infty} E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\} E\left\{w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\} \\
& \times E\left\{\mathbb{I}\left\{\boldsymbol{\theta}_{i}(\cdot)\right\}_{\left\{A_{\mathbf{j}_{1}, k_{0}}\right\}}\right\} E\left\{\mathbb{I}\left\{\boldsymbol{\theta}_{i 1}(\cdot)\right\}_{\left\{A_{\mathbf{j}_{2}, k_{0}}\right\}}\right\} \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid \beta\left(k_{0}, \boldsymbol{y}\right)\right)
\end{aligned}
\]

Since \(\left\{\boldsymbol{\theta}_{i}\right\}_{i \geq 1}\) are independent, then
\[
\lim _{l \rightarrow \infty} \operatorname{Cov}\left[F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right]
\]
\[
\begin{aligned}
= & \sum_{\substack{\mathbf{j}_{1} \in \mathcal{H}_{k_{0}, m} \\
\mathbf{j}_{2} \in \mathcal{H}_{k_{0}, m}}} \sum_{i=1}^{\infty} E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right\} E\left\{w_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right\} \operatorname{Cov}\left\{\mathbb{I}\left\{\boldsymbol{\theta}_{i}(\cdot)\right\}_{\left\{A_{\mathbf{j}_{1}, k_{0}}\right\}}, \mathbb{I}\left\{\boldsymbol{\theta}_{i}(\cdot)\right\}_{\left\{A_{\mathbf{j}_{2}, k_{0}}\right\}}\right\}, \\
& \times \bar{M}\left(\mathbf{j}_{1}, \mathbf{j}_{2} \mid k_{0}+m-1, \boldsymbol{y}\right) .
\end{aligned}
\]

Finally, by dominated convergence theorem, it follows that, for every \(\boldsymbol{y} \in \tilde{\Delta}_{m}\),
\[
\begin{aligned}
\lim _{l \rightarrow \infty} \operatorname{Cov}[ & \left.F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right] \\
= & E\left\{\lim _{l \rightarrow \infty} \operatorname{Cov}\left[F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right]\right\}, \\
& \quad+\operatorname{Cov}\left[\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right\}, \lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k(\cdot)\right\}\right], \\
= & \sum_{k_{0}=1}^{\infty} \lim _{l \rightarrow \infty} \operatorname{Cov}\left[F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right), F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right] P\left(\left\{\omega \in \Omega: k(\omega)=k_{0}\right\}\right), \\
& \quad+\operatorname{Cov}\left[E\left\{F\left(\boldsymbol{x}_{1}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k_{0}\right\}, E\left\{F\left(\boldsymbol{x}_{2}, \cdot\right)\left(B_{\boldsymbol{y}}\right) \mid k(\cdot)\right\}\right],
\end{aligned}
\]
which completes the proof of the theorem.

\section*{D.10. Proof of Theorem 10}

The proof of this theorem is an extension of the proof of Theorem 10 in Barrientos et al. (2017). For completeness we state the proof in what follows. We prove this theorem using the definition of correlation. Expectations are obtained by the law of total expectation, conditioning on the degree of the polynomial. Assume that \(F\) is a \(\operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \boldsymbol{\Psi}_{z}, \mathscr{V}, \mathscr{H}\right)\), a \(\theta \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{\eta}, \mathscr{V}, \boldsymbol{\Psi}_{\boldsymbol{\theta}}\right)\) or a \(w \operatorname{DMBPP}\left(\lambda, \boldsymbol{\Psi}_{v}, \boldsymbol{\Psi}_{z}, \mathscr{H}\right)\). By assumption, for every \(i \geq 1\), and every \(\left\{\left(\boldsymbol{x}_{1 l}\right.\right.\), \(\left.\left.\boldsymbol{x}_{2 l}\right)\right\}_{l=1}^{\infty}\), with \(\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right) \in \mathscr{X}^{2},\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \in \mathscr{X}^{2}\), such that
\[
\lim _{l \rightarrow \infty}\left(\boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}\right)=\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right),
\]
the joint processes \(\left(\eta_{i}\left(\boldsymbol{x}_{1 l}, \cdot\right), \eta_{i}\left(\boldsymbol{x}_{2 l}, \cdot\right)\right)\) and \(\left(\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1 l}, \cdot\right), \boldsymbol{z}_{i}\left(\boldsymbol{x}_{2 l}, \cdot\right)\right)\) converge in distribution to \(\left(\eta_{i}\left(\boldsymbol{x}_{1}, \cdot\right), \eta_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\), and \(\left(\boldsymbol{z}_{i}\left(\boldsymbol{x}_{1}, \cdot\right), \boldsymbol{z}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\), as \(l \rightarrow \infty\), respectively. By the same arguments used in the proof of the first part of Theorem 8, it follows that for every \(\boldsymbol{a}_{1} \in \tilde{\Delta}_{m}\) and \(\boldsymbol{a}_{2} \in \tilde{\Delta}_{m}, \lim _{l \rightarrow \infty} G_{0, \boldsymbol{x}_{1 l}, \boldsymbol{x}_{2 l}}\left(\boldsymbol{a}_{1}, \boldsymbol{a}_{2}\right)= G_{0, \boldsymbol{x}_{1}, \boldsymbol{x}_{2}}\left(\boldsymbol{a}_{1}, \boldsymbol{a}_{2}\right)\), where \(G_{0, \boldsymbol{x}_{1}, \boldsymbol{x}_{2}}\) denotes the joint distribution function of \(\left(\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{1}, \cdot\right)\right.\), \(\left.\boldsymbol{\theta}_{i}\left(\boldsymbol{x}_{2}, \cdot\right)\right)\), and mappings \(\boldsymbol{x} \mapsto E\left\{w_{i}(\boldsymbol{x}, \cdot)\right\}\) and
\[
\left(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}\right) \mapsto E\left\{w_{i}\left(\boldsymbol{x}_{1}, \cdot\right) w_{i_{1}}\left(\boldsymbol{x}_{2}, \cdot\right)\right\},
\]
\(i, i_{1} \in \mathbb{N}\), are continuous. By few applications of dominated convergence theorem, it follows that for \(m=1,2\),
\[
\begin{aligned}
\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{m l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\} & =E\left\{F\left(\boldsymbol{x}_{m}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\} \\
\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{m l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)^{2}\right\} & =E\left\{F\left(\boldsymbol{x}_{m}, \cdot\right)\left(B_{\boldsymbol{y}}\right)^{2}\right\}
\end{aligned}
\]
and
\[
\lim _{l \rightarrow \infty} E\left\{F\left(\boldsymbol{x}_{1 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{2 l}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\}=E\left\{F\left(\boldsymbol{x}_{1}, \cdot\right)\left(B_{\boldsymbol{y}}\right) F\left(\boldsymbol{x}_{2}, \cdot\right)\left(B_{\boldsymbol{y}}\right)\right\}
\]

Finally, for \(\boldsymbol{y} \in \tilde{\Delta}_{m}\) and by the definition of correlation, the proof of the theorem is completed.

\section*{D.11. Proof of Theorem 11}

The proof of this theorem is an extension of the proof of Theorem 11 in Barrientos et al. (2017). For completeness we state the proof in what follows. Let \(m(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) g_{\boldsymbol{x}}(\boldsymbol{y})\) be the random joint distribution for the response and predictors arising when \(\left\{g_{\boldsymbol{x}}(\boldsymbol{y}): \boldsymbol{x} \in \mathscr{X}\right\}\) is a DMBPP, \(w \mathrm{DMBPP}\) or \(\theta \mathrm{DMBPP}\). Since the KL divergence between \(m_{0}\) and the implied joint distribution \(m\) can be bounded by the supremum over the predictor space of KL divergences between the predictor-dependent probability measures,
\[
\begin{aligned}
\mathrm{KL}\left(m_{0}, m\right) & =\int_{\mathscr{X}} \int_{\Delta_{m}} m_{0}(\boldsymbol{y}, \boldsymbol{x}) \log \left(\frac{m_{0}(\boldsymbol{y}, \boldsymbol{x})}{m(\boldsymbol{y}, \boldsymbol{x})}\right) d \boldsymbol{y} d \boldsymbol{x} \\
& =\int_{\mathscr{X}} q(\boldsymbol{x}) \int_{\Delta_{m}} q_{0}(\boldsymbol{y} \mid \boldsymbol{x}) \log \left(\frac{q_{0}(\boldsymbol{y} \mid \boldsymbol{x})}{g_{\boldsymbol{x}}(\boldsymbol{y})}\right) d \boldsymbol{y} d \boldsymbol{x} \\
& \leq \sup _{\boldsymbol{x} \in \mathscr{X}} \int_{\Delta_{m}} q_{0}(\boldsymbol{y} \mid \boldsymbol{x}) \log \left(\frac{q_{0}(\boldsymbol{y} \mid \boldsymbol{x})}{g_{\boldsymbol{x}}(\boldsymbol{y})}\right) d \boldsymbol{y}
\end{aligned}
\]
when \(\boldsymbol{x}\) contains only continuous predictors, it follows that, for every \(\delta>0\),
\(\operatorname{Pr}\left\{\mathrm{KL}\left(m_{0}, m\right)<\delta\right\} \geq \operatorname{Pr}\left\{\sup _{\boldsymbol{x} \in \mathscr{X}} \int_{\Delta_{m}} q_{0}(\boldsymbol{y} \mid \boldsymbol{x}) \log \left(\frac{q_{0}(\boldsymbol{y} \mid \boldsymbol{x})}{g_{\boldsymbol{x}}(\boldsymbol{y})}\right) d \boldsymbol{y}<\delta\right\}>0\),
under the assumptions of Theorem 4. Thus, by Schwartz's theorem (Schwartz, 1965), it follows that the posterior distribution associated with the random joint distribution induced by any of the proposed models is weakly consistent, that is, the posterior measure of any weak neighborhood, of any joint distribution of the form \(m_{0}(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) q_{0}(\boldsymbol{y} \mid \boldsymbol{x})\), converges to one as the sample size goes to infinity.

\section*{D.12. Proof of Theorem 12}

The proof of this theorem is a direct extension and follows the same arguments of the proof of Theorems 12, in Barrientos et al. (2017). Here the sequence of sieves is given by
\[
\mathscr{F}_{n}=\left\{\tilde{m}: \tilde{m}(\boldsymbol{y}, \boldsymbol{x})=q(\boldsymbol{x}) \sum_{j=1}^{\infty} w_{j}(\boldsymbol{x}) \tilde{\mathrm{d}}_{j}(\boldsymbol{y}), \quad\left\{\tilde{\mathrm{d}}_{j}\right\}_{j=1}^{m_{n}} \in \tilde{\mathcal{B}}_{k_{n}}^{m_{n}}\right.
\]
\[
\left.\eta_{j} \in \mathcal{B}_{j, n}, j=1, \ldots, m_{n}, \quad \sup _{\boldsymbol{x} \in \mathscr{X}} \sum_{j=m_{n}+1}^{\infty} w_{j}(\boldsymbol{x})<\epsilon\right\}
\]
where
\[
\tilde{\mathcal{B}}_{k_{n}}=\left\{\operatorname{dir}\left(\cdot \mid \mathbf{j}, \bar{k}+m-\|\mathbf{j}\|_{1}\right): \mathbf{j} \in \mathcal{H}_{\bar{k}, m}^{0}, \bar{k}=1, \ldots, k_{n}\right\}
\]

Finally, we only need to note that the cardinality of \(\tilde{\mathcal{B}}_{k_{n}}\) is \(\sum_{k=1}^{k_{n}} \frac{(k+m-1)!}{m!(k-1)!}= \frac{k_{n}\left(k_{n}+m\right)!}{k_{n}!(m+1)!}\).

\section*{Appendix E: The MCMC sampling scheme}

In this section we provide further details on the Markov chain Monte Carlo (MCMC) sampling scheme of the posterior distributions involved in the model. We used the multivariate slice sampler, proposed by Neal (2003), to scan the posterior distribution of the parameters \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j, l}^{z}\). Convergence of the posterior samples (not shown here) were evaluated by standard test as implemented in the CODA R library (Plummer et al., 2006) and by looking at the trace plots.

As described in the main document, the model considers a truncated version of the stick breaking representation of the predictor dependent mixing measures to a level \(N\). As usual, for every \(\boldsymbol{x} \in \mathscr{X}\), we set \(v_{N}(\boldsymbol{x})=1\) to ensure the weights to add up to one. In what follows we provide expressions for the joint posterior distribution of \(\left\{\boldsymbol{\beta}_{j}^{\eta}\right\}_{j=1}^{N-1},\left\{\boldsymbol{\beta}_{j l}^{\boldsymbol{z}}\right\}_{j=1, l=1}^{N, m}, \beta_{0 j}^{\eta}\), and \(\beta_{0 j l}^{\boldsymbol{z}}\), and describe the steps of the MCMC updating scheme.

Given data set \(\boldsymbol{Y}=\left(\boldsymbol{y}_{1}, \ldots, \boldsymbol{y}_{n}\right)^{t}\), the likelihood can be written as \(L(\boldsymbol{Y} \mid \ldots)\), with
\[
L(\boldsymbol{Y} \mid \ldots)=\prod_{i=1}^{n}\left\{\sum_{j=1}^{N} w_{j}\left(\boldsymbol{x}_{i}\right) \operatorname{dir}\left(\boldsymbol{y}_{i} \mid\left\lceil k \boldsymbol{\theta}_{j}\left(\boldsymbol{x}_{i}\right)\right\rceil, k+m-\left\|\left\lceil k \boldsymbol{\theta}_{j}\left(\boldsymbol{x}_{i}\right)\right\rceil\right\|_{1}\right)\right\}
\]

The full conditional distributions of \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j l}^{z}\) are given by
\[
\begin{aligned}
\pi\left(\boldsymbol{\beta}_{j}^{\eta} \mid \ldots\right) & \propto L(\boldsymbol{Y} \mid \ldots) \times\left[N_{p}\left(\boldsymbol{\beta}_{j}^{\eta} \mid \mathbf{0}, \Sigma_{1}^{\eta}\right)^{1-\gamma^{\eta}} \times N_{p}\left(\boldsymbol{\beta}_{j}^{\eta} \mid \mathbf{0}, \Sigma_{2}^{\eta}\right)^{\gamma^{\eta}}\right] \\
\pi\left(\boldsymbol{\beta}_{j, l}^{z} \mid \ldots\right) & \propto L(\boldsymbol{Y} \mid \ldots) \times\left[N_{p}\left(\boldsymbol{\beta}_{j, l}^{z} \mid \mathbf{0}, \Sigma_{1}^{z}\right)^{1-\gamma^{z}} \times N_{p}\left(\boldsymbol{\beta}_{j, l}^{z} \mid \mathbf{0}, \Sigma_{2}^{z}\right)^{\gamma^{z}}\right]
\end{aligned}
\]

The full conditional distributions of the degree of the polynomial, \(k\), is given by
\[
\pi(k \mid \ldots) \propto L(\boldsymbol{Y} \mid \ldots) \times \operatorname{Poisson}(k \mid \lambda) \mathbb{I}_{\{k \geq 1\}}
\]

Parameters \(\gamma^{\eta}\) and \(\gamma^{z}\) can be sampled from its conjugate posterior distribution. Note that \(\left(\gamma^{\eta}, \gamma^{z}\right) \mid \ldots \sim \operatorname{Discrete}\left(w_{1}, w_{2}, w_{3}, w_{4}\right)\), where
\[
w_{l}=\tilde{w}_{l} / \sum_{j=1}^{4} \tilde{w}_{j}
\]
\(l=1, \ldots, 4\), with
\[
\begin{aligned}
& \tilde{w}_{1} \propto \prod_{j=1}^{N-1} N_{p}\left(\boldsymbol{\beta}_{j}^{\eta} \mid \mathbf{0}, \Sigma_{1}^{\eta}\right) \times \prod_{j=1}^{N} \prod_{l=1}^{d} N_{p}\left(\boldsymbol{\beta}_{j, l}^{z} \mid \mathbf{0}, \Sigma_{1}^{z}\right) \times \pi_{1} \\
& \tilde{w}_{2} \propto \prod_{j=1}^{N-1} N_{p}\left(\boldsymbol{\beta}_{j}^{\eta} \mid \mathbf{0}, \Sigma_{2}^{\eta}\right) \times \prod_{j=1}^{N} \prod_{l=1}^{d} N_{p}\left(\boldsymbol{\beta}_{j, l}^{z} \mid \mathbf{0}, \Sigma_{1}^{z}\right) \times \pi_{2} \\
& \tilde{w}_{3} \propto \prod_{j=1}^{N-1} N_{p}\left(\boldsymbol{\beta}_{j}^{\eta} \mid \mathbf{0}, \Sigma_{1}^{\eta}\right) \times \prod_{j=1}^{N} \prod_{l=1}^{d} N_{p}\left(\boldsymbol{\beta}_{j, l}^{z} \mid \mathbf{0}, \Sigma_{2}^{z}\right) \times \pi_{3} \\
& \tilde{w}_{4} \propto \prod_{j=1}^{N-1} N_{p}\left(\boldsymbol{\beta}_{j}^{\eta} \mid \mathbf{0}, \Sigma_{2}^{\eta}\right) \times \prod_{j=1}^{N} \prod_{l=1}^{d} N_{p}\left(\boldsymbol{\beta}_{j, l}^{z} \mid \mathbf{0}, \Sigma_{2}^{z}\right) \times \pi_{4}
\end{aligned}
\]
where \(\pi_{1}, \pi_{2}, \pi_{3}\), and \(\pi_{4}\) denote the a priori probability of the binary pairs \((1,1)\), \((0,1),(1,0)\), and \((0,0)\).

\section*{Appendix F: Model specification for the simulation study}

As mentioned in the main document, the selection of the parameters \(\tau_{l}^{\eta}\) and \(\tau_{l}^{z}\), for \(l=1,2\), play a key role when selecting the version of the model that best fits the data. Recall that
\[
\begin{aligned}
& \boldsymbol{\beta}_{j}^{\eta} \mid \gamma^{\eta} \stackrel{i i d}{\sim} N_{p}\left(\mathbf{0}, \tau_{1}^{\eta}\left(\mathbb{X}^{t} \mathbb{X}\right)^{-1}\right)^{1-\gamma^{\eta}} \times N_{p}\left(\mathbf{0}, \tau_{2}^{\eta}\left(\mathbb{X}^{t} \mathbb{X}\right)^{-1}\right)^{\gamma^{\eta}} \\
& \boldsymbol{\beta}_{j l}^{\boldsymbol{z}} \mid \gamma^{\boldsymbol{z}} \stackrel{i i d}{\sim} N_{p}\left(\mathbf{0}, \tau_{1}^{\boldsymbol{z}}\left(\mathbb{X}^{t} \mathbb{X}\right)^{-1}\right)^{1-\gamma^{\boldsymbol{z}}} \times \gamma^{\boldsymbol{z}} N_{p}\left(\mathbf{0}, \tau_{2}^{\boldsymbol{z}}\left(\mathbb{X}^{t} \mathbb{X}\right)^{-1}\right)^{\gamma^{\boldsymbol{z}}}
\end{aligned}
\]
for \(j \geq 1\) and \(l=1, \ldots, m\). The prior distributions on \(\boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{\beta}_{j l}^{z}\) induce prior distributions for the bounded stochastic processes defining the predictor dependent weights, \(e^{\boldsymbol{\beta}_{0 j}^{\eta}+\boldsymbol{x}^{\boldsymbol{t}} \boldsymbol{\beta}_{j}^{\eta}} /\left(1+e^{\boldsymbol{\beta}_{0 j}^{\eta}+\boldsymbol{x}^{t} \boldsymbol{\beta}_{j}^{\eta}}\right)\), and the predictor dependent atoms, \(\left(e^{\boldsymbol{\beta}_{0 j 1}^{\boldsymbol{z}}+\boldsymbol{x}^{t} \boldsymbol{\beta}_{j 1}^{\boldsymbol{z}}}, \ldots, e^{\boldsymbol{\beta}_{0 j m}^{\boldsymbol{z}}+\boldsymbol{x}^{t} \boldsymbol{\beta}_{j m}^{\boldsymbol{z}}}\right) /\left(1+\sum_{l=1}^{m} e^{\boldsymbol{\beta}_{0 j l}^{\boldsymbol{z}}+\boldsymbol{x}^{t} \boldsymbol{\beta}_{j l}^{\boldsymbol{z}}}\right)\). We aim to choose the parameters \(\tau_{1}^{\eta}\) and \(\tau_{1}^{\boldsymbol{z}}\) such that with high probability \(\boldsymbol{x}^{t} \boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{x}^{t} \beta_{j l}^{\boldsymbol{z}}\) are close to zero and \(\tau_{2}^{\eta}\) and \(\tau_{2}^{\boldsymbol{z}}\) such that with high probability \(\boldsymbol{x}^{t} \boldsymbol{\beta}_{j}^{\eta}\) and \(\boldsymbol{x}^{t} \beta_{j l}^{\boldsymbol{z}}\) are away from zero.

Let \(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{G}\) be a grid of values for the [ 0,1 ]-valued predictor that was considered in the simulation study. For each \(\boldsymbol{x}_{g}, g=1, \ldots, G\), we define a sequence of values \(\tilde{\tau}_{1}, \ldots, \tilde{\tau}_{L}\) such that with high probability \(\boldsymbol{x}_{g}^{t} \boldsymbol{\beta}_{j}^{\eta}\) lies between -4 and 4 . We choose \(\tau_{1 g}^{\eta}\) as the largest \(\tilde{\tau}_{l}\) such that with high probability \(\boldsymbol{x}_{g}^{t} \boldsymbol{\beta}_{j}^{\eta}\) is between -0.2 and 0.2 and we choose \(\tau_{2 g}^{\eta}\) as the smallest \(\tilde{\tau}_{l}>\tau_{1 g}^{\eta}\) such that with high probability \(\boldsymbol{x}_{g}^{t} \boldsymbol{\beta}_{j}^{\eta}\) is between -2.2 and 2.2. Finally, we set \(\tau_{1}^{\eta}=\min \left\{\tau_{11}^{\eta}, \ldots, \tau_{1 G}^{\eta}\right\}\) and \(\tau_{2}^{\eta}=\min \left\{\tau_{21}^{\eta}, \ldots, \tau_{2 G}^{\eta}\right\}\). We follow a very similar procedure to choose \(\tau_{1}^{z}\) and \(\tau_{2}^{z}\).

\section*{Appendix G: Additional results for the simulation study}

In this section we provide additional results for the simulation study. More specifically, we provide the results measuring the performance of the model based on the estimate to the \(L_{\infty}\) distance and the contour plots of the density estimates obtained under Prior II for the binary latent variables.

Table 4 shows the mean of the \(L_{\infty}\) estimates across replicates for each Scenario, sample size and prior distribution for the binary latent variables. As expected, the integrated \(L_{\infty}\) distances between the true density and estimates decrease as the sample size increases. Regarding this criteria, the model shows the best density estimation performance for Scenario III, the single-atom true model, while the worst performance is observed for Scenario I, the fully predictor dependent true model, for every sample size and for both prior distributions for the binary latent variables.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4
Mean (across Monte Carlo replicates) of the \(L_{\infty}\) distances between the truth and random measure estimates for each scenario, prior for the binary latent variables, and sample size.}
\begin{tabular}{cccccccc}
\hline \hline & \multicolumn{3}{c}{ Prior I } & & \multicolumn{3}{c}{ Prior II } \\
\cline { 2 - 4 } \cline { 6 - 8 } Scenario & \(n=250\) & \(n=500\) & \(n=1000\) & & \(n=250\) & \(n=500\) & \(n=1000\) \\
\hline I & 32.638 & 28.078 & 26.892 & & 32.709 & 28.153 & 27.047 \\
II & 19.863 & 21.504 & 24.936 & & 19.481 & 22.041 & 23.610 \\
III & 11.912 & 9.876 & 9.203 & & 11.900 & 9.919 & 9.148 \\
IV & 26.837 & 17.814 & 11.267 & & 26.825 & 17.838 & 11.089 \\
\hline \hline
\end{tabular}
\end{table}

Figures 6 to 9 display the contour plot of the conditional density estimates mean (across replicates) for each sample size, selected values of the predictor, and Prior II for \(\left(\gamma^{\eta}, \gamma^{\boldsymbol{z}}\right)\), for Scenarios I to IV, respectively.

\section*{Appendix H: Model specification for the application to solid waste data}

The selection of the parameters \(\tau_{l}^{\eta}\) and \(\tau_{l}^{\boldsymbol{z}}\), for \(l=1,2\) follow the same reasoning as for the simulation study. Here, the \([0,1]\)-valued grid for the predictor is replaced by the six possible values that the predictor, in its dummy representation, can take.

\section*{Appendix I: Additional results for application to solid waste data}

Figure 10 displays the conditional density estimates, as the posterior predictive mean, for the DMBPP model for each value of the categorical predictor "low-low", "low", "medium-low", "medium", "medium-high", and "high", under Prior II for the binary latent variables. The LPML and \(-n W A I C\) values for the DMBPP model were 778.043 and 778.413, respectively, under Prior II for \(\left(\gamma^{\eta}, \gamma^{z}\right)\).

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-53.jpg?height=1651&width=1101&top_left_y=416&top_left_x=605}
\captionsetup{labelformat=empty}
\caption{Fig 6. Simulation study - Scenario I - Prior 2: contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for \(\left(\gamma^{\eta}, \gamma^{z}\right)\). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-54.jpg?height=1643&width=1093&top_left_y=418&top_left_x=416}
\captionsetup{labelformat=empty}
\caption{Fig 7. Simulation study - Scenario II - Prior II: Contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for \(\left(\gamma^{\eta}, \gamma^{z}\right)\). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-55.jpg?height=1651&width=1105&top_left_y=416&top_left_x=601}
\captionsetup{labelformat=empty}
\caption{Fig 8. Simulation study - Scenario III - Prior II: Contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for \(\left(\gamma^{\eta}, \gamma^{z}\right)\). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-56.jpg?height=1643&width=1094&top_left_y=419&top_left_x=416}
\captionsetup{labelformat=empty}
\caption{Fig 9. Simulation study - Scenario IV - Prior II: Contour plots of the true density (first row) and mean across replicates of the posterior mean of the conditional density for \(n=250\) (second row), \(n=500\) (third row), and \(n=1000\) (fourth row). The results are shown under Prior I for \(\left(\gamma^{\eta}, \gamma^{z}\right)\). Results are displayed for \(x=0.25\) (first column), \(x=0.50\) (second column), and \(x=0.75\) (third column).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/2544e0b1-391c-4c88-803e-54a244ebc847-57.jpg?height=1404&width=912&top_left_y=537&top_left_x=693}
\captionsetup{labelformat=empty}
\caption{Fig 10. Application to Solid Waste Data - Prior II: Contour plot of conditional density estimates and data points for colorredDMBPP, under Prior II for \(\left(\gamma^{\eta}, \gamma^{z}\right)\), for each value of the discrete socioeconomic level predictor, low-low (panel (a)), low (panel (b)), medium-low (panel (c)), medium (panel (d)), medium-high (panel (e)), and high (panel (f)). The x-axis and \(y\)-axis denote the proportion of food and hygienic waste, respectively.}
\end{figure}

\section*{References}

Aitchison, J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society, Series B 44 139-160. MR0676206
Atchison, J. \& Shen, S. M. (1980). Logistic-normal distributions: Some properties and uses. Biometrika 67 261-272. MR0581723
Babu, G. J. \& Chaubey, Y. P. (2006). Smooth estimation of a distribution and density function on a hypercube using Bernstein polynomials for dependent random vectors. Statistics \& Probability Letters 76 959-969. MR2270097
Banerjee, S., Carlin, B. P. \& Gelfand, A. E. (2003). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC. MR3362184
Barndorff-Nielsen, O. (1973). On M-ancillarity. Biometrika 60 447-455. MR0345255
Barndorff-Nielsen, O. (1978). Information and exponential families in statistical theory. John Wiley \& Sons, Ltd., Chichester. Wiley Series in Probability and Mathematical Statistics. MR0489333
Barrientos, A. F., Jara, A. \& Quintana, F. A. (2012). On the support of MacEachern's dependent Drichlet processes and extensions. Bayesian Analysis 7 277- 310. MR2934952
Barrientos, A. F., Jara, A. \& Quintana, F. A. (2015). Bayesian density estimation for compositional data using random Bernstein polynomials. Journal of Statistical Planning and Inference 166 116-125. MR3390138
Barrientos, A. F., Jara, A. \& Quintana, F. A. (2017). Fully nonparametric regression for bounded data using dependent bernstein polynomials. Journal of the American Statistical Association 112 806-825. MR3671772
Di Marzio, M., Panzera, A. \& Venieri, C. (2015). Non-parametric regression for compositional data. Statistical Modelling 15 113-133. MR3325749
Epifani, I. \& Lijoi, A. (2010). Nonparametric priors for vectors of survival functions. Statistica Sinica 20 1455-1484. MR2777332
Florens, J.-P., Mouchart, M. \& Rolin, J.-M. (1990). Elements of Bayesian statistics, vol. 134 of Monographs and Textbooks in Pure and Applied Mathematics. Marcel Dekker, Inc., New York. MR1051656
Geisser, S. \& Eddy, W. (1979). A predictive approach to model selection. Journal of the American Statistical Association 74 153-160. MR0529531
Gelfand, A. E. \& Dey, D. (1994). Bayesian model choice: asymptotics and exact calculations. Journal of the Royal Statistical Society, Series B 56 501514. MR1278223

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. \& Rubin, D. B. (2013). Bayesian data analysis. CRC press. MR3235677
George, E. I. \& McCulloch, R. E. (1993). Variable selection via gibbs sampling. Journal of the American Statistical Association 88 881-889.
Gueorguieva, R., Rosenheck, R. \& Zelterman, D. (2008). Dirichlet component regression and its applications to psychiatric data. Computational Statistics \& Data Analysis 52 5344-5355. MR2526600
Hijazi, R. H. (2003). Analysis of compositional data using Dirichlet covariate models. Ph.D. thesis, American University. MR2704518

Hijazi, R. H. \& Jernigan, R. W. (2009). Modelling compositional data using Dirichlet regression models. Journal of Applied Probability \& Statistics 47791. MR2668780

Ishwaran, H. \& James, L. F. (2001). Gibbs sampling methods for stickbreaking priors. Journal of the American Statistical Association 96 161-173. MR1952729
Jara, A. \& Hanson, T. (2011). A class of mixtures of dependent tail-free processes. Biometrika 98 553-566. MR2836406
Karabatsos, G. \& Walker, S. G. (2012). Adaptive-modal Bayesian nonparametric regression. Electronic Journal of Statistics 6 2038-2068. MR3020256
Klinger, R., Olaya, J., Marmolejo, L. \& Madera, C. (2009). A sampling plan for residentially generated solid waste quantification at urban zones of middle sized cities. Revista Facultad de Ingeniería Universidad de Antioquia 4876-86.
Leisen, F. \& Lijoi, A. (2011). Vectors of two-parameter Poisson-Dirichlet processes. Journal of Multivariate Analysis 102 482-495. MR2755010
Lijoi, A., Nipoti, B. \& Prünster, I. (2014). Bayesian inference with dependent normalized completely random measures. Bernoulli 20 1260-1291. MR3217444
MacEachern, S. N. (1999). Dependent nonparametric processes. In ASA Proceedings of the Section on Bayesian Statistical Science, Alexandria, VA. American Statistical Association.
MacEachern, S. N. (2000). Dependent Dirichlet processes. Tech. rep., Department of Statistics, The Ohio State University.
Müller, P., Quintana, F. A., Jara, A. \& E, H. T. (2015). Bayesian Nonparametric Data Analysis. New York, USA: Springer. MR3309338
Neal, R. (2003). Slice sampling. The Annals of Statistics 31 705-767. MR1994729
Ouimet, F. (2021). Asymptotic properties of bernstein estimators on the simplex. Journal of Multivariate Analysis 185 104784. MR4287788
Parthasarathy, K. R. (1967). Probability Measures in Metric Spaces. Providence, USA: AMS Chelsea Publishing. MR2169627
Pati, D., Dunson, D. B. \& Tokdar, S. T. (2013). Posterior consistency in conditional distribution estimation. Journal of Multivariate Analysis 116 456-472. MR3049916
Petrone, S. (1999a). Bayesian density estimation using Bernstein polynomials. The Canadian Journal of Statistics 27 105-126. MR1703623
Petrone, S. (1999b). Random Bernstein polynomials. Scandinavian Journal of Statistics 26 373-393. MR1712051
Petrone, S. \& Wasserman, L. (2002). Consistency of Bernstein polynomial posterior. Journal of the Royal Statistical Society, Series B 6479-100. MR1881846
Plummer, M., Best, N., Cowles, K. \& Vines, K. (2006). CODA: convergence diagnosis and output analysis for MCMC. R news 67-11.
Quintana, F., Müller, P., Jara, A. \& MacEachern, S. (2022). The dependent Dirichlet process and related models. Statistical Science 37 24-41.

MR4371095
Resnick, S. (2019). A probability path. Springer. MR3135152
Schwartz, L. (1965). On Bayes procedures. Z. Wahrscheinlichkeitstheorie und Verw. Gebiete 4 10-26. MR0184378
Shimizu, T. K., Louzada, F., Suzuki, A. K. \& Ehlers, R. S. (2021). Modeling compositional regression with uncorrelated and correlated errors: A Bayesian approach. Journal of Data Science 16 221-250.
Smithson, M. \& Verkuilen, J. (2006). A better lemon squeezer? maximumlikelihood regression with beta-distributed dependent variables. Psychological methods 1154 .
Tenbusch, A. (1994). Two-dimensional Bernstein polynomial density estimators. Metrika 41 233--253. MR1293514
Tierney, L. (1994). Markov chains for exploring posterior distributions. The Annals of Statistics 22 1701-1762. MR1329166
Trippa, L., Müller, P. \& Johnson, W. (2011). The multivariate beta process and an extension of the Polya tree model. Biometrika 98 17-34. MR2804207
Tsagris, M., Alenazi, A. \& Stewart, C. (2020). Non-parametric regression models for compositional data. arXiv preprint arXiv:2002.05137. MR3789428
Van der Merwe, S. (2019). A method for Bayesian regression modelling of composition data. South African Statistical Journal 5355-64. MR3966355
Wang, H., Meng, J. \& Tenenhaus, M. (2010). Regression modelling analysis on compositional data. In Handbook of Partial Least Squares. Springer, 381406. MR2762258

Watanabe, S. \& Opper, M. (2010). Asymptotic equivalence of bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research 11. MR2756194
Zellner, A. (1983). Applications of Bayesian analysis in econometrics. Journal of the Royal Statistical Society, Series D 32 23-34. MR0926622
Zheng, Y., Zhu, J. \& Roy, A. (2010). Nonparametric bayesian inference for the spectral density function of a random field. Biometrika 97 238-245. MR2594432