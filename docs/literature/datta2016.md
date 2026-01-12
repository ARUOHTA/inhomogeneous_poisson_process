\title{
Hierarchical Nearest-Neighbor Gaussian Process Models for Large Geostatistical Datasets
}

\author{
Abhirup Datta, Sudipto Banerjee, Andrew O. Finley and Alan E. Gelfand
}

\begin{abstract}
Spatial process models for analyzing geostatistical data entail computations that become prohibitive as the number of spatial locations become large. This manuscript develops a class of highly scalable Nearest Neighbor Gaussian Process (NNGP) models to provide fully model-based inference for large geostatistical datasets. We establish that the NNGP is a well-defined spatial process providing legitimate finite-dimensional Gaussian densities with sparse precision matrices. We embed the NNGP as a sparsityinducing prior within a rich hierarchical modeling framework and outline how computationally efficient Markov chain Monte Carlo (MCMC) algorithms can be executed without storing or decomposing large matrices. The floating point operations (flops) per iteration of this algorithm is linear in the number of spatial locations, thereby rendering substantial scalability. We illustrate the computational and inferential benefits of the NNGP over competing methods using simulation studies and also analyze forest biomass from a massive United States Forest Inventory dataset at a scale that precludes alternative dimension-reducing methods.
\end{abstract}

Keywords: Bayesian modeling; hierarchical models; Gaussian process; Markov chain Monte Carlo; nearest neighbors; predictive process; reduced-rank models; sparse precision matrices; spatial cross-covariance functions.

\section*{1 Introduction}

With the growing capabilities of Geographical Information Systems (GIS) and user-friendly software, statisticians today routinely encounter geographically referenced datasets containing a large number of irregularly located observations on multiple variables. This has, in turn, fueled considerable interest in statistical modeling for location-referenced spatial data; see, for example, the books by Stein (1999), Moller and Waagepetersen (2003), Banerjee et al. (2014), Schabenberger and Gotway (2004), and Cressie and Wikle (2011) for a variety of methods and applications. Spatial process models introduce spatial dependence between observations using an underlying random field, \(\{w(\mathbf{s}): \mathbf{s} \in \mathcal{D}\}\), over a region of interest \(\mathcal{D}\), which is endowed with a probability law that specifies the joint distribution for any finite set of random variables. For example, a zero-centered Gaussian process ensures that \(\mathbf{w}=\left(w\left(\mathbf{s}_{1}\right), w\left(\mathbf{s}_{2}\right) \ldots, w\left(\mathbf{s}_{n}\right)\right)^{\prime} \sim N(\mathbf{O}, \mathbf{C}(\boldsymbol{\theta}))\), where \(\mathbf{C}(\boldsymbol{\theta})\) is a family of covariance matrices, indexed by an unknown set of parameters \(\boldsymbol{\theta}\). Such processes offer a rich modeling framework and are being widely deployed to help researchers comprehend complex spatial phenomena in the sciences. However, model fitting usually involves the inverse and determinant of \(\mathbf{C}(\boldsymbol{\theta})\), which typically require \(\sim n^{3}\) floating point operations (flops) and storage of the order of \(n^{2}\). These become prohibitive when \(n\) is large and \(\mathbf{C}(\boldsymbol{\theta})\) has no exploitable structure.

Broadly speaking, modeling large spatial datasets proceeds from either exploiting "lowrank" models or using sparsity. The former attempts to construct spatial processes on a lower-dimensional subspace (see, e.g., Higdon 2001; Kammann and Wand 2003; Stein |2007, 2008; Banerjee et al. 2008; Cressie and Johannesson|2008; Crainiceanu et al. 2008; Rasmussen and Williams 2005; Finley et al. 2009) by regressing the original (parent) process on its realizations over a smaller set of \(r \ll n\) locations ("knots" or "centers"). The algorithmic cost for model fitting typically decreases from \(O\left(n^{3}\right)\) to \(O\left(n r^{2}+r^{3}\right) \approx O\left(n r^{2}\right)\) flops since \(n \gg r\). However, when \(n\) is large, empirical investigations suggest that \(r\) must be fairly large to adequately approximate the parent process and the \(n r^{2}\) flops becomes exorbitant (see

Section 5.1. Furthermore, low rank models perform poorly when neighboring observations are strongly correlated and the spatial signal dominates the noise (Stein 2014). Although bias-adjusted low-rank models tend to perform better (Finley et al. 2009; Banerjee et al. 2010; Sang and Huang 2012), they increase the computational burden.

Sparse methods include covariance tapering (see, e.g., Furrer et al. 2006, Kaufman et al. 2008; Du et al. 2009; Shaby and Ruppert 2012), which introduces sparsity in \(\mathbf{C}(\boldsymbol{\theta})\) using compactly supported covariance functions. This is effective for parameter estimation and interpolation of the response ("kriging"), but it has not been fully developed or explored for more general inference on residual or latent processes. Introducing sparsity in \(\mathbf{C}(\boldsymbol{\theta})^{-1}\) is prevalent in approximating Gaussian process likelihoods using Markov random fields (e.g., \begin{tabular}{|l|l|l|l|l|}
\hline Rue and Held 2005), products of lower dimensional conditional distributions (Vecchia & 1988,
\end{tabular} 1992; Stein et al. 2004), or composite likelihoods (e.g., Bevilacqua and Gaetan 2014; Eidsvik et al. 2014). However, unlike low rank processes, these do not, necessarily, extend to new random variables at arbitrary locations. There may not be a corresponding process, which restricts inference to the estimation of spatial covariance parameters. Spatial prediction ("kriging") at arbitrary locations proceeds by imputing estimates into an interpolator derived from a different process model. This may not reflect accurate estimates of predictive uncertainty and is undesirable.

Our intended inferential contribution is to offer substantial scalability for fully processbased inference on underlying, perhaps completely unobserved, spatial processes. Moving from finite-dimensional sparse likelihoods to sparsity-inducing spatial processes can be complicated. We first introduce sparsity in finite-dimensional probability models using specified neighbor sets constructed from directed acyclic graphs. We use these sets to extend these finite-dimensional models to a valid spatial process over uncountable sets. We call this process a Nearest-Neighbor Gaussian Process (NNGP). Its finite-dimensional realizations have sparse precision matrices available in closed form. While sparsity has been effectively ex-
ploited by Vecchia (1988); Stein et al. (2004); Emory (2009); Gramacy and Apley (2014); Gramacy et al. (2014) and Stroud et al. (2014) for approximating expensive likelihoods cheaply, a fully process-based modeling and inferential framework has, hitherto, proven elusive. The NNGP fills this gap and enriches the inferential capabilities of existing methods by subsuming estimation of model parameters, prediction of outcomes and interpolation of underlying processes into one highly scalable unifying framework.

To demonstrate its full inferential capabilities, we deploy the NNGP as a sparsity-inducing prior for spatial processes in a Bayesian framework. Unlike low rank processes, the NNGP always specifies non-degenerate finite dimensional distributions making it a legitimate proper prior for random fields and is applicable to any class of distributions that support a spatial stochastic process. It can, therefore, model an underlying process that is never actually observed. The modeling provides structured dependence for random effects, e.g. intercepts or coefficients, at a second stage of specification where the first stage need not be Gaussian. We cast a multivariate NNGP within a versatile spatially-varying regression framework Gelfand et al. 2003; Banerjee et al. 2008) and conveniently obtain entire posteriors for all model parameters as well as for the spatial processes at both observed and unobserved locations. Using a forestry example, we show how the NNGP delivers process-based inference for spatiallyvarying regression models at a scale where even low-rank processes, let alone full Gaussian processes, are unimplementable even in high-performance computing environments.

Here is a brief outline. Section 2 formulates the NNGP using multivariate Gaussian processes. Section 3 outlines Bayesian estimation and prediction within a very flexible hierarchical modeling setup. Section 4 discusses alternative NNGP models and algorithms. Section 5 presents simulation studies to highlight the inferential benefits of the NNGP and also analyzes forest biomass from a massive USDA dataset. Finally, Section 6 concludes the manuscript with a brief summary and pointers toward future work.

\section*{2 Nearest-Neighbor Gaussian Process}

\subsection*{2.1 Gaussian density on sparse directed acyclic graphs}

We will consider a \(q\)-variate spatial process over \(\Re^{d}\). Let \(\mathbf{w}(\mathbf{s}) \sim G P(\mathbf{0}, \mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta}))\) denote a zero-centered \(q\)-variate Gaussian process, where \(\mathbf{w}(\mathbf{s}) \in \Re^{q}\) for all \(\mathbf{s} \in \mathcal{D} \subseteq \Re^{d}\). The process is completely specified by a valid cross-covariance function \(\mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta})\), which maps a pair of locations \(\mathbf{s}\) and \(\mathbf{t}\) in \(\mathcal{D} \times \mathcal{D}\) into a \(q \times q\) real valued matrix \(\mathbf{C}(\mathbf{s}, \mathbf{t})\) with entries \(\operatorname{cov}\left\{w_{i}(\mathbf{s}), w_{j}(\mathbf{t})\right\}\). Here, \(\boldsymbol{\theta}\) denotes the parameters associated with the cross-covariance function. Let \(\mathcal{S}= \left\{\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{k}\right\}\) be a fixed collection of distinct locations in \(\mathcal{D}\), which we call the reference set. So, \(\mathbf{w}_{\mathcal{S}} \sim N\left(\mathbf{0}, \mathbf{C}_{\mathcal{S}}(\boldsymbol{\theta})\right)\), where \(\mathbf{w}_{\mathcal{S}}=\left(\mathbf{w}\left(\mathbf{s}_{1}\right)^{\prime}, \mathbf{w}\left(\mathbf{s}_{2}\right)^{\prime}, \ldots, \mathbf{w}\left(\mathbf{s}_{k}\right)^{\prime}\right)^{\prime}\) and \(\mathbf{C}_{\mathcal{S}}(\boldsymbol{\theta})\) is a positive definite \(q k \times q k\) block matrix with \(\mathbf{C}\left(\mathbf{s}_{i}, \mathbf{s}_{j}\right)\) as its blocks. Henceforth, we write \(\mathbf{C}_{\mathcal{S}}(\boldsymbol{\theta})\) as \(\mathbf{C}_{\mathcal{S}}\), the dependence on \(\boldsymbol{\theta}\) being implicit, with similar notation for all spatial covariance matrices.

The reference set \(\mathcal{S}\) need not coincide with or be a part of the observed locations, so \(k\) need not equal \(n\), although we later show that the observed locations are a convenient practical choice for \(\mathcal{S}\). When \(k\) is large, parameter estimation becomes computationally cumbersome, perhaps even unfeasible, because it entails the inverse and determinant of \(\tilde{\mathbf{C}}_{\mathcal{S}}\). Here, we benefit from expressing the joint density of \(\mathbf{w}_{\mathcal{S}}\) as the product of conditional densities, i.e.,
\[
p\left(\mathbf{w}_{\mathcal{S}}\right)=p\left(\mathbf{w}\left(\mathbf{s}_{1}\right)\right) p\left(\mathbf{w}\left(\mathbf{s}_{2}\right) \mid \mathbf{w}\left(\mathbf{s}_{1}\right)\right) \ldots p\left(\mathbf{w}\left(\mathbf{s}_{k}\right) \mid \mathbf{w}\left(\mathbf{s}_{k-1}\right), \ldots, \mathbf{w}\left(\mathbf{s}_{1}\right)\right)
\]
and replacing the larger conditioning sets on the right hand side of (1) with smaller, carefully chosen, conditioning sets of size at most \(m\), where \(m \ll k\) (see, e.g., Vecchia 1988; Stein et al. 2004; Gramacy and Apley 2014; Gramacy et al. 2014). So, for every \(\mathbf{s}_{i} \in \mathcal{S}\), a smaller conditioning set \(N\left(\mathbf{s}_{i}\right) \subset \mathcal{S} \backslash\left\{\mathbf{s}_{i}\right\}\) is used to construct
\[
\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)=\prod_{i=1}^{k} p\left(\mathbf{w}\left(\mathbf{s}_{i}\right) \mid \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\right)
\]
where \(\mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\) is the vector formed by stacking the realizations of \(\mathbf{w}(\mathbf{s})\) over \(N\left(\mathbf{s}_{i}\right)\).
Let \(N_{\mathcal{S}}=\left\{N\left(\mathbf{s}_{i}\right) ; i=1,2, \ldots, k\right\}\) be the collection of all conditioning sets over \(\mathcal{S}\). We can view the pair \(\left\{\mathcal{S}, N_{\mathcal{S}}\right\}\) as a directed graph \(\mathcal{G}\) with \(\mathcal{S}=\left\{\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{k}\right\}\) being the set of nodes and \(N_{\mathcal{S}}\) the set of directed edges. For every two nodes \(\mathbf{s}_{i}\) and \(\mathbf{s}_{j}\), we say \(\mathbf{s}_{j}\) is a directed neighbor of \(\mathbf{s}_{i}\) if there is a directed edge from \(\mathbf{s}_{i}\) to \(\mathbf{s}_{j}\). So, \(N\left(\mathbf{s}_{i}\right)\) denotes the set of directed neighbors of \(\mathbf{s}_{i}\) and is, henceforth, referred to as the "neighbor set" for \(\mathbf{s}_{i}\). A "directed cycle" in a directed graph is a chain of nodes \(\mathbf{s}_{i_{1}}, \mathbf{s}_{i_{2}}, \ldots, \mathbf{s}_{i_{b}}\) such that \(\mathbf{s}_{i_{1}}=\mathbf{s}_{i_{b}}\) and there is a directed edge between \(\mathbf{s}_{i_{j}}\) and \(\mathbf{s}_{i_{j+1}}\) for every \(j=1,2, \ldots, b-1\). A directed graph with no directed cycles is known as a 'directed acyclic graph'.

If \(\mathcal{G}\) is a directed acyclic graph, then \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\), as defined above, is a proper multivariate joint density (see Appendix A or Lauritzen (1996) for a similar result). Starting from a joint multivariate density \(p\left(\mathbf{w}_{\mathcal{S}}\right)\), we derive a new density \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) using a directed acyclic graph \(\mathcal{G}\). While this holds for any original density \(p\left(\mathbf{w}_{\mathcal{S}}\right)\), it is especially useful in our context, where \(p\left(\mathbf{w}_{\mathcal{S}}\right)\) is a multivariate Gaussian density and \(\mathcal{G}\) is sufficiently sparse. To be precise, let \(\mathbf{C}_{N\left(\mathbf{s}_{i}\right)}\) be the covariance matrix of \(\mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\) and let \(\mathbf{C}_{\mathbf{s}_{i}, N\left(\mathbf{s}_{i}\right)}\) be the \(q \times m q\) cross-covariance matrix between the random vectors \(\mathbf{w}\left(\mathbf{s}_{i}\right)\) and \(\mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\). Standard distribution theory reveals
\[
\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)=\prod_{i=1}^{k} N\left(\mathbf{w}\left(\mathbf{s}_{i}\right) \mid \mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}, \mathbf{F}_{\mathbf{s}_{i}}\right)
\]
where \(\mathbf{B}_{\mathbf{s}_{i}}=\mathbf{C}_{\mathbf{s}_{i}, N\left(\mathbf{s}_{i}\right)} \mathbf{C}_{N\left(\mathbf{s}_{i}\right)}^{-1}\) and \(\mathbf{F}_{\mathbf{s}_{i}}=\mathbf{C}\left(\mathbf{s}_{i}, \mathbf{s}_{i}\right)-\mathbf{C}_{\mathbf{s}_{i}, N\left(\mathbf{s}_{i}\right)} \mathbf{C}_{N\left(\mathbf{s}_{i}\right)}^{-1} \mathbf{C}_{N\left(\mathbf{s}_{i}\right), \mathbf{s}_{i}}\). Appendix B shows that \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) in (3) is a multivariate Gaussian density with covariance matrix \(\tilde{\mathbf{C}}_{\mathcal{S}}\), which, obviously, is different from \(\mathbf{C}_{\mathcal{S}}\). Furthermore, if \(N\left(\mathbf{s}_{i}\right)\) has at most \(m\) members for each \(\mathbf{s}_{i}\) in \(\mathcal{S}\), where \(m \ll k\), then \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\) is sparse with at most \(k m(m+1) q^{2} / 2\) non-zero entries. Thus, for a very general class of neighboring sets, \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) defined in (2) is the joint density of a multivariate Gaussian distribution with a sparse precision matrix.

Turning to the neighbor sets, choosing \(N\left(\mathbf{s}_{i}\right)\) to be any subset of \(\left\{\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{i-1}\right\}\) ensures an acyclic \(\mathcal{G}\) and, hence, a valid probability density in (3). Several special cases exist in
likelihood approximation contexts. For example, Vecchia (1988) and Stroud et al. (2014) specified \(N\left(\mathbf{s}_{i}\right)\) to be the \(m\) nearest neighbors of \(\mathbf{s}_{i}\) among \(\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{i-1}\) with respect to Euclidean distance. Stein et al. (2004) considered nearest as well as farthest neighbors from \(\left\{\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{i-1}\right\}\). Gramacy and Apley (2014) offer greater flexibility in choosing \(N\left(\mathbf{s}_{i}\right)\), but may require several approximations to be efficient.

All of the above choices depend upon an ordering of the locations. Spatial locations are not ordered naturally, so one imposes order by, for example, ordering on one of the coordinates. Of course, any other function of the coordinates can be used to impose order. However, the aforementioned authors have cogently demonstrated that the choice of the ordering has no discernible impact on the approximation of (1) by (3). Our own simulation experiments (see Appendix C) concur with these findings; inference based upon \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) is extremely robust to the ordering of the locations. This is not entirely surprising. Clearly, whatever order we choose in (1), \(p\left(\mathbf{w}_{\mathcal{S}}\right)\) produces the full joint density. Note that we reduce (1) to (2) based upon neighbor sets constructed with respect to the specific ordering in (1). A different ordering in (1) will produce a different set of neighbors for (2). Since \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) ultimately relies upon the information borrowed from the neighbors, its effectiveness is often determined by the number of neighbors we specify and not the specific ordering.

In the following section, we will extend the density \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) to a legitimate spatial process. We remark that our subsequent development holds true for any choice of \(N\left(\mathbf{s}_{i}\right)\) that ensures an acyclic \(\mathcal{G}\). In general, identifying a "best subset" of \(m\) locations for obtaining optimal predictions for \(\mathbf{s}_{i}\) is a non-convex optimization problem, which is difficult to implement and defeats our purpose of using smaller conditioning sets to ease computations. Nevertheless, we have found Vecchia's choice of \(m\)-nearest neighbors from \(\left\{\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{i-1}\right\}\) to be simple and to perform extremely well for a wide range of simulation experiments. In what ensues, this will be our choice for \(N\left(\mathbf{s}_{i}\right)\) and the corresponding density \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) will be referred to as the 'nearest neighbor' density of \(\mathbf{w}_{\mathcal{S}}\).

\subsection*{2.2 Extension to a Gaussian Process}

Let \(\mathbf{u}\) be any location in \(\mathcal{D}\) outside \(\mathcal{S}\). Consistent with the definition of \(N\left(\mathbf{s}_{i}\right)\), let \(N(\mathbf{u})\) be the set of \(m\)-nearest neighbors of \(\mathbf{u}\) in \(\mathcal{S}\). Hence, for any finite set \(\mathcal{U}=\left\{\mathbf{u}_{1}, \mathbf{u}_{2}, \ldots, \mathbf{u}_{r}\right\}\) such that \(\mathcal{S} \cap \mathcal{U}\) is empty, we define the nearest neighbor density of \(\mathbf{w}_{\mathcal{U}}\) conditional on \(\mathbf{w}_{\mathcal{S}}\) as
\[
\tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right)=\prod_{i=1}^{r} p\left(\mathbf{w}\left(\mathbf{u}_{i}\right) \mid \mathbf{w}_{N\left(\mathbf{u}_{i}\right)}\right)
\]

This conditional density is akin to (2) except that all the neighbor sets are subsets of \(\mathcal{S}\). This ensures a proper conditional density. Indeed (2) and (4) are sufficient to describe the joint density of any finite set over the domain \(\mathcal{D}\). More precisely, if \(\mathcal{V}=\left\{\mathbf{v}_{1}, \mathbf{v}_{2}, \ldots, \mathbf{v}_{n}\right\}\) is any finite subset in \(\mathcal{D}\), then, using (4) we obtain the density of \(\mathbf{w}_{\mathcal{V}}\) as,
\[
\tilde{p}\left(\mathbf{w}_{\mathcal{V}}\right)=\int \tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right) \tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right) \prod_{\left\{\mathbf{s}_{i} \in \mathcal{S} \backslash \mathcal{V}\right\}} d\left(\mathbf{w}\left(\mathbf{s}_{i}\right)\right) \text { where } \mathcal{U}=\mathcal{V} \backslash \mathcal{S}
\]

If \(\mathcal{U}\) is empty, then (4) implies that \(\tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right)=1\) in (5). If \(\mathcal{S} \backslash \mathcal{V}\) is empty, then the integration in (5) is not needed.

These probability densities, defined on finite topologies, conform to Kolmogorov's consistency criteria and, hence, correspond to a valid spatial process over \(\mathcal{D}\) (Appendix D). So, given any original (parent) spatial process and any fixed reference set \(\mathcal{S}\), we can construct a new process over the domain \(\mathcal{D}\) using a collection of neighbor sets in \(\mathcal{S}\). We refer to this process as the 'nearest neighbor process' derived from the original parent process. If the parent process is \(G P(\mathbf{0}, \mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta}))\), then
\[
\tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right)=\prod_{i=1}^{r} N\left(\mathbf{w}\left(\mathbf{u}_{i}\right) \mid \mathbf{B}_{\mathbf{u}_{i}} \mathbf{w}_{N\left(\mathbf{u}_{i}\right)}, \mathbf{F}_{\mathbf{u}_{i}}\right)=N\left(\mathbf{B}_{\mathcal{U}} \mathbf{w}_{\mathcal{S}}, \mathbf{F}_{\mathcal{U}}\right)
\]
for any finite set \(\mathcal{U}=\left\{\mathbf{u}_{1}, \mathbf{u}_{2}, \ldots, \mathbf{u}_{r}\right\}\) in \(\mathcal{D}\) outside \(\mathcal{S}\), where \(\mathbf{B}_{\mathbf{u}_{i}}\) and \(\mathbf{F}_{\mathbf{u}_{i}}\) are defined analogous to (3) based on the neighbor sets \(N\left(\mathbf{u}_{i}\right), \mathbf{F}_{\mathcal{U}}=\operatorname{diag}\left(\mathbf{F}_{\mathbf{u}_{1}}, \mathbf{F}_{\mathbf{u}_{2}}, \ldots, \mathbf{F}_{\mathbf{u}_{r}}\right)\) and \(\mathbf{B}_{\mathcal{U}}\) is
a sparse \(n q \times k q\) matrix with each row having at most \(m q\) non-zero entries (see Appendix E).
For any finite set \(\mathcal{V}\) in \(\mathcal{D}, \tilde{p}\left(\mathbf{w}_{\mathcal{V}}\right)\) is the density of the realizations of a Gaussian Process over \(\mathcal{V}\) with cross covariance function
\[
\tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} ; \boldsymbol{\theta}\right)=\left\{\begin{array}{l}
\tilde{\mathbf{C}}_{\mathbf{s}_{i}, \mathbf{s}_{j}}, \quad \text { if } \mathbf{v}_{1}=\mathbf{s}_{i} \text { and } \mathbf{v}_{2}=\mathbf{s}_{j} \text { are both in } \mathcal{S}, \\
\mathbf{B}_{\mathbf{v}_{1}} \tilde{\mathbf{C}}_{N\left(\mathbf{v}_{2}\right), \mathbf{s}_{j}} \quad \text { if } \mathbf{v}_{1} \notin \mathcal{S} \text { and } \mathbf{v}_{2}=\mathbf{s}_{j} \in \mathcal{S}, \\
\mathbf{B}_{\mathbf{v}_{1}} \tilde{\mathbf{C}}_{N\left(\mathbf{v}_{1}\right), N\left(\mathbf{v}_{2}\right)} \mathbf{B}_{\mathbf{v}_{2}}^{\prime}+\delta_{\left(\mathbf{v}_{1}=\mathbf{v}_{2}\right)} \mathbf{F}_{\mathbf{v}_{1}} \quad \text { if } \mathbf{v}_{1} \text { and } \mathbf{v}_{2} \text { are not in } \mathcal{S}
\end{array}\right.
\]
where \(\mathbf{v}_{1}\) and \(\mathbf{v}_{2}\) are any two locations in \(\mathcal{D}, \tilde{\mathbf{C}}_{A, B}\) denotes submatrices of \(\tilde{\mathbf{C}}_{\mathcal{S}}\) indexed by the locations in the sets \(A\) and \(B\), and \(\delta_{\left(\mathbf{v}_{1}=\mathbf{v}_{2}\right)}\) is the Kronecker delta. Appendix E also shows that \(\tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} \mid \boldsymbol{\theta}\right)\) is continuous for all pairs ( \(\mathbf{v}_{1}, \mathbf{v}_{2}\) ) outside a set of Lebesgue measure zero.

This completes the construction of a well-defined Nearest Neighbor Gaussian Process, \(N N G P(\mathbf{0}, \tilde{\mathbf{C}}(\cdot, \cdot \mid \boldsymbol{\theta}))\), derived from a parent Gaussian process, \(\operatorname{GP}(\mathbf{0}, \mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta}))\). In the NNGP, the size of \(\mathcal{S}\), i.e., \(k\), can be as large, or even larger than the size of the dataset. The reduction in computational complexity is achieved through sparsity of the NNGP precision matrices. Unlike low-rank processes, the NNGP is not a degenerate process. It is a proper, sparsity-inducing Gaussian process, immediately available as a prior in hierarchical modeling, and, as we show in the next section, delivers massive computational benefits.

\section*{3 Bayesian estimation and implementation}

\subsection*{3.1 A hierarchical model}

Consider a vector of \(l\) dependent variables, say \(\mathbf{y}(\mathbf{t})\), at location \(\mathbf{t} \in \mathcal{D} \subseteq \Re^{d}\) in a spatiallyvarying regression model,
\[
\mathbf{y}(\mathbf{t})=\mathbf{X}(\mathbf{t})^{\prime} \boldsymbol{\beta}+\mathbf{Z}(\mathbf{t})^{\prime} \mathbf{w}(\mathbf{t})+\boldsymbol{\epsilon}(\mathbf{t})
\]
where \(\mathbf{X}(\mathbf{t})^{\prime}\) is the \(l \times p\) matrix of fixed spatially-referenced predictors, \(\mathbf{w}(\mathbf{t})\) is a \(q \times 1\) spatial process forming the coefficients of the \(l \times q\) fixed design matrix \(\mathbf{Z}(\mathbf{t})^{\prime}\), and \(\boldsymbol{\epsilon}(\mathbf{t}) \stackrel{i i d}{\sim} N(\mathbf{0}, \mathbf{D})\) is an \(l \times 1\) white noise process capturing measurement error or micro-scale variability with dispersion matrix \(\mathbf{D}\), which we assume is diagonal with entries \(\tau_{j}^{2}, j=1,2, \ldots, l\). The matrix \(\mathbf{X}(\mathbf{t})^{\prime}\) is block diagonal with \(p=\sum_{i=1}^{l} p_{i}\), where the \(1 \times p_{i}\) vector \(\mathbf{x}_{i}(\mathbf{t})^{\prime}\), including perhaps an intercept, is the \(i\)-th block for each \(i=1,2, \ldots, l\). The model in (8) subsumes several specific spatial models. For instance, letting \(q=l\) and \(\mathbf{Z}(\mathbf{t})^{\prime}=\mathbf{I}_{l \times l}\) leads to a multivariate spatial regression model where \(\mathbf{w}(\mathbf{t})\) acts as a spatially-varying intercept. On the other hand, we could envision all coefficients to be spatially-varying and set \(q=p\) with \(\mathbf{Z}(\mathbf{t})^{\prime}=\mathbf{X}(\mathbf{t})^{\prime}\).

For scalability, instead of a customary Gaussian process prior for \(\mathbf{w}(t)\) in (8), we assume \(\mathbf{w}(\mathbf{t}) \sim N N G P(\mathbf{0}, \tilde{\mathbf{C}}(\cdot, \cdot \mid \boldsymbol{\theta}))\) derived from the parent \(G P(\mathbf{0}, \mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta}))\). Any valid isotropic cross covariance function (see, e.g., Gelfand and Banerjee 2010) can be used to construct \(\mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta})\). To elucidate, let \(\mathcal{T}=\left\{\mathbf{t}_{1}, \mathbf{t}_{2}, \ldots, \mathbf{t}_{n}\right\}\) be the set of locations where the outcomes and predictors have been observed. This set may, but need not, intersect with the reference set \(\mathcal{S}=\left\{\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{k}\right\}\) for the NNGP. Without loss of generality, we split up \(\mathcal{T}\) into \(\mathcal{S}^{*}\) and \(\mathcal{U}\), where \(\mathcal{S}^{*}=\mathcal{S} \cap \mathcal{T}=\left\{\mathbf{s}_{i_{1}}, \mathbf{s}_{i_{2}}, \ldots, \mathbf{s}_{i_{r}}\right\}\) with \(\mathbf{s}_{i_{j}}=\mathbf{t}_{j}\) for \(j=1,2, \ldots, r\) and \(\mathcal{U}=\mathcal{T} \backslash \mathcal{S}= \left\{\mathbf{t}_{r+1}, \mathbf{t}_{r+2}, \ldots, \mathbf{t}_{n}\right\}\). Since \(\mathcal{S} \cup \mathcal{T}=\mathcal{S} \cup \mathcal{U}\), we can completely specify the realizations of the NNGP in terms of the realizations of the parent process over \(\mathcal{S}\) and \(\mathcal{U}\), hierarchically, as \(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}} \sim N\left(\mathbf{B}_{\mathcal{U}} \mathbf{w}_{\mathcal{S}}, \mathbf{F}_{\mathcal{U}}\right)\) and \(\mathbf{w}_{\mathcal{S}} \sim N\left(\mathbf{0}, \tilde{\mathbf{C}}_{\mathcal{S}}\right)\). For a full Bayesian specification, we further specify prior distributions on \(\boldsymbol{\beta}, \boldsymbol{\theta}\) and the \(\tau_{j}^{2}\) 's. For example, with customary prior specifications, we obtain the joint distribution
\[
\begin{array}{rl}
p(\boldsymbol{\theta}) \times \prod_{j=1}^{q} & I G\left(\tau_{j}^{2} \mid a_{\tau_{j}}, b_{\tau_{j}}\right) \times N\left(\boldsymbol{\beta} \mid \boldsymbol{\mu}_{\beta}, \mathbf{V}_{\beta}\right) \times N\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{B}_{\mathcal{U}} \mathbf{w}_{\mathcal{S}}, \mathbf{F}_{\mathcal{U}}\right) \\
& \times N\left(\mathbf{w}_{\mathcal{S}} \mid \mathbf{0}, \tilde{\mathbf{C}}_{\mathcal{S}}\right) \times \prod_{i=1}^{n} N\left(\mathbf{y}\left(\mathbf{t}_{i}\right) \mid \mathbf{X}\left(\mathbf{t}_{i}\right)^{\prime} \boldsymbol{\beta}+\mathbf{Z}\left(\mathbf{t}_{i}\right)^{\prime} \mathbf{w}\left(\mathbf{t}_{i}\right), \mathbf{D}\right)
\end{array}
\]
where \(p(\boldsymbol{\theta})\) is the prior on \(\boldsymbol{\theta}\) and \(I G\left(\tau_{j}^{2} \mid a_{\tau_{j}}, b_{\tau_{j}}\right)\) denotes the Inverse-Gamma density.

\subsection*{3.2 Estimation and prediction}

To describe a Gibbs sampler for estimating (9), we define \(\mathbf{y}=\left(\mathbf{y}\left(\mathbf{t}_{1}\right)^{\prime}, \mathbf{y}\left(\mathbf{t}_{2}\right)^{\prime}, \ldots, \mathbf{y}\left(\mathbf{t}_{n}\right)^{\prime}\right)^{\prime}\), and \(\mathbf{w}\) and \(\boldsymbol{\epsilon}\) similarly. Also, we introduce \(\mathbf{X}=\left[\mathbf{X}\left(\mathbf{t}_{1}\right): \mathbf{X}\left(\mathbf{t}_{2}\right): \ldots: \mathbf{X}\left(\mathbf{t}_{n}\right)\right]^{\prime}, \mathbf{Z}= \operatorname{diag}\left(\mathbf{Z}\left(\mathbf{t}_{1}\right)^{\prime}, \ldots, \mathbf{Z}\left(\mathbf{t}_{n}\right)^{\prime}\right)\), and \(\mathbf{D}_{n}=\operatorname{Cov}(\boldsymbol{\epsilon})=\operatorname{diag}(\mathbf{D}, \ldots, \mathbf{D})\). The full conditional distribution for \(\boldsymbol{\beta}\) is \(N\left(\mathbf{V}_{\beta}^{*} \boldsymbol{\mu}_{\beta}^{*}, \mathbf{V}_{\boldsymbol{\beta}}^{*}\right)\), where \(\mathbf{V}_{\beta}^{*}=\left(\mathbf{V}_{\beta}^{-1}+\mathbf{X}^{\prime} \mathbf{D}_{n}^{-1} \mathbf{X}\right)^{-1}, \boldsymbol{\mu}_{\beta}^{*}=\left(\mathbf{V}_{\beta}^{-1} \boldsymbol{\mu}_{\beta}+\mathbf{X}^{\prime} \mathbf{D}_{n}^{-1}(\mathbf{y}-\right. \mathbf{Z w})\) ). Inverse-Gamma priors for the \(\tau_{j}^{2}\) 's leads to conjugate full conditional distribution \(I G\left(a_{\tau_{j}}+\frac{n}{2}, b_{\tau_{j}}+\frac{1}{2}\left(\mathbf{y}_{* j}-\mathbf{X}_{* j} \boldsymbol{\beta}-\mathbf{Z}_{* j} \mathbf{w}\right)^{\prime}\left(\mathbf{y}_{* j}-\mathbf{X}_{* j} \boldsymbol{\beta}-\mathbf{Z}_{* j} \mathbf{w}\right)\right.\) where \(\mathbf{y}_{* j}\) refers to the \(n \times 1\) vector containing the \(j^{\text {th }}\) co-ordinates of the \(\mathbf{y}\left(\mathbf{t}_{i}\right)\) 's, \(\mathbf{X}_{* j}\) and \(\mathbf{Z}_{* j}\) are the corresponding fixed and spatial effect covariate matrices respectively. For updating \(\boldsymbol{\theta}\), we use a random walk Metropolis step with target density \(p(\boldsymbol{\theta}) \times N\left(\mathbf{w}_{\mathcal{S}} \mid \mathbf{0}, \tilde{\mathbf{C}}_{\mathcal{S}}\right) \times N\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{B}_{\mathcal{U}} \mathbf{w}_{\mathcal{S}}, \mathbf{F}_{\mathcal{U}}\right)\), where
\[
\begin{gathered}
N\left(\mathbf{w}_{\mathcal{S}} \mid \mathbf{0}, \tilde{\mathbf{C}}_{\mathcal{S}}\right)=\prod_{i=1}^{k} N\left(\mathbf{w}\left(\mathbf{s}_{i}\right) \mid \mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}, \mathbf{F}_{\mathbf{s}_{i}}\right) \text { and } \\
N\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{B}_{\mathcal{U}} \mathbf{w}_{\mathcal{S}}, \mathbf{F}_{\mathcal{U}}\right)=\prod_{i=r+1}^{n} N\left(\mathbf{w}\left(\mathbf{t}_{i}\right) \mid \mathbf{B}_{\mathbf{t}_{i}} \mathbf{w}_{N\left(\mathbf{t}_{i}\right)}, \mathbf{F}_{\mathbf{t}_{i}}\right)
\end{gathered}
\]

Each of the component densities under the product sign on the right hand side of (10) can be evaluated without any \(n\)-dimensional matrix operations rendering the NNGP suitable for efficient Metropolis (Hastings) block updates for \(\boldsymbol{\theta}\).

Since the components of \(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\) are independent, we can update \(\mathbf{w}\left(\mathbf{t}_{i}\right)\) from its full conditional \(N\left(\mathbf{V}_{\mathbf{t}_{i}} \boldsymbol{\mu}_{\mathbf{t}_{i}}, \mathbf{V}_{\mathbf{t}_{i}}\right)\) for \(i=r+1, r+2, \ldots, n\) where \(\mathbf{V}_{\mathbf{t}_{i}}=\left(\mathbf{Z}\left(\mathbf{t}_{i}\right) \mathbf{D}^{-1} \mathbf{Z}\left(\mathbf{t}_{i}\right)^{\prime}+\mathbf{F}_{\mathbf{t}_{i}}^{-1}\right)^{-1}\) and \(\boldsymbol{\mu}_{\mathbf{t}_{i}}=\mathbf{Z}\left(\mathbf{t}_{i}\right) \mathbf{D}^{-1}\left(\mathbf{y}\left(\mathbf{t}_{i}\right)-\mathbf{X}\left(\mathbf{t}_{i}\right)^{\prime} \boldsymbol{\beta}\right)+\mathbf{F}_{\mathbf{t}_{i}}^{-1} \mathbf{B}_{\mathbf{t}_{i}} \mathbf{w}_{N\left(\mathbf{t}_{i}\right)}\). Finally, we update the components of \(\mathbf{w}_{\mathcal{S}}\) individually. For any two locations \(\mathbf{s}\) and \(\mathbf{t}\) in \(\mathcal{D}\), if \(\mathbf{s} \in N(\mathbf{t})\) and is the \(l\)-th component of \(N(\mathbf{t})\), i.e., say \(\mathbf{s}=N(\mathbf{t})(l)\), then define \(\mathbf{B}_{\mathbf{t}, \mathbf{s}}\) as the \(l \times l\) submatrix formed by columns \((l-1) q+1,(l-1) q+2, \ldots, l q\) of \(\mathbf{B}_{\mathbf{t}}\). Let \(U\left(\mathbf{s}_{i}\right)=\left\{\mathbf{t} \in \mathcal{S} \cup \mathcal{T} \mid \mathbf{s}_{i} \in N(\mathbf{t})\right\}\) and for every \(\mathbf{t} \in U\left(\mathbf{s}_{i}\right)\) define, \(\mathbf{a}_{\mathbf{t}, \mathbf{s}_{i}}=\mathbf{w}(\mathbf{t})-\sum_{\mathbf{s} \in N(\mathbf{t}), \mathbf{s} \neq \mathbf{s}_{i}} \mathbf{B}_{\mathbf{t}, \mathbf{s}} \mathbf{w}(\mathbf{s})\). Then, for \(i=1,2, \ldots, k\), we have the full conditional \(\mathbf{w}_{\mathbf{s}_{i}} \mid \cdot \sim N\left(\mathbf{V}_{\mathbf{s}_{i}} \boldsymbol{\mu}_{\mathbf{s}_{i}}, \mathbf{V}_{\mathbf{s}_{i}}\right)\) where \(\mathbf{V}_{\mathbf{s}_{i}}=\left(\operatorname{In}\left(\mathbf{s}_{i} \in \mathcal{S}^{*}\right) \mathbf{Z}\left(\mathbf{s}_{i}\right) \mathbf{D}^{-1} \mathbf{Z}\left(\mathbf{s}_{i}\right)^{\prime}+\mathbf{F}_{\mathbf{s}_{i}}^{-1}+\right. \left.\sum_{\mathbf{t} \in U\left(\mathbf{s}_{i}\right)} \mathbf{B}_{\mathbf{t}, \mathbf{s}_{i}}^{\prime} \mathbf{F}_{\mathbf{t}}^{-1} \mathbf{B}_{\mathbf{t}, \mathbf{s}_{i}}\right)^{-1}, \boldsymbol{\mu}_{\mathbf{s}_{i}}=\operatorname{In}\left(\mathbf{s}_{i} \in \mathcal{S}^{*}\right) \mathbf{Z}\left(\mathbf{s}_{i}\right) \mathbf{D}^{-1}\left(\mathbf{y}\left(\mathbf{s}_{i}\right)-\mathbf{X}\left(\mathbf{s}_{i}\right)^{\prime} \boldsymbol{\beta}\right)+\mathbf{F}_{\mathbf{s}_{i}}^{-1} \mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}+ \sum_{\mathbf{t} \in U\left(\mathbf{s}_{i}\right)} \mathbf{B}_{\mathbf{t}, \mathbf{s}_{i}}^{\prime} \mathbf{F}_{\mathbf{t}}^{-1} \mathbf{a}_{\mathbf{t}, \mathbf{s}_{i}}\) and \(\operatorname{In}(\cdot)\) denotes the indicator function. Hence, the \(\mathbf{w}^{\prime}\) s can also be
updated without requiring storage or factorization of any \(n \times n\) matrices.
Turning to predictions, let \(\mathbf{t}\) be a new location where we intend to predict \(\mathbf{y}(\mathbf{t})\) given \(\mathbf{X}(\mathbf{t})\) and \(\mathbf{Z}(\mathbf{t})\). The Gibbs sampler for estimation also generates the posterior samples \(\mathbf{w}_{\mathcal{S}} \mid \mathbf{y}\). So, if \(\mathbf{t} \in \mathcal{S}\), then we simply get samples of \(\mathbf{y}(\mathbf{t}) \mid \mathbf{y}\) from \(N\left(\mathbf{X}(\mathbf{t})^{\prime} \boldsymbol{\beta}+\mathbf{Z}(\mathbf{t})^{\prime} \mathbf{w}(\mathbf{t}), \mathbf{D}\right)\). If \(\mathbf{t}\) is outside \(\mathcal{S}\), then we generate samples of \(\mathbf{w}(\mathbf{t})\) from its full conditional, \(N\left(\mathbf{V}_{\mathbf{t}} \boldsymbol{\mu}_{\mathbf{t}}, \mathbf{V}_{\mathbf{t}}\right)\), where \(\mathbf{V}_{\mathbf{t}}=\left(\mathbf{Z}(\mathbf{t}) \mathbf{D}^{-1} \mathbf{Z}(\mathbf{t})^{\prime}+\mathbf{F}_{\mathbf{t}}^{-1}\right)^{-1}\) and \(\boldsymbol{\mu}_{\mathbf{t}}=\mathbf{Z}(\mathbf{t}) \mathbf{D}^{-1}\left(\mathbf{y}(\mathbf{t})-\mathbf{X}(\mathbf{t})^{\prime} \boldsymbol{\beta}\right)+\mathbf{F}_{\mathbf{t}}^{-1} \mathbf{B}_{\mathbf{t}} \mathbf{w}_{N(\mathbf{t})}\), and subsequently generate posterior samples of \(\mathbf{y}(\mathbf{t}) \mid \mathbf{y}\) similar to the earlier case.

\subsection*{3.3 Computational complexity}

Implementing the NNGP model in Section 3.2 reveals that one entire pass of the Gibbs sampler can be completed without any large matrix operations. The only difference between (9) and a full geostatistical hierarchical model is that the spatial process is modeled as an NNGP prior as opposed to a standard GP. For comparisons, we offer rough estimates of the flop counts to generate \(\boldsymbol{\theta}\) and w per iteration of the sampler. We express the computational complexity only in terms of the sample size \(n\), size of the reference set \(k\) and the size of the neighbor sets \(m\) as other dimensions are assumed to be small. For all locations, \(\mathbf{t} \in \mathcal{S} \cup \mathcal{T}\), \(\mathbf{B}_{\mathbf{t}}\) and \(\mathbf{F}_{\mathbf{t}}\) can be calculated using \(O\left(m^{3}\right)\) flops. So, from (10) it is easy to see that \(p(\boldsymbol{\theta} \mid \cdot)\) can be calculated using \(O\left((n+k) m^{3}\right)\) flops. All subsequent calculations to generate a set of posterior samples for \(\mathbf{w}\) and \(\boldsymbol{\theta}\) require around \(O\left((n+k) m^{2}\right)\) flops.

So, the total flop counts is of the order \((n+k) m^{3}\) and is, therefore, linear in the total number of locations in \(\mathcal{S} \cup \mathcal{T}\). This ensures scalability of the NNGP to large datasets. Compare this with a full GP model with a dense correlation matrix, which requires \(O\left(n^{3}\right)\) flops for updating \(\mathbf{w}\) in each iteration. Simulation results in Section 5.1 and Appendix F indicate that NNGP models with usually very small values of \(m(\approx 10)\) provides inference almost indistinguishable to full geostatistical models. Therefore, for large \(n\), this linear flop count is drastically less. Also, linearity with respect to \(k\) ensures a feasible implementation even for
\(k \approx n\). This offers substantial improvement over low rank models where the computational cost is quadratic in the number of "knots," limiting the size of the set of knots. Also, both the full geostatistical and the predictive process models require storage of the \(n \times n\) distance matrix, which can potentially exhaust storage resources for large datasets. An NNGP model only requires the distance matrix between neighbors for every location, thereby storing \(n+k\) small matrices, each of order \(m \times m\). Hence, NNGP accrues substantial computational benefits over existing methods for very large spatial datasets and may be the only feasible option for fully model-based inference in certain cases, as seen in the forestry data example (Section 5.2).

\subsection*{3.4 Model comparison and choice of \(\mathcal{S}\) and \(m\)}

As elaborated in Section 2, given any parent Gaussian process and any fixed reference set of locations \(\mathcal{S}\), we can construct a valid NNGP. The resulting finite dimensional likelihoods of the NNGP depend upon the choice of the reference set \(\mathcal{S}\) and the size of each \(N\left(\mathbf{s}_{i}\right)\), i.e., \(m\). Choosing the reference set is similar to selecting the knots for a predictive process. Unlike the number of "knots" in low rank models, the number of points in \(\mathcal{S}\) do not thwart computational scalability. From Section 3.3, we observe that the flop count in an NNGP model only increases linearly with the size of \(\mathcal{S}\). Hence, the number of locations in \(\mathcal{S}\) can, in theory, be large and this provides a lot of flexibility in choosing \(\mathcal{S}\).

Points over a grid across the entire domain seem to be a plausible choice for \(\mathcal{S}\). For example, we can construct a large \(\mathcal{S}\) using a dense grid to improve performance without adversely affecting computational costs. Another, perhaps even simpler, option for large datasets is to simply fix \(\mathcal{S}=\mathcal{T}\), the set of observed locations. Since the NNGP is a legitimate process for any fixed \(\mathcal{S}\), this choice is legitimate and it reduces computational costs even further by avoiding additional sampling of \(\mathbf{w}_{\mathcal{U}}\) in the Gibbs sampler. Our empirical investigations (see Section 5.1) reveal that choosing \(\mathcal{S}=\mathcal{T}\) deliver inference almost indistinguishable from
choosing \(\mathcal{S}\) to be a grid over the domain for large datasets.
Stein et al. (2004) and Eidsvik et al. (2014) proposed using a sandwich variance estimator for evaluating the inferential abilities of neighbor-based pseudo-likelihoods. Shaby (2012) developed a post sampling sandwich variance adjustment for posterior credible intervals of the parameters for quasi-Bayesian approaches using pseudo-likelihoods. However, all these adjustments concede accrual of additional computational costs. Also, the asymptotic results used to obtain the sandwich variance estimators are based on assumptions which are hard to verify in spatial settings with irregularly placed data points. Moreover, we view the NNGP as an independent model for fitting the data and not as an approximation to the original GP. Hence, we refrain from such sandwich variance adjustments. Instead, we can simply use any standard model comparison metrics such as DIC (Spiegelhalter et al. 2002), GPD (Gelfand and Ghosh 1998) or RMSPE(RMSECV) Yeniay and Goktas 2002) to compare the performance of NNGP and any other candidate model. The same model comparison metrics are also used for selecting \(m\). However, as we illustrate later in Section 5.1, usually a small value of \(m\) between 10 to 15 produces performance at par with the full geostatistical model. While larger \(m\) may be beneficial for massive datasets, perhaps under a different design scheme, it is still going to be much smaller than the number of knots required in low rank models (see Section 5.1).

\section*{4 Alternate NNGP models and algorithms}

\subsection*{4.1 Block update of \(\mathrm{w}_{\mathcal{S}}\) using sparse Cholesky}

The Gibbs' sampling algorithm detailed in Section 3.2 is extremely efficient for large datasets with linear flop counts per iteration. However, it can sometimes experience slow convergence issues due to sequential updating of the elements in \(\mathbf{w}_{\mathcal{S}}\). An alternative to sequential updating is to perform block updates of \(\mathbf{w}_{\mathcal{S}}\). We choose \(\mathcal{S}=\mathcal{T}\) so that \(\mathbf{s}_{i}=\mathbf{t}_{i}\) for all \(i=1,2, \ldots, n\)
and we denote \(\mathbf{w}_{\mathcal{S}}=\mathbf{w}_{\mathcal{T}}\) by \(\mathbf{w}\). Then,
\[
\mathbf{w} \mid \cdot \sim N\left(\mathbf{V}_{\mathcal{S}} \mathbf{Z}^{\prime} \mathbf{D}_{n}^{-1}(\mathbf{y}-\mathbf{X} \boldsymbol{\beta}), \mathbf{V}_{\mathcal{S}}\right), \quad \text { where } \quad \mathbf{V}_{\mathcal{S}}=\left(\mathbf{Z}^{\prime} \mathbf{D}_{n}^{-1} \mathbf{Z}+\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\right)^{-1}
\]

Recall that \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\) is sparse. Since \(\mathbf{Z}\) and \(\mathbf{D}_{n}\) are block diagonal, \(\mathbf{V}_{\mathcal{S}}^{-1}\) retains the sparsity of \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\). So, a sparse Cholesky factorization of \(\mathbf{V}_{\mathcal{S}}^{-1}\) will efficiently produce the Cholesky factors of \(\mathbf{V}_{\mathcal{S}}\). This will facilitate block updating of \(\mathbf{w}\) in the Gibbs sampler.

\subsection*{4.2 NNGP models for the response}

Another possible approach involves NNGP models for the response \(\mathbf{y}(\mathbf{s})\). If \(\mathbf{w}(\mathbf{s})\) is a Gaussian Process, then so is \(\mathbf{y}(\mathbf{s})=\mathbf{Z}(\mathbf{s})^{\prime} \mathbf{w}(\mathbf{s})+\boldsymbol{\epsilon}\) (without loss of generality we assume \(\beta=\mathbf{0}\) ). One can directly use the NNGP specification for \(\mathbf{y}(\mathbf{s})\) instead of \(\mathbf{w}(\mathbf{s})\). That is, we derive \(\mathbf{y}(\mathbf{s}) \sim N N G P(\mathbf{0}, \tilde{\boldsymbol{\Sigma}}(\cdot, \cdot))\) from the parent Gaussian process \(\operatorname{GP}(\mathbf{0}, \boldsymbol{\Sigma}(\cdot, \cdot \mid \boldsymbol{\theta}))\). The Gibbs sampler analogous to Section 3 now enjoys the additional advantage of avoiding full conditionals for \(\mathbf{w}\). This results in a Bayesian analogue for Vecchia (1988) and Stein et al. (2004) but precludes inference on the spatial residual surface \(\mathbf{w}(\mathbf{s})\). Modeling w(s) provides additional insight into residual spatial contours and is often important in identifying lurking covariates or eliciting unexplained spatial patterns. Vecchia (1992) used the nearest neighbor approximation on a spatial model for observations ( \(\mathbf{y}\) ) with independent measurement error (nuggets) in addition to the usual spatial component (w). However, it may not be possible to recover \(\mathbf{w}\) using this approach. For example, a univariate stationary process \(\mathbf{y}(\mathbf{s})\) with a nugget effect can be decomposed as \(\mathbf{y}(\mathbf{s})=\mathbf{w}(\mathbf{s})+\boldsymbol{\epsilon}(\mathbf{s})\) (letting \(\boldsymbol{\beta}=\mathbf{0}\) ) for some \(\mathbf{w}(\mathbf{s}) \sim G P(\mathbf{0}, \mathbf{C}(\cdot, \cdot \mid \boldsymbol{\theta}))\) and white noise process \(\boldsymbol{\epsilon}(\mathbf{s})\). If \(\mathbf{y}=\mathbf{w}+\boldsymbol{\epsilon}\), where \(\mathbf{w} \sim N(\mathbf{0}, \mathbf{C}), \boldsymbol{\epsilon} \sim N\left(\mathbf{0}, \tau^{2} \mathbf{I}_{n}\right)\), then \(\operatorname{Cov}(\mathbf{y})= \mathbf{C}+\tau^{2} \mathbf{I}=\boldsymbol{\Sigma}\), all eigenvalues of \(\boldsymbol{\Sigma}\) are greater than \(\tau^{2}\) and \(\operatorname{Cov}(\mathbf{w} \mid \mathbf{y})=\tau^{2} \mathbf{I}_{n}-\tau^{4} \boldsymbol{\Sigma}^{-1}\). For \(\mathbf{y}(\mathbf{s}) \sim N N G P(\mathbf{0}, \tilde{\boldsymbol{\Sigma}}(\cdot, \cdot))\), however, the eigenvalues of \(\tilde{\boldsymbol{\Sigma}}\) may be less than \(\tau^{2}\), so \(\tau^{2} \mathbf{I}_{n}-\tau^{4} \tilde{\boldsymbol{\Sigma}}^{-1}\) need not be positive definite for every \(\tau^{2}>0\) and \(p(\mathbf{w} \mid \mathbf{y})\) is no longer well-defined.

A different model is obtained by using an NNGP prior for \(\mathbf{w}\), as in (9), and then in-
tegrating out \(\mathbf{w}\). The resulting likelihood is \(N\left(\mathbf{y} \mid \mathbf{X} \boldsymbol{\beta}, \boldsymbol{\Sigma}_{y}\right)\), where \(\boldsymbol{\Sigma}_{y}=\mathbf{Z} \tilde{\mathbf{C}}_{\mathcal{S}} \mathbf{Z}^{\prime}+\mathbf{D}_{n}\) and the Bayesian specification is completed using priors on \(\boldsymbol{\beta}, \tau_{j}^{2}\) 's and \(\boldsymbol{\theta}\) as in (9). This model drastically reduces the number of variables in the Gibbs sampler, while preserving the nugget effect in the parent model. We can generate the full conditionals for the parameters in the marginalized model as follows: \(\boldsymbol{\beta} \mid \mathbf{y}, \phi \sim N\left(\left(\mathbf{V}_{\beta}^{-1}+\mathbf{X}^{\prime} \boldsymbol{\Sigma}_{y}^{-1} \mathbf{X}\right)^{-1}\left(\mathbf{V}_{\beta}^{-1} \boldsymbol{\mu}_{\beta}+\mathbf{X}^{\prime} \boldsymbol{\Sigma}_{y}^{-1} \mathbf{y}\right),\left(\mathbf{V}_{\beta}^{-1}+\right.\right. \left.\mathbf{X}^{\prime} \boldsymbol{\Sigma}_{y}^{-1} \mathbf{X}\right)^{-1}\) ). It is difficult to factor out \(\tau_{j}^{2}\), s from \(\boldsymbol{\Sigma}_{y}^{-1}\), so conjugacy is lost with respect to any standard prior. Metropolis block updates for \(\boldsymbol{\theta}\) are feasible for any tractable prior \(p(\boldsymbol{\theta})\). This involves computing \(\mathbf{X}^{\prime} \boldsymbol{\Sigma}_{y}^{-1} \mathbf{X}, \mathbf{X}^{\prime} \boldsymbol{\Sigma}_{y}^{-1} \mathbf{y}\) and \((\mathbf{y}-\mathbf{X} \boldsymbol{\beta})^{\prime} \boldsymbol{\Sigma}_{y}^{-1}(\mathbf{y}-\mathbf{X} \boldsymbol{\beta})\). Since \(\boldsymbol{\Sigma}_{y}^{-1}=\mathbf{D}_{n}^{-1}-\mathbf{D}_{n}^{-1} \mathbf{Z}\left(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}+\mathbf{Z}^{\prime} \mathbf{D}_{n}^{-1} \mathbf{Z}\right)^{-1} \mathbf{Z}^{\prime} \mathbf{D}_{n}^{-1}=\mathbf{D}_{n}^{-1}-\mathbf{D}_{n}^{-1} \mathbf{Z} \mathbf{V}_{\mathcal{S}} \mathbf{Z}^{\prime} \mathbf{D}_{n}^{-1}\), where \(\mathbf{V}_{\mathcal{S}}\) is given by (11), a sparse Cholesky factorization of \(\mathbf{V}_{\mathcal{S}}^{-1}\) will be beneficial. We draw posterior samples for \(\mathbf{w}\) from \(p(\mathbf{w} \mid \mathbf{y})=\int p\left(\mathbf{w} \mid \boldsymbol{\theta}, \boldsymbol{\beta},\left\{\tau_{j}^{2}\right\}, \mathbf{y}\right) p\left(\boldsymbol{\theta}, \boldsymbol{\beta},\left\{\tau_{j}^{2}\right\} \mid \mathbf{y}\right)\) using composition sampling-we draw \(\mathbf{w}^{(g)}\) from \(p\left(\mathbf{w} \mid \boldsymbol{\theta}^{(g)}, \boldsymbol{\beta}^{(g)},\left\{\tau_{j}^{2(g)}\right\}, \mathbf{y}\right)\) one-for-one for each sampled parameter.

Using block updates for \(\mathbf{w}_{\mathcal{S}}\) in (9) and fitting the marginalized version of (9) both require an efficient sparse Cholesky solver for \(\mathbf{V}_{\mathcal{S}}^{-1}\). Note that computational expenses for most sparse Cholesky algorithms depend on the precise nature of the sparse structure (mostly on the bandwidth) of \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\) (see, e.g. Davis 2006). The number of flops required for Gibbs sampling and prediction in this marginalized model depends upon the sparse structure of \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\) and may, sometimes, heavily exceed the linear usage achieved by the unmarginalized model with individual updates for \(\mathbf{w}_{i}\). Therefore, a prudent choice of the precise fitting algorithms should be based on the sparsity structure of \(\mathbf{\tilde { \mathbf { C } }}_{\mathcal{S}}^{-1}\) for the given dataset.

\subsection*{4.3 Spatiotemporal and GLM versions}

In spatiotemporal settings where we seek spatial interpolation at discrete time-points (e.g., weekly, monthly or yearly data), we write the response (possibly vector-valued) as \(\mathbf{y}_{t}(\mathbf{s})\) and the random effects as \(\mathbf{w}_{t}(\mathbf{s})\). One could, for example, envision that the data arise as a time series of spatial processes, i.e., there is a time series at each location. An alternative scenario
is cross-sectional data being collected at a set of locations associated with each time point and these locations can differ from time point to time point. Desired inference includes spatial interpolation for each time point. Spatial dynamic models incorporating the NNGP are easily formulated as below:
\[
\begin{gathered}
\mathbf{y}_{t}(\mathbf{s})=\mathbf{X}_{t}(\mathbf{s})^{\prime} \boldsymbol{\beta}_{t}+\mathbf{u}_{t}(s)+\boldsymbol{\epsilon}_{t}(\mathbf{s}), \boldsymbol{\epsilon}_{t}(\mathbf{s}) \stackrel{i i d}{\sim} N(0, D) \\
\boldsymbol{\beta}_{t}=\boldsymbol{\beta}_{t-1}+\boldsymbol{\eta}_{t}, \boldsymbol{\eta}_{t} \stackrel{i i d}{\sim} N\left(0, \boldsymbol{\Sigma}_{\eta}\right), \boldsymbol{\beta}_{0} \sim N\left(\mathbf{m}_{0}, \boldsymbol{\Sigma}_{0}\right) \\
\mathbf{u}_{t}(\mathbf{s})=\mathbf{u}_{t-1}(\mathbf{s})+\mathbf{w}_{t}(\mathbf{s}), \mathbf{w}_{t}(\mathbf{s}) \stackrel{i n d}{\sim} N N G P\left(\mathbf{0}, \tilde{\mathbf{C}}\left(\cdot, \cdot \mid \boldsymbol{\theta}_{t}\right)\right) .
\end{gathered}
\]

Thus, one retains exactly the same structure of process-based spatial dynamic models, e.g., as in Gelfand et al. (2005), and simply replaces the independent Gaussian process priors for \(\mathbf{w}_{t}(\mathbf{s})\) with independent NNGP's to achieve computational tractability.

The above is illustrative of how attractive and extremely convenient the NNGP is for model building. One simply writes down the parent model and subsequently replaces the full GP with an NNGP. Being a well-defined process, the NNGP ensures a valid spatial dynamic model. Similarly NNGP versions of dynamic spatiotemporal Kalman-filtering Wikle and Cressie 1999, as, e.g., in) can be constructed.

Handling non-Gaussian (e.g., binary or count) data is also straightforward using spatial generalized linear models (GLM's) (Diggle et al. 1998; Lin et al. 2000; Kammann and Wand 2003; Banerjee et al. 2014). Here, the NNGP provides structured dependence for random effects at the second stage. First, we replace \(\mathrm{E}[\mathbf{y}(\mathbf{t})]\) in (8) with \(g(E(\mathbf{y}(\mathbf{t})))\) where \(g(\cdot)\) is a suitable link function such that \(\boldsymbol{\eta}(\mathbf{t})=g(E(\mathbf{y}(\mathbf{t})))=\mathbf{X}(\mathbf{t})^{\prime} \boldsymbol{\beta}+\mathbf{Z}(\mathbf{t})^{\prime} \mathbf{w}(\mathbf{t})\). In the second stage, we model the \(\mathbf{w}(\mathbf{t})\) as an NNGP. The benefits of the algorithms in Sections 3.2 and 3.3 still hold, but some of the alternative algorithms in Section 4 may not apply. For example, we do obtain tractable marginalized likelihoods by integrating out the spatial effects.

\section*{5 Illustrations}

We conduct simulation experiments and analyze a large forestry dataset. Additional simulation experiments are detailed in Appendices C through I. Posterior inference for subsequent analysis were based upon three chains of 25000 iterations (with a burn-in of 5000 iterations). All the samplers were programmed in C++ and leveraged Intels Math Kernel Library's (MKL) threaded BLAS and LAPACK routines for matrix computations on a Linux workstation with 384 GB of RAM and two Intel Nehalem quad-Xeon processors.

\subsection*{5.1 Simulation experiment}

We generated observations using 2500 locations within a unit square domain from the model (8) with \(q=l=1\) (univariate outcome), \(p=2, \mathbf{Z}(\mathbf{t})^{\prime}=1\) (scalar), the spatial covariance matrix \(\mathbf{C}(\boldsymbol{\theta})=\sigma^{2} \mathbf{R}(\boldsymbol{\phi})\), where \(\mathbf{R}(\boldsymbol{\phi})\) is a \(n \times n\) correlation matrix, and \(\mathbf{D}=\tau^{2}\) (scalar). The model included an intercept and a covariate \(\mathbf{x}_{1}\) drawn from \(N(0,1)\). The \((i, j)\) th element of \(\mathbf{R}(\phi)\) was calculated using the Matérn function
\[
\rho\left(\mathbf{t}_{i}, \mathbf{t}_{j} ; \boldsymbol{\phi}\right)=\frac{1}{2^{\nu-1} \Gamma(\nu)}\left(\left\|\mathbf{t}_{i}-\mathbf{t}_{j}\right\| \phi\right)^{\nu} \mathcal{K}_{\nu}\left(\left\|\mathbf{t}_{i}-\mathbf{t}_{j}\right\| \phi\right) ; \phi>0, \nu>0,
\]
where \(\left\|\mathbf{t}_{i}-\mathbf{t}_{j}\right\|\) is the Euclidean distance between locations \(\mathbf{t}_{i}\) and \(\mathbf{t}_{j}, \boldsymbol{\phi}=(\phi, \nu)\) with \(\phi\) controlling the decay in spatial correlation and \(\nu\) controlling the process smoothness, \(\Gamma\) is the usual Gamma function while \(\mathcal{K}_{\nu}\) is a modified Bessel function of the second kind with order \(\nu\) (Stein 1999) Evaluating the Gamma function for each matrix element within each iteration requires substantial computing time and can obscure differences in sampler run times; hence, we fixed \(\nu\) at 0.5 which reduces (13) to the exponential correlation function. The first column in Table 1 gives the true values used to generate the responses. Figure 2(a) illustrates the \(w(\mathbf{t})\) surface interpolated over the domain.

We then estimated the following models from the full data: \(i)\) the full Gaussian Process

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1: Univariate synthetic data analysis parameter estimates and computing time in minutes for NNGP and full GP models. Parameter posterior summary 50 (2.5, 97.5) percentiles.}
\begin{tabular}{|l|l|l|l|l|l|}
\hline \multirow{2}{*}{} & \multicolumn{3}{|c|}{NNGP ( \(\mathcal{S} \neq \mathcal{T}\) )} & \multicolumn{2}{|c|}{NNGP \((\mathcal{S}=\mathcal{T})\)} \\
\hline & True & \(m=10, k=2000\) & \(m=20, k=2000\) & \(m=10\) & \(m=20\) \\
\hline \(\beta_{0}\) & 1 & 0.99 (0.71, 1.48) & 1.02 (0.73, 1.49) & 1.00 (0.62, 1.31) & 1.03 (0.65, 1.34) \\
\hline \(\beta_{1}\) & 5 & 5.00 (4.98, 5.03) & 5.01 (4.98, 5.03) & 5.01 (4.99, 5.03) & 5.01 (4.99, 5.03) \\
\hline \(\sigma^{2}\) & 1 & 1.09 (0.89, 1.49) & 1.04 (0.85, 1.40) & 0.96 (0.78, 1.23) & 0.94 (0.77, 1.20) \\
\hline \(\tau^{2}\) & 0.1 & 0.07 (0.04, 0.10) & 0.07 (0.04, 0.10) & 0.10 (0.08, 0.13) & 0.10 (0.08, 0.13) \\
\hline \(\phi\) & 12 & 11.81 (8.18, 15.02) & 12.21 (8.83, 15.62) & 12.93 (9.70, 16.77) & 13.36 (9.99, 17.15) \\
\hline \(\mathrm{p}_{D}\) & - & 1491.08 & 1478.61 & 1243.32 & 1249.57 \\
\hline DIC & - & 1856.85 & 1901.57 & 2390.65 & 2377.51 \\
\hline G & - & 33.67 & 35.68 & 77.84 & 76.40 \\
\hline P & - & 253.03 & 259.13 & 340.40 & 337.88 \\
\hline D & - & 286.70 & 294.82 & 418.24 & 414.28 \\
\hline RMSPE & - & 1.22 & 1.22 & 1.2 & 1.2 \\
\hline 95\% CI cover \% & - & 97.2 & 97.2 & 97.6 & 97.6 \\
\hline 95\% CI width & - & 2.19 & 2.18 & 2.13 & 2.12 \\
\hline Time & - & 14.2 & 47.08 & 9.98 & 33.5 \\
\hline
\end{tabular}
\end{table}

\begin{tabular}{|l|l|l|l|}
\hline \multirow{2}{*}{} & \multirow[b]{2}{*}{True} & Predictive Process & Full \\
\hline & & 64 knots & Gaussian Process \\
\hline \(\beta_{0}\) & 1 & 1.30 (0.54, 2.03) & 1.03 (0.69, 1.34) \\
\hline \(\beta_{1}\) & 5 & 5.03 (4.99, 5.06) & 5.01 (4.99, 5.03) \\
\hline \(\sigma^{2}\) & 1 & 1.29 (0.96, 2.00) & 0.94 (0.76, 1.23) \\
\hline \(\tau^{2}\) & 0.1 & 0.08 (0.04, 0.13) & 0.10 (0.08, 0.12) \\
\hline \(\phi\) & 12 & 5.61 (3.48, 8.09) & 13.52 (9.92, 17.50) \\
\hline \(\mathrm{P}_{D}\) & - & 1258.27 & 1260.68 \\
\hline DIC & - & 13677.97 & 2364.80 \\
\hline G & - & 1075.63 & 74.80 \\
\hline P & - & 200.39 & 333.27 \\
\hline D & - & 1276.03 & 408.08 \\
\hline RMSPE & - & 1.68 & 1.2 \\
\hline 95\% CI cover \% & - & 95.6 & 97.6 \\
\hline 95\% CI width & - & 2.97 & 2.12 \\
\hline Time & - & 43.36 & 560.31 \\
\hline
\end{tabular}
(Full \(G P\) ); \(i i\) ) the NNGP with \(m=\{1,2, \ldots, 25\}\) for \(\mathcal{S} \neq \mathcal{T}\) and \(\mathcal{S}=\mathcal{T}\), and; \(i i i\) ) a Gaussian Predictive Process (GPP) model (Banerjee et al. 2008) with 64 knots placed on a grid over the domain. For the NNGP with \(\mathcal{S} \neq \mathcal{T}\) we considered 2000 randomly placed reference locations within the domain. The 64 knot GPP was chosen because its computing time was comparable to that of NNGP models. We used an efficient marginalized sampling algorithm for the Full GP and GPP models as implemented in the spBayes package in R (Finley et al. 2013). All the models were trained using 2000 of the 2500 observed locations, while the remaining 500 observations were withheld to assess predictive performance.

For all models, the intercept and slope regression parameters, \(\beta_{0}\) and \(\beta_{1}\), were given flat prior distributions. The variance components \(\sigma^{2}\) and \(\tau^{2}\) were assigned inverse-Gamma \(I G(2,1)\) and \(I G(2,0.1)\) priors, respectively, and the spatial decay \(\phi\) received a uniform prior \(U(3,30)\), which corresponds to a spatial range between approximately 0.1 and 1 units.

Parameter estimates and performance metrics for the NNGP (with \(m=10\) and \(m=20\) ), GPP, and the Full GP models are provided in Table 1. All model specifications produce

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-20.jpg?height=596&width=1467&top_left_y=399&top_left_x=315}
\captionsetup{labelformat=empty}
\caption{Figure 1: Choice of \(m\) in NNGP models: Out-of-sample Root Mean Squared Prediction Error (RMSPE) and mean width between the upper and lower \(95 \%\) posterior predictive credible intervals for a range of \(m\) for the univariate synthetic data analysis}
\end{figure}
similar posterior median and \(95 \%\) credible intervals estimates, with the exception of \(\phi\) in the 64 knot GPP model. Larger values of DIC and D suggest that the GPP model does not fit the data as well as the NNGP and Full GP models. The NNGP \(\mathcal{S}=\mathcal{T}\) models provide DIC, GPD scores that are comparable to those of the Full GP model. These fit metrics suggest the NNGP \(\mathcal{S} \neq \mathcal{T}\) models provide better fit to the data than that achieved by the full GP model which is probably due to overfitting caused by a very large reference set \(\mathcal{S}\). The last row in Table 1 shows computing times in minutes for one chain of 25000 iterations reflecting on the enormous computational gains of NNGP models over full GP model.

Turning to out-of-sample predictions, the Full model's RMSPE and mean width between the upper and lower \(95 \%\) posterior predictive credible interval is 1.2 and 2.12 , respectively. As seen in Figure 1, comparable RMSPE and mean interval width for the NNGP \(\mathcal{S}=\mathcal{T}\) model is achieved within \(m \approx 10\). There are negligible difference between the predictive performances of the NNGP \(\mathcal{S} \neq \mathcal{T}\) and \(\mathcal{S}=\mathcal{T}\) models. Both the NNGP and Full GP model have better predictive performance than the Predictive Process models when the number of knots is small, e.g., 64 . All models showed appropriate \(95 \%\) credible interval coverage rates.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-21.jpg?height=1092&width=1503&top_left_y=365&top_left_x=296}
\captionsetup{labelformat=empty}
\caption{Figure 2: Univariate synthetic data analysis: Interpolated surfaces of the true spatial random effects and posterior median estimates for different models}
\end{figure}

Figures 2(b-f) illustrate the posterior median estimates of the spatial random effects from the Full GP, NNGP ( \(\mathcal{S}=\mathcal{T}\) ) with \(m=10\) and \(m=20\), NNGP ( \(\mathcal{S} \neq \mathcal{T}\) ) with \(m=10\) and GPP models. These surfaces can be compared to the true surface depicted in Figure 2(a). This comparison shows: \(i)\) the NNGP models closely approximates the true surface and that estimated by the Full GP model, and; \(i i)\) the reduced rank predictive process model based on 64 knots greatly smooths over small-scale patterns. This last observation highlights one of the major criticisms of reduced rank models Stein (2014) and illustrates why these models often provide compromised predictive performance when the true surface has fine spatial resolution details. Overall, we see the clear computational advantage of the NNGP over the

Full GP model, and both inferential and computational advantage over the GPP model.

\subsection*{5.2 Forest biomass data analysis}

Information about the spatial distribution of forest biomass is needed to support global, regional, and local scale decisions, including assessment of current carbon stock and flux, bio-feedstock for emerging bio-economies, and impact of deforestation. In the United States, the Forest Inventory and Analysis (FIA) program of the USDA Forest Service collects the data needed to support these assessments. The program has established field plot centers in permanent locations using a sampling design that produces an equal probability sample (Bechtold and Patterson 2005). Field crews recorded stem measurements for all trees with diameter at breast height (DBH); 1.37 m above the forest floor) of 12.7 cm or greater. Given these data, established allometric equations were used to estimate each plot's forest biomass. For the subsequent analysis, plot biomass was scaled to metric tons per ha then square root transformed. The transformation ensures that back transformation of subsequent predicted values have support greater than zero and helps to meet basic regression models assumptions.

Figure 3(a) illustrates the georeferenced forest inventory data consisting of 114, 371 forested FIA plots measured between 1999 and 2006 across the conterminous United States. The two blocks of missing observations in the Western and Southwestern United States correspond to Wyoming and New Mexico, which have not yet released FIA data. Figure 3(b) shows a deterministic interpolation of forest biomass observed on the FIA plots. Dark blue indicates high forest biomass, which is primarily seen in the Pacific Northwest, Western Coastal ranges, Eastern Appalachian Mountains, and in portions of New England. In contrast, dark red indicates regions where climate or land use limit vegetation growth.

A July 2006 Normalized Difference Vegetation Index (NDVI) image from the MODerateresolution Imaging Spectroradiometer (MODIS); http://glcf.umd.edu/data/ndvi) sensor was used as a single predictor. NDVI is calculated from the visible and near-infrared light
reflected by vegetation, and can be viewed as a measure of greenness. In this image, Figure 3(c), dark green corresponds to dense vegetation whereas brown identifies regions of sparse or no vegetation, e.g., in the Southwest. NDVI is commonly used as a covariate in forest biomass regression models, see, for e.g., Zhang and Kondraguanta (2006). Results from these and similar studies show a positive linear relationship between forest biomass and NDVI. The strength of this relationship, however, varies by forest tree species composition, age, canopy structure, and level of reflectance. We expect a space-varying relationship between biomass and NDVI, given tree species composition and disturbance regimes generally exhibit strong spatial dependence across forested landscapes.

The \(\sim 38\) gigabytes of memory in our workstation was insufficient for storage of distance matrices required to fit a Full GP or GPP model. Subsequently, we explore the relationship between forest biomass and NDVI using a non-spatial model, a NNGP space-varying intercept (SVI) model (i.e., \(q=l=1\) and \(\mathbf{Z}(\mathbf{t})=1\) ) in (8), and a NNGP spatially-varying coefficients (SVC) regression model with \(l=1, q=p=2\) and \(\mathbf{Z}(\mathbf{t})=\mathbf{X}(\mathbf{t})\) in (8). The reference sets for the NNGP models were again the observed locations and \(m\) was chosen to be 5 or 10 . The parent process \(\mathbf{w}(\mathbf{t})\) is a bivariate Gaussian process with a isotropic cross-covariance specification \(\mathbf{C}\left(\mathbf{t}_{i}, \mathbf{t}_{j} \mid \boldsymbol{\theta}\right)=\mathbf{A} \boldsymbol{\Gamma}(\boldsymbol{\phi}) \mathbf{A}^{\prime}\), where \(\mathbf{A}\) is \(2 \times 2\) lower-triangular with positive diagonal elements, \(\boldsymbol{\Gamma}\) is \(2 \times 2\) diagonal with \(\rho\left(\mathbf{t}_{i}, \mathbf{t}_{j} ; \boldsymbol{\phi}_{b}\right)\) (defined in (13)) as the \(b^{t h}\) diagonal entry, \(b=1,2\) and \(\boldsymbol{\phi}_{b}=\left(\phi_{b}, \nu_{b}\right)^{\prime}\) (see, e.g., Gelfand and Banerjee 2010).

For all models, the intercept and slope regression parameters were given flat prior distributions. The variance components \(\tau^{2}\) and \(\sigma^{2}\) were assigned inverse-Gamma \(\operatorname{IG}(2,1)\) priors, the SVC model cross-covariance matrix \(\mathbf{A A}^{\prime}\) followed an inverse-Wishart \(I W(3,0.1)\), and the Matérn spatial decay and smoothness parameters received uniform prior supports \(U(0.01,3)\) and \(U(0.1,2)\), respectively. These prior distributions on \(\phi\) and \(\nu\) correspond to support between approximately 0.5 and 537 km . Candidate models are assessed using the metrics described in Section 3.4, inference drawn from mapped estimates of the regression

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2: Forest biomass data analysis parameter estimates and computing time in hours for candidate models. Parameter posterior summary 50 (2.5, 97.5) percentiles.}
\begin{tabular}{|l|l|l|l|}
\hline & Non-spatial & NNGP Space-varying intercept & NNGP Space-varying coefficients \\
\hline \(\beta_{0}\) & 1.043 (1.02, 1.065) & 1.44 (1.39, 1.48) & 1.23 (1.20, 1.26) \\
\hline \(\beta_{\text {NDVI }}\) & 0.0093 (0.009, 0.0095) & 0.0061 (0.0059, 0.0062) & 0.0072 (0.0071, 0.0074) \\
\hline \(\sigma^{2}\) & - & 0.16 (0.15, 0.17) & - \\
\hline \(\mathbf{A A}_{1,1}^{\prime}\) & - & - & 0.24 (0.23, 0.24) \\
\hline \(\mathbf{A A}_{2,1}^{\prime}\) & - & - & -0.00088 (-0.00093, -0.00083) \\
\hline \(\mathbf{A A}_{2,2}^{\prime}\) & - & - & 0.0000052 (0.0000047, 0.0000056) \\
\hline \(\tau^{2}\) & 0.52 (0.51, 0.52) & 0.39 (0.39, 0.40) & 0.39 (0.38, 0.40) \\
\hline \(\phi_{1}\) & - & 0.016 (0.015, 0.016) & 0.022 (0.021, 0.023) \\
\hline \(\phi_{2}\) & - & - & 0.030 (0.029, 0.031) \\
\hline \(\nu_{1}\) & - & 0.66 (0.64, 0.67) & 0.92 (0.90, 0.93) \\
\hline \(\nu_{2}\) & - & - & 0.92 (0.89, 0.93) \\
\hline \(\mathrm{p}_{D}\) & 2.94 & 6526.95 & 4976.13 \\
\hline DIC & 250137 & 224484.2 & 222845.1 \\
\hline G & 59765.30 & 42551.08 & 43117.37 \\
\hline P & 59667.15 & 47603.47 & 46946.49 \\
\hline D & 119432.45 & 90154.55 & 90063.86 \\
\hline Time & - & 14.53 & 41.35 \\
\hline
\end{tabular}
\end{table}
coefficients, and out-of-sample prediction.
Parameter estimates and performance metrics for NNGP with \(m=5\) are shown in Table 2. The corresponding numbers for \(m=10\) were similar. Relative to the spatial models, the non-spatial model has higher values of DIC and D which suggests NDVI alone does not adequately capture the spatial structure of forest biomass. This observation is corroborated using a variogram fit to the non-spatial model's residuals, Figure 3(d). The variogram shows a nugget of \(\sim 0.42\), partial sill of \(\sim 0.05\), and range of \(\sim 150 \mathrm{~km}\). This residual spatial dependence is apparent when we map the SVI model spatial random effects as shown in Figure 3(e). This map, and the estimate of a non-negligible spatial variance \(\sigma^{2}\) in Table 2 , suggests the addition of a spatial random effect was warranted and helps satisfy the model assumption of uncorrelated residuals.

The values of the SVC model's goodness of fit metrics suggest that allowing the NDVI

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-25.jpg?height=1754&width=1578&top_left_y=388&top_left_x=264}
\captionsetup{labelformat=empty}
\caption{Figure 3: Forest biomass data analysis: (a) locations of observed biomass, (b) interpolated biomass response variable, (c) NDVI regression covariate, (d) variogram of non-spatial model residuals, and (e) surface of the SVI model random spatial effects posterior medians. Following our FIA data sharing agreement, plot locations depicted in (a) have been "fuzzed" to hide the true coordinates.}
\end{figure}
regression coefficient to vary spatially improves model fit over that achieved by the SVI model. Figures 4(a) and 4(b) show maps of posterior estimates for the spatially varying intercept and NDVI, respectively. The clear regional patterns seen in Figure 4(b) suggest the relationship between NDVI and biomass does vary spatially - with stronger positive regression coefficients in the Pacific Northwest and northern California areas. Forest in the Pacific Northwest and northern California is dominated by conifers and support the greatest range in biomass per unit area within the entire conterminous United States. The other strong regional pattern seen in Figure 4(b) is across western New England, where near zero regression coefficients suggest that NDVI is not as effective at discerning differences in forest biomass. This result is not surprising. For deciduous forests, NDVI can explain variability in low to moderate vegetation density. However, in high biomass deciduous forests, like those found across western New England, NDVI saturates and is no longer sensitive to changes in vegetation structure (Wang et al. 2005). Hence, we see a higher intercept in this region but lower slope coefficient on NDVI. Figures 4(c) and 4(d) map each location's posterior predictive median and the range between the upper and lower \(95 \%\) credible interval, respectively, from the SVC model. Figure 4(c) shows strong correspondence with the deterministic interpolation of biomass in Figure 3(b). The prediction uncertainty in Figure 4(d) provides a realistic depiction of the model's ability to quantify forest biomass across the United States.

We also used prediction mean squared error (PMSE) to assess predictive performance. We fit the candidate models using 100,000 observations and withheld 14,371 for validation. PMSE for the non-spatial, SVI, and SVC models was \(0.52,0.41\), and 0.42 respectively. Lower PMSE for the spatial models, versus the non-spatial model, corroborates the results from the model fit metrics and further supports the need for spatial random effects in the analysis.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-27.jpg?height=1137&width=1606&top_left_y=380&top_left_x=259}
\captionsetup{labelformat=empty}
\caption{Figure 4: Forest biomass data analysis using SVC model: (a) Posterior medians of the intercept, (b) NDVI regression coefficients, (c) median of biomass posterior predictive distribution, and (d) range between the upper and lower \(95 \%\) percentiles of the posterior predictive distribution.}
\end{figure}

\section*{6 Summary and conclusions}

We regard the NNGP as a highly scalable model, rather than a likelihood approximation, for large geostatistical datasets. It significantly outperforms competing low-rank processes such as the GPP, in terms of inferential capabilities as well as scalability. A reference set \(\mathcal{S}\) and the resulting neighbor sets (of size \(m\) ) define the NNGP. Larger \(m\) 's would increase costs, but there is no apparent benefit to increasing \(m\) for larger datasets (see Appendix F). Selecting \(\mathcal{S}\) is akin to choosing the "knots" or "centers" in low-rank methods. While some sensitivity to
\(m\) and the choice of points in \(\mathcal{S}\) is expected, our results indicate that inference is very robust with respect to \(\mathcal{S}\) and very modest values of \(m(\ll 20)\) typically suffice. Larger reference sets may be needed for larger datasets, but its size does not thwart computations. In fact, we observed that a very convenient choice for the reference set is the observed locations.

A potential concern with this choice is that if the observed locations have large gaps, then the resulting NNGP may be a poor approximation of the full Gaussian Process. This arises from the fact that observations at locations outside the reference set are correlated via their respective neighbor sets and large gaps may imply two very near points have very different neighbor sets leading to low correlation. Our simulations in Appendix G indeed reveal that in such a situation, the NNGP covariance field is very flat at points in the gap. However, even with this choice of \(\mathcal{S}\) the NNGP model performs at par with the full GP model as the latter also fails to provide strong information about observations located in large gaps. Of course, one can always choose a grid over the entire domain as \(\mathcal{S}\) to construct a NNGP with covariance function similar to the full GP (see Figure 9). Another choice for \(\mathcal{S}\) could be based upon configurations for treed Gaussian processes (Gramacy and Lee 2008). .

Our simulation experiments revealed that estimation and kriging based on NNGP models closely emulate those from the true Matérn GP models, even for slow decaying covariances (see Appendix H). The Matérn covariance function is monotonically decreasing with distance and satisfies theoretical screening conditions, i.e. the ability to predict accurately based on a few neighbors (Stein 2002). This, perhaps, explains the excellent performance of NNGP models with Matérn covariances. We also investigated the performance of NNGP models using a wave covariance function, which does not satisfy the screening conditions, in a setting where a significant proportion of nearest neighbors had negative correlation with the corresponding locations. The NNGP estimates were still close to the true model parameters and the kriged surface closely resembled the true surface (see Appendix II).

Most wave covariance functions (like the damped cosine or the cardinal sine function)
produce covariance matrices with several small eigenvalues. The full GP model cannot be implemented for such models because the matrix inversion is numerically unstable. The NNGP model involves much smaller matrix inversions and can be implemented in some cases (e.g. for the damped cosine model). However, for the cardinal sine covariance, the NNGP also faces numerical issues as even the small \(m \times m\) covariance matrices are numerically unstable. Bias-adjusted low-rank GPs (Finley et al. 2009) possess a certain advantage in this aspect as the covariance matrix is guaranteed to have eigen values bounded away from zero. However, computations involving low-rank processes with numerically unstable covariance functions cannot be carried out with the efficient Sherman-Woodbury-Morrison type matrix identities and more expensive full Cholesky decompositions will be needed.

Apart from being easily extensible to multivariate and spatiotemporal settings with discretized time, the NNGP can fuel interest in process-based modeling over graphs. Examples include networks, where data arising from nodes are posited to be similar to neighboring nodes. It also offers new modeling avenues and alternatives to the highly pervasive Markov random field models for analyzing regionally aggregated spatial data. Also, there is scope for innovation when space and time are jointly modeled as processes using spatiotemporal covariance functions. One will need to construct neighbor sets both in space and time and effective strategies, in terms of scalability and inference, will need to be explored. Comparisons with alternate approaches (see, e.g., Katzfuss and Cressie 2012) will also need to be made. Finally, a more comprehensive study on the alternate algorithms, including direct methods for executing sparse Cholesky factorizations, in Section 4 is being undertaken. More immediately, we plan to migrate our lower-level C++ code to the existing spBayes package (Finley et al. | 2013) in the \(R\) statistical environment (http://cran.r-project.org/web/packages/spBayes) to facilitate wider user accessibility to NNGP models.

Acknowledgments: We express our gratitude to Professors Michael Stein and Noel Cressie for discussions which helped to enrich this work. The work of the first three authors was
supported by federal grants NSF/DMS 1106609 and NIH/NIGMS RC1-GM092400-01.

\section*{Appendix}

\section*{A Densities on directed acyclic graphs}

We will show that if \(\mathcal{G}=\left(\mathcal{S}, N_{\mathcal{S}}\right)\) is acyclic then \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) defined in (3) corresponds to a true density over \(\mathcal{S}\). For any directed acyclic graph, there exists a node with zero in-degree i.e. no directed edge pointing towards it. We denote this node by \(\mathbf{s}_{\pi(1)}\) This means \(\mathbf{s}_{\pi(1)}\) does not belong to the neighbor set of any other location in \(\mathcal{S}\). The only term where it appears on the right hand side of (2) is \(p\left(\mathbf{w}\left(\mathbf{s}_{\pi(1)} \mid \mathbf{w}_{N\left(\mathbf{s}_{\pi(1)}\right)}\right)\right.\) which integrates out to one with respect to \(d \mathbf{w}\left(\mathbf{s}_{\pi(1)}\right)\). We now have a new acyclic directed graph \(\mathcal{G}_{1}\) obtained by removing vertex \(\mathbf{s}_{\pi(1)}\) and its directed edges from \(\mathcal{G}\). Now we can find a new vertex \(\mathbf{s}_{\pi(2)}\) with zero out-degree in \(\mathcal{G}_{1}\) and continue as before to get a permutation \(\pi(1), \pi(2), \ldots, \pi(k)\) of \(1,2, \ldots, k\) such that
\[
\int \prod_{i=1}^{k} p\left(\mathbf{w}\left(\mathbf{s}_{i}\right) \mid \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\right) d \mathbf{w}\left(\mathbf{s}_{\pi(1)}\right) d \mathbf{w}\left(\mathbf{s}_{\pi(2)}\right) \ldots d \mathbf{w}\left(\mathbf{s}_{\pi(k)}\right)=1
\]

An easy application of Fubini's theorem now ensures that this is a proper joint density.

\section*{B Properties of \(\tilde{\mathrm{C}}_{\mathcal{S}}^{-1}\)}

If \(p\left(\mathbf{w}_{\mathcal{S}}\right)=N\left(\mathbf{w}_{\mathcal{S}} \mid \mathbf{0}, \mathbf{C}_{\mathcal{S}}\right)\), then \(\mathbf{w}\left(\mathbf{s}_{i}\right) \mid \mathbf{w}_{N\left(\mathbf{s}_{i}\right)} \sim N\left(\mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}, \mathbf{F}_{\mathbf{s}_{i}}\right)\), where \(\mathbf{B}_{\mathbf{s}_{i}}\) and \(\mathbf{F}_{\mathbf{s}_{i}}\) are defined in (3). So, the likelihood in (2) is proportional to
\[
\frac{1}{\prod_{i=1}^{k} \sqrt{\operatorname{det}\left(\mathbf{F}_{\mathbf{s}_{i}}\right)}} \exp \left(-\frac{1}{2} \sum_{i=1}^{k}\left(\mathbf{w}\left(\mathbf{s}_{i}\right)-\mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\right)^{\prime} \mathbf{F}_{\mathbf{s}_{i}}^{-1}\left(\mathbf{w}\left(\mathbf{s}_{i}\right)-\mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\right)\right)
\]

For any matrix \(\mathbf{A}\), let \(\mathbf{A}\left[, j: j^{\prime}\right]\) denote the submatrix formed using columns \(j\) to \(j^{\prime}\) where \(j<j^{\prime}\). For \(j=1,2, \ldots, k\), we define \(q \times q\) blocks \(\mathbf{B}_{\mathbf{s}_{i}, j}\) as
\[
\mathbf{B}_{\mathbf{s}_{i}, j}=\left\{\begin{array}{l}
\mathbf{I}_{q} \text { if } j=i \\
-\mathbf{B}_{\mathbf{s}_{i}}[,(l-1) q+1: l q] \text { if } \mathbf{s}_{j}=N\left(\mathbf{s}_{i}\right)(l) \text { for some } l \\
\mathbf{O} \text { otherwise }
\end{array}\right.
\]
where, for any location \(\mathbf{s}, N(\mathbf{s})(l)\) is the \(l\)-th neighbor of \(\mathbf{s}\). So, \(\mathbf{w}_{\mathbf{s}_{i}}-\mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}=\mathbf{B}_{\mathbf{s}_{i}}^{*} \mathbf{w}_{\mathcal{S}}\), where \(\mathbf{B}_{\mathbf{s}_{i}}^{*}=\left[\mathbf{B}_{\mathbf{s}_{i}, 1}, \mathbf{B}_{\mathbf{s}_{i}, 2}, \ldots, \mathbf{B}_{\mathbf{s}_{i}, k}\right]\) is \(q \times k q\) and sparse with at most \(m+1\) non-zero blocks. Then,
\[
\sum_{i=1}^{k}\left(\mathbf{w}\left(\mathbf{s}_{i}\right)-\mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\right)^{\prime} \mathbf{F}_{\mathbf{s}_{i}}^{-1}\left(\mathbf{w}\left(\mathbf{s}_{i}\right)-\mathbf{B}_{\mathbf{s}_{i}} \mathbf{w}_{N\left(\mathbf{s}_{i}\right)}\right)=\sum_{i=1}^{k} \mathbf{w}_{\mathcal{S}}^{\prime}\left(\mathbf{B}_{\mathbf{s}_{i}}^{*}\right)^{\prime} \mathbf{F}_{\mathbf{s}_{i}}^{-1} \mathbf{B}_{\mathbf{s}_{i}}^{*} \mathbf{w}_{\mathcal{S}}=\mathbf{w}_{\mathcal{S}}^{\prime} \mathbf{B}_{\mathcal{S}}^{\prime} \mathbf{F}_{\mathcal{S}}^{-1} \mathbf{B}_{\mathcal{S}} \mathbf{w}_{\mathcal{S}}
\]
where \(\mathbf{F}=\operatorname{diag}\left(\mathbf{F}_{\mathbf{s}_{1}}, \mathbf{F}_{\mathbf{s}_{2}}, \ldots, \mathbf{F}_{\mathbf{s}_{k}}\right)\) and \(\mathbf{B}_{\mathcal{S}}=\left(\left(\mathbf{B}_{\mathbf{s}_{1}}^{*}\right)^{\prime},\left(\mathbf{B}_{\mathbf{s}_{2}}^{*}\right)^{\prime}, \ldots,\left(\mathbf{B}_{\mathbf{s}_{k}}^{*}\right)^{\prime}\right)^{\prime}\). So, we have:
\[
\left(\tilde{\mathbf{C}}_{\mathcal{S}}\right)^{-1}=\mathbf{B}_{\mathcal{S}}^{\prime} \mathbf{F}_{\mathcal{S}}^{-1} \mathbf{B}_{\mathcal{S}}
\]

From the form of \(\mathbf{B}_{\mathbf{s}_{i}, j}\), it is clear that \(\mathbf{B}_{\mathcal{S}}\) is sparse and lower triangular with ones on the diagonals. So, \(\operatorname{det}\left(\mathbf{B}_{\mathcal{S}}\right)=1, \operatorname{det}\left(\left(\mathbf{B}_{\mathcal{S}}^{\prime} \mathbf{F}_{\mathcal{S}}^{-1} \mathbf{B}_{\mathcal{S}}\right)^{-1}\right)=\prod \operatorname{det}\left(\mathbf{F}_{\mathbf{s}_{i}}\right)\) and (2) simplifies to \(N\left(\mathbf{w}_{\mathcal{S}} \mid \mathbf{0}, \tilde{\mathbf{C}}_{\mathcal{S}}\right)\).

Let \(\tilde{\mathbf{C}}_{\mathcal{S}}^{i j}\) denote the \((i, j)^{\text {th }}\) block of \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\). Then from equation (14) we see that for \(i<j\), \(\tilde{\mathbf{C}}_{\mathcal{S}}^{i j}=\sum_{l=j}^{k}\left(\mathbf{B}_{\mathbf{s}_{l}, i}^{*}\right)^{\prime} \mathbf{F}_{\mathbf{s}_{l}}^{-1} \mathbf{B}_{\mathbf{s}_{l}, j}^{*}\). So, \(\tilde{\mathbf{C}}_{\mathcal{S}}^{i j}\) is non-zero only if there exists at least one location \(\mathbf{s}_{l}\) such that \(\mathbf{s}_{i} \in N\left(\mathbf{s}_{l}\right)\) and \(\mathbf{s}_{j}\) is either equal to \(\mathbf{s}_{l}\) or is in \(N\left(\mathbf{s}_{l}\right)\). Since every neighbor set has at most \(m\) elements, there are at most \(k m(m+1) / 2\) such pairs \((i, j)\). This demonstrates the sparsity of \(\tilde{\mathbf{C}}_{\mathcal{S}}^{-1}\) for \(m \ll k\).

\section*{C Simulation Experiment: Robustness of NNGP to ordering of locations}

We conduct a simulation experiment demonstrating the robustness of NNGP to the ordering of the locations. We generate the data for \(n=2500\) locations using the model in Section 5.1. However instead of a square domain we choose a long skinny domain (see Figure 5(a)) which can bring out possible sensitivity to ordering due to scale disparity between the \(x\) and \(y\) axes. We use three different orderings for the locations: ordering by \(x\)-coordinates, by \(y\)-coordinates and by the function \(f(x, y)=x+y\).

Table 3 demonstrates that the point estimates and the \(95 \%\) credible intervals for the process parameters from all three NNGP models are extremely consistent with the estimates from the full Gaussian process model.

Posterior estimates of the spatial residual surface from the different models are shown in Figure 5. Again, the impact of the different ordering is negligible. As one of the reviewers suggested, we also plotted the difference between the posterior estimates of the random effects of the true GP and NNGP for all 3 orderings in Figure 6. It was seen that this difference was negligible compared to the difference between the true spatial random effects and full GP estimates. This shows the inference obtained from the NNGP (using any ordering) closely emulates the corresponding full GP inference.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3: Univariate synthetic data analysis parameter estimates and computing time in minutes for NNGP \(m=10\) and full GP models. Parameter posterior summary \(50(2.5,97.5)\) percentiles.}
\begin{tabular}{|l|l|l|l|l|l|}
\hline & True & Full Gaussian Process & Order by \(y\)-coordinates & NNGP ( \(\mathcal{S}=\mathcal{T}\) ) Order by \(x\)-coordinates & Order by \(x+y\)-coordinates \\
\hline \(\sigma^{2}\) & 1 & 0.640 (0.414, 1.297) & 0.712 (0.449, 1.530) & 0.757 (0.479, 1.501) & 0.718 (0.464, 1.436) \\
\hline \(\tau^{2}\) & 0.1 & 0.107 (0.098, 0.117) & 0.106 (0.097, 0.114) & 0.107 (0.099, 0.117) & 0.107 (0.098, 0.115) \\
\hline \(\phi\) & 6 & 8.257 (4.056, 13.408) & 8.294 (3.564, 12.884) & 7.130 (3.405, 11.273) & 7.497 (3.600, 11.911) \\
\hline
\end{tabular}
\end{table}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-33.jpg?height=1655&width=1533&top_left_y=384&top_left_x=264}
\captionsetup{labelformat=empty}
\caption{Figure 5: Robustness of NNGP to ordering: Figures (a) and (b) show interpolated surfaces of the true spatial random effects and posterior median estimates for full geostatistical model respectively. Figures (c), (d), and (e) show interpolated surfaces of the posterior median estimates for NNGP model with \(\mathcal{S}=\mathcal{T}, m=10\), and alternative coordinate ordering. Corresponding true and estimated process parameters are given in Table 3.}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-34.jpg?height=2033&width=1611&top_left_y=268&top_left_x=236}
\captionsetup{labelformat=empty}
\caption{Figure 6: Difference between Full GP and NNGP estimates of spatial effects: Figure (a) shows the difference between the true spatial random effects and the full GP posterior median estimates. Figures (b), (c) and (d) plots the difference between posterior median estimates of full GP and NNGP ordered by \(x, y\) and \(\mathfrak{B} 4^{+} y\) co-ordinates respectively. All the figures are in the same color scale.}
\end{figure}

\section*{D Kolmogorov Consistency for NNGP}

Let \(\{\mathbf{w}(\mathbf{s}) \mid \mathbf{s} \in \mathcal{D}\}\) be a random process over some domain \(\mathcal{D}\) with density \(p\) and let \(\tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\) be a probability density for observations over a fixed finite set \(\mathcal{S} \subset \mathcal{D}\). The conditional density \(\tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right)\) for any finite set \(\mathcal{U} \subset \mathcal{D}\) outside of \(\mathcal{S}\) is defined in (4).

We will first show that for every finite set \(\mathcal{V}=\left\{\mathbf{v}_{1}, \mathbf{v}_{2}, \ldots, \mathbf{v}_{n}\right\}\) in \(\mathcal{D}, n \in\{1,2, \ldots\}\) and for every permutation \(\pi(1), \pi(2), \ldots, \pi(n)\) of \(1,2, \ldots, n\) we have,
\[
\tilde{p}\left(\mathbf{w}\left(\mathbf{v}_{1}\right), \mathbf{w}\left(\mathbf{v}_{2}\right), \ldots, \mathbf{w}\left(\mathbf{v}_{n}\right)\right)=\tilde{p}\left(\mathbf{w}\left(\mathbf{v}_{\pi(1)}\right), \mathbf{w}\left(\mathbf{v}_{\pi(2)}\right), \ldots, \mathbf{w}\left(\mathbf{v}_{\pi(n)}\right)\right)
\]
. We begin by showing that for any finite set \(\mathcal{V}\), the expression given in (5) is a proper density. Let \(\mathcal{U}=\mathcal{V} \backslash \mathcal{S}\). Since \(\mathcal{V} \cup(\mathcal{S} \backslash \mathcal{V})=\mathcal{S} \cup \mathcal{U}\), we obtain
\[
\begin{aligned}
& \int \tilde{p}\left(\mathbf{w}_{\mathcal{V}}\right) \prod_{\mathbf{v}_{i} \in \mathcal{V}} d\left(\mathbf{w}\left(\mathbf{v}_{i}\right)\right)=\int \tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right) \tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right) \prod_{\mathbf{v}_{i} \in \mathcal{U}} d\left(\mathbf{w}\left(\mathbf{v}_{i}\right)\right) \prod_{\mathbf{s}_{i} \in \mathcal{S}} d\left(\mathbf{w}\left(\mathbf{s}_{i}\right)\right) \\
= & \int \tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right)\left(\int \tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right) \prod_{\mathbf{v}_{i} \in \mathcal{U}} d\left(\mathbf{w}\left(\mathbf{v}_{i}\right)\right)\right) \prod_{\mathbf{s}_{i} \in \mathcal{S}} d\left(\mathbf{w}\left(\mathbf{s}_{i}\right)\right)=\int \tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right) \prod_{\mathbf{s}_{i} \in \mathcal{S}} d\left(\mathbf{w}\left(\mathbf{s}_{i}\right)\right)=1
\end{aligned}
\]

Note that \(\mathcal{S}\) is fixed. Therefore, the expression for the joint density of \(\mathbf{w}_{\mathcal{V}}\) depends only on the the neighbor sets \(N\left(\mathbf{v}_{i}\right)\) for \(\mathbf{v}_{i} \in \mathcal{U}\). So the NNGP density for \(\mathcal{V}\) is invariant under any permutation of locations inside \(\mathcal{V}\).

We now prove that for every location \(\mathbf{v}_{0} \in \mathcal{D}\), we have, \(\tilde{p}\left(\mathbf{w}_{\mathcal{V}}\right)=\int \tilde{p}\left(\mathbf{w}_{\mathcal{V} \cup\left\{\mathbf{v}_{0}\right\}}\right) d\left(\mathbf{w}\left(\mathbf{v}_{0}\right)\right)\). let \(\mathcal{V}_{1}=\mathcal{V} \cup\left\{\mathbf{v}_{0}\right\}\). We split the proof into two cases. If \(\mathbf{v}_{0} \in \mathcal{S}\), then using the fact \(\mathcal{V}_{1} \backslash \mathcal{S}=\mathcal{V} \backslash \mathcal{S}=\mathcal{U}\), we obtain
\[
\begin{aligned}
\int \tilde{p}\left(\mathbf{w}_{\mathcal{V}_{1}}\right) d\left(\mathbf{w}\left(\mathbf{v}_{0}\right)\right) & =\int \tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right) \tilde{p}\left(\mathbf{w}_{\mathcal{V}_{1} \backslash \mathcal{S}} \mid \mathbf{w}_{\mathcal{S}}\right) \prod_{\mathbf{s}_{i} \in \mathcal{S} \backslash \mathcal{V}_{1}} d\left(\mathbf{w}\left(\mathbf{s}_{i}\right)\right) d\left(\mathbf{w}\left(\mathbf{v}_{0}\right)\right. \\
& =\int \tilde{p}\left(\mathbf{w}_{\mathcal{S}}\right) \tilde{p}\left(\mathbf{w}_{\mathcal{V} \backslash \mathcal{S}} \mid \mathbf{w}_{\mathcal{S}}\right) \prod_{\mathbf{s}_{i} \in \mathcal{S} \backslash \mathcal{V}} d\left(\mathbf{w}\left(\mathbf{s}_{i}\right)\right)=\tilde{p}\left(\mathbf{w}_{\mathcal{U}}\right)
\end{aligned}
\]

If \(\mathbf{v}_{0} \notin \mathcal{S}\), then \(\mathbf{w}\left(\mathbf{v}_{0}\right)\) does not appear in the neighborhood set of any other term. So, \(p\left(\mathbf{w}\left(\mathbf{v}_{0}\right) \mid \mathbf{w}_{\mathcal{S}}\right)\) integrates to one with respect to \(d\left(\mathbf{w}\left(\mathbf{v}_{0}\right)\right)\). The result now follows from \(\int p\left(\mathbf{w}_{\mathcal{V}_{1}} \mid \mathbf{w}_{\mathcal{S}}\right) d\left(\mathbf{w}\left(\mathbf{v}_{0}\right)\right)=p\left(\mathbf{w}_{\mathcal{V}} \mid \mathbf{w}_{\mathcal{S}}\right)\).

\section*{E Properties of NNGP}

Standard Gaussian conditional distribution facts reveal that the conditional distribution \(\mathbf{w}\left(\mathbf{u}_{i}\right) \mid \mathbf{w}_{\mathcal{S}} \sim N\left(\mathbf{B}_{\mathbf{u}_{i}} \mathbf{w}_{N\left(\mathbf{u}_{i}\right)}, \mathbf{F}_{\mathbf{u}_{i}}\right)\) where \(\mathbf{B}_{\mathbf{u}_{i}}\) and \(\mathbf{F}_{\mathbf{u}_{i}}\) be defined analogous to (3) based on the neighbor sets \(N\left(\mathbf{u}_{i}\right)\). From (4), we see that
\(\tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right)=\frac{1}{\prod_{i=1}^{r} \sqrt{\operatorname{det}\left(\mathbf{F}_{\mathbf{u}_{i}}\right)}} \exp \left(-\frac{1}{2} \sum_{i=1}^{r}\left(\mathbf{w}\left(\mathbf{u}_{i}\right)-\mathbf{B}_{\mathbf{u}_{i}} \mathbf{w}_{N\left(\mathbf{u}_{i}\right)}\right)^{\prime} \mathbf{F}_{\mathbf{u}_{i}}^{-1}\left(\mathbf{w}\left(\mathbf{u}_{i}\right)-\mathbf{B}_{\mathbf{u}_{i}} \mathbf{w}_{N\left(\mathbf{u}_{i}\right)}\right)\right)\)
It then follows that \(\tilde{p}\left(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\right) \sim N\left(\mathbf{B}_{\mathcal{U}} \mathbf{w}_{\mathcal{S}}, \mathbf{F}_{\mathcal{U}}\right)\) where \(\mathbf{B}_{\mathcal{U}}=\left(\mathbf{B}_{\mathbf{u}_{1}}^{\prime}, \mathbf{B}_{\mathbf{u}_{2}}^{\prime}, \ldots, \mathbf{B}_{\mathbf{u}_{r}}^{\prime}\right)^{\prime}\) and \(\mathbf{F}_{\mathcal{U}}= \operatorname{diag}\left(\mathbf{F}_{\mathbf{u}_{1}}, \mathbf{F}_{\mathbf{u}_{2}}, \ldots, \mathbf{F}_{\mathbf{u}_{r}}\right)\). Since each row of \(\mathbf{B}_{\mathcal{U}}\) has at most \(m\) non-zero entries, \(\mathbf{B}_{\mathcal{U}}\) is sparse for \(m \ll k\).

As the nearest neighbor densities of \(\mathbf{w}_{\mathcal{S}}\) and \(\mathbf{w}_{\mathcal{U}} \mid \mathbf{w}_{\mathcal{S}}\) for every finite \(\mathcal{U}\) outside \(\mathcal{S}\) are Gaussian, all finite dimensional realizations of an NNGP process will be Gaussian. Let \(\mathbf{v}_{1}\) and \(\mathbf{v}_{2}\) be any two locations in \(\mathcal{D}\) and let \(\widetilde{E}\) and \(\widetilde{C o v}\) denote, respectively, the expectation and covariance operator for a NNGP. Then, if \(\mathbf{v}_{1}=\mathbf{s}_{i}\) and \(\mathbf{v}_{2}=\mathbf{s}_{j}\) are both in \(\mathcal{S}\) then we obviously have \(\widetilde{\operatorname{Cov}}\left(\mathbf{w}\left(\mathbf{v}_{1}\right), \mathbf{w}\left(\mathbf{v}_{2}\right) \mid \boldsymbol{\theta}\right)=\tilde{\mathbf{C}}_{\mathbf{s}_{i}, \mathbf{s}_{j}}\). If \(\mathbf{v}_{1}\) is outside \(\mathcal{S}\) and \(\mathbf{v}_{2}=\mathbf{s}_{j}\), then
\[
\begin{aligned}
\widetilde{\operatorname{Cov}}\left(\mathbf{w}\left(\mathbf{v}_{1}\right), \mathbf{w}\left(\mathbf{v}_{2}\right) \mid \boldsymbol{\theta}\right) & \left.=\widetilde{E}\left(\widetilde{\operatorname{Cov}}\left(\mathbf{w}\left(\mathbf{v}_{1}\right), \mathbf{w}\left(\mathbf{v}_{2}\right) \mid \mathbf{w}_{\mathcal{S}}, \boldsymbol{\theta}\right)\right)+\widetilde{\operatorname{Cov}}\left(\widetilde{E}\left(\mathbf{w}\left(\mathbf{v}_{1}\right)\right), \widetilde{E}\left(\mathbf{w}\left(\mathbf{v}_{2}\right)\right) \mid \mathbf{w}_{\mathcal{S}}, \boldsymbol{\theta}\right)\right) \\
\therefore \tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} \mid \boldsymbol{\theta}\right) & =0+\widetilde{\operatorname{Cov}}\left(\mathbf{B}_{\mathbf{v}_{1}} \mathbf{w}_{N\left(\mathbf{v}_{1}\right)}, \mathbf{w}\left(\mathbf{s}_{j}\right) \mid \boldsymbol{\theta}\right)=\mathbf{B}_{\mathbf{v}_{1}} \tilde{\mathbf{C}}_{N\left(\mathbf{v}_{1}\right), \mathbf{s}_{j}}
\end{aligned}
\]

If both \(\mathbf{v}_{1}\) and \(\mathbf{v}_{2}\) are outside \(\mathcal{S}\), then \(\tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} \mid \boldsymbol{\theta}\right)=\delta\left(\mathbf{v}_{1}=\mathbf{v}_{2}\right) \mathbf{F}_{\mathbf{v}_{1}}+\mathbf{B}_{\mathbf{v}_{1}} \tilde{\mathbf{C}}_{N\left(\mathbf{v}_{1}\right), N\left(\mathbf{v}_{2}\right)} \mathbf{B}_{\mathbf{v}_{2}}^{\prime}\), which yields (7).

For any two set of locations \(A\) and \(B\), let \(\|A, B\|\) denote the pairwise Euclidean distance matrix. Let \(\mathcal{Z}_{1}\) denote set of all points \(\mathbf{v}\) such that \(\mathbf{v}\) is equidistant from any two points in
\(\mathcal{S}\). Since \(\mathcal{S}\) is finite, the set \(\mathcal{Z}_{2}=\left(\mathcal{Z}_{1} \times \mathcal{Z}_{1}\right) \cup\{(\mathbf{v}, \mathbf{v}) \mid \mathbf{v} \in \mathcal{D}\}\) has Lebesgue measure zero in the Euclidean domain \(\Re^{d} \times \Re^{d}\). We will show that \(\tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} \mid \boldsymbol{\theta}\right)\) is continuous for any pair \(\left(\mathbf{v}_{1}, \mathbf{v}_{2}\right)\) in \(\mathcal{D} \times \mathcal{D} \backslash \mathcal{Z}_{2}\). Observe that for any pair of points ( \(\mathbf{v}_{1}, \mathbf{v}_{2}\) ) in \(\mathcal{D} \times \mathcal{D} \backslash \mathcal{Z}_{2}\), it is easy to verify that \(\lim _{\mathbf{h}_{i} \rightarrow 0} \|\left(\mathbf{v}_{i}+\mathbf{h}_{i}, N\left(\mathbf{v}_{i}+\mathbf{h}_{i}\right)\|\rightarrow\| \mathbf{v}_{i}, N\left(\mathbf{v}_{i}\right) \|\right.\), for \(i=1,2\), and \(\lim _{\mathbf{h}_{1} \rightarrow 0, \mathbf{h}_{2} \rightarrow 0} \| N\left(\mathbf{v}_{1}+\right. \left.\mathbf{h}_{1}\right), N\left(\mathbf{v}_{2}+\mathbf{h}_{2}\right)\|\rightarrow\| N\left(\mathbf{v}_{1}\right), N\left(\mathbf{v}_{2}\right) \|\). We prove the continuity of \(\tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} \mid \boldsymbol{\theta}\right)\) for the case when \(\mathbf{v}_{1}\) is outside \(\mathcal{S}\) and \(\mathbf{v}_{2}=\mathbf{s}_{j}\). The other cases are proved similarly. We assume that the covariance function for the original GP is isotropic and continuous. The two distance results yield \(\mathbf{B}_{\mathbf{v}_{1}+\mathbf{h}_{1}}=\mathbf{C}_{\mathbf{v}_{1}+\mathbf{h}_{1}, N\left(\mathbf{v}_{1}+\mathbf{h}_{1}\right)} \mathbf{C}_{N\left(\mathbf{v}_{1}+\mathbf{h}_{1}\right)}^{-1} \rightarrow \mathbf{C}_{\mathbf{v}_{1}, N\left(\mathbf{v}_{1}\right)} \mathbf{C}_{N\left(\mathbf{v}_{1}\right)}^{-1}=\mathbf{B}_{\mathbf{v}_{1}}\). Also, as \(\mathbf{v}_{2}+\mathbf{h}_{2} \rightarrow \mathbf{v}_{2}=\mathbf{s}_{j}\), then \(\mathbf{s}_{j} \in N\left(\mathbf{v}_{2}+\mathbf{h}_{2}\right)\) for small enough \(\mathbf{h}_{2}\). Let \(\mathbf{s}_{j}=N\left(\mathbf{v}_{2}+\mathbf{h}_{2}\right)\) (1) and, hence, \(\mathbf{C}_{\mathbf{v}_{2}+\mathbf{h}_{2}, N\left(\mathbf{v}_{2}+\mathbf{h}_{2}\right)} \mathbf{C}_{N\left(\mathbf{v}_{2}+\mathbf{h}_{2}\right)}^{-1} \rightarrow \mathbf{e}_{1}\) where \(\mathbf{e}_{1}=(1,0, \ldots, 0)_{m \times 1}\). Therefore,
\[
\begin{aligned}
\lim _{\mathbf{h}_{1} \rightarrow 0, \mathbf{h}_{2} \rightarrow 0} & \tilde{\mathbf{C}}\left(\mathbf{v}_{1}+\mathbf{h}_{1}, \mathbf{v}_{2}+\mathbf{h}_{2} \mid \boldsymbol{\theta}\right)=\mathbf{B}_{\mathbf{v}_{1}} \lim _{\mathbf{h}_{1} \rightarrow \mathbf{0}, \mathbf{h}_{2} \rightarrow \mathbf{0}} \widetilde{\operatorname{Cov}}\left(\mathbf{w}_{N\left(\mathbf{v}_{1}+\mathbf{h}_{1}\right)}, \mathbf{w}_{N\left(\mathbf{v}_{2}+\mathbf{h}_{2}\right)} \mid \boldsymbol{\theta}\right) \mathbf{e}_{1} \\
& =\mathbf{B}_{\mathbf{v}_{1}} \lim _{\mathbf{h}_{1} \rightarrow \mathbf{0}} \widetilde{\operatorname{Cov}}\left(\mathbf{w}_{N\left(\mathbf{v}_{1}+\mathbf{h}_{1}\right)}, \mathbf{w}\left(\mathbf{s}_{j}\right) \mid \boldsymbol{\theta}\right)=\mathbf{B}_{\mathbf{v}_{1}} \tilde{\mathbf{C}}_{N\left(\mathbf{v}_{1}\right), \mathbf{s}_{j}}=\tilde{\mathbf{C}}\left(\mathbf{v}_{1}, \mathbf{v}_{2} \mid \boldsymbol{\theta}\right)
\end{aligned}
\]

\section*{F Simulation experiment: NNGP credible intervals as function of \(m\)}

From a classical viewpoint, NNGP can be regarded as a computationally convenient approximation to the full GP model. The accuracy of the approximation is expected to improve with increase in \(m\) as NNGP model becomes identical to the full model when \(m\) equals the sample size. However, we construct the NNGP as an independent model and found that inference from this model closely emulates that from the full GP model. Figure 1 demonstrates how root mean square predictive error and parameter CI width vary with choice of \(m\). We conduct another simulation experiment to investigate how the parameter estimation of the hierarchical NNGP model depends on \(m\).

We generated a dataset of size 1000 using the model described in Section 5.1 for 4 combination of values of \(\phi\) and \(\sigma^{2}\). Other parameters and prior choices were similar to those in section H. Figure 7 gives true values of \(\sigma^{2}\) and effective range ( \(3 / \phi\) ) alongwith the posterior medians and credible intervals for the full GP, NNGP with \(m=10\) and \(m=100\). We see that the CI's for NNGP \(m=10\) and \(m=100\) are almost identical and are very close to the CI for full GP. This suggests that even for small values of \(m\) NNGP, parameter CI's closely resemble full GP parameter CI's.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-39.jpg?height=1402&width=1514&top_left_y=689&top_left_x=244}
\captionsetup{labelformat=empty}
\caption{Figure 7: NNGP credible intervals for small and large values of \(m\)}
\end{figure}

\section*{G Simulation experiment: Data with gaps}

One possible area of concern for NNGP is that if the data have large gaps and the NNGP is constructed using the data locations as the reference set \(\mathcal{S}\), then NNGP covariance function may be a poor approximation of the full GP covariance function. This arises from the fact that if the reference set has large gaps then two very close locations outside \(\mathcal{S}\) can have very different neighbor sets. Since, in a NNGP, locations outside \(\mathcal{S}\) are correlated through their neighbors sets this may lead to little correlation among very close points in certain regions of the domain.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-40.jpg?height=1120&width=1647&top_left_y=1050&top_left_x=236}
\captionsetup{labelformat=empty}
\caption{Figure 8: Full GP and NNGP ( \(m=10\) ) covariance function for data with gaps}
\end{figure}

Figure 8 demonstrates this issue. We generate a set \(\mathcal{T}\) of 100 locations (topleft) on the domain \([0,3] \times[0,1]\) with half the locations in \([0,1] \times[0,1]\) and the remaining half in \([2,3] \times[0,1]\). This creates a large gap in the middle where there are no datapoints. The topright panel shows the heatmap of the full GP covariance function with \(\sigma^{2}=1\) and \(\phi=2\) (so that the effective range is 1.5 ). The NNGP is a non-stationary process and the covariance function depends on the locations. We evaluate this covariance at two points (red dots in the topleft figure \()-(0.5,0.5)\) (which is surrounded by many points in \(\mathcal{S}\) ) and ( \(1.5,0.5\) ) (which is at the middle of the gap and equidistant from the two sets of locations in \(\mathcal{S}\) ).

The NNGP field at ( \(0.5,0.5\) ) (bottomleft) closely resembles the GP field. This is because the neighbors of ( \(0.5,0.5\) ) are close to the point and provides strong information about the true GP at that point. The NNGP field at (1.5, 0.5) (bottomright) is almost non-existent with near zero correlations even at very small distances. This is an expected consequence of the way NNGP is constructed. Any two points outside \(\mathcal{S}\) are correlated via their neigbhor sets only. The neighbors for \((1.5,0.5)\) are far away from the point it provides weak information about the point as it is in the middle of the gap.

This suggests that a NNGP constructed using a reference set with large gaps is a poor approximation to the full GP as a process in certain regions of the domain. If the data locations do have large gaps, perhaps a NNGP with \(\mathcal{S}\) as a grid over the domain provides a much better approximation to the full GP. To observe this we used a \(14 \times 7\) grid over the domain \([0,3] \times[0,1]\) as \(\mathcal{S}\). So the size of this new \(\mathcal{S}\) was similar to the original sample size of 100. Figure 9 demonstrates the NNGP covariance function at the two points using this new \(\mathcal{S}\). We see that using the grid \(\mathcal{S}\), the NNGP covariance function at the two points are very similar and closely resemble the true GP covariance function. This suggests that in order for the NNGP to resemble full GP, the reference set needs to have points uniformly distributed over the domain.

However, from a kriging perspective, if the data have large gaps, inference from a NNGP

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-42.jpg?height=662&width=1408&top_left_y=320&top_left_x=423}
\captionsetup{labelformat=empty}
\caption{Figure 9: NNGP covariance function using a grid \(\mathcal{S}\)}
\end{figure}
with \(\mathcal{S}=\mathcal{T}\) may not differ a lot from the full GP inference. Even when one uses the full GP, kriging is usually done one point at a time and thereby ignores the covariances between points outside the data locations and assumes conditional independence. Figure 10 plots the

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-42.jpg?height=690&width=1445&top_left_y=1452&top_left_x=386}
\captionsetup{labelformat=empty}
\caption{Figure 10: Kriging means and variances for full GP and NNGP ( \(\mathcal{S}=\) data locations)}
\end{figure}
kriging mean and variances over the entire domain for the full GP and the NNGP. They are

\begin{table}
\begin{tabular}{|l|l|l|l|l|}
\hline & True & Full GP & NNGP \(\mathrm{m}=10\) & NNGP \(\mathrm{m}=20\) \\
\hline \(\beta\) & 1 & 0.72 (0.00, 1.32) & 0.65 (-0.14, 1.30) & 0.69 (0.02, 1.16) \\
\hline \(\tau\) & 0.01 & 0.03 (0.01, 0.05) & 0.03 (0.01, 0.06) & 0.03 (0.01, 0.06) \\
\hline \(\sigma^{2}\) & 1 & 0.63 (0.38, 1.31) & 0.65 (0.39, 1.29) & 0.62 (0.38, 1.27) \\
\hline \(\phi\) & 2 & 2.94 (1.27, 5.19) & 2.76 (1.27, 5.25) & 2.91 (1.34, 5.20) \\
\hline RMSPE & - & 0.58 (ind) & 0.57 & 0.57 \\
\hline & - & 0.58 (joint) & - & - \\
\hline 95\% CI cover & - & 94.00 (ind) & 95.66 & 95.33 \\
\hline & - & 95.33 (joint) & - & - \\
\hline Mean 95\% CI width & - & 2.12 (ind) & 2.12 & 2.13 \\
\hline & - & 2.11 (joint) & - & - \\
\hline
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 4: Data analysis for locations with gaps}
\end{table}
very close. This suggests even for data with gaps the kriging performance of NNGP and GP are similar.

We also generated a dataset over \(\mathcal{T}\) and fitted the full GP and NNGP ( \(\mathcal{S}=\mathcal{T}\) ) model to compare parameter estimation and kriging performance. In addition to the conventional independent kriging, we also used the computationally expensive joint kriging for the full GP to see if it improves kriging quality at locations in the gap. Table 4 provide the parameter estimates and model fitting metrics. Figures 11 and 12 gives the posterior median and the variance surface over the domain. We see that the the NNGP and full GP produce very similar parameter estimates and kriging. Hence, for data with large gaps both the full GP and NNGP \((\mathcal{S}=\mathcal{T})\) doesn't provide enough information for locations inside the gaps. So even if NNGP ( \(\mathcal{S}=\mathcal{T}\) ) poorly approximates the full GP as a process, in terms of model fitting, their performances are very similar.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-44.jpg?height=708&width=1619&top_left_y=408&top_left_x=259}
\captionsetup{labelformat=empty}
\caption{Figure 11: Posterior median surface for data with gaps}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-44.jpg?height=712&width=1621&top_left_y=1516&top_left_x=249}
\captionsetup{labelformat=empty}
\caption{Figure 12: Posterior variance surface for data with gaps}
\end{figure}

\section*{H Simulation experiment: Slow decaying covariance functions}

We note in Section 2.1 that several valid choices of neighbor sets can be used to construct a NNGP. However, our choice of using \(m\)-nearest neighbors to construct neighbor sets performed extremely well for all the data analysis in Section 5. Since, our design of NNGP just includes \(m\)-nearest neighbors it is natural to be skeptical of the performance of NNGP when the data arises from a Gaussian process with very flat tailed covariance function. Such a covariance function implies that even distant observations are significantly correlated with the given observation and \(m\)-nearest neighbors may fail to capture all the information about the covariance parameters.

We generate datasets of size 2500 in a unit domain using the model described in Section 5.1 for a wide range of values for the parameters \(\sigma^{2}\) and \(\phi\). The marginal variance \(\sigma^{2}\) was varied over \((0.05,0.1,0.2,0.5)\) and the 'true effective range' \(3 / \phi\) phi was varied over \((0.1,0.2, \ldots, 1)\). Larger values of the 'true effective range' indicate higher correlation between points at large distances. The nugget variance \(\tau^{2}\) was held constant at 0.1 . The prior on \(\phi\) was \(\mathrm{U}(3,300)\) or 0.01 to 1 distance units. Also both \(\tau^{2}\) and \(\sigma^{2}\) were given Inverse \(\operatorname{Gamma}(2,0.1)\) priors in all cases.

Figure 13 gives the results for NNGP and full GP CIs. We see that for all choices of parameters, the posterior samples from the NNGP and full GP look identical. This strongly suggests that the NNGP model deliver inference similar to that of a full GP even for slow decaying covariance functions and justifies the choice of the neighbor sets.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-46.jpg?height=1670&width=1569&top_left_y=483&top_left_x=249}
\captionsetup{labelformat=empty}
\caption{Figure 13: Univariate synthetic data analysis: true versus posterior \(50 \%(2.5 \%, 97.5 \%)\) percentiles for the effective spatial range simulated for various values of \(\sigma^{2}\) and \(\tau^{2}=0.1\). NNGP model fit with \(\mathcal{S}=\mathcal{T}\) and \(m=10\).}
\end{figure}

\section*{I Simulation experiment: Wave covariance function}

We have restricted most of our simulation experiments to Matérn (or in particular exponential) covariance functions. Matérn covariance functions like many other covariance functions decrease monotonically with distance and hence nearest neighbors of a location have highest correlation with that location. We wanted to investigate the performance of NNGP for covariance functions which do not monotonically decrease with distance. We use the two-dimensional damped cosine covariance function given by:
\[
C(d)=\exp (-d / a) \cos (\phi d), a \leq 1 / \phi
\]

First, we generated the Kullback-Leibler (KL) divergence numbers for the NNGP model with respect to the full GP model using damped cosine covariance. In addition to the default neighbor selection scheme, we also used an alternate scheme described by Stein et al. (2004). This scheme includes \(m^{\prime}=\lceil 0.75 m\rceil\) nearest neighbors and \(m-m^{\prime}\) neighbors whose ranked distances from the \(i^{\text {th }}\) location equal \(m+\left\lfloor l(i-m-1) /\left(m-m^{\prime}\right)\right\rfloor\) for \(l=1,2, \ldots, m-m^{\prime}\). Stein et al. (2004) suggested that this scheme choice often improves parameter estimation. The two schemes are referred to as NNGP and NNGP (alt) respectively. We used \(\phi=10\), \(a=.099\), sample sizes of 100,200 and 500 and varied \(m\) from 5 to 50 in increments of 5 .

Figure 14 plots the KL divergence numbers (in log-scale) for varying \(m, n\) and neighbor selection schemes. We see that larger sample size implies higher KL divergence numbers which is expected as with increasing sample size the size of the neighbor set \(m\) becomes smaller in proportion. Also, we see that KL numbers for the alternate neighbor selection scheme are always higher indicating that nearest neighbors perform better even for such wave covariance functions. In general we observed that the KL numbers are quite small for \(m \geq 25\) for all \(n\) and neighbor selection schemes indicating that the NNGP models closely approximate the true damped cosine GP.
damped cosine: phi= 10, a=0.099

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-48.jpg?height=645&width=1595&top_left_y=363&top_left_x=242}
\captionsetup{labelformat=empty}
\caption{Figure 14: NNGP KL divergence numbers \({ }^{\text {m }}\) (log scale) for damped cosine covariance}
\end{figure}

Next, we conducted a data analysis using the wave covariance function. We choose \(n= 500, m=10,20\). The two values of \(m\) yielded around \(3.4 \%\) and \(18.7 \%\) nearest neighbors which were negatively correlated with the corresponding locations. Table 5 gives the parameter estimates for the NNGP model. Figure 15 demonstrates how the NNGP approximates the wave covariance function while figure 16 plots the true and fitted random effect surface. We observe that NNGP provides an excellent approximation of the the true wave GP in terms of model parameter estimation and kriging.

We could not fit the full GP model due to computation instability of the large wave covariance matrix. NNGP does not involve inverting large matrices and hence we could use it for model fitting.

\begin{table}
\begin{tabular}{cccc} 
& True & \(\mathrm{m}=10\) & \(\mathrm{~m}=20\) \\
\hline\(\beta_{0}\) & 1 & \(1.03(0.65,1.34)\) & \(1.06(0.70,1.32)\) \\
\(\beta_{1}\) & 5 & \(5.00(4.95,5.06)\) & \(5.00(4.95,5.06)\) \\
\(\tau^{2}\) & 0.1 & \(0.06(0.02,0.12)\) & \(0.05(0.03,0.11)\) \\
\(\sigma^{2}\) & 1 & \(1.13(0.90,1.57)\) & \(1.14(0.90,1.57)\) \\
\(\phi\) & 10 & \(7.41(1.63,11.59)\) & \(6.31(1.61,10.50)\) \\
\(a\) & 0.099 & \(0.093(0.067,0.135)\) & \(0.09(0.07,0.14)\)
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 5: Damped cosine GP data analysis using NNGP}
\end{table}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-49.jpg?height=669&width=1557&top_left_y=539&top_left_x=259}
\captionsetup{labelformat=empty}
\caption{Figure 15: Wave covariance function estimates using NNGP}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/d2e826ff-0bff-40ec-bbd4-bb6a2eaa87b8-49.jpg?height=513&width=1509&top_left_y=1677&top_left_x=292}
\captionsetup{labelformat=empty}
\caption{Figure 16: True and estimated (posterior median) random effect surface of the damped cosine GP}
\end{figure}

\section*{References}

Banerjee, S., Carlin, B. P., and Gelfand, A. E. (2014), Hierarchical Modeling and Analysis for Spatial Data, Boca Raton, FL: Chapman \& Hall/CRC, 2nd ed.

Banerjee, S., Finley, A. O., Waldmann, P., and Ericcson, T. (2010), "Hierarchical Spatial Process Models for Multiple Traits in Large Genetic Trials," Journal of the American Statistical Association, 105, 506-521.

Banerjee, S., Gelfand, A. E., Finley, A. O., and Sang, H. (2008), "Gaussian Predictive Process Models for Large Spatial Datasets," Journal of the Royal Statistical Society, Series B, 70, 825-848.

Bechtold, W. A. and Patterson, P. L. (2005), "The Enhanced Forest Inventory and Analysis National Sample Design and Estimation Procedures," SRS-80, U.S. Department of Agriculture, Forest Service, Southern Research Station: Asheville, NC.

Bevilacqua, M. and Gaetan, C. (2014), "Comparing Composite Likelihood Methods Based on Pairs for Spatial Gaussian Random Fields," Statistics and Computing, 1-16.

Crainiceanu, C. M., Diggle, P. J., and Rowlingson, B. (2008), "Bivariate Binomial Spatial Modeling of Loa Loa Prevalence in Tropical Africa," Journal of the American Statistical Association, 103, 21-37.

Cressie, N. A. C. and Johannesson, G. (2008), "Fixed Rank Kriging for Very Large Data Sets," Journal of the Royal Statistical Society, Series B, 70, 209-226.

Cressie, N. A. C. and Wikle, C. K. (2011), Statistics for spatio-temporal data, Hoboken, NJ: Wiley, Wiley Series in Probability and Statistics.

Davis, T. A. (2006), Direct Methods for Sparse Linear Systems, Philadelphia, PA: Society for Industrial and Applied Mathematics.

Diggle, P. J., Tawn, J. A., and Moyeed, R. A. (1998), "Model-based Geostatistics (with discussion)," Applied Statistics, 47, 299-350.

Du, J., Zhang, H., and Mandrekar, V. S. (2009), "Fixed-domain Asymptotic Properties of Tapered Maximum Likelihood Estimators," Annals of Statistics, 37, 3330-3361.

Eidsvik, J., Shaby, B. A., Reich, B. J., Wheeler, M., and Niemi, J. (2014), "Estimation and Prediction in Spatial Models with Block Composite Likelihoods," Journal of Computational and Graphical Statistics, 23, 295-315.

Emory, X. (2009), "The Kriging Update Equations and Their Application to the Selection of Neighboring Data," Computational Geosciences, 13, 269-280.

Finley, A. O., Banerjee, S., and Gelfand, A. E. (2013), "spBayes for Large Univariate and Multivariate Point-referenced Spatio-temporal Data Models," Journal of Statistical Software, 0, In press.

Finley, A. O., Banerjee, S., and McRoberts, R. E. (2009), "Hierarchical Spatial Models for Predicting Tree Species Assemblages across Large Domains," The Annals of Applied Statistics, 0, 1-32.

Furrer, R., Genton, M. G., and Nychka, D. W. (2006), "Covariance Tapering for Interpolation of Large Spatial Datasets," Journal of Computational and Graphical Statistics, 15, 503523.

Gelfand, A. E. and Banerjee, S. (2010), "Multivariate Spatial Process Models," in Handbook of Spatial Statistics, eds. Gelfand, A. E., Diggle, P. J., Fuentes, M., and Guttorp, P., Boca Raton, FL: Chapman \& Hall/CRC, pp. 495-516.

Gelfand, A. E., Banerjee, S., and Gamerman, D. (2005), "Spatial Process Modelling for Univariate and Multivariate Dynamic Spatial Data," Environmetrics, 16, 465-479.

Gelfand, A. E. and Ghosh, S. K. (1998), "Model Choice: A Minimum Posterior Predictive Loss Approach," Biometrika, 85, 1-11.

Gelfand, A. E., Kim, H.-J., Sirmans, C., and Banerjee, S. (2003), "Spatial Modeling with Spatially Varying Coefficient Processes," Journal of the American Statistical Association, 98, 387-396.

Gramacy, R. B. and Apley, D. W. (2014), "Local Gaussian Process Approximation for Large Computer Experiments," http://arxiv.org/abs/1303.0383.

Gramacy, R. B. and Lee, H. (2008), "Bayesian Treed Gaussian Process Models with an Application to Computer Experiments," Journal of the American Statistical Association, 103, 1119-1130.

Gramacy, R. B., Niemi, J., and Weiss, R. M. (2014), "Massively Parallel Approximate Gaussian Process Regression," http://arxiv.org/abs/1310.5182.

Higdon, D. (2001), "Space and Space Time Modeling using Process Convolutions," Technical Report, Institute of Statistics and Decision Sciences, Duke University, Durham.

Kammann, E. E. and Wand, M. P. (2003), "Geoadditive Models," Applied Statistics, 52, 1-18.

Katzfuss, M. and Cressie, N. (2012), "Bayesian Hierarchical Spatio-temporal Smoothing for Very Large Datasets," Environmetrics, 23, 94-107.

Kaufman, C. G., Scheverish, M. J., and Nychka, D. W. (2008), "Covariance Tapering for Likelihood-based Estimation in Large Spatial Data Sets," Journal of the American Statistical Association, 103, 1545-1555.

Lauritzen, S. L. (1996), Graphical Models, Oxford, United Kingdom: Clarendon Press.

Lin, X., Wahba, G., Xiang, D., Gao, F., Klein, R., and Klein, B. (2000), "Smoothing Spline ANOVA Models for Large Data Sets with Bernoulli Observations and the Randomized GACV," Annals of Statistics, 28, 1570-1600.

Moller, J. and Waagepetersen, R. P. (2003), Statistical Inference and Simulation for Spatial Point Processes, Boca Raton, FL: Chapman \& Hall/CRC, 1st ed.

Rasmussen, C. E. and Williams, C. K. I. (2005), Gaussian Processes for Machine Learning, Cambridge, MA: The MIT Press, 1st ed.

Rue, H. and Held, L. (2005), Gaussian Markov Random Fields : Theory and Applications, Boca Raton, FL: Chapman \& Hall/CRC, Monographs on Statistics and Applied Probability.

Sang, H. and Huang, J. Z. (2012), "A Full Scale Approximation of Covariance Functions for Large Spatial Data Sets," Journal of the Royal Statistical Society, Series B, 74, 111-132.

Schabenberger, O. and Gotway, C. A. (2004), Statistical Methods for Spatial Data Analysis, Boca Raton, FL: Chapman \& Hall/CRC, 1st ed.

Shaby, B. A. (2012), "The Open-faced Sandwich Adjustment for MCMC using Estimating Functions," http://arxiv.org/abs/1204.3687.

Shaby, B. A. and Ruppert, D. (2012), "Tapered Covariance: Bayesian Estimation and Asymptotics," Journal of Computational and Graphical Statistics, 21, 433-452.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P., and van der Linde, A. (2002), "Bayesian Measures of Model Complexity and Fit," Journal of the Royal Statistical Society, Series B, 64, 583-639.

Stein, M. L. (1999), Interpolation of Spatial Data: Some Theory for Kriging, New York, NY: Springer, 1st ed.
- (2002), "The Screening Effect in Kriging," The Annals of Statistics, 30, 298-323.
- (2007), "Spatial Variation of Total Column Ozone on a Global Scale," Annals of Applied Statistics, 1, 191-210.
- (2008), "A Modeling Approach for Large Spatial Datasets," Journal of the Korean Statistical Society, 37, 3-10.
- (2014), "Limitations on Low Rank Approximations for Covariance Matrices of Spatial Data," Spatial Statistics, 8, 1-19.

Stein, M. L., Chi, Z., and Welty, L. J. (2004), "Approximating Likelihoods for Large Spatial Data Sets," Journal of the Royal Statistical Society, Series B, 66, 275-296.

Stroud, J. R., Stein, M. L., and Lysen, S. (2014), "Bayesian and Maximum Likelihood Estimation for Gaussian Processes on an Incomplete Lattice," http://arxiv.org/abs/1402.4281.

Vecchia, A. V. (1988), "Estimation and Model Identification for Continuous Spatial Processes," Journal of the Royal Statistical Society, Series B, 50, 297-312.
- (1992), "A New Method of Prediction for Spatial Regression Models with Correlated Errors," Journal of the Royal Statistical Society, Series B, 54, 813-830.

Wang, Q., Adiku, S., Tenhunen, J., et al. (2005), "On the Relationship of NDVI with Leaf Area Index in a Deciduous Forest Site," Remote Sensing of Environment, 94, 244-255.

Wikle, C. and Cressie, N. A. C. (1999), "A Dimension-reduced Approach to Space-time Kalman Fltering," Biometrika, 86, 815-829.

Yeniay, O. and Goktas, A. (2002), "A Comparison of Partial Least Squares Regression with Other Prediction Methods," Hacettepe Journal of Mathematics and Statistics, 31, 99-111.

Zhang, X. and Kondraguanta, S. (2006), "Estimating Forest Biomass in the USA Using Generalized Allometric Models and MODIS Land Products," Geophysical Research Letters, 33, L09402.