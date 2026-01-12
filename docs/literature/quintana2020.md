\title{
The Dependent Dirichlet Process and Related Models
}

\author{
Fernand A. Quintana \({ }^{*, \dagger}\) Peter Müller \({ }^{\ddagger}\), Alejandro Jara \({ }^{*, \dagger}\) and Steven N. MacEachern \({ }^{\S}\) \\ Pontificia Universidad Católica de Chile * and Millennium Nucleus Center for the Discovery of Structures in Complex Data \({ }^{\dagger}\) and The University of Texas at Austin \({ }^{\ddagger}\) and The Ohio State University \({ }^{§}\)
}

\begin{abstract}
Standard regression approaches assume that some finite number of the response distribution characteristics, such as location and scale, change as a (parametric or nonparametric) function of predictors. However, it is not always appropriate to assume a location/scale representation, where the error distribution has unchanging shape over the predictor space. In fact, it often happens in applied research that the distribution of responses under study changes with predictors in ways that cannot be reasonably represented by a finite dimensional functional form. This can seriously affect the answers to the scientific questions of interest, and therefore more general approaches are indeed needed. This gives rise to the study of fully nonparametric regression models. We review some of the main Bayesian approaches that have been employed to define probability models where the complete response distribution may vary flexibly with predictors. We focus on developments based on modifications of the Dirichlet process, historically termed dependent Dirichlet processes, and some of the extensions that have been proposed to tackle this general problem using nonparametric approaches.
\end{abstract}

Key words and phrases: Related random probability distributions, Bayesian nonparametrics, Nonparametric regression, Quantile regression.

Departamento de Estadística,Santiago, Chile (e-mail: quintana@mat.uc.cl) (e-mail: atjara@uc.cl)
Department of Statistics and Data Science, Austin, Texas (e-mail: pmueller@math.utexas.edu)
Department of Statistics, Columbus, Ohio (e-mail: snm@stat.osu.edu)

\section*{1. INTRODUCTION}

We review the popular class of dependent Dirichlet process (DDP) models. These define a widely used fully nonparametric Bayesian regression for a response \(\boldsymbol{y} \in \mathscr{Y}\), based on a set of predictors \(\boldsymbol{x} \in \mathscr{X} \subseteq \mathbb{R}^{p}\). Despite a barrage of related literature over the past 25 years, to date there is no good review of such models. This paper fills this gap.

Fully nonparametric regression can be seen as an extension of traditional regression models, where, starting from some elements in \(\mathscr{X}\) and a corresponding set of responses in \(\mathscr{Y}\), the goal is to model the distribution of \(\boldsymbol{y}\) given \(\boldsymbol{x}\). Standard linear regression models proceed under the assumption of a Gaussian distribution for \(\boldsymbol{y} \mid \boldsymbol{x}\) with a mean modeled as a linear combination of \(\boldsymbol{x}\). Further extensions of this idea to exponential families gave rise to the popular class of generalized linear models, where a transformation of the mean response is modeled as a linear combination of \(\boldsymbol{x}\). Many other similar extensions are available. We focus on a nonparametric version of this idea, which involves going beyond the notion that the effect of predictors is restricted to change some particular functional of the response distribution, such as the mean, a quantile, or the parameters in a generalized linear model.

The fully nonparametric regression problem that we focus on arises when we assume that \(\boldsymbol{y}_{i} \mid F_{\boldsymbol{x}_{i}} \stackrel{i n d .}{\sim} F_{\boldsymbol{x}_{i}}, i=1, \ldots, n\). The parameter of interest is the complete set of predictor-dependent random probability measures \(\mathscr{F}=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\), where \(F_{\boldsymbol{x}}\) is a probability measure defined on the response sample space \(\mathscr{Y}\), whose elements can flexibly change with the values of the predictors \(\boldsymbol{x}\), i.e. the entire shape of the distribution can change with \(\boldsymbol{x}\). From a Bayesian point of view, the fully nonparametric regression model is completed by defining a prior distribution for \(\mathscr{F}\), which is taken to be the probability law of a probability measure-valued stochastic process with index \(\boldsymbol{x}\). At the risk of abusing notation, we use from now on the same symbol to refer to the probability measure and its cumulative distribution function (CDF). The distinction should be clear from the context.

Several popular approaches have been developed to formalize Bayesian inference for such nonparametric regression. These include additive random tree models like the BART (Chipman, George and McCulloch, 2010), approaches based on basis expansions such as wavelet regression and more. Also, there is of course extensive literature on non-Bayesian approaches to nonparametric regression. Many of these approaches are based on a model of the form
\[
\boldsymbol{y}_{i}=f\left(\boldsymbol{x}_{i}\right)+\epsilon_{i}, \quad i=1, \ldots, n,
\]
with \(E\left(\epsilon_{i}\right)=0\), and are concerned with finding a function \(f: \mathscr{X} \rightarrow \mathscr{Y}\) such that \(\left\|y_{i}-f\left(\boldsymbol{x}_{i}\right)\right\|\) is small, for \(f\) in some some class, often represented as being spanned by some basis functions. Such
methods include the following: Under local averaging \(f(\boldsymbol{x})\) is estimated from those \(\boldsymbol{y}_{i}\) 's such that \(\boldsymbol{x}_{i}\) is "close" to \(\boldsymbol{x}\); local modeling estimates \(f(\boldsymbol{x})\) by locally fitting some function or kernel such as a Gaussian function or a polynomial; global modeling or least squares estimation finds \(f\) that minimizes \(\frac{1}{n} \sum_{i=1}\left\|\boldsymbol{y}_{i}-f\left(\boldsymbol{x}_{i}\right)\right\|^{2}\) in the class; and penalized modeling is based on finding \(f\) that minimizes \(\frac{1}{n} \sum_{i=1}\left\|\boldsymbol{y}_{i}-f\left(\boldsymbol{x}_{i}\right)\right\|^{2}+J_{n}(\boldsymbol{y})\) in the class, where \(J_{n}(f)\) is a penalization term, such as \(J_{n}(f)=\lambda_{n} \int_{\mathscr{X}}\left|f^{\prime \prime}(t)\right|^{2} d t\). See, for example, Györfi et al. (2002); Klemelä (2014); Faraway (2016) and references within. Many of these classical frequentist approaches could be construed to imply nonparametric Bayesian models, but they are not usually cast as prior probability models for a family \(\mathscr{F}\) of random probability measures indexed by covariates.

In the Bayesian nonparametric (BNP) literature, the problem of defining priors over related random probability distributions has received increasing attention over the past few years. To date, most of the BNP priors to account for the dependence of a set of probability distributions on predictors are generalizations and extensions of the celebrated Dirichlet process (DP) (Ferguson, 1973, 1974) and Dirichlet process mixture (DPM) models (Lo, 1984). A DPM model defines a random probability measure as
\[
f(\boldsymbol{y} \mid G)=\int_{\Theta} \psi(\boldsymbol{y}, \boldsymbol{\theta}) G(d \boldsymbol{\theta}), \quad \boldsymbol{y} \in \mathscr{Y}
\]
where \(\psi(\bullet, \boldsymbol{\theta})\) is a continuous density function, for every \(\boldsymbol{\theta} \in \Theta\), and \(G\) is a discrete random probability measure with a DP prior. If \(G\) is DP with parameters \(\left(M, G_{0}\right)\), where \(M \in \mathbb{R}_{0}^{+}\)and \(G_{0}\) is a probability measure on \(\Theta\), written as \(G \mid M, G_{0} \sim \operatorname{DP}\left(M G_{0}\right)\), then the trajectories of the process can be a.s. represented by the stick-breaking representation (Sethuraman, 1994):
\[
G(B)=\sum_{h=1}^{\infty} w_{h} \delta_{\boldsymbol{\theta}_{h}}(B)
\]
where \(B\) is any measurable set, \(\delta_{\boldsymbol{\theta}}(\cdot)\) is the Dirac measure at \(\boldsymbol{\theta}, w_{h}=V_{h} \prod_{\ell<h}\left(1-V_{\ell}\right)\), with \(V_{h} \mid M \stackrel{i i d}{\sim} \operatorname{Be}(1, M), \boldsymbol{\theta}_{h} \mid G_{0} \stackrel{i i d}{\sim} G_{0}\), and the \(\left\{w_{h}\right\}\) and \(\left\{\boldsymbol{\theta}_{h}\right\}\) collections are independent. Discussion of properties and applications of DPs can be found, for instance, in Müller et al. (2015). Many BNP priors for nonparametric regressions \(\mathscr{F}=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) are based on extensions of model (1). They incorporate dependence on predictors via the mixing distribution in (1), by replacing \(G\) with \(G_{\boldsymbol{x}}\), and the prior specification problem is related to the modeling of the collection of predictor-dependent mixing probability measures \(\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\).

Consider first the simplest case, where a finite number of dependent RPMs \(\mathcal{G}=\left\{G_{j}, j=\right. 1, \ldots, J\}\) are judged to be exchangeable so that the prior model \(p(\mathcal{G})\) should accordingly be invariant with respect to all permutations of the indices. Consider, for example, an application to borrowing
strength across \(J\) related clincal studies. This can be achieved, for example, through joint modeling of study-specific effects distributions \(G_{j}\) for \(j=1, \ldots, J\). A main aim here is that subjects under study \(j_{1}\) should inform inference about subjects enrolled in a different but related study \(j_{2} \neq j_{1}\). Two extreme modeling choices would be (i) to pool all patients and assume one common effects distribution, or (ii) to assume \(J\) distinct distributions with independent priors. Formally, the earlier choice assumes \(G_{j} \equiv G, j=1, \ldots, J\), with a prior \(p(G)\), such as \(G \sim D P\left(M, G_{0}\right)\). The latter assumes \(G_{j} \sim p\left(G_{j}\right)\), independently, \(j=1, \ldots, J\). We refer to the two choices as extremes since the first choice implies maximum borrowing of strength, and the other choice implies no borrowing of strength. In most applications, the desired level of borrowing strength is somewhere in-between these two extremes.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0877b326-0bf3-4110-817e-01cbc5995ace-04.jpg?height=331&width=1306&top_left_y=898&top_left_x=403}
\captionsetup{labelformat=empty}
\caption{Fig 1. One common RPM \(G\) (panel a) versus distinct RPMs \(G_{j}\), independent across studies (panel b). Here \(\eta\) is a fixed hyperparameter.}
\end{figure}

Figure 1 illustrates the two modeling approaches. Note that in Figure 1 we added a hyperparameter \(\eta\) to index the prior model \(p\left(G_{j} \mid \eta\right)\) and \(p(G \mid \eta)\), which was implicitly assumed fixed. The use of a random hyperparameter \(\eta\) allows for some borrowing of strength even in the case of conditionally independent \(p\left(G_{j} \mid \eta\right)\). Learning across studies can happen through learning about the hyperparameter \(\eta\). However, the nature of the learning across studies is determined by the parametric form of \(\eta\). This is illustrated in Figure 2. Assume \(G_{j} \sim \mathrm{DP}\left(M, G_{\eta}^{\star}\right)\), independently, \(j=1,2\), and a base measure \(G_{\eta}^{\star}=\mathrm{N}(m, B)\) with unknown hyperparameter \(\eta=(m, B)\). In this case, prediction for a future study \(G_{3}\) can not possibly learn about the multimodality of \(G_{1}\) and \(G_{2}\), beyond general location and orientation.

The previous simple example illustrates the need to develop classes of models with the ability to relate collections of nonparametric distributions in more complex fashions. When this collection is indexed by a set of predictors \(\boldsymbol{x} \in \mathscr{X}\), the nonparametric regression approach mentioned earlier arises, and the definition of a prior on this collection enables one to borrow information across the distributions for responses, \(F_{\boldsymbol{x}}\). For modeling, one important property is the notion of distributions changing smoothly with respect to \(\boldsymbol{x} \in \mathscr{X}\), just as is the case of generalized linear models in the scale of the transformed mean. The smoothness could be expressed as continuity of \(F_{\boldsymbol{x}}\) (with respect

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0877b326-0bf3-4110-817e-01cbc5995ace-05.jpg?height=568&width=1531&top_left_y=236&top_left_x=317}
\captionsetup{labelformat=empty}
\caption{FIG 2. \(G_{j} \sim D P\left(M, G^{\star}\right)\) with common \(G^{\star}=N(m, B)\). Learning across studies is restricted to the parametric form of \(\eta\). The obvious common structure of \(G_{1}\) and \(G_{2}\) as defining three well separated clusters can not be learned by the model, which is restricted to learning through the common hyperparameters \(\eta\).}
\end{figure}
to some conveniently chosen topology) or as the notion that \(F_{\boldsymbol{x}}\) "approaches" \(F_{\boldsymbol{x}_{0}}\) as \(\boldsymbol{x} \rightarrow \boldsymbol{x}_{0}\), for instance, \(\operatorname{Corr}\left\{F_{\boldsymbol{x}}(A), F_{\boldsymbol{x}_{0}}(A)\right\} \rightarrow 1\) as \(\boldsymbol{x} \rightarrow \boldsymbol{x}_{0}\) for any event \(A\). Many of the models to be discussed later satisfy some version of this property.

An early reference on predictor-dependent DP models is Cifarelli and Regazzini (1978), who defined a model for related probability measures by introducing a regression model in the centering measure of a collection of independent DP random measures. This approach is used, for example, by Muliere and Petrone (1993), who considered a linear regression model for the centering distribution of the form \(G_{\boldsymbol{x}}^{0} \equiv N\left(\boldsymbol{x}^{\prime} \beta, \sigma^{2}\right)\), where \(\beta \in \mathbb{R}^{p}\) is a vector of regression coefficients, and \(N\left(\mu, \sigma^{2}\right)\) stands for a normal distribution with mean \(\mu\) and variance \(\sigma^{2}\). This is the type of construction illustrated in Figure 2. Similar models were discussed by Mira and Petrone (1996) and Giudici, Mezzetti and Muliere (2003). Linking the related nonparametric models through a regression on the baseline parameters of nonparametric models, however, limits the nature of the trajectories and the type of dependent processes that can be thus generated. Indeed, realizations of the resulting process \(\mathcal{G}=\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) are not continuous as a function of the predictors. The very limited type of association structure motivated the development of alternative extensions of the DP model to a prior for \(\mathcal{G}\). In this paper, we provide an overview of the main constructions of such predictordependent extensions of DP priors and and their main properties. The discussion centers on different ways of constructing the nonparametric component of models. A few of the many successful types of applications that have been proposed are mentioned. In reviewing the various models to be presented, we discuss some of the main corresponding works without attempting to provide a complete catalog of references. We include a brief discussion of other popular constructions of
dependent DP random measures, without the explicit notion of a conditioning covariate \(x\).
While we focus on DP-based constructions, we note that several interesting alternatives to develop predictor-driven random probability measures have been considered in the recent literature. Tokdar, Zhu and Ghosh (2010) develop a logistic Gaussian process that allows for smoothly varying dependence on conditioning variables. Still using Gaussian process priors, but starting from a rather different construction, Jara and Hanson (2011) proposed another alternative, putting the Gaussian process prior on the (logit transformation) of the branch probabilities in a Polya tree prior (Lavine, 1992). Another covariate-dependent extension of the Polya tree model was introduced in Trippa, Müller and Johnson (2011) who define a dependent multivariate process for the branch probabilities based on a simple gamma process construction.

Finally, although issues pertaining to the implementation of posterior simulation are relevant for practical application of these methods, our discussion does not focus on computational aspects.

In Section 2 we describe MacEachern's dependent Dirichlet process (DDP) and its basic properties. In Section 3 we discuss the main variations and alternative constructions to MacEachern's DDP. In Section 4 we discuss approaches to handle endogenous predictors. In Section 5 we discuss the implied partition structure of DDP models. In Section 6 we illustrate the main approaches. A final discussion in Section 7 concludes the article, including some thoughts on future research directions.

\section*{2. DEPENDENT DIRICHLET PROCESS (DDP)}

We start our discussion with the general definition of DDP and then give details for popular special cases.

\subsection*{2.1 General definition}

MacEachern \((1999,2000)\) introduced the DDP model as a flexible class of predictor-dependent random probability distributions. The key idea behind the DDP construction is to define of a set of random measures that are marginally (i.e. for every possible predictor value \(\boldsymbol{x} \in \mathscr{X}\) ) DPdistributed random measures. In this framework, dependence is introduced through a modification of the stick-breaking representation of each element in the set,
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty} \underbrace{\left\{V_{h}(\boldsymbol{x}) \prod_{\ell<h}\left[1-V_{\ell}(\boldsymbol{x})\right]\right\}}_{w_{h}(\boldsymbol{x})} \delta_{\boldsymbol{\theta}_{h}(\boldsymbol{x})}(\bullet)
\]
where \(V_{h}(\boldsymbol{x}), h \in \mathbb{N}\), are [ 0,1 ]-valued independent stochastic processes with index set \(\mathscr{X}\) and \(\operatorname{Be}\left(1, M_{\boldsymbol{x}}\right)\) marginal distributions, and \(\boldsymbol{\theta}_{h}(\boldsymbol{x}), h \in \mathbb{N}\), are independent stochastic processes with
index set \(\mathscr{X}\) and \(G_{\boldsymbol{x}}^{0}\) marginal distributions. The processes associated to the weights and atoms are independent. From an intuitive viewpoint, the constructed DDP can be thought of as taking an ordinary DP and modifying some of its components (i.e. weights and atoms) according to the type of desired indexing or functional dependence of predictors \(\boldsymbol{x} \in \mathscr{X}\). Conditions on the \(V_{h}(\boldsymbol{x})\) and \(\boldsymbol{\theta}_{h}(\boldsymbol{x})\) processes can be established to ensure smoothness of the resulting random measures \(G_{\boldsymbol{x}}(\bullet)\) when \(\boldsymbol{x}\) ranges over \(\mathscr{X}\).

Canonical DDP construction. MacEachern \((1999,2000)\) defined and provided a canonical construction of the DDP by using transformations of two independent sets of stochastic processes, \(Z_{\mathscr{X}}^{V_{h}}=\left\{Z_{h}^{V}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\), and \(Z_{\mathscr{X}}^{\boldsymbol{\theta}_{h}}=\left\{Z_{h}^{\boldsymbol{\theta}}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\), for \(h \geq 1\), the former used for defining \(\left\{V_{h}(\boldsymbol{x})\right\}\), and the latter for defining \(\left\{\boldsymbol{\theta}_{h}(\boldsymbol{x})\right\}\). To induce the desired marginal distributions for \(\left\{V_{h}(\boldsymbol{x})\right\}\) and \(\left\{\boldsymbol{\theta}_{h}(\boldsymbol{x})\right\}\), MacEachern resorted to the well-known inverse transformation method (see, e.g. Devroye, 1986). For instance, let \(Z(x)\) denote a zero-mean Gaussian process on \(\mathscr{X}=\mathbb{R}\) having constant variance \(\sigma^{2}\). Let \(\Phi(\cdot)\) and \(B(\cdot)\) denote the cumulative distribution functions of the \(N(0,1)\) and \(\operatorname{Be}(1, M)\) distributions, respectively. Then \(V(x)=B^{-1}\left(\Phi\left(\sigma^{-1} Z(x)\right)\right)\) is a stochastic process on \(\mathscr{X}\) that satisfies \(V(x) \sim \operatorname{Be}(1, M)\) for all \(x \in \mathscr{X}\). The same type of transformation can be applied to construct suitable atom processes \(\left\{\theta_{h}(x), h \geq 1\right\}\) such that \(\theta_{h}(x) \sim G_{0}\) for all \(x \in \mathscr{X}\) and \(h \geq 1\).

Practical application of this general model requires specification of its various components, which has traditionally motivated the adoption of some specific forms. The most commonly used DDPs assume that covariate dependence is introduced either in the atoms or weights, leaving the other as a collection of random variables exhibiting no covariate indexing, so that the basic DP definition is partially modified but the distributional properties retained. We review these forms in the next section.

Support and an alternative definition. One particularity of MacEachern's DDP definition is that given the sets of stochastic processes, \(Z_{\mathscr{X}}^{V_{h}}=\left\{Z_{h}^{V}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\) and \(Z_{\mathscr{X}}^{\boldsymbol{\theta}_{h}}=\left\{Z_{h}^{\boldsymbol{\theta}}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\), and all other parameters involved in the transformations described above, the collection of dependent probability distributions given in (3) are not random: they are just deterministic functions of these quantities. To facilitate the study of theoretical properties of the DDP, Barrientos, Jara and Quintana (2012) gave an alternative definition. This alternative definition exploits the connection between copulas and stochastic processes. Since under certain regularity conditions a stochastic process is completely characterized by its finite-dimensional distributions, it is possible -and usefulto define stochastic processes with given marginal distributions via copulas. The basic idea is to specify the collection of finite dimensional distributions of a process through a collection of copulas and marginal distributions.

Copulas are functions that are useful for describing and understanding the dependence structure between random variables. If \(H\) is a \(d\)-variate CDF with marginal CDF's given by \(F_{1}, \ldots, F_{d}\), then by Sklar's theorem (Sklar, 1959), there exists a copula function \(C:[0,1]^{d} \longrightarrow[0,1]\) such that \(H\left(t_{1}, \ldots, t_{d}\right)=C\left(F_{1}\left(t_{1}\right), \ldots, F_{d}\left(t_{d}\right)\right)\), for all \(t_{1}, \ldots, t_{d} \in \mathbb{R}\), and this representation is unique if the marginal distributions are absolutely continuous. Thus by the probability integral transform, a copula function is a \(d\)-variate CDF on \([0,1]^{d}\) with uniform marginals on \([0,1]\), which fully captures the dependence among the associated random variables, irrespective of the marginal distributions.

Let \(\mathcal{C}_{\mathcal{X}}^{V}=\left\{C_{\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{d}}^{V}: \boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{d} \in \mathcal{X}, d>1\right\}\) and \(\mathcal{C}_{\mathscr{X}}^{\theta}=\left\{C_{\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{d}}^{\theta}: \boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{d} \in \mathcal{X}, d>1\right\}\) be two sets of copulas satisfying Kolmogorov's consistency conditions. In Barrientos, Jara and Quintana (2012)'s definition, \(V_{h}(\boldsymbol{x}), h \in \mathbb{N}\), are [ 0,1 ]-valued independent stochastic processes with index set \(\mathscr{X}\), with common finite dimensional distributions determined by the set of copulas \(\mathcal{C}_{\mathscr{X}}^{V}\), and \(\operatorname{Be}\left(1, M_{\boldsymbol{x}}\right)\) marginal distributions. Similarly, \(\boldsymbol{\theta}_{h}(\boldsymbol{x}), h \in \mathbb{N}\), are independent stochastic processes with index set \(\mathscr{X}\), with common finite dimensional distributions determined by the set of copulas \(\mathcal{C}_{\mathscr{X}}^{\theta}\), and \(G_{\boldsymbol{x}}^{0}\) marginal distributions. This alternative construction produces a definition of the DDP exactly as in (3), and in particular, the interpretation of the DDP obtained as modifying a basic DP persists. Furthermore, based on this alternative definition, Barrientos, Jara and Quintana (2012) established basic properties of MacEachern's DDP and other dependent-stick breaking processes. Specifically, they provided sufficient conditions for the full weak support of different versions of the process and also to ensure smoothness of trajectories of \(G_{\boldsymbol{x}}(\bullet)\) as \(\boldsymbol{x}\) ranges over \(\mathscr{X}\). In addition, they also characterized the Hellinger and Kullback-Leibler support of mixtures induced by different versions of the DDP and extended the results to the general class of dependent stick-breaking processes.

\subsection*{2.2 The single-weights DDP}

MacEachern considered the case of common weights across the values of \(\boldsymbol{x}\), also referred to as "single-weights" DDP model, defined as
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty} \underbrace{\left\{V_{h} \prod_{\ell<h}\left[1-V_{\ell}\right]\right\}}_{w_{h}} \delta_{\boldsymbol{\theta}_{h}(\boldsymbol{x})}(\bullet)=\sum_{h=1}^{\infty} w_{h} \delta_{\boldsymbol{\theta}_{h}(\boldsymbol{x})}(\bullet),
\]
where the \(V_{h}\) 's are iid \(\operatorname{Be}(1, M)\) random variables, which are common across all levels of \(\boldsymbol{x}\). The \(\boldsymbol{\theta}_{h}(\boldsymbol{x})\) 's are independent stochastic processes with index set \(\mathcal{X}\) and marginal distributions \(G_{\boldsymbol{x}}^{0}\). In the literature, to this day, this is the most popular form of DDP, mainly due to the fact that posterior simulation can be implemented using the same type of sampling algorithms available for the case of the DP.
2.2.1 The ANOVA-DDP and linear DDP models One of the earliest versions of DDP models was the ANOVA-DDP of De Iorio et al. (2004). Let \(\boldsymbol{y}=\left(y_{1}, \ldots, y_{n}\right)\) be a vector of responses (possibly vector-valued) for each of \(n\) subjects, and suppose that \(\boldsymbol{x}=\left(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{n}\right)\) is a corresponding set of covariates. Assume each \(\boldsymbol{x}_{i}\) is in turn a vector of \(c\) categorical covariates, \(\boldsymbol{x}_{i}=\left(x_{i 1}, \ldots, x_{i c}\right)\). Interpret \(\boldsymbol{x}_{i}\) as factors in an ANOVA model, and let \(d_{i}\) denote corresponding design vectors. Assume then that \(\boldsymbol{x}_{i}\) contains all the desired main effects and interactions, as well as desired identifiability constraints. Note that the covariate space \(\mathscr{X}\) in this setup is discrete, and so we have a finite number of RPMs. The idea of the ANOVA-DDP models is to encode the covariate dependence in the form of simple linear regressions for the atom processes \(\left\{\theta_{h}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\). Specifically, this approach uses \(\theta_{h}(\boldsymbol{x})=\lambda_{h}^{\prime} d_{\boldsymbol{x}}\) for \(h \geq 1\) where \(\left\{\lambda_{h}: h \geq 1\right\}\) is a sequence of iid random vectors with distribution \(G_{0}\) and \(d_{\boldsymbol{x}}\) is the design vector that corresponds to a generic combination of observed factorial covariates \(\boldsymbol{x}\). The model just described implies that (4) becomes
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty} w_{h} \delta_{\lambda_{h}^{\prime} d_{\boldsymbol{x}}}(\bullet),
\]
i.e., a DP mixture of linear models \(\lambda_{h}^{\prime} d_{\boldsymbol{x}}\). Each element of the collection \(\mathcal{G}=\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) has a DP prior distribution with atoms given by \(\left\{\lambda_{h}^{\prime} d_{\boldsymbol{x}}: h \geq 1\right\}\). The elements of \(\mathcal{G}\) are correlated because they share a common set of weights and the atoms are originated as linear combinations computed from a single set of parameters, namely \(\left\{w_{h}: h \geq 1\right\}\) and \(\left\{\lambda_{h}: h \geq 1\right\}\).

To accommodate a continous response, De Iorio et al. (2004) extended the above construction through a convolution with a continuous kernel, e.g., a normal kernel, leading to
\[
y_{i} \mid G_{\boldsymbol{x}_{i}} \stackrel{i n d}{\sim} \int N\left(y_{i} \mid m \mu, \phi\right) d G_{\boldsymbol{x}_{i}}(m \mu)=\int N\left(y_{i} \mid \lambda^{\prime} d_{\boldsymbol{x}_{i}}, \phi\right) d G(\lambda)
\]

The model can be restated by breaking the mixture with the introduction of latent parameters:
\[
y_{i}\left|\lambda_{i}, \phi \sim N\left(\lambda_{i}^{\prime} d_{i}, \phi\right), \quad \lambda_{1}, \ldots, \lambda_{n}\right| G \stackrel{\mathrm{iid}}{\sim} G, \quad G \sim D P\left(M, G_{0}\right) .
\]

The last expression highlights the nature of the model as just a DP mixture of, in this case, normal linear models. The same simplification is possible whenever the atoms \(\left\{\theta_{h}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\) are indexed by a finite-dimensional parameter vector, like the linear model \(\theta_{h}(\boldsymbol{x})=\lambda_{h}^{\prime} d_{\boldsymbol{x}}\) in this case. The model in (5) is completed with a suitable prior for the precision parameter \(\phi\), e.g. \(\phi \sim G a(a, b)\) if a scalar, or \(\phi \sim \operatorname{Wishart}(\nu, S)\) if a matrix. The above model can be easily modified to mix over scale parameters as well. An immediate consequence of (5) is that the induced marginal distribution for a single response \(y\) with design vector \(d_{\boldsymbol{x}}\) then becomes a flexible infinite mixture model:
\[
y \sim \sum_{h=1}^{\infty} w_{h} N\left(y \mid \lambda_{h}^{\prime} d_{\boldsymbol{x}}, \phi\right)
\]

We remark here that the hierarchical structure leading to (5) reflects a common practice in the use and application of the DDP. Since marginally each element of the \(\mathcal{G}\) family is almost surely discrete (because it is drawn from a DP), models for discrete outcomes are frequently built on convolving the DPs with a continuous kernel, thus yielding a mixture of continuous distributions, which is itself a continuous distribution. In the ANOVA-DDP model of De Iorio et al. (2004), the normal kernel plays precisely this role.

De la Cruz-Mesía, Quintana and Müller (2007) applied the ANOVA-DDP construction to model random effects for longitudinal hormone profiles of pregnant women, where the dependence was on a normal/abnormal pregnancy indicator. This setting was particularly useful for classification purposes. More recently, Gutiérrez et al. (2019) use the ANOVA-DDP framework to propose a multiple testing procedure for comparing several treatments against a control. A further extension of the ANOVA-DDP construction was given in De Iorio et al. (2009), who considered the modeling of nonproportional hazards for survival analysis. They considered a cancer clinical trial, where interest centered on whether high doses of a treatment are more effective than lower doses. The data included additional discrete and continuous covariates, so the model was under the extended ANCOVA-style framework that adds linear combinations of continuous covariates to the ANOVA factorial design.

This same idea can be extended to linear combinations of any given set of covariates, giving rise to the linear DDP (LDDP). Specifically, such models involve a linear combination of a set of covariates, as in, e.g. general linear models, and so the infinite mixture on the right-hand side of (6) becomes \(\sum_{h=1}^{\infty} w_{h} N\left(y \mid \lambda_{h}^{\prime} \boldsymbol{x}, \phi\right)\), where \(\boldsymbol{x}\) is now the generic value of the (typically vector-valued) covariate. As earlier, the weights \(\left\{w_{h}\right\}\) follow a DP-style stick-breaking specification. An analogous expression for a more general kernel function \(k\) can be immediately derived. The same type of construction was explored in Jara et al. (2010) in the context of doubly censored outcomes. Their model involves an interval-valued response, corresponding to the observed onset and event times (cavities in the teeth of children from Flanders, Belgium, in their example). Associated with each such response is a latent bivariate vector of true onset and event times, and these are modeled (in the logarithmic scale) using a linear DDP defined in terms of covariates that include deciduous second molars health status and the age at which children started brushing.
2.2.2 Spatial DDP Gelfand, Kottas and MacEachern (2005) define what can be interpreted as a spatial case of a common weight DDP (4) for \(G_{s}\), with \(s \in D \subset \mathbb{R}^{d}\) being spatial locations and \(\theta_{h}(s)\) generated by a baseline GP, as in the common-weight DDP. However, the focus is not on \(G_{s}\) as in (4), but instead on \(\boldsymbol{\theta}_{D} \sim \sum w_{h} \delta_{\boldsymbol{\theta}_{h, D}}\), where \(\boldsymbol{\theta}_{h, D}=\left\{\theta_{h}(s), s \in D\right\}\). Let \(\boldsymbol{s}=\left(s_{1}, \ldots, s_{n}\right)\)
denote a set of \(n\) locations at which observations \(\boldsymbol{y}=\left(y_{1}, \ldots, y_{n}\right)\) are made. They consider repeat observations \(\boldsymbol{y}_{t}, t=1, \ldots, T\), with occasion-specific covariates \(\boldsymbol{x}_{t}\). Writing a mixture with respect to a DP random measure as a hierarchical model, they assume
\[
\boldsymbol{y}_{t}\left|\boldsymbol{\theta}_{t}, \boldsymbol{\beta}, \tau^{2} \stackrel{\mathrm{ind}}{\sim} N\left(\boldsymbol{x}_{t}^{\prime} \boldsymbol{\beta}+\boldsymbol{\theta}_{t}, \tau^{2} \boldsymbol{I}\right), \quad \boldsymbol{\theta}_{t}\right| G^{\eta} \stackrel{\mathrm{iid}}{\sim} G^{\eta}, \quad G^{\eta} \sim D P\left(M, G_{0}^{\eta}\right),
\]
where \(G_{0}^{\eta} \equiv N\left(\mathbf{0}, \sigma^{2} \boldsymbol{H}(\eta)\right)\) and \(\boldsymbol{H}(\eta)\) is a suitable covariance function depending on hyperparameters \(\eta\).

Dunson and Herring (2006) considered a model for a collection of random functions based on a finite set of latent trajectories described by Gaussian processes. The observations are thus seen as arising from the convolution of a smooth latent trajectory and a noisy Gaussian process. Their motivation came from the study of the relationship between disinfection by-products in the water in early pregnancy and later outcomes. Specifically, denoting by \(g_{i}\) the stochastic process, i.e. \(\left\{g_{i}(t): t>0\right\}\), associated with subject \(1 \leq i \leq n\), Dunson and Herring (2006) assume that
\[
g_{i}=\gamma_{i}+\epsilon_{i}, \quad \gamma_{i} \stackrel{\mathrm{iid}}{\sim} G, \quad \epsilon_{i} \stackrel{\mathrm{iid}}{\sim} G P(\boldsymbol{H}(\eta)),
\]
where \(\gamma_{i}\) is the latent trajectory, and \(G P(\boldsymbol{H}(\eta))\) denotes a Gaussian process with covariance function \(\boldsymbol{H}(\eta)\). Their approach specifies the RPM \(G\) as \(G(\cdot)=\sum_{h=1}^{k} p_{h} \delta_{\Theta_{h}}(\cdot)\) with \(\Theta_{h} \sim G P\left(\boldsymbol{H}\left(\eta_{\kappa_{h}}\right)\right)\), i.e., a finite mixture of atoms given by Gaussian processes with suitable covariance functions. By choosing \(\kappa_{h}=\kappa\) for all \(h\) and \(\left(p_{1}, \ldots, p_{k}\right) \sim \operatorname{Dir}(M / k, \ldots, M / k)\), the resulting RPM \(G\) approaches \(G(\cdot)=\sum_{h=1}^{\infty} w_{h} \delta_{\Theta_{h}}(\cdot)\) as \(k \rightarrow \infty\) with DP-style weights (see, e.g. Green and Richardson, 2001).
2.2.3 Dynamic DDP The DDP framework has also been used to model dynamic phenomena, by means of a sequence of random distributions that evolve in time. Caron et al. (2008) considered a dynamic linear model formulation to solve this problem, where the state and observation noise distributions where modeled as DP mixtures using two independent DPs so that the mean of the underlying processes is allowed to change in time.

Rodríguez and ter Horst (2008) considered a related model, based on a DDP formulation, where now the atoms in the infinite mixture are allowed to change in time. Letting \(y_{i t}\) denote the \(i\) th observation at time \(1 \leq t \leq T\), they proposed the model
\[
y_{i t} \mid G_{t} \sim \int N\left(\boldsymbol{F}_{i t}^{\prime} \boldsymbol{\theta}_{t}, \sigma^{2}\right) d G_{t}\left(\boldsymbol{\theta}_{t}, \sigma^{2}\right), \quad G_{t}(\cdot)=\sum_{h=1}^{\infty} w_{h} \delta_{\left(\boldsymbol{\theta}_{h t}^{*}, \sigma_{h}^{* 2}\right)}(\cdot), \quad \boldsymbol{\theta}_{h t}^{*} \sim N\left(\boldsymbol{H}_{t} \boldsymbol{\theta}_{h, t-1}^{*}, \sigma_{h}^{* 2} \boldsymbol{W}_{t}\right)
\]
completed with conjugate priors for \(\sigma_{h}^{* 2}\) and \(\boldsymbol{\theta}_{h, 0}^{*}\). Matrices \(\boldsymbol{F}_{i t}, \boldsymbol{H}_{t}\) and \(\boldsymbol{W}_{t}\) are assumed known and can be used to represent many patterns such as trends, periodicity, etc. The resulting model
for \(\mathcal{G}=\left\{G_{t}: 1 \leq t \leq T\right\}\) is thus a DDP, where the components of the atoms controlling the distribution means evolve in time in an autoregressive fashion.

Di Lucca et al. (2012) considered a model for a sequence of random variables \(\left\{y_{t}: t \geq 1\right\}\) featuring a general autoregressive formulation by means of \(y_{t} \mid\left(y_{t-1}, \ldots, y_{t-p}\right)=\boldsymbol{y} \sim G_{\boldsymbol{y}}\) and the problem of defining a prior for \(\mathcal{G}=\left\{G_{\boldsymbol{y}}: \boldsymbol{y} \in \mathscr{Y}\right\}\). They discussed a general prior DDP model of the form \(G_{\boldsymbol{y}}(\cdot)=\sum_{h=1}^{\infty} w_{h}(\boldsymbol{y}) \delta_{\boldsymbol{y}}(\cdot)\). Lau and So (2008) considered similar types of model, where each atom can be expressed as an infinite mixture of autoregressions of order \(p\). Di Lucca et al. (2012) focused on the particular single-weights case and an order \(p=1\) process where the atom processes are expressed as simple linear autoregression: \(\theta_{h}(\boldsymbol{y})=\beta_{h}+\alpha_{h} y\). The full model in this case can be expressed as
\[
y_{t}\left|y_{t-1}=y, \alpha_{t}, \beta_{t}, \sigma^{2} \sim N\left(\beta_{t}+\alpha_{t} y, \sigma^{2}\right), \quad\left(\beta_{t}, \alpha_{t}\right)\right| G \stackrel{\mathrm{iid}}{\sim} G, \quad G \sim D P\left(M, G_{0}\right)
\]

However, they also considered the case when atoms are defined as \(\theta_{h}(y)=b+a_{h} y+O U\left(\rho, \tau^{2}\right)\), where \(O U\left(\rho, \tau^{2}\right)\) denotes the Ornstein-Uhlenbeck process, a particular Gaussian process with covariance function of the form \(\operatorname{Cov}[\theta(s), \theta(t)]=\tau^{2} \rho^{|s-t|}\). Di Lucca et al. (2012) extended this approach for sequences of binary outcomes defined in terms of an autoregressive process \(Z_{t}\) with a flexible DDP prior distribution, where dependence is on the previous \(p\) binary responses.

An interesting variation of a dynamic DDP construction is proposed by Ascolani, Lijoi and Ruggiero (2020) who define a family \(\mathcal{G}=\left\{G_{t}, t \geq 0\right\}\) of dependent random probability measures indexed by time. Their construction is motivated by a Fleming-Viot process. The random probability measures \(G_{t}\) share some, but not all atoms. The set \(D_{t}\) of atoms in the original \(G_{0}\) which are shared in \(G_{t}\) is defined as a pure death process over time. Importantly, each \(G_{t}\) marginally remains a DP random measure. They refer to the model as the Fleming-Viot-DDP. In Prünster and Ruggiero (2013) this construction is applied to model market shares over time. Mena and Ruggiero (2016) construct another common-atoms DDP over time by setting up a Wrights-Fisher diffusion on the fractions \(v_{t, \ell}\) in the stick-breaking construction of the marginal DP prior for \(G_{t}\).

\subsection*{2.3 The single-atoms DDP}

A parallel construction to the common weights DDP in the previous section considers a set of common atoms across all values of \(\boldsymbol{x}\). This is the so called "single-atoms" DDP model, for which (3) takes the form
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty} \underbrace{\left\{V_{h}(\boldsymbol{x}) \prod_{\ell<h}\left[1-V_{\ell}(\boldsymbol{x})\right]\right\}}_{w_{h}(\boldsymbol{x})} \delta_{\boldsymbol{\theta}_{h}}(\bullet)
\]
where \(V_{h}(\boldsymbol{x}), h \in \mathbb{N}\), are \([0,1]\)-valued independent stochastic processes with index set \(\mathscr{X}\) and marginal distributions \(\operatorname{Be}\left(1, M_{\boldsymbol{x}}\right)\). The locations \(\boldsymbol{\theta}_{h}, h \in \mathbb{N}\), are independent with marginal distributions \(G^{0}\); and the \(\left\{V_{h}(\boldsymbol{x})\right\}\) and \(\left\{\boldsymbol{\theta}_{h}\right\}\) collections are mutually independent.

Under the single-atoms model, all the covariate-dependence is expressed through the weights of the stick-breaking representation. One advantage of doing so is that, unlike the single-weights case, the implied prior probability model on partitions changes with the values of \(x \in \mathscr{X}\). This is important when the implied partition is of interest. Another important feature is that problems related to extrapolation of \(\boldsymbol{\theta}_{h}(\boldsymbol{x})\) are avoided, which could otherwise arise for inference for a new value of \(\boldsymbol{x}\) beyond the range of the observed data. This is the case because under the single-atoms DDP all atoms are linked with observed data, in contrast to the single-weights DDP which includes atoms for new covariate values that are not linked with any observed data.

Duan, Guindani and Gelfand (2007) describe a model motivated by the analysis of spatially varying responses. Let \(\{y(s): s \in D\}\) be a stochastic process indexed by locations in a set \(D \subset \mathbb{R}^{d}\), and let \(s_{1}, \ldots, s_{n}\) the locations at which observations are collected. Their general construction involves a RPM \(G\) over the space of surfaces of \(D\) having finite-dimensionals adopting the following form: for any \(s_{1}, \ldots, s_{n} \in D\) and \(A_{1}, \ldots, A_{n}\) Borel-measurable sets in \(\mathbb{R}\),
\[
P\left(y\left(s_{1}\right) \in A_{1}, \ldots, y\left(s_{n}\right) \in A_{n}\right)=\sum_{i_{1}=1}^{\infty} \cdots \sum_{i_{n}=1}^{\infty} p_{i\left(s_{1}\right), \ldots, i\left(s_{n}\right)} \delta_{\theta_{i\left(s_{1}\right)}}\left(s_{1}\right) \cdots \delta_{\theta_{i\left(s_{n}\right)}}\left(s_{n}\right)
\]
where the \(\theta_{j}\) 's are iid from \(G_{0}\) and the weights \(\left\{p_{i\left(s_{1}\right), \ldots, i\left(s_{n}\right)}\right\}\) determine the site-specific joint selection probabilities. Conditions can be given so that the above specification follows a DP at any given location.

Always in the spatial context, specifically of modeling for hurricane surface wind fields, Reich and Fuentes (2007) propose a general framework that includes the single-atoms DDP as a special case. Their model is specially designed for spatial dependence as well, so that the covariates are geographical coordinates. Letting \(s\) denote such coordinates, their construction involves weights computed as \(w_{1}(\boldsymbol{s})=V_{1}(\boldsymbol{s})\) and \(w_{h}(\boldsymbol{s})=V_{h}(\boldsymbol{s}) \prod_{\ell=1}^{h-1}\left(1-V_{\ell}(\boldsymbol{s})\right)\) for \(h>1\), where \(V_{h}(\boldsymbol{s})=\omega_{h}(\boldsymbol{s}) V_{h}\), and \(V_{h} \stackrel{\text { iid }}{\sim} \operatorname{Beta}(a, b)\). The function \(\omega_{h}(\boldsymbol{s})\) is centered at knot \(\boldsymbol{\psi}_{h}=\left(\psi_{h 1}, \psi_{h 2}\right)\), and the spread is controlled by parameters \(\boldsymbol{e}_{h}=\left(e_{h 1}, e_{h 2}\right)\). Reich and Fuentes (2007) discuss several possible choices for the \(\omega_{h}\) functions and related parameters.

Griffin and Steel (2006) define another interesting variation of the basic DDP by keeping both sets of parameters, locations and the fractions \(\left(V_{h}\right)\), unchanged across \(\boldsymbol{x}\). They use instead permutations of how the weights are matched with locations. The permutations change with \(\boldsymbol{x}\). One advantage of such models is the fact that the support of \(G_{\boldsymbol{x}}\) remains constant over \(\boldsymbol{x}\), a feature that can be
important for extrapolation beyond the observed data. A modification of this idea was explored by Griffin and Steel (2010) to generate what they called the DP regression smoother. The construction is centered over a class of regression models, and dependence is on the weights. More recently, similar ideas are used by Griffin and Steel (2011) to construct a family of prior distributions for a sequence of time dependent general RPMs that include the DDP setting as a special case. Another simple sequence of time-dependent DDPs was proposed by Gutiérrez, Mena and Ruggiero (2016), with a Markov chain structure for the sequence of time-varying sticks, and with application to the analysis of air quality data.

\section*{3. VARIATIONS OF MACEACHERN'S DDP}

In this section we discuss a variety of models extending the original definition (3). Many of these extensions are based on constructing independent weights and atoms processes indexed by covariates, but that do not necessarily produce a DP-distributed random measure. From an intuitive viewpoint, these classes of models can be seen as taking the basic DP construction and altering some of their basic components in terms of predictors \(\boldsymbol{x} \in \mathscr{X}\) to a form that may differ from the initial distributional properties. While this typically modifies the marginal DP property, the extra flexibility allows one to tailor the properties of the model to fit specific applications.

\subsection*{3.1 Weighted mixture of DPs (WMDP)}

Dunson, Pillai and Park (2007) proposed a data-based prior using the observed predictors \(\boldsymbol{x}_{1}, \ldots, \boldsymbol{x}_{n}\). For every \(\boldsymbol{x} \in \mathscr{X} \subset \mathbb{R}^{p}\), they considered the following construction
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{j=1}^{n}\left(\frac{\gamma_{j} K\left(\boldsymbol{x}, \boldsymbol{x}_{j}\right)}{\sum_{\ell=1}^{n} \gamma_{\ell} K\left(\boldsymbol{x}, \boldsymbol{x}_{\ell}\right)}\right) G_{j}(\bullet),
\]
with
\[
\gamma_{j}\left|\kappa \stackrel{i i d}{\sim} \Gamma(\kappa, n \kappa), \quad G_{j}\right| M, G_{0} \stackrel{i i d}{\sim} D P\left(M, G_{0}\right),
\]
where \(K: \mathscr{X} \times \mathscr{X} \longrightarrow \mathbb{R}^{+}\)is a bounded kernel function. The choice of \(K\) impacts the degree of borrowing of information from the neighbors in estimating the distribution at any particular predictor value \(\boldsymbol{x}\). Some choices are discussed in the original technical report. In the paper, they considered
\[
K\left(\boldsymbol{x}, \boldsymbol{x}^{\prime}\right)=\exp \left\{\psi\left\|\boldsymbol{x}-\boldsymbol{x}^{\prime}\right\|^{2}\right\}, \quad \psi \mid \mu_{\psi}, \sigma_{\psi}^{2} \sim L N\left(\mu_{\psi}, \sigma_{\psi}^{2}\right),
\]
where \(L N(a, b)\) denotes the log-normal distribution with parameters \(a \in \mathscr{R}\) and \(b>0\). With this choice, the resulting model for a given \(\boldsymbol{x}\) borrows more heavily from those \(G_{j}\) 's for which
the corresponding \(\boldsymbol{x}_{j}\) is close to \(\boldsymbol{x}\). One primary application of this particular construction is in the context of density regression i.e. in measuring how a probability distribution on the space of responses \(\mathscr{Y}\) changes according to predictors \(\boldsymbol{x} \in \mathscr{X}\).

\subsection*{3.2 Kernel stick-breaking}

The kernel stick-breaking process (KSBP) was introduced by Dunson and Park (2008). For all \(\boldsymbol{x} \in \mathscr{X} \subset \mathbb{R}^{p}\), the KSBP is defined as follows
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty}\left\{W\left(\boldsymbol{x} ; \boldsymbol{\Gamma}_{h}, V_{h}\right) \prod_{\ell<h}\left(1-W\left(\boldsymbol{x} ; \Gamma_{\ell}, V_{\ell}\right)\right)\right\} G_{h}(\bullet),
\]
where \(W\left(\boldsymbol{x} ; \Gamma_{h}, V_{h}\right)=V_{h} K\left(\boldsymbol{x}, \boldsymbol{\Gamma}_{h}\right)\), with \(K: \mathscr{X} \times \mathscr{X} \longrightarrow[0,1]\), e.g. as given in (9), \(V_{h} \mid a_{h}, b_{h} \stackrel{\text { ind }}{\sim}\). \(\operatorname{Be}\left(a_{h}, b_{h}\right), \boldsymbol{\Gamma}_{h} \mid H \stackrel{i i d}{\sim} H\) (random kernel locations), and \(G_{h} \mid \mathcal{G} \stackrel{i i d}{\sim} \mathcal{G}\) (random probability measures). The KSBP thus begins with an infinite sequence of basis random distributions \(\left\{G_{h}\right\}\) and then constructs covariate-dependent random measures by mixing according to distance from the random locations \(\Gamma_{h}\), with stick-breaking probabilities that are defined as a kernel multiplied by Betadistributed weights. It is also possible to simplify the definition of KSBP, adopting the particular form
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty}\left\{W\left(\boldsymbol{x} ; \boldsymbol{\Gamma}_{h}, V_{h}\right) \prod_{\ell<h}\left(\left(1-W\left(\boldsymbol{x} ; \Gamma_{\ell}, V_{\ell}\right)\right)\right\} \delta_{\boldsymbol{\theta}_{h}}(\bullet)\right.
\]
where \(W\left(\boldsymbol{x} ; \Gamma_{h}, V_{h}\right)=V_{h} K\left(\boldsymbol{x}, \boldsymbol{\Gamma}_{h}\right)\), with \(K: \mathscr{X} \times \mathscr{X} \longrightarrow[0,1], V_{h}\left|M \stackrel{i i d}{\sim} \operatorname{Be}(1, M), \boldsymbol{\Gamma}_{h}\right| H \stackrel{i i d}{\sim} H\) (random kernel locations), and \(\boldsymbol{\theta}_{h} \mid G_{0} \stackrel{i i d}{\sim} G_{0}\). This amounts to replacing the random measure \(G_{h}(\bullet)\) defined in (10) by just a single atom \(\boldsymbol{\theta}_{h}\). Compared to the former, this latter version of KSBP greatly reduces model complexity while still retaining some flexibility.

\subsection*{3.3 Probit and logit stick-breaking}

Chung and Dunson (2009) introduced a modification of the stick-breaking representation for DPs where the Beta random variables are replaced by normally distributed random variables transformed using the standard normal CDF. They refer to the resulting measure as the probit-stick breaking (PSB) process. The PSB is defined by
\[
G(\bullet)=\sum_{h=1}^{\infty}\left\{\Phi\left(\eta_{h}\right) \prod_{\ell<h}\left(1-\Phi\left(\eta_{\ell}\right)\right)\right\} \delta_{\boldsymbol{\theta}_{h}}(\bullet)
\]
where \(\eta_{h} \mid \mu \stackrel{i i d}{\sim} N(\mu, 1)\) and \(\boldsymbol{\theta}_{h} \mid G_{0} \stackrel{i i d}{\sim} G_{0}\). If \(\mu=0\), (11) reduces to a regular DP with \(M=1\), i.e. uniformly distributed sticks. Chung and Dunson (2009) also consider a covariate-dependent version of the PSB to model sets of related probability distributions. This is done by replacing the \(\eta_{h}\)
variables with suitable stochastic processes or regression functions. For instance, if \(\left\{\eta_{h}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\) denote independent Gaussian processes with unit variance, a dependent PSB can be defined as
\[
G_{\boldsymbol{x}}(\bullet)=\sum_{h=1}^{\infty}\left\{\Phi\left(\eta_{h}(\boldsymbol{x})\right) \prod_{\ell<h}\left[1-\Phi\left(\eta_{\ell}(\boldsymbol{x})\right)\right]\right\} \delta_{\boldsymbol{\theta}_{h}}(\bullet)
\]

A similar modification can be obtained by taking \(\eta_{h}(\boldsymbol{x})=\boldsymbol{x}^{T} \gamma_{h}\). More generally, let \(\eta_{h}(\boldsymbol{x})= \alpha_{h}+f_{h}(\boldsymbol{x})\) with \(\alpha_{h} \sim N(\mu, 1)\) and \(f_{h}: \mathbb{R}^{p} \rightarrow \mathbb{R}\) an unknown regression function, characterized by finitely many parameters \(\phi_{h}\), with \(\phi_{h} \sim \boldsymbol{H}\). Denote this model as \(\operatorname{PSBP}\left(\mu, \boldsymbol{H}, G_{0}\right)\). One main focus of the proposal in Chung and Dunson (2009) was variable selection. To that end, they assume the model
\[
y \mid \boldsymbol{x} \sim f(y \mid \boldsymbol{x})=\int N\left(y \mid \boldsymbol{x}^{\prime} \boldsymbol{\beta}, \tau^{-1}\right) d P_{\mathscr{X}}(\boldsymbol{\beta}, \tau), \quad P_{\mathscr{X}}=\left\{P_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\} \sim \operatorname{PSBP}\left(\mu, \boldsymbol{H}, G_{0}\right)
\]
where the variable selection structure is here introduced in \(\boldsymbol{H}\) and in \(G_{0}\), and by considering inclusion/exclusion indicators at the level of the atoms in (12). See further discussion on PSBP in Rodríguez and Dunson (2011). A related construction, termed the logit-stick breaking process was proposed in Ren et al. (2011), which essentially replaces the probit by a logit link in (11). Applications of logit-stick breaking processes to density regression can be found in Rigon and Durante (2020).

\subsection*{3.4 Hierarchical mixture of DP}

Consider again the case \(\mathscr{X}=\{1, \ldots, J\}\), as in the example presented in Section 1, and let \(\mathcal{G}=\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}=\left\{G_{1}, \ldots, G_{J}\right\}\). Motivated by the need to borrow strength across related studies (a situation also arising in applications of meta-analysis), Müller, Quintana and Rosner (2004) proposed a hierarchical DP model. In this construction, the probability distribution for group \(j\) is a weighted mixture of independent random measures. Specifically, the probability model for a group is defined as a mixture of a common distribution \(H_{0}\), shared by all groups, and an idiosyncratic component \(H_{j}\), which is specific to each group,
\[
G_{j}(\bullet)=\epsilon H_{0}(\bullet)+(1-\epsilon) H_{j}(\bullet),
\]
where \(\epsilon \in[0,1]\) controls the level of dependence in the set \(\mathcal{G}\), and \(H_{0}, H_{1}, \ldots, H_{J}\) are assumed to be independent DPs. The two extreme cases depicted in Figure 1 correspond to \(\epsilon=1\) for panel (a), i.e. a single common measure, and \(\epsilon=0\) for panel (b), i.e. independent model and no borrowing of strength. Model (13) represents then a trade-off between these two extreme options, allowing one to borrow strength through the common part, while retaining flexibility for the study-specific
part of the model. More recently, Wang and Rosner (2019) used this construction to propose a propensity score-based mixture model to combine subject-level information from randomized and registry studies, their goal being inference on a causal treatment effect.

Extending (13) to the case of continuous predictors can be easily accomplished by combining a study index, \(j\), continuous predictors \(\boldsymbol{z}\), and setting up
\[
G_{j, z}(\bullet)=\epsilon H_{0, z}(\bullet)+(1-\epsilon) H_{j, z}(\bullet),
\]
where \(H_{0, \boldsymbol{z}}, H_{1, \boldsymbol{z}}, \ldots, H_{J, \boldsymbol{z}}\) are now independent MacEachern's DDPs based on the continuous predictors \(\boldsymbol{z}\), incorporating dependence on predictors as in the LDDP or ANCOVA-DDP of Section 2.2.1, according to the available covariates types. The construction is easily modified to allow for study-specific variation in the weight assigned to the idiosyncratic component \(H_{j}\) by replacing \(\epsilon\) with \(\epsilon_{j}\).

A clever variation of this construction is introduced in Kolossiatis, Griffin and Steel (2013) who chose the weight \(\epsilon\) to ensure that \(G_{j}\) remains marginally a DP again. A more general version of the same construction appears in Camerlenghi et al. (2019).

\subsection*{3.5 Hierarchical DP of Teh et al. (2006)}

In the context of \(\mathscr{X}=\{1, \ldots, J\}\), Teh et al. (2006) proposed a model that induces an ANOVA type of dependence. In their construction, referred to as the hierarchical DP (HDP), the random probability measure for the \(j\) th group \(G_{j}, j=1, \ldots, J\), is a DP conditional on a common measure \(G\), which in turn is also a DP,
\[
G_{j}\left|M_{j}, G \stackrel{i n d}{\sim} D P\left(M_{j}, G\right), \quad j=1, \ldots, J, \quad G\right| M, G_{0} \sim D P\left(M, G_{0}\right)
\]

A main motivation behind the particular form adopted in (14) was to provide a model that allows for sharing clusters among related subpopulations. Teh et al. (2006) consider the analysis of text, where a primary goal was to share clusters among various documents within a cluster, and also to share clusters among various corpora. The HDP facilitates the construction of clusters at various levels, due to its hierarchical formulation. In fact, this clustering structure can be described in terms of a Chinese restaurant franchise, where at each of a collection of restaurants customers sit at tables organized by dishes, and dishes can be ordered from a global menu available to all restaurants. This construction, if restricted to a single restaurant, reduces to the usual Chinese restaurant process (Aldous, 1985) that is colloquially used to describe the DP.

\subsection*{3.6 The nested DP}

Also in the context of \(\mathscr{X}=\{1, \ldots, J\}\), Rodríguez, Dunson and Gelfand (2008) proposed an alternative model, referred to as the nested DP. In their construction the law of the random probability measure for the \(j\) th group \(G_{j}, j=1, \ldots, J\), is an infinite mixture of trajectories of DPs,
\[
G_{j} \stackrel{i n d}{\sim} \sum_{h=1}^{\infty} \pi_{h} \delta_{G_{h}^{*}}(\bullet), \quad j=1, \ldots, J, \quad G_{h}^{*} \mid M_{2}, H \stackrel{i . i . d .}{\sim} D P\left(M_{2}, H\right)
\]
where \(\pi_{h}=V_{h} \prod_{\ell<h}\left(1-V_{\ell}\right)\), with \(V_{h} \mid M_{1} \stackrel{\text { i.i.d. }}{\sim} \operatorname{Be}\left(1, M_{1}\right)\), for \(h=1,2, \ldots\). The main motivation behind (15) was to construct a clustering of individuals across the different groups, e.g. patients within different medical centers. The NDP model aims to simultaneously cluster patients within centers, borrowing information across centers for which similar clusters are detected, and to cluster different centers. This is then a type of multilevel clustering.

By way of comparison, it can be noted that in the HDP of Teh et al. (2006), the random measures in \(\mathcal{G}=\left\{G_{1}, \ldots, G_{J}\right\}\) share the same atoms but assign them different weights, while in the NDP two distributions \(G_{j_{1}}\) and \(G_{j_{2}}\) either share both atoms and weights (i.e. they are identical) or share nothing at all. Thus, the NDP allows for clusters at the level of the responses and also at the level of distributions, while the HDP allows for clusters only at the level of observations.

One of the limitations of the NDP is that for any two random measures \(G_{j_{1}}, G_{j_{2}}\) it supports only the two extreme cases of either all atoms and weights shared, i.e., \(G_{j_{1}}=G_{j_{2}}\), or no atoms shared, but does not allow any intermediate configuration with some atoms being shared. As a consequence, whenever there are ties of atoms between \(G_{j_{1}}\) and \(G_{j_{2}}\), the nested structure forces the two random distributions to be identical. For a discussion of this problem see Camerlenghi et al. (2019) who introduce the latent nested process as a more general hierarchical prior for random probability measures that avoids this restriction. More recently, Beraha, Guglielmi and Quintana (2020) propose the semi-hierarchical DP as an alternative solution to the limitations inherent to latent nested processes, with the added benefit of computationally efficient implementations to the comparison and clustering of potentially many subpopulations.

Like any discrete random probability measure, the NDP can be used to define random partitions. Model (15) could be written in short as \(G_{j} \sim \operatorname{DP}\left\{M_{1}, \mathrm{DP}\left(M_{2}, H\right)\right\}\). The outer DP, with total mass \(M_{1}\) gives rise to a partition of \(\mathscr{X}\). Consider now samples \(y_{j i} \sim G_{j}, i=1, \ldots, n_{j}\). The inner DP gives rise to random partitions of \(\mathscr{Y}_{j}=\left\{1, \ldots, n_{i}\right\}\), i.e., the NDP defines a nested partition of \(\mathscr{X}\) and \(\mathscr{Y}_{j}\), with the prior for the random partitions for \(\mathscr{Y}_{j}\) and \(\mathscr{Y}_{j^{\prime}}\) being equal in distribution when \(G_{j}=G_{j^{\prime}}\). Curiously, exactly the same random nested partition on \(\mathscr{X}\) and \(\mathscr{Y}_{j}\) is implied by the enriched DP (EDP) defined in Wade, Mongelluzzo and Petrone (2011). The EDP defines a random
probability measure for pairs \(\left(x_{i}, y_{i}\right)\) as \(P_{X}\left(x_{i}\right) P_{Y \mid X}\left(y_{i} \mid x_{i}\right)\), which, as discrete random probability measures, gives rise to the same random nested partition.

\subsection*{3.7 The product of independent DPs}

Alternatively, Gelfand and Kottas (2001) proposed an approach based on the product of independent random measures. In this construction the distribution for the \(j\) th group \(G_{j}, j=1, \ldots, J\), is given by
\[
G_{j}(\bullet) \equiv H_{j}(\bullet) \prod_{\ell<j} H_{\ell}(\bullet), \quad j=1, \ldots, J
\]
where
\[
H_{j} \mid M_{j}, H_{0 j} \stackrel{i n d .}{\sim} D P\left(M_{j}, H_{0 j}\right), \quad j=1, \ldots, J .
\]

The motivation for this construction arises from the need to define models that induce stochastic ordering for the random group specific distributions \(G_{j}\). The ordering holds with probability 1 in the prior and so is also satisfied a posteriori.

\subsection*{3.8 Other constructions}

Chung and Dunson (2009) proposed a similar construction, referred to as the local DP, where the stick-breaking weights selected to define the probability weights depend on a set of random locations and their distances to a given predicted value. In this construction, the support points also depend on predictors.

Fuentes-García, Mena and Walker (2009) considered a dependent variation of geometric-weights stick-breaking processes (Mena, Ruggiero and Walker, 2011). In this construction, the stick-breaking weights are replaced by their expected value, thus reducing the number of parameters.

Dependent neutral to the right processes and correlated two-parameter Poisson-Dirichlet processes have been proposed by Epifani and Lijoi (2010) and Leisen and Lijoi (2011), respectively, by considering suitable Lévy copulas. The general class of dependent normalized completely random measures has been discussed, for instance, by Lijoi, Nipoti and Prünster (2014).

Another type of construction stems from the fact that the Dirichlet process is also a special case of a normalized random measure with independent increments (NRMI), as described in Regazzini, Lijoi and Prünster (2003). This means that if \(F\) has a DP distribution, then it can be expressed in the form
\[
F(\bullet)=\frac{\mu(\bullet)}{\mu(\Omega)}
\]
where \(\Omega\) is the space where the DP is defined, and \(\mu\) is a completely random measure on \((\Omega, \mathcal{B}(\Omega))\), that is, for any collection of disjoint sets \(A_{1}, A_{2}, \ldots\) in \(\mathcal{B}(\Omega)\), the Borel \(\sigma\)-field in \(\Omega\), the random
variables \(\mu\left(A_{1}\right), \mu\left(A_{2}\right), \ldots\) are independent, and \(\mu\left(\cup_{j=1}^{\infty} A_{j}\right)=\sum_{j=1}^{\infty} \mu\left(A_{j}\right)\) holds true a.s. See, e.g. James, Lijoi and Prünster (2009). As shown in Ferguson (1973), the Dirichlet process arises as the normalized version of a Gamma process. Barrios et al. (2013), Favaro and Teh (2013) and Argiento, Guglielmi and Pievatolo (2010) discuss modeling with mixtures of NRMIs, and in particular discuss practical implementation of posterior simulation for such models. See additional MCMC implementation details in Building on related ideas, Epifani and Lijoi (2010) and Leisen and Lijoi (2011), proposed dependent neutral to the right processes and correlated two-parameter Poisson-Dirichlet processes, respectively, by considering suitable Lévy copulas. A more general class of dependent normalized completely random measures has been discussed, for instance, by Lijoi, Nipoti and Prünster (2014). This construction has also motivated work on defining DDPs by way of introducing dependence in NRMIs. Lin, Grimson and Fisher (2010) used this idea to propose a Markov chain of Dirichlet processes, and other extensions to normalized random measured are described in Chen, Ding and Buntine (2012) and in Chen et al. (2013).

\section*{4. THE INDUCED CONDITIONAL DENSITY APPROACH}

The approaches described so far yield valid inferences when the set of predictors \(\boldsymbol{x}\) are fixed by design or are random but exogenous. Notice that the exogeneity assumption permits us to focus on the problem of conditional density estimation, regardless of the data generating mechanism of the predictors, that is, if they are randomly generated or fixed by design (see, e.g., Barndorff-Nielsen, 1973, 1978). Under the presence of endogenous predictors, both the response and the predictors should be modeled jointly.

In the context of continuous responses and predictors, Müller, Erkanli and West (1996) proposed a DPM of multivariate Gaussian distributions for the complete data \(\boldsymbol{d}_{i}=\left(y_{i}, \boldsymbol{x}_{i}\right)^{\prime}, i=1, \ldots, n\), and looked at the induced conditional distributions. Although Müller, Erkanli and West (1996) focused on the mean function only, \(m(\boldsymbol{x})=E(y \mid \boldsymbol{x})\), their method can be easily extended to provide inferences for the conditional density at covariate level \(\boldsymbol{x}\). The model is given by
\[
\boldsymbol{d}_{i} \mid G \stackrel{i i d}{\sim} \int N_{k}\left(\boldsymbol{d}_{i} \mid \boldsymbol{\mu}, \boldsymbol{\Sigma}\right) d G(\boldsymbol{\mu}, \boldsymbol{\Sigma}),
\]
and
\[
G \mid M, G_{0} \sim D P\left(M, G_{0}\right)
\]
where \(k=p+1\) is the dimension of the complete data vector \(\boldsymbol{d}_{i}\), and the baseline distribution \(G_{0}\) is the conjugate normal-inverted-Wishart (IW) distribution \(G_{0} \equiv N_{k}\left(\boldsymbol{\mu} \mid \boldsymbol{m}_{1}, \kappa_{0}^{-1} \boldsymbol{\Sigma}\right) \times I W_{k}\left(\boldsymbol{\Sigma} \mid \nu_{1}, \boldsymbol{\Psi}_{1}\right)\). The model is completed with conditionally conjugate priors and hyperpriors on
\(m_{1}, \kappa_{0}\) and \(\boldsymbol{\Psi}\), and, if desired, a gamma hyperprior on \(M\). The model induces a weight-dependent mixture model for the regression,
\[
f_{\boldsymbol{x}}(y)=\sum_{h=1}^{\infty} \omega_{h}(\boldsymbol{x}) N\left(y \mid \beta_{0 h}+\boldsymbol{x}^{\prime} \boldsymbol{\beta}_{h}, \sigma_{h}^{2}\right),
\]
where
\[
\omega_{h}(\boldsymbol{x})=\frac{w_{h} N_{p}\left(\boldsymbol{x} \mid \boldsymbol{\mu}_{2 h}, \boldsymbol{\Sigma}_{22 h}\right)}{\sum_{\ell=1}^{\infty} w_{\ell} N_{p}\left(\boldsymbol{x} \mid \boldsymbol{\mu}_{2 \ell}, \boldsymbol{\Sigma}_{22 \ell}\right)}, \quad h=1,2, \ldots
\]
\(\beta_{0 h}=\mu_{1 h}-\boldsymbol{\Sigma}_{12 h} \boldsymbol{\Sigma}_{22 h}^{-1} \boldsymbol{\mu}_{2 h}, \boldsymbol{\beta}_{h}=\boldsymbol{\Sigma}_{12 h} \boldsymbol{\Sigma}_{22 h}^{-1}\), and \(\sigma_{h}^{2}=\sigma_{11 h}^{2}-\boldsymbol{\Sigma}_{12 h} \boldsymbol{\Sigma}_{22 h}^{-1} \boldsymbol{\Sigma}_{21 h}\). Here, the weights \(w_{h}\) follow the usual DP stick-breaking construction, and the remaining elements arise from the standard partition of the vectors of means and (co)variance matrices given by
\[
\boldsymbol{\mu}_{h}=\binom{\mu_{1 h}}{\boldsymbol{\mu}_{2 h}} \quad \text { and } \quad \boldsymbol{\Sigma}_{h}=\left(\begin{array}{cc}
\sigma_{11 h}^{2} & \boldsymbol{\Sigma}_{12 h} \\
\boldsymbol{\Sigma}_{21 h} & \boldsymbol{\Sigma}_{22 h}
\end{array}\right)
\]
respectively.
The induced conditional density approach of Müller, Erkanli and West (1996) can be easily extended to handle mixed continuous, \(\boldsymbol{x}_{C}\), and discrete predictors, \(\boldsymbol{x}_{D}\), by considering a DPM model of product of appropriate kernels for discrete \(k_{D}\) and continuous \(k_{D}\) variables,
\[
\boldsymbol{d}_{i} \mid G \stackrel{i i d}{\sim} \int k_{D}\left(\boldsymbol{x}_{i D} \mid \boldsymbol{\theta}_{1}\right) k_{C}\left(y_{i}, \boldsymbol{x}_{C} \mid \boldsymbol{\theta}_{2}\right) d G\left(\boldsymbol{\theta}_{1}, \boldsymbol{\theta}_{2}\right)
\]
i.e., assuming a multiplicative structure in the joint model for ( \(y, \boldsymbol{x}_{D}, \boldsymbol{x}_{C}\) ) that mimics conditional independence of \(\left(y, \boldsymbol{x}_{C}\right)\) and \(\boldsymbol{x}_{D}\) given suitable parameter vectors \(\boldsymbol{\theta}_{1}\) and \(\boldsymbol{\theta}_{2}\). Similar types of models, but looking only at the induced partition structures, are discussed in Müller and Quintana (2010). In particular, Müller, Quintana and Rosner (2011) proposed a version of (17) that may be viewed as integrating out the random measure \(G\) in (17), retaining only the random partition model, while still allowing for covariate dependence in the prior. This approach exploits the connection between the DP and product partition models. See, e.g., Quintana and Iglesias (2003).

We introduced the conditional density regression approach assuming endogenous predictors, when the construction of a joint probability model for ( \(y_{i}, \boldsymbol{x}_{i}\) ) is natural. However, the same construction can be used to achieve the desired smooth locally weighted mixture of linear regressions even when the \(\boldsymbol{x}_{i}\) are exogenous, or even if they are not random at all. The choice of model depends largely on properties of the model and ease of prior specification, tempered by computational concerns.

\section*{5. IMPLIED RANDOM PARTITIONS AND OTHER USES OF THE DDP MODEL}

One of the common applications of the DP mixture model (1) is to define a random partition and allow statistical inference on such partitions. Consider an equivalent statement of i.i.d. sampling from (1) as a hierarchical model
\[
y_{i} \mid \theta_{i} \sim p\left(y_{i} \mid \theta_{i}\right) \quad \text { and } \quad \theta_{i} \sim G,
\]
\(i=1, \ldots, n\). The discrete nature of the DP random measure \(G\) implies positive probabilities of ties among the \(\theta_{i}\) with \(K \leq N\) unique values \(\left\{\theta_{1}^{\star}, \ldots, \theta_{K}^{\star}\right\}\). Defining \(S_{j}=\left\{i: \theta_{i}=\theta_{j}^{\star}\right\}\) defines a partition \(\{1, \ldots, n\}=\bigodot_{j} S_{j}\). A common application of the DP mixture model is to derive inference on such partitions \(\rho=\left\{S_{1}, \ldots, S_{K}\right\}\), and interpret the partitioning subsets as meaningful subpopulations of the experimental units (e.g., patient subpopulations). In anticipation of the upcoming generalization to the DDP, we introduce a slightly different but equivalent definition of the clusters \(S_{j}\). Recall the representation (2) of DP random measure, \(G=\sum w_{h} \delta_{\widetilde{\theta}_{h}}\) Then the non-empty sets \(R_{h}=\left\{i: \theta_{i}=\right. \left.\widetilde{\theta}_{h}\right\}\) describe the same partition \(\rho\). We switched from indexing clusters by their common unique \(\theta_{i}\) values to identifying clusters by the matching atoms in \(G\). Similarly we can set up a model for independent sampling using a DDP prior. Specifically, consider
\[
y_{i} \mid \theta_{i} \sim p\left(y_{i} \mid \theta_{i}\right) \quad \text { and } \quad \theta_{i} \mid x_{i}=x \sim G_{x},
\]
\(i=1, \ldots, n\), with a DDP prior on \(\mathcal{G}=\left\{G_{x}, x \in X\right\}\). For the moment assume a categorical covariate \(x_{i} \in\left\{1, \ldots, n_{x}\right\}\), and let \(G_{x}, x=1, \ldots, n_{x}\) denote the (marginal) random measures, and let \(I_{x}=\left\{i: x_{i}=x\right\}\) denote the subpopulation with covariate \(x\). First, by the earlier argument the model implies a random partition \(\rho_{x}\) of \(I_{x}\), marginally, for each \(x\). Indexing clusters by the corresponding atom in \(G_{x}\) implicitly defines a joint prior on \(\left\{\rho_{x}, x \in X\right\}\), or, alternatively, defines a partition of \(\{1, \ldots, n\}\) with clusters \(S_{j}\) that cut across \(I_{x}\). In particular, the model implies a joint prior on \(\left(\rho_{x}, \rho_{x^{\prime}}\right)\) for any \(x \neq x^{\prime}\), and it allows for shared clusters across subpopulations. Different assumptions on various model aspects, such as dispersion in the baseline distribution, or total mass parameter, would have a practical effect on the this joint prior. Curiously, in contrast to the DP mixture model, the DDP model is not commonly used for inference on these implied random partition(s).

Another feature of the DDP model is inference about distributional homogeneity. To be specific, consider again the context of independent sampling in (18) with a categorical covariate \(x \in\left\{1, \ldots, n_{x}\right\}\) and let \(f_{x}(y)=\int p(y \mid \theta) d G_{x}(\theta)\) denote the implied marginal distribution of \(y_{i} \mid x_{i}=x\). In many applications investigators might be interested in the event \(f_{x}=f_{x^{\prime}}\) for \(x \neq x^{\prime}\).

While the DDP prior, short of a pathological special case, implies zero prior probability for exact equality, posterior inference includes meaningful posterior probabilities for \(\left\{d\left(f_{x}, f_{x^{\prime}}\right)>\epsilon\right\}\) for any well defined distance of the two distributions. Specifics would depend on particular applications. Related summaries, for example, by displaying posterior means for \(f_{x}\) over \(x\) are shown in some papers using DDP priors for density regression. See, e.g. Gutiérrez et al. (2019).

\section*{6. APPLICATION TO AUTOREGRESSIVE MODELS}

We illustrate some of the discussed DDP-based nonparametric regression models. We implement inference under the ANOVA-DDP or LDDP model of (5) and conditional density regression as in (16) to model (auto-)regression on \(x_{t}=y_{t-1}\) in time series data, using the LDDP model
\[
y_{t}\left|y_{t-1}=y, \beta_{t 0}, \beta_{t 1}, \sigma_{t}^{2} \sim N\left(\beta_{t 0}+\beta_{t 1} y, \sigma_{t}^{2}\right), \quad\left(\beta_{t 0}, \beta_{t 1}, \sigma_{t}^{2}\right)\right| G \stackrel{\mathrm{iid}}{\sim} G, \quad G \sim D P\left(M, G_{0}(\cdot \mid \boldsymbol{\eta})\right),
\]
where \(t=2, \ldots, n\), i.e. we mix over the linear coefficients and the variance. The dependence in (20) is conveyed through linear functions of the first lagged response in the atoms, keeping common weights. Here, \(G_{0}(\cdot \mid \boldsymbol{\eta})\) is the centering measure with hyperparameters \(\boldsymbol{\eta}\). Following Jara et al. (2011) we use \(G_{0} \equiv N_{2}\left(\boldsymbol{\beta} \mid \boldsymbol{\mu}_{b}, \boldsymbol{S}_{b}\right) \Gamma\left(\sigma^{-2} \mid \tau_{1} / 2, \tau_{2} / 2\right)\), and complete the prior specification as
\[
\begin{aligned}
M \mid a_{0}, b_{0} \sim \Gamma\left(a_{0}, b_{0}\right), & \tau_{2} \mid \tau_{s_{1}}, \tau_{s_{2}} \sim \Gamma\left(\tau_{s_{1}} / 2, \tau_{s_{2}} / 2\right), \\
\boldsymbol{\mu}_{b} \mid \boldsymbol{m}_{0}, \boldsymbol{S}_{0} \sim N_{p}\left(\boldsymbol{m}_{0}, \boldsymbol{S}_{0}\right), & \boldsymbol{S}_{b} \mid \nu, \boldsymbol{\Psi} \sim I W_{p}(\nu, \boldsymbol{\Psi}) .
\end{aligned}
\]

For this illustration, we consider two specific datasets:
Data set D1 are the Old Faithful geyser data (Härdle, 1991), available as part of the datasets library available in R , consisting of \(n=272\) observations on eruption times (in minutes) and waiting times to the next eruption (also in minutes). Data set D2 is a time series of the Standard \& Poor's 500 index, from February 9, 1993 through February 9, 2015. It is available in the R package pdfetch (Reinhart, 2019), using the command
```
pdfetch_YAHOO("SPY",fields = "adjclose",
    from = as.Date("1993-02-09"), to = as.Date("2015-02-09"))
```


In the following results we compare inference under model (20) with inference under density regression, as in (16), again using \(x_{t}=y_{t-1}\). Recall that a conditional density approach is based on a DPM model for \(\left\{\left(y_{t}, y_{t-1}\right): t=2, \ldots, n\right\}\).

In all cases, we used hyperparameters as in Jara et al. (2011). Results for D1 are shown in Figure 3. In particular, we show a comparison posterior inference for \(G_{x}\) for (a) \(y_{t-1}=58\), (b)

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0877b326-0bf3-4110-817e-01cbc5995ace-24.jpg?height=783&width=1604&top_left_y=255&top_left_x=242}
\captionsetup{labelformat=empty}
\caption{FIG 3. Old Faithful Geyser data: Posterior estimated \(G_{x}\), i.e., posterior predictive densities (mean and point-wise 95\% HPD intervals) for the waiting times at lagged times (a) \(y_{t-1}=58\), (b) \(y_{t-1}=76\) and (c) \(y_{t-1}=82\), including \(95 \%\) HPD credibility bands. The red curve shows inference under the LDDP model. The green curve shows inference under the conditional density approach.}
\end{figure}
\(y_{t-1}=76\) and (c) \(y_{t-1}=82\). While there are some model-specific differences in the estimated distributions \(G_{x}\), they both largely agree on the bimodal nature.

Figure 4 shows \(G_{x}\) over a grid of lagged values \(x_{t}=y_{t-1}\), under the conditional density approach. In this figure, the bimodality is also seen in the data (red dots). The solid blue curve in Figure 4 shows the posterior mean \(E\left(y_{t} \mid y_{t-1}\right)\) with \(95 \%\) credibility bands (blue dashed curves).

In contrast, similar inference for the LDDP (not shown) shows a straight line for the mean process \(E\left(G_{x} \mid \boldsymbol{y}\right)\), as a function of \(x\), as is implied by the linear structure of \(\theta_{h}(x)\) under the LDDP. See also Figure 5 below

Figure 5 shows the same results for the S\&P500 data, using the same models as above. As before, the data are shown as red dots, and blue lines show the posterior predictive means (solid line) together with a \(95 \%\) HPD interval (dashed line). Interestingly, for the LDDP model the HPD lines fall outside the plotting region.

\section*{7. CONCLUDING REMARKS}

DDPs have come a long way since they were originally proposed. By its very definition, a DDP has the potential to incorporate covariate indexing (dependence) either in the atoms or the weights or both. The results in Barrientos, Jara and Quintana (2012) show that under full support of the stochastic processes that are used to convey covariate dependence, the resulting DDP has full sup-

\section*{Conditional DDP}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0877b326-0bf3-4110-817e-01cbc5995ace-25.jpg?height=1030&width=1247&top_left_y=760&top_left_x=474}
\captionsetup{labelformat=empty}
\caption{Fig 4. Old Faithful Geyser data: Posterior estimated densities \(G_{x}\) for a grid of lagged values \(x_{t}=y_{t-1}\). The blue curve shows conditional mean process (solid line), with 95\% credible intervals (dashed line).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0877b326-0bf3-4110-817e-01cbc5995ace-26.jpg?height=1035&width=1526&top_left_y=766&top_left_x=189}
\captionsetup{labelformat=empty}
\caption{Fig 5. S\&P500 data: Posterior estimated densities \(G_{x}\) for a grid of \(x=y_{t-1}\). The blue lines show the conditional mean process, with dashed lines for 95\% HDP intervals. Panel (a) shows inference under the LDDP model, and panel (b) shows inference under the conditional density approach.}
\end{figure}
port in the space \(\mathscr{F}=\left\{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\). This holds true for all of the basic DDP constructions: single-atoms, single-weights, and with dependence in both. A natural question is then: which DDP version is the best? There is no final answer to this question, although DDPs with dependence in both atoms and weights are less commonly found, mostly due to the computational complexity of implementation entails. An exception to this is the conditional approach described in Section 4. In broad terms, the single-weights models are typically easier to fit, as the standard algorithms designed to implement posterior simulation in the context of DPs can be applied with minor adjustments. See, e.g., the computational aspects in De Iorio et al. (2004). The same applies for the LDDP.

On the other hand, the single-atoms models are typically less attractive from a computational viewpoint, mainly due to how covariate dependence is encoded in the definition of the weight processes \(\left\{w_{h}(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\right\}\). However, the single-atoms DDP allows for the prior probability distribution on the partitions to change with \(\boldsymbol{x}\), a feature that is not supported by the singleweights DDP. For a formal description of this feature, let \(\mathcal{G}=\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) denote the family of random probability measures with DDP prior, as before. Let \(\rho_{\boldsymbol{x}}\) denote the partition of \(\{1, \ldots, n\}\) that is implied by a hypothetical sample from \(G_{\boldsymbol{x}}\), of size \(n\). Under the single-weights DDP, \(p\left(\rho_{\boldsymbol{x}} \mid \mathcal{G}\right)\) is invariant across \(\boldsymbol{x}\); but not so under the single-atoms DDP. This is the case since the prior on the random partition \(\rho_{\boldsymbol{x}}\) is determined by the weights in \(G_{\boldsymbol{x}}\).

Models for dependent probability distributions do not easily allow for the incorporation of existing prior information about arbitrary functionals. A modeler is unlikely to have prior knowledge about all aspects of a collection of probability measures, but could have real historical prior information about specific functionals (such as the mean or quantile functions). For example, such information could be obtained as the product of applying parametric or (classical) nonparametric approaches to previous data. Furthermore, even in models for single (non-dependent) probability measures, the derivation of the induced distribution for arbitrary functionals is challenging and, thus, usually not exploited. This makes the prior elicitation process difficult. We refer the reader to Lijoi and Prünster (2009) for an exhaustive summary of existing results concerning distributional properties of functional of single and discrete random probability measures.

In the context of a single probability measure, Kessler, Hoff and Dunson (2015) proposed a clever construction of a BNP model with a given distribution on a finite set of functionals. Their approach is based on the conditional distribution of a standard BNP prior, given the functionals of interest. A Metropolis-Hastings MCMC algorithm is proposed to explore the posterior distribution under the marginally specified BNP model, where the standard BNP model is used as a candidate
generating model, and that is closely related to the well-known importance-sampling approach for assessing prior sensitivity. Their MCMC algorithm is developed for DP-based models and relies on the marginalization of the random probability measure. Thus, a Monte Carlo approximation of the functionals of interest is employed at any step of the MCMC algorithm to obtain approximated posterior samples of the functionals of interest. The study of extensions of the approach proposed by Kessler, Hoff and Dunson (2015) to the context of sets of predictor-dependent probability measures is a topic of interest for future research.

An interesting topic has been recently brought up by Campbell et al. (2019). They introduced a relaxed version of the notion of exchangeability, local exchangeability, which considers bounded changes in total variation norm of the distribution of observations under permutations of data having nearby covariate values. This notion generalizes that of exchangeability and partial exchangeability. The work by Campbell et al. (2019) discusses conditions under which a version of de Finetti's Theorem holds in such a way that a DDP is the corresponding de Finetti measure, i.e. conditional independence of the observations under a DDP is still true. The study of extensions and applications of these and related results is another topic of interest for future research.

The bulk of work on the DDP and related methods focuses on the family of conditional distributions \(\mathcal{G}=\left\{G_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\right\}\) and models where an observation \(y\) is associated with a single value of the covariate \(\boldsymbol{x}\). When data are longitudinal, spatial or functional, the observations may be considered to have dependence that cannot be captured by the marginal distributions \(G_{\boldsymbol{x}}\). See, for example, Xu , MacEachern and Xu (2015) who separate dependence in financial data series from the marginal distributions. Many open questions remain in this direction.

Finally, the idea of introducing dependence through normalization, e.g. as mentioned earlier in Section 3.8 can be further exploited and extended to more general cases, including going beyond the context of DDPs.

\section*{ACKNOWLEDGEMENTS}
A. Jara's and F. Quintana's research is supported by Millennium Science Initiative of the Ministry of Economy, Development, and Tourism, grant "Millennium Nucleus Center for the Discovery of Structures in Complex Data". A. Jara is also supported by Fondecyt grant 1180640, F. Quintana is also supported by Fondecyt grant 1180034. P. Müller acknowledges partial support from grant NSF/DMS 1952679 from the National Science Foundation, and under R01 CA132897 from the U.S. National Cancer Institute.

\section*{REFERENCES}

Aldous, D. J. (1985). Exchangeability and related topics. In École d'été de probabilités de Saint-Flour, XIII—1983. Lecture Notes in Math. 1117 1-198. Springer, Berlin.
Argiento, R., Guglielmi, A. and Pievatolo, A. (2010). Bayesian density estimation and model selection using nonparametric hierarchical mixtures. Computational Statistics \& Data Analysis 54 816-832.
Ascolani, F., Lijoi, A. and Ruggiero, M. (2020). Predictive inference with FlemingViot-driven dependent Dirichlet processes. Bayesian Analysis. Advance publication.
Barndorff-Nielsen, O. (1973). On M-ancillarity. Biometrika 60 447-455.
Barndorff-Nielsen, O. (1978). Information and exponential families in statistical theory. John Wiley \& Sons, Ltd., Chichester Wiley Series in Probability and Mathematical Statistics.
Barrientos, A. F., Jara, A. and Quintana, F. A. (2012). On the support of MacEachern's dependent Drichlet processes and extensions. Bayesian Analysis 7277-310.
Barrios, E., Lijoi, A., Nieto-Barajas, L. E. and Prünster, I. (2013). Modeling with normalized random measure mixture models. Statistical Science 28 313-334.
Beraha, M., Guglielmi, A. and Quintana, F. A. (2020). The semi-hierarchical Dirichlet Process and its application to clustering homogeneous distributions.
Camerlenghi, F., Dunson, D. B., Lijoi, A., Prünster, I. and Rodríguez, A. (2019). Latent Nested Nonparametric Priors (with Discussion). Bayesian Anal. 14 1303-1356.
Campbell, T., Syed, S., Yang, C.-Y., Jordan, M. I. and Broderick, T. (2019). Local Exchangeability.
Caron, F., Davy, M., Doucet, A., Duflos, E. and Vanheeghe, P. (2008). Bayesian Inference for Linear Dynamic Models with Dirichlet Process Mixtures. IEEE Transactions on Signal Processing 56 71-84.
Chen, C., Ding, N. and Buntine, W. (2012). Dependent Hierarchical Normalized Random Measures for Dynamic Topic Modeling. In Proceedings of the 29th International Conference on Machine Learning (ICML-12) (J. Langford and J. Pineau, eds.). ICML '12 895-902. Omnipress, New York, NY, USA.
Chen, C., Rao, V., Buntine, W. and Teh, Y. W. (2013). Dependent Normalized Random Measures. In Proceedings of the 30th International Conference on Machine Learning (S. Dasgupta and D. McAllester, eds.). Proceedings of Machine Learning Research 28 969-977. PMLR, Atlanta, Georgia, USA.
Chipman, H. A., George, E. I. and McCulloch, R. E. (2010). BART: Bayesian additive regression trees. Ann. Appl. Stat. 4 266-298.
Chung, Y. and Dunson, D. B. (2009). Nonparametric Bayes conditional distribution modeling with variable selection. Journal of the American Statistical Association 104 1646-1660.
Cifarelli, D. and Regazzini, E. (1978). Problemi statistici non parametrici in condizioni di scambialbilita parziale e impiego di medie associative Technical Report, Quaderni Istituto Matematica Finanziaria, Torino.
De Iorio, M., Müller, P., Rosner, G. L. and MacEachern, S. N. (2004). An ANOVA model for dependent random measures. Journal of the American Statistical Association 99 205-215.
De Iorio, M., Johnson, W. O., Müller, P. and Rosner, G. L. (2009). Bayesian nonparametric non-proportional hazards survival modelling. Biometrics 65 762-771.
De la Cruz-Mesía, R., Quintana, F. A. and Müller, P. (2007). Semiparametric Bayesian classification with longitudinal markers. Journal of the Royal Statistical Society. Series C. Applied Statistics 56 119-137.
Devroye, L. (1986). Non-Uniform Random Variate Generation(originally published with. Springer-Verlag.
Di Lucca, M. A., Guglielmi, A., Müller, P. and Quintana, F. A. (2012). A simple class of Bayesian nonpara-
metric autoregression models. Bayesian Analysis 8 63-88.
Duan, J. A., Guindani, M. and Gelfand, A. E. (2007). Generalized spatial Dirichlet process models. Biometrika 94 809-825.

Dunson, D. B. and Herring, A. H. (2006). Semiparametric Bayesian latent trajectory models Technical Report, ISDS Discussion Paper 16, Duke University, Durham, NC, USA.
Dunson, D. B. and Park, J. H. (2008). Kernel stick-breaking processes. Biometrika 95 307-323.
Dunson, D. B., Pillai, N. and Park, J. H. (2007). Bayesian density regression. Journal of the Royal Statistical Society, Series B69163-183.
Epifani, I. and Lijoi, A. (2010). Nonparametric priors for vectors of survival functions. Statistica Sinica 20 14551484.

Faraway, J. J. (2016). Extending the linear model with R. Chapman \& Hall/CRC Texts in Statistical Science Series. CRC Press, Boca Raton, FL Generalized linear, mixed effects and nonparametric regression models, Second edition [of MR2192856].
Favaro, S. and Teh, Y. W. (2013). MCMC for Normalized Random Measure Mixture Models. Statistical Science 28 335-359.

Ferguson, T. S. (1973). A Bayesian analysis of some nonparametric problems. Annals of Statistics 1 209-230.
Ferguson, T. S. (1974). Prior distribution on the spaces of probability measures. Annals of Statistics 2615-629.
Fuentes-García, R., Mena, R. and Walker, S. G. (2009). A nonparametric dependent process for Bayesian regression. Statistics and Probability Letters 79 1112-1119.
Gelfand, A. E. and Kottas, A. (2001). Nonparametric Bayesian modeling for stochastic order. Annals of the Institute of Statistical Mathematics 53 865-876.
Gelfand, A. E., Kottas, A. and MacEachern, S. N. (2005). Bayesian nonparametric spatial modeling with Dirichlet process mixing. Journal of the American Statistical Association 100 1021-1035.
Giudici, P., Mezzetti, M. and Muliere, P. (2003). Mixtures of Dirichlet process priors for variable selection in survival analysis. Journal of Statistical Planning and Inference 111 101-115.
Green, P. J. and Richardson, S. (2001). Modelling heterogeneity with and without the Dirichlet process. Scandinavian Journal of Statistics. Theory and Applications 28 355-375.
Griffin, J. E. and Steel, M. F. J. (2006). Order-based dependent Dirichlet processes. Journal of the American Statistical Association 101 179-194.
Griffin, J. E. and Steel, M. F. J. (2010). Bayesian nonparametric modelling with the Dirichlet process regression smoother. Statistica Sinica 20 1507-1527.
Griffin, J. E. and Steel, M. F. J. (2011). Stick-breaking autoregressive processes. Journal of Econometrics 162 383-396.

Gutiérrez, L., Mena, R. H. and Ruggiero, M. (2016). A time dependent Bayesian nonparametric model for air quality analysis. Computational Statistics \& Data Analysis 95 161-175.
Gutiérrez, L., Barrientos, A. F., González, J. and Taylor-Rodríguez, D. (2019). A Bayesian Nonparametric Multiple Testing Procedure for Comparing Several Treatments Against a Control. Bayesian Analysis 14 649-675.
Györfi, L., Kohler, M., Krzy* zak, A. and Walk, H. (2002). A distribution-free theory of nonparametric regression. Springer Series in Statistics. Springer-Verlag, New York.
Härdle, W. (1991). Smoothing techniques. Springer Series in Statistics. Springer-Verlag, New York With implementation in S .

James, L. F., Lijoi, A. and Prünster, I. (2009). Posterior Analysis for Normalized Random Measures with Independent Increments. Scandinavian Journal of Statistics 36 76-97.
Jara, A. and Hanson, T. (2011). A class of mixtures of dependent tail-free processes. Biometrika 98 553-566.
Jara, A., Lesaffre, E., De Iorio, M. and Quintana, F. A. (2010). Bayesian semiparametric inference for multivariate doubly-interval-censored data. The Annals of Applied Statistics 4 2126-2149.
Jara, A., Hanson, T., Quintana, F., Müller, P. and Rosner, G. L. (2011). DPpackage: Bayesian Semi- and Nonparametric Modeling in R. Journal of Statistical Software 40 1-30.
Kessler, D., Hoff, P. and Dunson, D. (2015). Marginally specified priors for non-parametric Bayesian estimation. Journal of the Royal Statistical Society, Series B 7735-58.
Klemelä, J. (2014). Multivariate nonparametric regression and visualization. Wiley Series in Computational Statistics. John Wiley \& Sons, Inc., Hoboken, NJ With R and applications to finance.
Kolossiatis, M., Griffin, J. E. and Steel, M. F. J. (2013). On Bayesian nonparametric modelling of two correlated distributions. Statistics and Computing 23 1-15. MR3018346
Lau, J. W. and So, M. K. P. (2008). Bayesian mixture of autoregressive models. Computational Statistics and Data Analysis 53 38-60.
Lavine, M. (1992). Some aspects of Polya tree distributions for statistical modeling. The Annals of Statistics 20 1222-1235.

Leisen, F. and Lijoi, A. (2011). Vectors of two-parameter Poisson-Dirichlet processes. Journal of Multivariate Analysis 102 482-495.
Lijoi, A., Nipoti, B. and Prünster, I. (2014). Bayesian inference with dependent normalized completely random measures. Bernoulli 20 1260-1291.
Lijoi, A. and Prünster, I. (2009). Distributional properties of means of random probability measures. Statistical Surveys 34795.
Lin, D., Grimson, E. and Fisher, J. W. (2010). Construction of Dependent Dirichlet Processes based on Poisson Processes. In Advances in Neural Information Processing Systems 23 (J. D. Lafferty, C. K. I. Williams, J. ShaweTaylor, R. S. Zemel and A. Culotta, eds.) 1396-1404. Curran Associates, Inc.
Lo, A. Y. (1984). On a class of Bayesian nonparametric estimates I: Density estimates. The Annals of Statistics 12 351-357.

MacEachern, S. N. (1999). Dependent nonparametric processes. In ASA Proceedings of the Section on Bayesian Statistical Science, Alexandria, VA. American Statistical Association.
MacEachern, S. N. (2000). Dependent Dirichlet processes Technical Report, Department of Statistics, The Ohio State University.
Mena, R. H., Ruggiero, M. and Walker, S. G. (2011). Geometric stick-breaking processes for continuous-time Bayesian nonparametric modeling. Journal of Statistical Planning and Inference 141 3217-3230.
Mena, R. H. and Ruggiero, M. (2016). Dynamic density estimation with diffusive Dirichlet mixtures. Bernoulli 22 901-926.

Mira, A. and Petrone, S. (1996). Bayesian hierarchical nonparametric inference for change-point problems. In Bayesian Statistics 5 (J. M. Bernardo, J. O. Berger, A. P. Dawid and A. F. M. Smith, eds.). Oxford University Press.
Muliere, P. and Petrone, S. (1993). A Bayesian predictive approach to sequential search for an optimal dose: parametric and nonparametric models. Journal of the Italian Statistical Society 2 349-364.

Müller, P., Erkanli, A. and West, M. (1996). Bayesian curve fitting using multivariate normal mixtures. Biometrika 83 67-79.
Müller, P., Quintana, F. A. and Rosner, G. (2004). A method for combining inference across related nonparametric Bayesian models. Journal of the Royal Statistical Society, Series B 66735-749.
Müller, P. and Quintana, F. A. (2010). Random partition models with regression on covariates. Journal of Statistical Planning and Inference 140 2801-2808.
Müller, P., Quintana, F. A. and Rosner, G. L. (2011). A product partition model with regression on covariates. Journal of Computational and Graphical Statistics 20 260-278.
Müller, P., Quintana, F. A., Jara, A. and E, H. T. (2015). Bayesian Nonparametric Data Analysis. Springer, New York, USA.
Prünster, I. and Ruggiero, M. (2013). A Bayesian nonparametric approach to modeling market share dynamics. Bernoulli 19 64-92.

Quintana, F. A. and Iglesias, P. L. (2003). Bayesian clustering and product partition models. Journal of The Royal Statistical Society Series B 65 557-574.
Regazzini, E., Lijoi, A. and Prünster, I. (2003). Distributional results for means of normalized random measures with independent increments. The Annals of Statistics 31 560-585.
Reich, B. J. and Fuentes, M. (2007). A multivariate semiparametric Bayesian spatial modeling framework for hurricane surface wind fields. The Annals of Applied Statistics 1 249-264.
Reinhart, A. (2019). pdfetch: Fetch Economic and Financial Time Series Data from Public Sources R package version 0.2.4.

Ren, L., Du, L., Carin, L. and Dunson, D. B. (2011). Logistic stick-breaking process. Journal of Machine Learning Research 12 203-239.
Rigon, T. and Durante, D. (2020). Tractable Bayesian density regression via logit stick-breaking priors. Journal of Statistical Planning and Inference (To appear).
Rodríguez, A., Dunson, D. B. and Gelfand, A. (2008). The nested Dirichlet process. Journal of the American Statistical Association 103 1131-1154.
Rodríguez, A. and Dunson, D. B. (2011). Nonparametric Bayesian models through probit stick-breaking processes. Bayesian Analysis 6 145-178.
Rodríguez, A. and ter Horst, E. (2008). Bayesian dynamic density estimation. Bayesian Anal. 3 339-365.
Sethuraman, J. (1994). A constructive definition of Dirichlet prior. Statistica Sinica 2 639-650.
Sklar, A. (1959). Fonctions de répartition à n dimensions et leurs marges. Publications de l'Institut de Statistique de L'Université d e Paris 8 229-231.
Teh, Y. W., Jordan, M. I., Beal, M. J. and Blei, D. M. (2006). Hierarchical Dirichlet processes. Journal of the American Statistical Association 101 1566-1581.
Tokdar, S. T., Zhu, Y. M. and Ghosh, J. K. (2010). Bayesian density regression with logistic Gaussian process and subspace projection. Bayesian Analysis 5 1-26.
Trippa, L., Müller, P. and Johnson, W. (2011). The multivariate beta process and an extension of the Polya tree model. Biometrika 98 17-34.
Wade, S., Mongelluzzo, S. and Petrone, S. (2011). An enriched conjugate prior for Bayesian nonparametric inference. Bayesian Anal. 6 359-385.
Wang, C. and Rosner, G. L. (2019). A Bayesian nonparametric causal inference model for synthesizing randomized
clinical trial and real-world evidence. Statistics in Medicine 38 2573-2588.
Xu, Z., MacEachern, S. N. and Xu, X. (2015). Modeling non-Gaussian time series with nonparametric Bayesian model. IEEE Transactions on Pattern Analysis and Machine Intelligence 37 372-382.