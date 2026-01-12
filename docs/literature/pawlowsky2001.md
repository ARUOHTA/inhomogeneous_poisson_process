\title{
Geometric approach to statistical analysis on the simplex
}

\author{
V. Pawlowsky-Glahn, J. J. Egozcue
}

\begin{abstract}
The geometric interpretation of the expected value and the variance in real Euclidean space is used as a starting point to introduce metric counterparts on an arbitrary finite dimensional Hilbert space. This approach allows us to define general reasonable properties for estimators of parameters, like metric unbiasedness and minimum metric variance, resulting in a useful tool to better understand the logratio approach to the statistical analysis of compositional data, who's natural sample space is the simplex.
\end{abstract}

Key words: Aitchison geometry, compositional data, Euclidean space, finite dimensional Hilbert space, metric center, metric variance.

\section*{1 \\ Introduction}

The logratio approach to the statistical analysis of compositional data proposed by Aitchison (1982) has been the source of many discussions over the last decades. This is due to the enormous importance compositional data have in practice, as measurements in proportions of some whole, like percentages, ppm, etc. are extremely common in applied sciences. This approach makes it possible to perform classical statistical analysis on transformed data and to back transform the results, which is a clear advantage due to the large amount of methods available for multivariate normally distributed phenomena and the robustness of those. But there has been a certain reluctance in using the new approach by practitioners, which, besides the usual resistance to new theories, is due to the lack of classical properties of backtransformed estimators and models, like unbiasedness and minimum variance.

\footnotetext{
V. Pawlowsky-Glahn ( ® )

Dept. d'Informàtica i Matemàtica Aplicada, Universitat de Girona,
Campus Montilivi - P-1, E-17071 Girona, Spain
e-mail: vera.pawlowsky@udg.es

\section*{J. J. Egozcue}

Dept. de Matemàtica Aplicada III, ETSECCPB, Universitat Politècnica de Catalunya, Barcelona, Spain
}

The authors want to thank H. Burger and another anonymous referee for their comments, which helped to greatly improve the paper.
This research has been supported by the Dirección General de Enseñanza Superior (DGES) of the Spanish Ministry for Education and Culture through the project PB96-0501-C02-01.

In a recent paper, we have given a partial answer to these problems, based on concepts of metric center and metric variance related to the geometric structure of the simplex (Pawlowsky-Glahn and Egozcue, 2001). Here it is shown that the concepts of metric center and metric variance make sense for random vectors with sample space an arbitrary finite dimensional real Hilbert space. Using this approach, it is easy to proof essential properties for statistical inference not only on the simplex, which is the natural sample space for compositional data, but also in other sample spaces. Obviously, the same reasoning can be applied to complex spaces and there is no need to constrain it to elementary concepts and properties. But precisely elementary concepts and properties are useful to convince of the appropriateness and naturality of the definitions, showing that interpretation of real phenomena are much easier if we work on an appropriate sample space using the appropriate measures of central tendency and variability.

Throughout this work, we use the term finite dimensional real Hilbert space, instead of Euclidean space, for spaces with the appropriate structure that are different from \(m\)-dimensional real space \(\mathbb{R}^{m}\). Although mathematically equivalent, we think that speaking about an Euclidean space, whether we refer to \(\mathbb{R}^{m}\), or to its positive orthant \(\mathbb{R}_{+}^{m}\), or to the interval \((0,1)\), or to the simplex \(\mathscr{S}_{c}^{d}\), can be easily misleading in this presentation.

The rationale behind the definitions and properties is related to that of Fréchet (1948) in his paper on random elements of arbitrary nature in a metric space. But Fréchet was primarily concerned with general, non-numerical spaces, while our interest lies in subsets of \(\mathbb{R}^{m}\) with an appropriate structure. Given his approach, Fréchet was naturally interested in probabilistic problems, whereas we emphasize the estimation of parameters.

To illustrate our procedere, let us start recalling basic definitions and properties related to random variables in real space: given a continuous random variable \(X\), the center or expected value is introduced as \(\mathrm{E}[X]=\int_{-\infty}^{+\infty} x \mathrm{~d} F_{X}(x)\), where \(F_{X}(x)\) stands for the distribution function of \(X\), and the variance as \(\operatorname{Var}[X]=\mathrm{E}[(X- \left.\mathrm{E}[X])^{2}\right]\). The geometric interpretation of these concepts is well known, and is often given either as a motivation or as an illustration. Nevertheless, the center can be defined as that value \(\mu\) which minimizes the expected squared Euclidean distance \(\mathrm{E}\left[d_{e}(X, \xi)^{2}\right]\), and the variance \(\sigma^{2}\) can be defined as the expected value of the squared Euclidean distance around \(\mu, \sigma^{2}=\mathrm{E}\left[d_{e}(X, \mu)^{2}\right]\). Obviously, \(\mu=\mathrm{E}[X]\) and \(\sigma^{2}=\operatorname{Var}[X]\). To our understanding, this geometric approach gives its real meaning to the center as a measure of central tendency and to the variance as a measure of dispersion. Fréchet (1948) uses this philosophy to introduce the center and variance of a random vector with support an arbitrary metric space which is not necessarily a vector space. A similar reasoning lead Aitchison (2001) to justify the closed geometric mean as the natural center of a random composition and, later, Pawlowsky-Glahn and Egozcue (2001) to define the concept of metric variance for a random composition. Here, we extend this approach first to an arbitrary finite dimensional real Hilbert space, then we give some simple examples to illustrate its general interest, and finally we particularize on the simplex.

\section*{2 \\ Notation and basic concepts}

Given an \(m\)-dimensional real Hilbert space \(\mathscr{E}\) ( \(m\)-Hilbert space for short), with internal operation \(\oplus\), external operation \(\otimes\), and inner product \(\langle.,\).\(\rangle , denote the\) associated norm by \(\|\cdot\|\) and the associated distance by \(d(.,\).\() . We will use \ominus\) and \(\oslash\) whenever needed for the corresponding operations on the inverses and denote by \(\mathbf{e}\) the neutral element with respect to the internal operation ⊕. This notation

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/66d0948b-d7ad-43aa-b9bc-af1d3d570c13-03.jpg?height=451&width=828&top_left_y=169&top_left_x=204}
\captionsetup{labelformat=empty}
\caption{Fig. 1. Schematic representation of the relationship of a random vector with sample space \(\varepsilon\) and a random vector with sample space \(\mathbb{R}^{m}\) through an isometric transformation \(h\)}
\end{figure}
has been chosen in order to relate and, at the same time, distinguish, operations in \(\mathscr{E}\) from their counterparts in \(\mathbb{R}^{m}\). For convenience, we will identify the scalar field with \((\mathbb{R},+, \cdot)\). Recall that, in an \(m\)-Hilbert space, the distance is invariant with respect to the internal operation, as well as to the external operation,
\[
d(\mathbf{x} \oplus \mathbf{z}, \mathbf{y} \oplus \mathbf{z})=d(\mathbf{x}, \mathbf{y}) ; \quad d(\alpha \otimes \mathbf{x}, \alpha \otimes \mathbf{y})=|\alpha| \cdot d(\mathbf{x}, \mathbf{y})
\]
and that an \(m\)-Hilbert space is always isometric to a real Euclidean space, both having the same dimension \(m\). This property is essential for subsequent developments. Thus, if we denote by \(h\) such an isometry, \(h\) is an isomorphism such that for any \(\mathbf{x}, \mathbf{y} \in \mathscr{E}\),
\[
\langle\mathbf{x}, \mathbf{y}\rangle=\langle h(\mathbf{x}), h(\mathbf{y})\rangle_{e} ; \quad d(\mathbf{x}, \mathbf{y})=d_{e}(h(\mathbf{x}), h(\mathbf{y})),
\]
where \(\langle., .\rangle_{e}\) stands for the Euclidean inner product and \(d_{e}(.,\).\() for the Euclidean\) distance in \(\mathbb{R}^{m}\). A detailed account of these and related properties can be found in (Berberian, 1961).

In the sequel, following a standard approach as i.e. in (Ash, 1972), we consider random vectors \(\mathbf{X}\), defined as measurable functions from a probability space ( \(\Omega\), \(\mathscr{F}, P_{\Omega}\) ) onto a sample space \(\mathscr{E}\). Here \(\Omega\) denotes an arbitrary set, \(\mathscr{F}\) a \(\sigma\)-field of subsets of \(\Omega\), and \(P_{\Omega}\) a probability measure on \(\mathscr{F}\). A consistent definition of a random vector X requires a \(\sigma\)-field of subsets of \(\mathscr{E}\). An easy and natural way to define it consists in taking \(h^{-1}\left(\mathscr{B}\left(\mathbb{R}^{m}\right)\right)\), being \(\mathscr{B}\left(\mathbb{R}^{m}\right)\) the class of Borel sets of \(\mathbb{R}^{m}\) (see Fig. 1 for a schematic illustration). With this definition \(h(\mathbf{X})\) is a measurable function (i.e., a random vector) that goes from ( \(\Omega, \mathscr{F}, P_{\Omega}\) ) to \(\mathbb{R}^{m} . \mathbf{X}\) induces a probability measure \(P\) on the \(\sigma\)-field \(h^{-1}\left(\mathscr{B}\left(\mathbb{R}^{m}\right)\right)\) of \(\mathscr{E}\) and \(h(\mathbf{X})\) induces a probability measure \(P_{e}\) on the \(\sigma\)-field \(\mathscr{B}\left(\mathbb{R}^{m}\right)\) of \(\mathbb{R}^{m}\). Note that for any set \(A \in h^{-1}\left(\mathscr{B}\left(\mathbb{R}^{m}\right)\right)\) we have \(A=h^{-1}(B)\), for some \(B \in \mathscr{B}\left(\mathbb{R}^{m}\right)\). Thus,
\[
\begin{aligned}
P[A] & =P\left[h^{-1}(B)\right]=P\left[\left\{\omega \mid \mathbf{X}(\omega) \in h^{-1}(B)\right\}\right]=P[\{\omega \mid h(\mathbf{X}(\omega)) \in B\}] \\
& =P_{e}[B]=P_{e}[h(A)]
\end{aligned}
\]

With these probability measures we can define expectation in both spaces,
\[
\mathrm{E}[\mathbf{X}]=\int_{x \in \mathscr{E}} x \mathrm{~d} P, \quad \mathrm{E}_{e}[h(\mathbf{X})]=\int_{h(x) \in \mathbb{R}^{m}} h(x) \mathrm{d} P_{e}
\]
satisfying \(h(\mathrm{E}[\mathbf{X}])=\mathrm{E}_{e}[h(\mathbf{X})]\). From now on, we will use the symbol \(\mathrm{E}[\cdot]\) for both expectations. This concept of expectation is extended to functions \(g\) of \(\mathbf{X} \in \mathscr{E}\) and \(g_{e}\) of \(h(\mathbf{X}) \in \mathbb{R}^{m}\) as usual, resulting in
\[
\mathrm{E}[g(\mathbf{X})]=\int_{x \in \mathscr{E}} g(x) \mathrm{d} P, \quad \mathrm{E}_{e}\left[g_{e}(h(\mathbf{X}))\right]=\int_{h(x) \in \mathbb{R}^{m}} g_{e}(h(x)) \mathrm{d} P_{e}
\]

In particular, let \(\xi \in \mathscr{E}\) be a fixed element and consider the function \(d^{2}(\mathbf{X}, \xi)\), where \(d(\ldots, \ldots)\) denotes the distance in \(\mathscr{E}\) (Fig. 1). The expectation of such a function is then
\[
\begin{aligned}
\mathrm{E}\left[d^{2}(\mathbf{X}, \xi)\right] & =\int_{x \in \mathscr{E}} d^{2}(x, \xi) \mathrm{d} P \\
& =\int_{h(x) \in \mathbb{R}^{m}} d_{e}^{2}(h(x), h(\xi)) \mathrm{d} P_{e}=\mathrm{E}\left[d_{e}^{2}(h(\mathbf{X}), h(\xi))\right]
\end{aligned}
\]

Note that the latter expectation can also be defined using the corresponding probability measure induced in \(\mathbb{R}_{+}\). In fact, \(d^{2}(\mathbf{X}, \xi)\) is a univariate random variable with sample space \(\mathbb{R}_{+}\)and its probability can be described by a univariate distribution function.

Now, following the rationale described in the introduction, let us introduce metric counterparts in \(\mathscr{E}\) to usual measures of dispersion and central tendency in real Euclidean space.

\section*{3 \\ Metric center and metric variance}

Definition 1 The dispersion or metric variance around \(\xi \in \mathscr{E}\) is the expected value of the squared distance between \(\mathbf{X}\) and \(\xi: \operatorname{Mvar}[\mathbf{X}, \xi]=\mathrm{E}\left[d^{2}(\mathbf{X}, \xi)\right]\), provided that the last expectation exists.

Note that the metric variance is well defined, given that the squared distance is a real function. Assuming the metric variance of \(\mathbf{X}\) exists, we can introduce now the metric center of its distribution as follows.

Definition 2 The metric center of the distribution of \(\mathbf{X}\) is that element \(\xi \in \mathscr{E}\) which minimizes \(\operatorname{Mvar}[\mathbf{X}, \xi]\). It is called metric center of \(\mathbf{X}\) and is denoted by Mcen [X] for short.

Following our strategy to paraphrase standard statistical concepts, to call metric variance the metric variance around Mcen \([\mathbf{X}]\) and metric standard deviation its square root is only natural. We state this as a definition for easy of reference.

Definition 3 The metric variance around the metric center Mcen \([\mathbf{X}]\) of the distribution of \(\mathbf{X}\) is given by \(\operatorname{Mvar}[\mathbf{X}, \operatorname{Mcen}[\mathbf{X}]]=\mathrm{E}\left[d^{2}(\mathbf{X}, \operatorname{Mcen}[\mathbf{X}])\right]\). It is called metric variance and is denoted by \(\operatorname{Mvar}[\mathbf{X}]\) for short. The square root of the metric variance of a random composition is called metric standard deviation and is denoted by \(\operatorname{Mstd}[\mathbf{X}]\).

Given the existence of an isometry \(h\) between the \(m\)-Hilbert space \(\mathscr{E}\) and real \(m\)-Euclidean space, it is clear that we can transfer directly properties derived from
the geometric structure between them. To do so, the following two propositions are essential.

Proposition 1: If \(h: \mathscr{E} \rightarrow \mathbb{R}^{m}\) is an isometry, then \(\operatorname{Mcen}[\mathbf{X}]=h^{-1}(\mathrm{E}[h(\mathbf{X})])\).

Proof: If \(h: \mathscr{E} \rightarrow \mathbb{R}^{m}\) is an isometry, then it holds that \(\mathrm{E}\left[d^{2}(\mathbf{X}, \xi)\right]=\mathrm{E}\left[d_{e}^{2}(h(\mathbf{X}), h(\xi))\right]\) (see Eq. (2)), and the vector that minimizes \(\mathrm{E}\left[d_{e}^{2}(h(\mathbf{X}), h(\xi))\right]\) is \(h(\xi)=\mathrm{E}[h(\mathbf{X})]\). Consequently, the vector that minimizes \(\mathrm{E}\left[d^{2}(\mathbf{X}, \xi)\right]\) is \(\xi=h^{-1}(\mathrm{E}[h(\mathbf{X})])\).

Proposition 2: If \(h: \mathscr{E} \rightarrow \mathbb{R}^{m}\) is an isometry and \(h(\mathbf{X})=\mathbf{Y} \in \mathbb{R}^{m}\), then
\(\operatorname{Mvar}[\mathbf{X}]=\sum_{i=1}^{m} \operatorname{Var}\left[Y_{i}\right]\).

Proof: Being \(h\) an isometry and taking into account Proposition 1,
\(\operatorname{Mvar}[\mathbf{X}]=\mathrm{E}\left[d^{2}(\mathbf{X}, \operatorname{Mcen}[\mathbf{X}])\right]=\mathrm{E}\left[d_{e}^{2}(h(\mathbf{X}), \mathrm{E}[h(\mathbf{X})])\right]\)
holds. Writing \(h(\mathbf{X})=\mathbf{Y}\) and using the definition of Euclidean distance the equality is obtained. \(\square\)

Note the assumption \(h\) is an isometry is actually too strong for Proposition 1 to hold, although it simplifies the proof. In fact, if \(h\) is an isomorphism, Proposition 1 holds too. Nevertheless, the isometry assumption is necessary for Proposition 2 to hold.

Given the linearity of \(h\) and \(\mathrm{E}[\cdot]\) in Propositions 1 and 2, classical properties of center and variance in real Euclidean space related to linearity and translation invariance, as well as unicity, are transfered to the metric center and metric variance in \(\mathscr{E}\). In particular, denoting by X, Y,... random vectors with sample space \(\mathscr{E}\) and assuming that the expectations involved exist, the following properties, which are stated without proof, hold:

Proposition 3: For all \(\mathbf{b} \in \mathscr{E}\) and all \(\alpha \in \mathbb{R}\),
\(\operatorname{Mcen}[(\alpha \otimes \mathbf{X}) \oplus \mathbf{b}]=(\alpha \otimes \operatorname{Mcen}[\mathbf{X}]) \oplus \mathbf{b}\).

Proposition 4: \(\operatorname{Mcen}[\mathbf{X} \ominus \operatorname{Mcen}[\mathbf{X}]]=\mathbf{e}\).
Proposition 5: Mcen \([\mathbf{X} \oplus \mathbf{Y}]=\operatorname{Mcen}[\mathbf{X}] \oplus \operatorname{Mcen}[\mathbf{Y}]\), and, in general,
\(\operatorname{Mcen}\left[\mathbf{X}_{1} \oplus \mathbf{X}_{2} \oplus \cdots \oplus \mathbf{X}_{N}\right]=\operatorname{Mcen}\left[\mathbf{X}_{1}\right] \oplus \operatorname{Mcen}\left[\mathbf{X}_{2}\right] \oplus \cdots \oplus \operatorname{Mcen}\left[\mathbf{X}_{N}\right]\).
For simplicity, in what follows we will use the notation
\(\oplus_{n=1}^{N} \mathbf{X}_{n}=\mathbf{X}_{1} \oplus \mathbf{X}_{2} \oplus \cdots \oplus \mathbf{X}_{N}\).

Proposition 6: For all \(\mathbf{b} \in \mathscr{E}\) and all \(\alpha \in \mathbb{R}\),
\(\operatorname{Mvar}[(\alpha \otimes \mathbf{X}) \oplus \mathbf{b}]=\alpha^{2} \operatorname{Mvar}[\mathbf{X}]\).

Proposition 7: For all \(\mathbf{b} \in \mathscr{E}\)
\(\operatorname{Mvar}[\mathbf{X}]=\operatorname{Mvar}[\mathbf{X}, \mathbf{b}]-d^{2}(\operatorname{Mcen}[\mathbf{X}], \mathbf{b})\).

Remark 1 Setting in the last proposition \(\mathbf{b}=\mathbf{e}\), the result is analogous to the standard relationship for univariate real random variables, namely
\[
\begin{aligned}
\operatorname{Var}[X] & =\mathrm{E}\left[X^{2}\right]-(\mathrm{E}[X])^{2}=\mathrm{E}\left[(X-0)^{2}\right]-(\mathrm{E}[X]-0)^{2} \\
& =\mathrm{E}\left[d_{e}^{2}(X, 0)\right]-d_{e}^{2}(\mathrm{E}[X], 0)
\end{aligned}
\]

In what follows, it is understood that independence between random vectors which sample space is \(\mathscr{E}\) is defined in an analogous manner to the standard way (Ash, 1972: p. 213).

Proposition 8: For \(\mathrm{X}, \mathrm{Y}\) independent,
\(\operatorname{Mvar}[\mathbf{X} \oplus \mathbf{Y}]=\operatorname{Mvar}[\mathbf{X}]+\operatorname{Mvar}[\mathbf{Y}]\),
and, in general, for \(\mathrm{X}_{1}, \mathrm{X}_{2}, \ldots, \mathrm{X}_{N}\) jointly independent,
\(\operatorname{Mvar}\left[\bigoplus_{n=1}^{N} \mathbf{X}_{n}\right]=\sum_{n=1}^{N} \operatorname{Mvar}\left[\mathbf{X}_{n}\right]\).
A simple property, but of primary importance, is the Chebyshev inequality that holds for the metric center and metric variance. It gives us a further interpretation of the metric variance, or its squared root, as a dispersion measure around the metric center, independently of the distribution of the random vector.

Proposition 9: Chebyshev inequality. For any \(k>0\),
\(\mathrm{P}[d(\mathrm{X}, \operatorname{Mcen}[\mathrm{X}]) \geq k \operatorname{Mstd}[\mathrm{X}]] \leq \frac{1}{k^{2}}\).
Proof: The standard proof of the Chebyshev inequality follows. Define the set
\(A=\{\mathbf{x} \in \mathscr{E}: d(\mathbf{x}, \operatorname{Mcen}[\mathbf{X}]) \geq k \operatorname{Mstd}[\mathbf{X}]\}\).
Then,
\(\operatorname{Mvar}[\mathbf{X}] \geq \mathrm{E}\left[d^{2}(\mathbf{X}, \operatorname{Mcen}[\mathbf{X}]) \mid A\right] \geq k^{2} \operatorname{Mvar}[X] P[\mathbf{X} \in A]\),
from which the statement holds.
Thus, the metric variance, defined as the expected value of the distance to the metric center, allows us a geometric understanding of the measures of central tendency and dispersion of a random vector in a given \(m\)-Hilbert space.

Before we proceed, it is worthwhile to note that Rao (1982) introduced the concept of quadratic entropy, which is equivalent to our metric variance, and that

Cuadras et al. (1997) introduced equivalent concepts with the purpose of classification. The difference is that they do not insist explicitly in the distinction between random vectors which sample space is not real Euclidean space, specially when it comes to estimation problems. The importance of this distinction will be seen in the next section.

\section*{4}

4

\section*{Estimation}

Consider a random vector \(\mathbf{X}\) which sample space is \(\mathscr{E}\), and a random sample of size \(N, \mathbf{X}_{1}, \mathbf{X}_{2}, \ldots, \mathbf{X}_{N}\). Assume that the probability measure of \(\mathbf{X}, P_{\mathscr{E}}\), depends on an unknown parameter, or vector of parameters, \(\theta \in \Theta\), where \(\Theta\) stands for the parameter space. \(\Theta\) is also assumed to have an \(m_{\theta}\)-dimensional real Hilbert space structure. The purpose of this section is to define basic desirable properties of estimators \(\theta\) of \(\theta\). Since these definitions, being parallel to ordinary ones, depend on the metric defined in \(\Theta\), we introduce an ' \(\mathrm{M}_{\theta}\) ' preceeding the concept in order to distinguish them; the subindex \(\theta\) is introduced to make clear that the parameter space and the sample space are not necessarily the same. Therefore, whenever we refer to the metric in the sample space \(\mathscr{E}\), we will use the notation ' \(\mathrm{M}_{\mathscr{E}}\) '. Although tedious, we insist in this notation to state clearly that the essential idea of metric properties is to consider them in the appropriate space. We also point out that we are going to use \(\operatorname{Mcen}_{\theta}\left[\hat{\theta}\left(\mathbf{X}_{1}, \ldots, \mathbf{X}_{N}\right)\right]\), where \(\hat{\theta}\) is a function of the random sample. This metric center implicitly requires the definition of expectation and the corresponding probability measure. If \(P\) is the probability measure induced by \(\mathbf{X}\) in \(\mathscr{E}\), the joint probability measure of the random sample is obtained as the direct product of \(P\) as many times as \(N\). The expectation is then taken as an integral of \(\hat{\theta}\) with respect to this probability measure. Note that the subscript in \(\operatorname{Mcen}_{\theta}[\hat{\theta}]\) is related to \(\Theta\), the sample space of the function \(\hat{\theta}\), and not to the definition of the probability measure.

At this point, the way of reasoning may differ from standard approaches because the \(m_{\theta}\)-Hilbert space structure of \(\Theta\) is normally not assumed to be different from \(\mathbb{R}^{m}\), being \(\Theta \subset \mathbb{R}^{m}\). Whenever \(\Theta\) can be identified with \(\mathbb{R}^{k}\) for some integer \(k\), the following approach coincides with the standard approach. But, if \(\Theta\) is not a linear subspace of \(\mathbb{R}^{m}\), then the present approach claims for a new \(m_{\theta}\)-Hilbert space structure in \(\Theta\) and the equivalence no longer holds.

Definition \(4 \hat{\theta}\) is an \(\mathrm{M}_{\theta}\)-centered or \(\mathrm{M}_{\theta}\)-unbiased estimator of \(\theta\) if, and only if, the metric center of \(\hat{\theta}\), with respect to the metric defined in the parameter space \(\Theta\), is the unknown parameter \(\theta: \operatorname{Mcen}_{\theta}[\hat{\theta}]=\theta\).

Using Proposition 4 we obtain the equivalent property:
Proposition 10: \(\operatorname{Mcen}_{\theta}[\hat{\theta}]=\theta\) is equivalent to \(\operatorname{Mcen}_{\theta}\left[\hat{\theta} \oplus_{\theta} \theta^{-1}\right]=\mathbf{e}_{\theta}\), the neutral element of the internal operation \(\oplus_{\theta}\) on \(\Theta\).

Definition \(5 \operatorname{Mcen}_{\theta}\left[\hat{\theta} \oplus_{\theta} \theta^{-1}\right]\) is called the \(\mathrm{M}_{\theta}\)-bias of \(\hat{\theta}\).
Using now Proposition 3 we obtain:
Proposition 11: \(\operatorname{Mcen}_{\theta}\left[\hat{\theta} \oplus_{\theta} \theta^{-1}\right]=\operatorname{Mcen}_{\theta}[\hat{\theta}] \oplus_{\theta} \theta^{-1}=\operatorname{Mcen}_{\theta}[\hat{\theta}] \ominus_{\theta} \theta\).
In order to compare the \(\mathrm{M}_{\theta}\)-bias from different estimators of the same parameter \(\theta\), the distance to \(\mathbf{e}_{\theta}\) can be used, as they belong to the same space. Thus, the adequate measure to be used is
\(d_{\theta}\left(\operatorname{Mcen}_{\theta}[\hat{\theta}] \oplus_{\theta} \theta^{-1}, \mathbf{e}_{\theta}\right)=d_{\theta}\left(\operatorname{Mcen}_{\theta}[\hat{\theta}], \theta\right)\),
where the equality is derived from Eq. (1).
The mean quadratic error is another important criterion in estimation of standard parameters. The analogous in our case is the following.

Definition \(6 \operatorname{Mvar}_{\theta}[\hat{\theta}, \theta]\) is called the \(\mathrm{M}_{\theta}\)-quadratic error of \(\hat{\theta}\).
This definition of \(\mathrm{M}_{\theta}\)-quadratic error has similar properties to standard quadratic error. Particularly, applying property 7, it is related with the \(\mathrm{M}_{\theta}\)-bias and the metric variance of the estimator in the same way the standard estimators are: quadratic error equals squared bias plus variance of the estimator.

Proposition 12: \(\operatorname{Mvar}_{\theta}[\hat{\theta}, \theta]=\operatorname{Mvar}_{\theta}[\hat{\theta}]+d_{\theta}^{2}\left(\operatorname{Mcen}_{\theta}[\hat{\theta}], \theta\right)\).
After these definitions, general concepts on estimation of standard parameters can be easily extended to metric counterparts (e.g. asymptotically unbiased estimators, consistency in mean quadratic error). In the present context we are specially interested in the following definitions, related to the so called best linear unbiased estimators (BLUE).

Definition 7 Given two estimators \(\hat{\theta}_{1}\) and \(\hat{\theta}_{2}\) of \(\theta \in \Theta, \hat{\theta}_{1}\) is said to be more \(\mathrm{M}_{\theta}\) efficient than \(\hat{\theta}_{2}\) with respect to the distance defined in \(\Theta\) if, and only if, \(\operatorname{Mvar}_{\theta}\left[\hat{\theta}_{1}, \theta\right]<\operatorname{Mvar}_{\theta}\left[\hat{\theta}_{2}, \theta\right]\).

Definition 8 Given a class \(\hat{\Theta} \subset \Theta\) of estimators of \(\theta, \hat{\theta} \in \hat{\Theta}\) is said to be \(\mathrm{M}_{\theta}\)-best within the class \(\hat{\Theta}\) if, and only if, it is \(\mathrm{M}_{\theta}\)-centered and \(\operatorname{Mvar}_{\theta}[\hat{\theta}]<\operatorname{Mvar}_{\theta}\left[\hat{\theta}_{i}\right]\) for all \(\hat{\theta}_{i} \in \hat{\Theta}\); i.e. it is the most \(\mathrm{M}_{\theta}\)-efficient among the \(\mathrm{M}_{\theta}\)-centered estimators in \(\hat{\Theta}\).

Obviously, other standard characterizations of estimators, usual in the context of random variables with support the real line, can be given, simply by substituting the Euclidean distance by the appropriate distance defined in the parameter space, and the expected value by the \(\mathrm{M}_{\theta}\)-center, but it goes beyond the purpose of this paper. Therefore, let us proceed to define an \(\mathrm{M}_{\theta}\)-best linear estimator of a parameter \(\theta\) within the class of linear estimators of \(\theta\), where linear is understood in the following sense:

Definition 9 Given a function \(g(\cdot)\) from \(\mathscr{E}\) onto \(\Theta\),
\(\hat{\theta}=\bigoplus_{n=1}^{N}\left(\alpha_{n} \otimes_{\theta} g\left(\mathbf{X}_{n}\right)\right)\),
is said to be an \(\mathrm{M}_{\theta}\)-linear \(g\)-function of the sample \(\mathbf{X}_{1}, \mathbf{X}_{2}, \ldots, \mathbf{X}_{N}\). Whenever \(g(\cdot)\) is suitable for estimation of \(\theta\), then \(\hat{\theta}\) is said to be an \(\mathrm{M}_{\theta}\)-linear \(g\)-estimator of \(\theta\).

Now, for \(g(\cdot)\) a function from \(\mathscr{E}\) onto \(\Theta\), the following propositions can be set forth.

Proposition 13: If, for any \(n=1, \ldots, N, \operatorname{Mcen}_{\theta}\left[g\left(\mathbf{X}_{n}\right)\right]=\theta\), then
\(\hat{\theta}=\bigoplus_{n=1}^{N}\left(\frac{1}{N} \otimes_{\theta} g\left(\mathbf{X}_{n}\right)\right)\)
is an \(\mathrm{M}_{\theta}\)-linear and \(\mathrm{M}_{\theta}\)-unbiased \(g\)-estimator of \(\theta\).

Proof: \(\hat{\theta}\) is an \(\mathrm{M}_{\theta}\)-linear function by definition. Taking metric centers in the definition of the \(g\)-estimator and using Propositions 5 and 3, the fact that \(\operatorname{Mcen}_{\theta}\left[g\left(\mathbf{X}_{n}\right)\right]=\theta\), and standard properties of vector spaces, the statement holds. Note that, as usual, independence of the sample is not a requirement for the proof. \(\square\)

Proposition 14: Given a function \(g: \mathscr{E} \rightarrow \Theta\) such that \(\operatorname{Mcen}_{\theta}\left[g\left(\mathbf{X}_{n}\right)\right]=\theta\), \(\hat{\theta}=\bigoplus_{n=1}^{N}\left(\frac{1}{N} \otimes_{\theta} g\left(\mathbf{X}_{n}\right)\right)\)
is the \(\mathrm{M}_{\theta}\)-best \(g\)-estimator of \(\theta\) within the class of \(\mathrm{M}_{\theta}\)-linear \(\mathrm{M}_{\theta}\)-unbiased \(g\) estimators of \(\theta\). Moreover,
\(\operatorname{Mvar}_{\theta}[\hat{\theta}]=\frac{\operatorname{Mvar}_{\theta}[g(\mathbf{X})]}{N}\).

Proof: Consider a general \(\mathrm{M}_{\theta}\)-linear \(\mathrm{M}_{\theta}\)-unbiased \(g\)-estimator of \(\theta\)
\(\tilde{\theta}=\bigoplus_{n=1}^{N}\left(\alpha_{n} \otimes_{\theta} g\left(\mathbf{X}_{n}\right)\right)\).
If \(\tilde{\theta}\) is \(\mathrm{M}_{\theta}\)-unbiased, then, using again Propositions 5 and 3 as well as standard operations in a vector space, we obtain \(\sum_{n=1}^{N} \alpha_{n}=1\). To see that it is \(\mathrm{M}_{\theta}\)-best we have to see that the metric variance reaches a minimum when \(\alpha_{n}=1 / N\). Given that by definition of random sample we have \(\operatorname{Mvar}_{\theta}\left[g\left(\mathbf{X}_{n}\right)\right]= \mathrm{Mvar}_{\theta}[g(\mathrm{X})]\) and independence of the sample, the metric variance of \(\theta\) can be expressed, using Propositions 8 and 6, as
\(\operatorname{Mvar}_{\theta}[\tilde{\theta}]=\operatorname{Mvar}_{\theta}[g(\mathbf{X})]\left(\sum_{n=1}^{N} \alpha_{n}^{2}\right)\),
which is minimum when, for \(n=1, \ldots, N, \alpha_{n}=1 / N\). Therefore, for minimum metric variance \(\tilde{\theta}=\hat{\theta}\) and \(\operatorname{Mvar}_{\theta}[\hat{\theta}]=\operatorname{Mvar}_{\theta}[g(\mathbf{X})] / N\). \(\square\)

In Proposition 14 the \(\mathrm{M}_{\theta}\)-linear, \(\mathrm{M}_{\theta}\)-unbiased, \(\mathrm{M}_{\theta}\)-best \(g\)-estimator has been identified, but the function \(g: \mathscr{E} \rightarrow \Theta\) was beforehand given. A natural extension may be to make \(g\) free within a given class of functions. Although a detailed discussion of such a case is out of the scope of this presentation, the proof of Proposition 14 points out that minimization of \(\operatorname{Mvar}_{\theta}[\hat{\theta}]\) can be decomposed into two steps: search for the best \(g: \mathscr{E} \rightarrow \Theta\) and optimization for the \(\alpha_{n}\) 's; and, then, the optimum value of \(\alpha_{n}\) is still \(1 / N\).

Proposition 15: Given a function \(g: \mathscr{E} \rightarrow \Theta\) such that \(\operatorname{Mcen}_{\theta}\left[g\left(\mathbf{X}_{n}\right)\right]=\theta\), \(\hat{\theta}=\bigoplus_{n=1}^{N}\left(\frac{1}{N} \otimes_{\theta} g\left(\mathbf{X}_{n}\right)\right)\) satisfies the weak law of large numbers given by \(P\left[d_{\theta}(\hat{\theta}, \theta) \geq \frac{\operatorname{Mstd}_{\theta}[g(\mathbf{X})]}{\sqrt{\varepsilon N}}\right] \leq \varepsilon\),
for any \(\varepsilon>0\). Consequently, \(\hat{\theta}\) converges in probability to \(\operatorname{Mcen}_{\theta}[g(\mathbf{X})]\) for \(N \rightarrow \infty\).

Proof: This standard result is obtained by applying the Chebyshev inequality stated in Proposition 9 to the random function \(\theta\), which metric center and metric variance are, respectively, \(\theta\) and \(\operatorname{Mvar}_{\theta}[g(\mathbf{X})] / N\). Setting \(1 / h^{2}=\varepsilon\), the desired result is obtained. \(\square\)

Propositions 13-15 clearly establish that \(\hat{\theta}=\bigoplus_{n=1}^{N}\left(\frac{1}{N} \otimes_{\theta} g\left(\mathbf{X}_{n}\right)\right)\) is an \(\mathrm{M}_{\theta}\)-linear \(\mathrm{M}_{\theta}\)-unbiased \(g\)-estimator of \(\theta\), which is \(\mathrm{M}_{\theta}\)-best with respect to the corresponding distance defined in the parameter space \(\Theta\).

Example 1: Univariate random variables with sample space the real line.
For univariate random variables with sample space \(\mathbb{R}=\mathscr{E}\) (i.e. support the real line), it is straightforward to obtain all the standard results. In fact, \(\mathbb{R}\) is a real 1Hilbert space with the usual operations: addition for the internal or Abelian group operation and product for the external operation. The inner product is the usual one (i.e. the product), and the distance is the Euclidean distance. The metric center is then nothing else but the usual expected value, and the best, linear, unbiased estimator associated to the Euclidean distance is the average or arithmetic mean of the sample.

Example 2: Univariate random variables with sample space the positive real line.
A particularly interesting case is that of an univariate random variable \(X\), which sample space (i.e. support) is the positive real line \(\mathbb{R}_{+}=\mathscr{E}\). This sample space is a 1 -Hilbert space with the following operations and definitions. Let \(x, y \in \mathbb{R}_{+}\)and \(\alpha \in \mathbb{R}\); then we have
1. Internal, Abelian group operation: \(x \oplus_{+} y=x \cdot y\).
2. External operation: \(\alpha \otimes_{+} x=x^{\alpha}\).
3. Inner product: \(\langle x, y\rangle_{+}=\ln x \cdot \ln y\).
4. Distance: \(d_{+}(x, y)=|\ln x-\ln y|\).
5. Norm: \(\|x\|_{+}=|\ln x|\).
6. Isometry \(h: \mathbb{R}_{+} \rightarrow \mathbb{R}\) such that \(h(x)=\ln x\), with inverse \(h^{-1}(y)=\exp (y)\).

The metric center in \(\mathbb{R}_{+}\)with this structure is given in Proposition 1 as
\(\gamma=\operatorname{Mcen}_{+}[X]=\exp (\mathrm{E}[\ln X])\).
\(\gamma\) is again a value in \(\mathbb{R}_{+}\)and, therefore, the sample space of \(\gamma\) is \(\mathbb{R}_{+}\), which 1Hilbert space structure has been previously defined. Given a random sample \(X_{1}, \ldots, X_{N}\), to estimate \(\gamma\) we can take in Definition 9 the function \(g=\mathrm{id}\), the identity in \(\mathbb{R}^{+}\). Then, for any \(n\), Mcen \(+\left[g\left(X_{n}\right)\right]=\operatorname{Mcen}_{+}[X]=\gamma\), and Propositions 13 and 14 state that \(\hat{\gamma}=\Pi_{n=1}^{N} X_{n}^{1 / N}\) is the \(\mathrm{M}_{+}\)-best id-estimator of \(\gamma\) within the class of \(\mathrm{M}_{+}\)-linear, \(\mathrm{M}_{+}\)-unbiased id-estimators of \(\gamma\). Note that the estimator obtained is the geometric mean of the sample and the estimated parameter is \(\gamma=\exp (\mathrm{E}[\ln X])\), also known as the theoretical geometric mean of a random variable. These facts recall us the standard treatment of lognormal variates and state that, in terms of the geometric structure of \(\mathbb{R}_{+}\), the natural measure of central tendency is the theoretical geometric mean and the best estimator is the geometric mean of the sample.

Note that this reasoning can be applied to any measure of difference, which sample space is by definition \(\mathbb{R}_{+}\). As a result, we obtain an estimator for the
metric variance which is different from the one obtained by the method of moments in real Euclidean space.

Example 3: Univariate random variable with sample space \(I=(0,1)\).
Let \(X\) be an univariate random variable which sample space is \(I=(0,1)\). This sample space is a 1 -Hilbert space with the following operations and definitions. Let \(x, y \in I\) and \(\alpha \in \mathbb{R}\), then we have:
1. Internal, Abelian group operation:
\[
x \oplus_{I} y=\frac{x y}{(1-x)(1-y)+x y}
\]
2. External operation:
\[
\alpha \otimes_{I} x=\frac{x^{\alpha}}{(1-x)^{\alpha}+x^{\alpha}}
\]
3. Inner product:
\[
\langle x, y\rangle_{I}=\ln \frac{x}{1-x} \cdot \ln \frac{y}{1-y}
\]
4. Distance:
\[
d_{I}(x, y)=\left|\ln \frac{x(1-y)}{y(1-x)}\right|
\]
5. Norm:
\[
\|x\|_{I}=\left|\ln \frac{x}{1-x}\right|
\]
6. Isometry (logit transformation):
\(h: I \rightarrow \mathbb{R}\) such that \(h(x)=\ln \frac{x}{1-x}, \quad\) with inverse \(h^{-1}(y)=\frac{\exp (y)}{1+\exp (y)}\).

According to Proposition 1, the metric center in \(I\) is given by
\(\gamma=\operatorname{Mcen}_{I}[X]=\frac{\exp (\mathrm{E}[\ln (X /(1-X))])}{1+\exp (\mathrm{E}[\ln (X /(1-X))])}\),
and \(\gamma\) is again a value in \(I\). Thus, the parameter space of \(\gamma\) is \(I\), which 1-Hilbert space structure has been just defined. Given a random sample \(X_{1}, \ldots, X_{N}\), to estimate \(\gamma\) we can take in Definition 9 the function \(g=\mathrm{id}\), the identity in \(I\). Then, for any \(n, \operatorname{Mcen}_{I}\left[g\left(X_{n}\right)\right]=\operatorname{Mcen}_{I}[X]=\gamma\), and Propositions 13 and 14 state that
\(\hat{\gamma}=\frac{\prod_{n=1}^{N} X_{n}^{1 / N}}{\prod_{n=1}^{N}\left(1-X_{n}\right)^{1 / N}+\prod_{n=1}^{N} X_{n}^{1 / N}}\)
is the \(\mathrm{M}_{I}\)-best id-estimator of \(\gamma\) within the class of \(\mathrm{M}_{I}\)-linear, \(\mathrm{M}_{I}\)-unbiased idestimators of \(\gamma\).

\section*{5}

\section*{Estimation on the simplex}

Recall that \(\mathbf{x}=\left(x_{1}, \ldots, x_{d}\right)^{\prime}\) is by definition a \(d\)-part composition if, and only if, all its components are strictly positive real numbers and their sum is a constant \(c\). The constant \(c\) is 1 if measurements are made in parts per unit, or 100 if measurements are made in percent. The sample space of \(d\)-part compositional data with constant sum \(c\) is thus the simplex
\(\mathscr{S}_{c}^{d}=\left\{\mathbf{x}=\left(x_{1}, \ldots, x_{d}\right)^{\prime} \mid x_{i}>0, i=1, \ldots, d ; \sum_{i=1}^{d} x_{i}=c\right\}\),
where the prime stands for transpose. Although mathematically less comfortable, we keep the constant \(c\) in the definition and in the notation, to avoid confusion arising from the fact that in geology it is more common to use \(c=100\) than \(c=1\). But, to simplify the mathematical developments, we include the constant in the closure operation as stated below.

Basic operations on the simplex have been introduced by Aitchison (1986). They are the perturbation operation, defined for any two vectors \(\mathbf{x}, \mathbf{y} \in \mathscr{S}_{c}^{d}\) as
\(\mathbf{x} \circ \mathbf{y}=\mathscr{C}\left(x_{1} y_{1}, \ldots, x_{d} y_{d}\right)^{\prime}\),
and the power transformation, defined for a vector \(\mathbf{x} \in \mathscr{S}_{c}^{d}\) and a scalar \(\alpha \in \mathbb{R}\) as
\(\alpha \diamond \mathbf{x}=\mathscr{C}\left(x_{1}^{\alpha}, \ldots, x_{d}^{\alpha}\right)^{\prime}\),
where the \(\mathscr{C}\) denotes the closure operation defined for a vector \(\mathbf{z}=\left(z_{1}, \ldots, z_{d}\right)^{\prime}\) as
\(\mathscr{C}(\mathbf{z})=\mathscr{C}\left(\begin{array}{c}z_{1} \\ z_{2} \\ \vdots \\ z_{d}\end{array}\right)=\left(\begin{array}{c}\frac{c \cdot z_{1}}{z_{1}+z_{2} \cdot \cdot \cdots+z_{d}} \\ \frac{z_{1}+z_{2}+\cdots+z_{d}}{\vdots} \\ \vdots \\ \frac{c \cdot z_{d}}{z_{1}+z_{2}+\cdots+z_{d}}\end{array}\right)\).
Perturbation and power transformation induce a vector space structure in the simplex. Then, to obtain a \((d-1)\)-Hilbert space structure on \(\mathscr{S}_{c}^{d}\), the following inner product and associated norm and distance can be used (Aitchison, 2001; Pawlowsky-Glahn and Egozcue, 2001):
\(\langle\mathbf{x}, \mathbf{y}\rangle_{a}=\frac{1}{d} \sum_{i<j} \ln \frac{x_{i}}{x_{j}} \ln \frac{y_{i}}{y_{j}} ;\)
\(\|\mathbf{x}\|_{a}=\sqrt{\frac{1}{d} \sum_{i<j}\left(\ln \frac{x_{i}}{x_{j}}\right)^{2}} ;\)
\(d_{a}(\mathbf{x}, \mathbf{y})=\sqrt{\frac{1}{d} \sum_{i<j}\left(\ln \frac{x_{i}}{x_{j}}-\ln \frac{y_{i}}{y_{j}}\right)^{2}}\).

To refer to the properties of \(\left(\mathscr{S}_{c}^{d}, \circ, \diamond\right)\) as a \((d-1)\)-Hilbert space, we shall talk globally about the Aitchison geometry on the simplex, and in particular about the Aitchison distance, norm and inner product.
\(\mathscr{S}_{c}^{d}\) is a \((d-1)\)-Hilbert space and thus there exists an isometry between \(\mathscr{S}_{c}^{d}\) and \(\mathbb{R}^{d-1}\) which could be used to analize results from previous sections in the particular case of random compositions. Nevertheless, to simplify presentation, we will use the clr transformation defined by Aitchison (1986), which is an isometry between the simplex with the Aitchison geometry and the \((d-1)\)-hyperplane going through the origin and parallel to the simplex in \(\mathbb{R}^{d}\) and equipped with the usual Euclidean geometry in \(\mathbb{R}^{d-1}\) projected from the real Euclidean space \(\mathbb{R}^{d}\). Recall that the clr transformation is defined as
\(\operatorname{clr}(\mathbf{x})=\left(\ln \frac{x_{1}}{g(\mathbf{x})}, \ln \frac{x_{2}}{g(\mathbf{x})}, \ldots, \ln \frac{x_{d}}{g(\mathbf{x})}\right)\),
where \(g(\mathbf{x})=\left(\prod_{i=1}^{d} x_{i}\right)^{1 / d}\), the geometric mean of \(\mathbf{x}\). The inverse is obtained by taking first exponentials and then applying the closure operation, which cancels out multiplicative constants.

With these elements at hand we can proceed to analyze the metric variance and metric center of the distribution of a random vector \(\mathbf{X}\) with sample space \(\mathscr{S}_{c}^{d}\).

Proposition 16: The metric center \(\operatorname{Mcen}_{a}[\mathbf{X}]\) of \(\mathbf{X}\) is the center or closed geometric mean of \(\mathbf{X}\),
\(\operatorname{Mcen}_{a}[\mathbf{X}]=\mathscr{C}\left(\exp \left\{\mathrm{E}\left[\ln \left(X_{1}\right)\right]\right\}, \ldots, \exp \left\{\mathrm{E}\left[\ln \left(X_{d}\right)\right]\right\}\right)^{\prime}\).
This result is actually the original definition of center of a random composition given by Aitchison (1997). It is obtained by applying Proposition 1 with the clr transformation.

Proposition 17: The metric variance can be expressed as,
\(\operatorname{Mvar}_{a}[\mathbf{X}]=\frac{1}{d} \sum_{i<j} \operatorname{Var}\left[\ln \frac{X_{i}}{X_{j}}\right]=\sum_{i=1}^{d} \operatorname{Var}\left[\ln \frac{X_{i}}{g(\mathbf{X})}\right]\).
The first equality states that the metric variance with respect to the Aitchison distance is identical to the total variance defined by Aitchison (1997). It is derived directly from the definition of metric variance and of the Aitchison distance. The second equality is obtained using Proposition 2 and the clr transformation.

Note that Proposition 16 implies that, writing \(\operatorname{Mcen}_{a}[\mathbf{X}]=\gamma\), for \(i, j=1, \ldots, d\),
\(\mathrm{E}\left[\ln \frac{X_{i}}{X_{j}}\right]=\ln \frac{\gamma_{i}}{\gamma_{j}}\).
As a result of these statements, we can say that the closed geometric mean of random compositions minimizes the metric variance on the simplex with respect to the Aitchison distance. We can also say that the total variance defined by Aitchison (1997) is an appropriate measure of compositional variability within
the simplex, as it coincides with the expected value of the squared Aitchison distance to the metric center of the distribution.

Looking at other properties of random compositions derived from propositions stated in Sect. 3, we see that perturbation of a random composition affects the metric center in that it leads to a perturbed metric center (Proposition 3), whereas it has no effect on the metric variance (Proposition 6). As a consequence, we can center the random composition by perturbing it with the inverse of the metric center (proposition 4), thus giving theoretical support to the approach presented in (Buccianti et al., 1999; Martín-Fernández et al., 1999; Eynatten et al., 2001). Furthermore, the metric center is linear with respect to perturbation (Proposition 5), whereas this property holds for the metric variance only in case of independence of the random compositions involved (Proposition 8).

Proposition 7 applied to random compositions on the simplex tells us, that the metric variance can be expressed as the metric variance around an arbitrary point b in the simplex minus the squared Aitchison distance between the metric center and the same point. Substituting \(\mathbf{b}\) by the baricenter or neutral element of perturbation, e , gives the following result:
\(\operatorname{Mvar}_{a}[\mathbf{X}]=\mathrm{E}\left[d_{a}^{2}(\mathbf{X}, \mathbf{e})\right]-d_{a}^{2}\left(\operatorname{Mcen}_{a}[\mathbf{X}], \mathbf{e}\right)\),
suggesting an analogy to central and non-central moments in real space.
Another interesting feature is related to the power transformation. The power transformation of a random composition multiplies the metric center (Proposition 3) and its square multiplies the metric variance (Proposition 6). Thus, we can introduce an equivalent concept to standardized random vectors by using perturbation with the inverse of the metric center and power transformation with the inverse of the metric standard deviation to obtain random compositions centered at the baricenter \(\mathbf{e}\) and with unit variance:
\(U=\frac{1}{\operatorname{Mstd}_{a}[\mathbf{X}]} \diamond\left(\mathbf{X} \circ\left(\operatorname{Mcen}_{a}[\mathbf{X}]\right)^{-1}\right)\).
Finally, the Chebyshev inequality stated in Proposition 9 gives us a way to obtain regions within the simplex where we have a probability smaller or equal to \(1 / k^{2}\) that the random composition is at a distance from the metric center larger then \(k\) times the metric standard deviation.

Concerning the estimation of the center, we can say that, taking in Proposition 14 the identity function \(g\left(\mathbf{X}_{n}\right)=\operatorname{id}\left(\mathbf{X}_{n}\right)=\mathbf{X}_{n}\) we obtain that
\(\overline{\mathbf{X}}_{a}=\bigcup_{n=1}^{N}\left(\frac{1}{N} \diamond \mathbf{X}_{n}\right)\)
is the \(\mathrm{M}_{a}\)-best id-estimator of \(\operatorname{Mcen}_{a}[\mathbf{X}]\) within the class of \(\mathrm{M}_{a}\)-linear \(\mathrm{M}_{a}\)-unbiased id-estimators of \(\operatorname{Mcen}_{a}[\mathbf{X}]\). Moreover,
\(\operatorname{Mvar}_{a}\left[\overline{\mathbf{X}}_{a}\right]=\frac{\operatorname{Mvar}_{a}[\mathbf{X}]}{N}\).
These are only the basic properties of the metric center, but it is clear that the same rationale would lead us to transfer whatsoever properties based on Euclidean reasoning from real space into the simplex. This approach assures us that we will obtain properties of optimality in the simplex, completely equivalent to those in real space.

\section*{6 \\ Conclusions}

The existence of an appropriate \(m\)-Hilbert space structure in the simplex suggests a different approach to the statistical analysis of compositional data based on geometric reasoning. Based on this approach, which is completely parallel to the usual one in Euclidean space, it is straightforward to define reasonable properties for estimators of compositional parameters. It assures us that we will obtain properties of optimality, completely equivalent to those in real space, in the simplex. In particular, the closed geometric mean is a linear, unbiased estimator that minimizes the metric variance with respect to the Aitchison geometry on the simplex.

But even more important is, that the same methodology allows us to study properties of probability measures on any sample space with an appropriate finite dimensional real Hilbert space structure, thus opening up a geometric approach to the study of statistical properties in general. Furthermore, it has been shown that this approach is also valid for estimators of unknown parameters of a probability measure and/or characteristics of a random vector, like the metric center. As a consequence, care has to be applied in analyzing the structure of the parameter space to assure the appropriateness of applied methods.

\section*{References}

Aitchison J (1982) The statistical analysis of compositional data (with discussion). J. Royal Stat. Soc., Series B (Statistical Methodology) 44(2): 139-177
Aitchison J (1986) The Statistical Analysis of Compositional Data. Monographs on Statistics and Applied Probability, 416 p. Chapman \& Hall Ltd., London
Aitchison J (1997) The one-hour course in compositional data analysis or compositional data analysis is simple. In: Pawlowsky-Glahn V (ed.) Proceedings of IAMG'97 - The Third Annual Conference of the International Association for Mathematical Geology, vols. I, II and addendum, pp. 3-35. International Center for Numerical Methods in Engineering (CIMNE), Barcelona (E)
Aitchison J (2001) Simplicial inference. In: Viana M, Richards D (eds) Algebraic Methods in Statistics, Contemporary Mathematics Series. American Mathematical Society, New York, NY (in press)
Ash RB (1972) Real Analysis and Probability, 476 p. Academic Press, Inc., New York, NY Berberian SK (1961) Introduction to Hilbert Space, 206 p. Oxford University Press, New York, NY
Buccianti A, Pawlowsky-Glahn V, Barceló-Vidal C, Jarauta-Bragulat E (1999) Visualization and modeling of natural trends in ternary diagrams: a geochemical case study. In: Lippard SJ, Næss A, Sinding-Larsen R (eds) Proceedings of IAMG '99 - The Fifth Annual Conference of the International Association for Mathematical Geology, vols. I and II, pp. 139-144. Tapir, Trondheim (N)
Cuadras CM, Fortiana J, Oliva F (1997) The proximity of an individual to a population with applications in discriminant analysis. J. Classification 14: 117-136
Eynatten Hv, Pawlowsky-Glahn V, Egozcue JJ (2001) Understanding perturbation on the simplex: a simple method to better visualise and interpret compositional data in ternary diagrams. Accepted for publication in Mathematical Geology
Fréchet \(\mathbf{M}\) (1948) Les éléments aléatoires de nature quelconque dans une espace distancié. Annales de L'Institut Henri Poincaré 10(4): 215-308
Martín-Fernández JA, Bren M, Barceló-Vidal C, Pawlowsky-Glahn V (1999) A measure of difference for compositional data based on measures of divergence. In: Lippard SJ, Nss A, Sinding-Larsen R (eds) Proceedings of IAMG '99 - The Fifth Annual Conference of the International Association for Mathematical Geology, vols. I and II, pp. 211-216. Tapir, Trondheim (N)
Pawlowsky-Glahn V, Egozcue JJ (2001) About BLU estimators and compositional data. Accepted for publication in Mathematical Geology
Rao CR (1982) Diversity: its measurement, decomposition, apportionment and analysis. Sankhya, Indian J. Stat., Series A 44: 1-22