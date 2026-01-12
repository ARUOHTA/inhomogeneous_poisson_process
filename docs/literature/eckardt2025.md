\title{
On spatial point processes with composition-valued marks
}

\author{
Matthias Eckardt \({ }^{a}\), Mari Myllymäki \({ }^{b}\) and Sonja Greven \({ }^{a}\) \\ \({ }^{\mathrm{a}}\) Chair of Statistics, Humboldt-Universität zu Berlin, Berlin, Germany \\ \({ }^{\mathrm{b}}\) Natural Resources Institute Finland (Luke), Helsinki, Finland
}

\begin{abstract}
Methods for marked spatial point processes with scalar marks have seen extensive development in recent years. While the impressive progress in data collection and storage capacities has yielded an immense increase in spatial point process data with highly challenging non-scalar marks, methods for their analysis are not equally well developed. In particular, there are no methods for compositionvalued marks, i.e. vector-valued marks with a sum-to-constant constrain (typically 1 or 100 ). Prompted by the need for a suitable methodological framework, we extend existing methods to spatial point processes with composition-valued marks and adapt common mark characteristics to this context. The proposed methods are applied to analyse spatial correlations in data on tree crown-to-base and business sector compositions.
\end{abstract}

Keywords: business sector composition; compositional data analysis; crown-to-base ratios; mark correlation function; marked spatial point processes; mark variogram

\section*{1 Introduction}

In recent years, there has been substantial interest in the analysis of marked spatial point processes. This type of data involves random spatial positions of \(n\) events with marks providing specific additional details about each event. Although spatial point process theory is mature, and various summaries for marked points have been created (see Stoyan and Stoyan, 1994, Illian et al., 2008, Chiu et al., 2013), there is still a gap in addressing non-scalar marks. Particularly, spatial point processes with compositional marks-where marks are made up of \(D\) components that together sum to a constant-are not fully examined. Such instances include pest infestation ratios in wood samples, proportions of business sectors locally, and parts of tree biomass. Here,
marks present relative information with interdependent \(D\) components; an increase in one results in a decrease in others, due to the fixed sum. The constrained nature of this data demands analysis techniques beyond traditional Euclidean methods. Individual component analysis within the current framework disregards these constraints, risking biased or spurious outcomes (Chayes, 1960; Pawlowsky, 1984). Thus, a methodological framework for analyzing such composition-valued marks remains crucial in modern applications. This paper seeks to address this by introducing a new class of spatial point processes with composition-valued marks and expanding the toolbox for mark summary characteristics. By applying compositional data analysis principles Aitchison, 1986) to marked spatial point processes, this work goes beyond existing methods, establishing novel summary characteristics for complex mark scenarios. To our knowledge, composition-valued marked point processes are examined here for the first time.

Various summary characteristics are standard tools to characterize both the properties of points and the interrelations between marks (and marks and points) (Stoyan and Stoyan, 1994; Illian et al., 2008; Chiu et al., 2013). Of all established methods, functional summary characteristics, in which the characteristic is defined as a function of the distance between points \(r\), play a particularly important role. They are used to investigate the marked point pattern, to determine a suitable model and to evaluate the goodness-of-fit of such models (e.g. Illian et al., 2008; Myllymäki et al., 2017).

Despite the growing availability of highly demanding spatial marked point process scenarios, the literature focused, so far, almost exclusively on the analysis of scalarvalued attributes (see Baddeley, 2010, for a general treatment). For integer-valued, i.e. qualitative, marks, where the points are assigned to exactly one out of \(k \geqslant 2\) distinct types, so-called cross- and dot-type extensions of classic nearest-neighbour and pairwise distance based summary characteristics exist (Lotwick and Silverman, 1982; Harkness and Isham, 1983; Diggle, 1986; Van Lieshout and Baddeley, 1999). For real-valued marks, the aim has mainly been to quantify the average association or variation among the mark values for pairs of distinct points at a distance \(r\) apart. This includes characteristics which highlight the spatial mark-to-mark (Stoyan, 1984b, 1987; Stoyan and Stoyan, 1994; Stoyan and Wälder, 2000; Schlather, 2001) and point-to-mark associations (Schlather et al., 2004; Guan, 2006; Guan and Afshartous, 2007; Ho and Stoyan, 2008, Zhang and Zhuang, 2014), (real-valued) mark weighted versions of classic point summary characteristics (Penttinen et al., 1992; van Lieshout, 2006) and frequency domain approaches (Eckardt and Mateu, 2019a|b). Further, Stoyan (1987) considered the analysis of two distinct real-valued marks, and Penttinen et al. (1992), Wiegand and Moloney (2013) and Eckardt and Mateu (2019ab) discussed methods for mixtures of integer- and real-valued marks. Further, function-valued marks Comas et al., 2008, 2011, 2013, Ghorbani et al., 2021) and generalizations to multivariate function-valued marks (Eckardt et al., 2024) have been investigated (see Eckardt and Moradi, 2024a b, for a general review). However, there are no methods available for the joint analysis of (constrained) vector-valued marks.

We extend the considered mark quantities. This paper aims to characterise the spatial association for pairs of points with \(D\)-part composition-valued marks. Instead of investigating the component parts separately, the underlying idea is to treat the observed marks as part of a total, and to transfer the main principles of compositional data analysis to marked spatial point processes. Composition-valued geostatistical and areal data has already received much attention, including various structural analysis and kriging approaches (Pawlowsky, 1986; Pawlowsky-Glahn and Ricardo, 2004; Tolosana Delgado, 2006) and areal regression models (Leininger et al., 2013; Huang et al., 2021). However, no extensions to marked point processes with composition-valued marks exist so far. As our main contribution, we introduce a general framework of compositionvalued marks and define appropriate mark characteristics for them.

The remainder of the paper is structured as follows. Section 2 summarizes the present methodological toolbox for real-valued marks we build on. Section 3 introduces composition-valued marks, suitable mark transformations and first-order tools (Section 3.1, 3.3), different mark summary characteristics (Sections 3.4, 3.5), extensions to composition-valued marks with total information (Section 3.6) and introduces estimation for the new methods (Section 3.7). Section 4 considers testing of the basic random labeling hypothesis for composition-valued marks. Applications of the proposed methods to a Finnish tree pattern and Spanish business sector data are provided in Section 5. The paper closes with a conclusion in Section 6.

\section*{2 Summary characteristics for real-valued point attributes}

As compositions can be seen as a constrained vector of real-valued marks, the following discussion mostly restricts to methods for the real-valued marks scenario, which are to be extended for composition-valued marks in Section 3.

\subsection*{2.1 Preliminaries}

Let \(X=\left\{x_{i}, m\left(x_{i}\right)\right\}_{i=1, \ldots, n}\) denote a marked spatial point process on \(\mathbb{R}^{2} \times \mathbb{M}\) with points \(x_{i}\) in a two-dimensional Euclidean space and associated marks \(m\left(x_{i}\right)\) living in a mark space \(\mathbb{M}\). The observed point pattern and the related unmarked point process will be denoted by \(\mathbf{x}\) and \(\breve{X}\), respectively. In what follows, \(X\) is assumed to be simple, where simplicity means that multiple coincident points do not occur. In general, we assume IM to be a Polish space equipped with a \(\sigma\)-algebra \(\mathcal{M}\) and an appropriate reference measure \(\varpi\). The mark distribution will be denoted by \(M\). The Borel \(\sigma\)-algebra of \(\mathbb{R}^{2}\) is denoted by \(\mathcal{B}\). For \(X\), the expected number \(N(\cdot)\) of points in \(B \in \mathcal{B}\) with marks in \(L \in \mathcal{M}\) is
\[
\Lambda(B \times L)=\mathbb{E}(N(B \times L))
\]

Further, \(X\) is said to be stationary if \(\left\{x_{i}, m\left(x_{i}\right)\right\} \stackrel{d}{=}\left\{x_{i}+s, m\left(x_{i}\right)\right\}\) for all \(s \in \mathbb{R}^{2}\) and isotropic if \(\left\{x_{i}, m\left(x_{i}\right)\right\} \stackrel{d}{=}\left\{\mathfrak{r} x_{i}, m\left(x_{i}\right)\right\}\) for any rotation \(\mathfrak{r} \in \mathbb{R}^{2}, \stackrel{d}{=}\) denoting equality in distribution. If \(X\) is stationary and isotropic, \(X\) is called motion-invariant. For stationary \(X, \Lambda(B \times L)\) simplifies to
\[
\Lambda(B \times L)=\lambda \nu(B) M(L)
\]
where \(\lambda\) is the intensity of \(\breve{X}\), and \(\nu(\cdot)\) is the Lebesgue measure. If \(\mathbb{M}=\mathbb{R}, M\) is completely determined by the mark distribution function \(F_{M}(m)=M((-\infty, m])\) for \(-\infty \leqslant m \leqslant \infty\), where
\[
F_{M}(m)=\int_{-\infty}^{m} f_{M}(\tilde{m}) \mathrm{d} \tilde{m}
\]
with \(f_{M}(\cdot)\) denoting the mark density function if it exists. Here, \(f_{M}(L):=\int_{L} f_{M}(l) \mathrm{d} l\) can be interpreted as the probability that at an arbitrarily chosen point the mark is in \(L\). Useful quantities of the mark distribution are the mean mark \(\mu_{M}\) and the mark variance \(\sigma_{M}^{2}\). Finally, the second-order factorial moment measure of \(X\) plays an important role in the second-order summary characteristics. It is defined as
\[
\begin{aligned}
\alpha^{(2)}\left(B_{1} \times B_{2}\right)= & \mathbb{E}\left[\sum_{x_{1}, x_{2} \in \breve{X}}^{\neq} \mathbb{1}_{B_{1}}\left(x_{1}\right) \mathbb{1}_{B_{2}}\left(x_{2}\right)\right] \\
& =\int_{B_{1}} \int_{B_{2}} \varrho^{(2)}\left(x_{1}, x_{2}\right) \mathrm{d} x_{1} \mathrm{~d} x_{2}
\end{aligned}
\]
where the symbol \(\sum^{\neq}\)denotes the sum over distinct pairs of points, \(\mathbb{1}_{B}(x)\) an indicator function of whether \(x \in B\) and \(\varrho^{(2)}\) is the second-order product density. Heuristically, for any \(x_{1}, x_{2} \in \mathbb{R}^{2}, \varrho^{(2)}\left(x_{1}, x_{2}\right) \mathrm{d} x_{1} \mathrm{~d} x_{2}\) can be interpreted as the probability of observing exactly one point in each of the infinitesimal areas \(\mathrm{d} x_{1}\) and \(\mathrm{d} x_{2}\).

\subsection*{2.2 Functional summary characteristics for real-valued marks}

Within the last decades various mark characteristics were introduced. These characteristics either describe the average pairwise variation or association of the marks at two point locations as a function of the interpoint distance \(r\). Prominent cases include Stoyan's mark correlation function (Stoyan and Stoyan, 1994), the mark variogram (Cressie, 1993), the mark covariance function (Stoyan, 1984a), Isham's mark correlation function (Isham, 1985), and Schlather's (Schlather et al., 2004) and Shimatani's (Shimatani, 2002) I functions. All these characteristics are defined exclusively for stationary point processes. These characteristics are conditional quantities (in a Palm sense, see Chiu et al., 2013), i.e. conditional on that there are indeed points at location \(\circ\) and \(\mathbf{r}\) in \(X\). They are commonly constructed by taking the conditional expectation \(\mathbb{E}_{\circ, r}\) of a so-called test function \(\mathfrak{t}_{f}: \mathbb{M} \times \mathbb{M} \rightarrow \mathbb{R}^{+}\), which takes the marks \(m(\circ)\) and \(m(\mathbf{r})\) at the origin \(\circ\) and any alternative points at distance \(\|\mathbf{r}\|=r>0\) from \(\circ\) as its
arguments (Penttinen and Stoyan, 1989). Without imposing any invariance assumptions, the expectation is formally defined with respect to the joint distribution of the marks \(m\left(x_{1}\right)\) and \(m\left(x_{2}\right)\) for any two points \(x_{1}, x_{2} \in \mathbb{R}^{2}\) under the condition that there are indeed points at locations \(x_{1}\) and \(x_{2}\) in \(\breve{X}\), i.e. the so-called two-point mark distribution \(M_{x_{1}, x_{2}}\left(\mathrm{~d} m\left(x_{1}\right) \mathrm{d} m\left(x_{2}\right)\right)\). For motion-invariant \(X, M_{x_{1}, x_{2}}\) depends on the points only through the Euclidean distance \(\left\|x_{1}-x_{2}\right\|=r\) and can be written as \(M_{r}\). For \(L_{1}, L_{2} \in \mathcal{M}, M_{r}\left(L_{1} \times L_{2}\right)\) corresponds to the probability of having \(m(\circ) \in L_{1}\) and \(m(\mathbf{r}) \in L_{2}\) under the condition that there are indeed points at \(\circ\) and \(\mathbf{r}\) at a distance \(r\) from each other in \(\breve{X}\). Note that \(M_{r}\left(L_{1} \times L_{2}\right)=\varrho^{(2)}\left(r, L_{1}, L_{2}\right) / \varrho^{(2)}(r)\) simplifies to \(\left(\varrho^{(2)}(r) M\left(L_{1}\right) M\left(L_{2}\right)\right) / \varrho^{(2)}(r)=M\left(L_{1}\right) M\left(L_{2}\right)\) under independent marks for all \(L_{1}\) and \(L_{2}\) in \(\mathcal{M}\), where \(\varrho^{(2)}\left(r, L_{1}, L_{2}\right)\) is the second-order product density of
\[
\alpha_{m}^{(2)}\left(B_{1} \times L_{1} \times B_{2} \times L_{2}\right)=\int_{B_{1} \times B_{2}} M_{r}\left(L_{1} \times L_{2}\right) \alpha^{(2)}\left\{d\left(x_{1}, x_{2}\right)\right\}
\]
and \(\varrho^{(2)}(r)\) is the second-order product density as in (1) for \(x_{1}, x_{2}\) at distance \(r\) (Penttinen and Stoyan, 1989).

Using the above notation, functional mark characteristics are constructed as follows. Denoting by \(\nabla_{\mathfrak{t}_{f}}(r)=\mathbb{E}_{\circ, r}\left[\mathfrak{t}_{f}(m(\circ), m(\mathbf{r}))\right]\) and writing \(\nabla_{\mathfrak{t}_{f}}=\nabla_{\mathfrak{t}_{f}}(\infty)\) for the expected value of the chosen test function \(\mathfrak{t}_{f}\) at very large distances, i.e. when the marks are expected to be independent, the specific mark characteristic itself is determined by the specific choice of the test function \(\mathfrak{t}_{f}\) (see Illian et al., 2008). Formally, \(\nabla_{\mathfrak{t}_{f}}(r)\) can be expressed as the ratio of two product density functions,
\[
\nabla_{\mathfrak{t}_{f}}(r)=\varrho_{\mathfrak{t}_{f}}^{(2)}(r) / \varrho^{(2)}(r)
\]
i.e. the densities of the \(\mathfrak{t}_{f}\)-factorial moment measure \(\alpha_{\mathfrak{t}_{f}}^{(2)}\left(B_{1} \times B_{2}\right)\),
\[
\alpha_{\mathfrak{t}_{f}}^{(2)}\left(B_{1} \times B_{2}\right)=\mathbb{E}\left[\sum_{\left(x_{2}, m\left(x_{2}\right)\right) \in X}^{\neq} \mathfrak{t}_{f}\left(m\left(x_{1}\right), m\left(x_{2}\right)\right) \mathbb{1}_{B_{1}}\left(x_{1}\right) \mathbb{1}_{B_{2}}\left(x_{2}\right)\right]
\]
and of the factorial moment measures \(\alpha^{(2)}\left(B_{1} \times B_{2}\right)\) of (1), respectively, where \(B_{1}, B_{2}\) are sets in \(\mathcal{B}\). When the distance \(r\) tends to infinity, the marks are assumed to be independent, \(\nabla_{\mathfrak{t}_{f}}(r)\) becomes independent of the points and simplifies to
\[
\nabla_{\mathfrak{t}_{f}}=\int_{\mathbb{M}} \int_{\mathbb{M}} \mathfrak{t}_{f}\left(m_{1}, m_{2}\right) F_{M}\left(\mathrm{~d} m_{1}\right) F_{M}\left(\mathrm{~d} m_{2}\right)
\]

An overview of the most prominent specifications for \(\mathfrak{t}_{f}\) is presented in Table 1 .
Evaluating the test functions of Table 1 immediately yields different functional mark summary characteristics, which we briefly discuss next. With the exception of the last

\begin{table}
\begin{tabular}{|l|l|l|l|l|}
\hline Name for \(\mathfrak{t}_{f}\) & Test function \(\mathfrak{t}_{f}\) & Normalising factor \(\nabla_{\mathrm{t}_{f}}\) & Notation for \(\nabla_{\mathfrak{t}_{f}}(r)\) & Notation for \(\kappa_{\mathrm{t}_{f}}(r)\) \\
\hline \(\mathfrak{t}_{1}\) & \(m(\circ) m(\mathbf{r})\) & \(\mu_{M}^{2}\) & \(\tau_{m m}(r)\) & \(\kappa_{\text {mm }}(r)\) \\
\hline \(\mathfrak{t}_{2}\) & \(m(\circ)\) & \(\mu_{M}\) & \(\tau_{m \bullet}(r)\) & \(\kappa_{m \bullet}(r)\) \\
\hline \(\mathfrak{t}_{3}\) & \(m(\mathbf{r})\) & \(\mu_{M}\) & \(\tau_{\bullet m}(r)\) & \(\kappa_{\bullet m}(r)\) \\
\hline \(\mathfrak{t}_{4}\) & \(0.5(m(\circ)-m(\mathbf{r}))^{2}\) & \(\sigma_{M}^{2}\) & \(\gamma_{m m}(r)\) & \(\gamma_{m m}^{\mathrm{n}}(r)\) \\
\hline \(\mathfrak{t}_{5}\) & \(\left(m(\circ)-\mu_{M}\right)\left(m(\mathbf{r})-\mu_{M}\right)\) & \(\sigma_{M}^{2}\) & \(l_{m m}^{\text {Shi }}(r)\) & \(I_{m m}^{\mathrm{Shi}}(r)\) \\
\hline \(\mathfrak{t}_{6}\) & \(\left(m(\circ)-\mu_{M}(\mathbf{r})\right)\left(m(\mathbf{r})-\mu_{M}(\mathbf{r})\right)\) & \(\sigma_{M}^{2}\) & \(l_{m m}^{\text {Sch }}(r)\) & \(I_{m m}^{\text {Sch }}(r)\) \\
\hline
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 1: Overview of prominent test function specifications with \(\mu_{M}\) and \(\sigma_{M}^{2}\) denoting the unconditional mark mean and mark variance, \(\mu_{M}(\mathbf{r})\) the mean for the second point at location \(\mathbf{r}\) under the condition that there are points at ∘ and \(\mathbf{r}\).}
\end{table}
two test function, the normalising factor \(\nabla_{\mathfrak{t}_{f}}\) (3rd column) corresponds to the expected valued of the test function if \(r \rightarrow \infty\) and follows directly from the solution of (4) with \(\mathfrak{t}_{f}\left(m_{1}, m_{2}\right)\) replaced by the specific test function (as given in the 2nd column of Table 1) (see Illian et al., 2008). For \(\mathfrak{t}_{5}\) and \(\mathfrak{t}_{6}\), however, \(\nabla_{\mathfrak{t}_{f}}\) is set to \(\sigma_{M}^{2}\) in close analogy to Moran's I (Moran, 1950). Instead of \(\nabla_{\mathfrak{t}_{f}}(r)\), it is sometimes preferable to compute the \(\mathfrak{t}_{f}\)-correlation function (Penttinen and Stoyan, 1989),
\[
\kappa_{\mathrm{t}_{f}}(r)=\frac{\nabla_{\mathrm{t}_{f}}(r)}{\nabla_{\mathrm{t}_{f}}}=\frac{\varrho_{\mathrm{t}_{f}}^{(2)}(r)}{\varrho^{(2)}(r)} / \nabla_{\mathrm{t}_{f}},
\]
which normalizes the conditional expectation of the test function \(\nabla_{\mathfrak{t}_{f}}(r)\) by its expectation \(\nabla_{\mathfrak{t}_{f}}\) for \(r \rightarrow \infty\) such that \(\kappa_{\mathfrak{t}_{f}}(r)=1\) for all \(r\) under mark independence (or by \(\sigma_{M}^{2}\) in the last two cases). For completeness, Table 1 covers both the unnormalized \(\left(\nabla_{\mathfrak{t}_{f}}(r)\right)\) and related normalized ( \(\kappa_{\mathfrak{t}_{f}}(r)\) ) expressions.

A classic summary is the conditional mean product of marks of two points being a distance \(r\) apart, \(\tau_{m m}(r)\), and its scaled version \(\kappa_{m m}(r)\) often termed the mark correlation function (Stoyan and Stoyan, 1994). If large (resp. small) marks systematically co-occur at interpoint distance \(r\), their pairwise products will also be large (resp. small) and deviate from the independent mark assumption, i.e. the mark mean squared. Both \(\mathbf{r}\)-mark functions \(\tau_{m \bullet}\) and \(\tau_{\bullet m}\) and the related \(\mathbf{r}\)-mark correlation functions \(\kappa_{m}\). and \(\kappa_{\bullet m}\) can be interpreted as the conditional expectation of the mark of a point given that there is another point with distance \(r\). This expectation often deviates from \(\mu_{M}\), the expectation under independent marks, when the marks are dependent on the existence of other points (Schlather et al., 2004; Myllymäki et al., 2015). The mark variogram \(\gamma_{m m}(r)\) is a measure of the average dispersion (Cressie, 1993). It helps to detect situations where the marks of points close together tend to be more similar (or different) than expected under mark independence. Shimantani's (Shimatani, 2002) and Schlather's (Schlather, 2001 \(I\)-functions \(\iota_{m m}^{\mathrm{Shi}}\) and \(\iota_{m m}^{\mathrm{Sch}}\) and their normalized versions \(I_{m m}^{\mathrm{Shi}}\) and \(I_{m m}^{\mathrm{Shi}}\) can be seen
as adaptations of Moran's I (Moran, 1950) to spatial point processes. They are helpful for identifying potential spatial autocorrelation among the marks. Although these test functions are similar in spirit, Schlather applies a centering not by the unconditional mean but by the conditional mean \(\mathbb{E}_{\circ, r}[m(\mathbf{r})]=\mu_{M}(\mathbf{r})\) of the second point.

We note that the above construction principle through differently specified test functions can also be applied to define so-called nearest-neighbour correlation indices, i.e. numerical mark summary characteristics, where instead of \(m(\circ)\) and \(m(\mathbf{r})\) the test functions only consider the marks at the origin and its nearest-neighbouring point(s) (see Stoyan and Stoyan, 1994, for a detailed discussion).

\subsection*{2.3 Mark-weighted summary characteristics for real-valued marks}

Apart from the above summary characteristics, some authors propose to use the test function as a weight for classic second-order spatial point process characteristics. Initially proposed by Penttinen et al. (1992), the mark-weighted \(K_{\mathfrak{t}_{f}}\) function, which is a generalisation of Ripley's \(K\)-function (Ripley, 1976) to real-valued marks, is defined as
\[
K_{\mathfrak{t}_{f}}(r)=\frac{\mathbb{E}_{\circ, r}\left[\sum_{\left(x_{i}, m\left(x_{i}\right)\right) \in X} \mathfrak{t}_{f}\left(m(\circ), m\left(x_{i}\right)\right) \mathbb{1}_{b(\circ, r)}\left\{x_{i}\right\}\right]}{\lambda \nabla_{\mathfrak{t}_{f}}}
\]

Here, \(b(\circ, r)\) is a disc of radius \(r\) centered at the origin, \(\lambda\) is the intensity of \(\breve{X}\), and \(\mathfrak{t}_{f}(\cdot)\) is any test function as presented in Table 1. Again, the precise interpretation of the \(K_{\mathfrak{t}_{f}}\) function depends on the specific test function under study. For \(\mathfrak{t}_{f}=\mathfrak{t}_{1}, \nabla_{\mathfrak{t}_{f}}\) equals \(\mu_{M}^{2}\) and \(K_{\mathfrak{t}_{f}}\) is denoted by \(K_{m m}\),
\[
K_{m m}(r)=\frac{\mathbb{E}_{\circ, r}\left[\sum_{\left(x_{i}, m\left(x_{i}\right)\right) \in X} m(\circ) \cdot m\left(x_{i}\right) \mathbb{1}_{b(\circ, r)}\left\{x_{i}\right\}\right]}{\lambda \mu_{M}^{2}}
\]

Recalling the definition of Ripley's \(K\) function, \(K_{m m}\) can be interpreted as the expected number of further points within a distance \(r\) weighted by the pairwise product of marks. If the average product of marks coincides with the mark mean squared, \(\nabla_{\mathfrak{t}_{1}}(r)\) equals \(\mu_{M}^{2}\) and \(K_{m m}\) reduces to \(K\). However, for dependent marks, i.e. when \(\nabla_{\mathfrak{t}_{1}}(r)>\mu_{M}^{2}\) (resp. \(\nabla_{\mathfrak{t}_{1}}(r)<\mu_{M}^{2}\) ) for some \(r, K_{m m}(r)>K(r)\) (resp. \(K_{m m}(r)<K(r)\) ). We note that, for nicer visualization, it is often preferable to use \(L_{\mathfrak{t}_{f}}(r)=\sqrt{K_{\mathfrak{t}_{f}}(r) / \pi}\), a variance stabilising and centered version of \(K_{\mathfrak{t}_{f}}\) instead of \(K_{\mathfrak{t}_{f}}\) to control for the strict monotonic behaviour of the \(K\)-function with respect to the distance \(r\).

\section*{3 Composition-valued marked point processes}

\subsection*{3.1 Composition-valued marks}

To extend spatial point processes to composition-valued marks, let \(\left\{\left(x_{i}, \mathbf{c}\left(x_{i}\right)\right)\right\}_{i=1}^{n}\) denote a set of \(n\) points \(x_{i} \in \mathbb{W} \subset \mathbb{R}^{2}\) with associated marks \(\mathbf{c}\left(x_{i}\right)=\left(c_{1}\left(x_{i}\right), \ldots, c_{D}\left(x_{i}\right)\right)^{\top}\) living in a \(D\)-part simplex \(\mathbb{S}^{D} \subset \mathbb{R}^{D}\), where for some \(\mathfrak{w} \in \mathbb{R}\)
\[
\mathbb{S}^{D}=\left\{\mathbf{c}=\left(c_{1}, c_{2}, \ldots, c_{D}\right)^{\top} \in \mathbb{R}^{D} \mid c_{j} \geqslant 0, j=1,2, \ldots, D ; \sum_{j=1}^{D} c_{j}=\mathfrak{w}\right\} .
\]

Common choices of \(\mathfrak{w}\) include \(\mathfrak{w}=1\) and \(\mathfrak{w}=100\) (per cent) depending on the marks at hand. Thus, the mark at each location is a composition of \(D\) non-negative parts summing to a constant and we call any such mark composition-valued. That is, compositionvalued marked spatial point processes are a particular type of spatial compositional data (Aitchison, 1986), where each mark \(\mathbf{c}\left(x_{i}\right)\) quantitatively describes the relative importance of the individual parts with respect to a given total. We note that compositionvalued marks could have two potential origins: they might (i) directly arise from the data collection or (ii) be constructed from non-simple spatial point process scenarios with integer-valued marks, by summarizing the absolute numbers \(\tilde{c}_{j}, j=1, \ldots, D\), at a point location a-posteriori into relative numbers for distinct categories. Generally, any such absolute information can always be transformed into a composition-valued mark by dividing each part by the sum over all components, i.e. by applying the closure operation \(\mathbf{c}=\operatorname{cls}(\tilde{\mathbf{c}})=\left(\tilde{c}_{1} / \sum_{j=1}^{D} \tilde{c}_{j}, \ldots, \tilde{c}_{D} / \sum_{j=1}^{D} \tilde{c}_{j}\right)^{\top}\). Consequently, both the absolute and the relative, composition-valued marks can be analysed depending on the question of interest, providing two different views of the marks as discussed below in Section 3.6.

The \(\mathbb{S}^{D}\) space of compositions can be equipped with a finite ( \(D-1\) ) dimensional Euclidean vector space structure, i.e. the Aitchison geometry, with the perturbation \(\mathbf{c} \oplus \mathbf{c}^{\prime}=\operatorname{cls}\left(c_{1} c_{1}^{\prime}, c_{2} c_{2}^{\prime}, \ldots, c_{D} c_{D}^{\prime}\right)\), a commutative group operation on the simplex with neutral element \(\mathbf{n}=\operatorname{cls}(1,1, \ldots, 1)\) and inverse operation \(\mathbf{c} \ominus \mathbf{c}^{\prime}=\mathbf{c} \oplus\left((-1) \odot \mathbf{c}^{\prime}\right)\), and the powering operation \(\xi \odot \mathbf{c}=\operatorname{cls}\left(c_{1}^{\xi}, c_{2}^{\xi}, \ldots, c_{D}^{\xi}\right)\), where \(\mathbf{c}, \mathbf{c}^{\prime} \in \mathbb{S}^{D}\) and \(\xi \in \mathbb{R}\) Aitchison, 2001). Additionally, the inner product is defined as
\[
\left\langle\mathbf{c}, \mathbf{c}^{\prime}\right\rangle_{A}=\frac{1}{2 D} \sum_{l=1}^{D} \sum_{j=1}^{D} \log \left(\frac{c_{l}}{c_{j}}\right) \log \left(\frac{c_{l}^{\prime}}{c_{j}^{\prime}}\right)
\]
yielding the norm \(\|\mathbf{c}\|_{A}=\sqrt{\langle\mathbf{c}, \mathbf{c}\rangle_{A}}\) and the associated distance \(d_{A}\left(\mathbf{c}, \mathbf{c}^{\prime}\right)=\left\|\mathbf{c} \ominus \mathbf{c}^{\prime}\right\|_{A}\). Note that the composition-valued marks can be transformed to real-valued coordinates through a map function \(\psi: \mathbb{S}^{D} \rightarrow \mathbb{R}^{\tilde{D}}, \mathbf{c} \mapsto \psi(\mathbf{c})\) where \(\tilde{D}\) is determined by the particular choice of \(\psi\). For particular useful choices, an isometric isomorphism exists between the \(\mathbb{S}^{D}\) and the \(\mathbb{R}^{\tilde{D}}\) (Billheimer et al., 2001; Pawlowsky-Glahn and Egozcue, 2001),
which allows expressing the composition-valued marks in real coordinates such that their relations are correspondent to those in the Aitchison geometry. Below in Section 3.2, we give alternatives for the transformation \(\psi\). After transformation, statistical analysis methods can be performed in \(\mathbb{R}^{\tilde{D}}\) (Mateu-Figueras et al., 2011).

\subsection*{3.2 Transformations}

Writing \(\psi_{j}(\mathbf{c})\) for the \(j\)-th element of \(\psi(\mathbf{c})=\left(\psi_{1}(\mathbf{c}), \ldots, \psi_{\tilde{D}}(\mathbf{c})\right)^{\top}\), different coordinate representations can be derived (Pawlowsky-Glahn and Buccianti, 2011). Apart from the log-ratio (lr) transformation early specifications of \(\psi\) include the additive log-ratio (alr) transformation Aitchinson and Shen, 1980) where \(\psi_{j}(\mathbf{c})=\log \left(c_{j} / c_{D}\right), j=1, \ldots, \tilde{D}= D-1\) i.e. the log-ratios relative to the \(D\)-th component, and the centered log-ratio (clr) (Aitchinson, 1983) transformation \(\psi_{j}(\mathbf{c})=\log \left(c_{j} / g(\mathbf{c})\right), j=1, \ldots, \tilde{D}=D\), i.e. the \(\log\) ratios relative to the geometric mean \(g(\mathbf{c})=\left(\prod_{j=1}^{D} c_{j}\right)^{1 / D}\). While the alr transformation yields an isomorphic but non-isometric relation between the two spaces, i.e. it does not preserve distances, the clr transformation results in an isometric isomorphism by mapping the marks from the simplex to a hyperplane \(H \subset \mathbb{R}^{D}\) that is orthogonal to the vector of ones. However, the clr imposes a sum-to-zero constraint on the transformed marks, which can make analysis more difficult e.g. due to degenerated distributions and singular covariance matrices.

We observe that neither the alr nor the clr transformation can be directly linked to an orthonormal coordinate system on the simplex. However, an orthonormal basis \(\left(\mathbf{e}_{1}, \mathbf{e}_{2}, \ldots, \mathbf{e}_{D-1}\right)\) on the simplex \(\mathbb{S}^{D}\) with respect to the inner product can be derived using the Gram-Schmidt procedure giving
\[
\mathbf{c}=\bigoplus_{j=1}^{D-1}\left\langle\mathbf{c}, \mathbf{e}_{j}\right\rangle_{A} \odot \mathbf{e}_{j}
\]

This coordinate representation corresponds to the isometric log-ratio (ilr) transformation, a class of orthonormal coordinate representations, defined by
\[
\mathrm{i} \operatorname{lr}(\mathbf{c})=\left(\left\langle\mathbf{c}, \mathbf{e}_{1}\right\rangle_{A},\left\langle\mathbf{c}, \mathbf{e}_{2}\right\rangle_{A}, \ldots,\left\langle\mathbf{c}, \mathbf{e}_{D-1}\right\rangle_{A}\right)
\]
and establishes an isometric isomorphism through the map between \(\mathbb{S}^{D}\) and \(\mathbb{R}^{D-1}\). Further, the ilr transformation is related to the clr and log transformations through \(\operatorname{ilr}(\mathbf{c})=\operatorname{clr}(\mathbf{c}) \mathbf{H}_{\mathrm{D}}^{\top}=\log (\mathbf{c}) \mathbf{H}_{\mathrm{D}}^{\top}\) where \(\mathbf{H}_{\mathrm{D}}\) denotes a \(((D-1) \times D)\)-dimensional Helmert matrix with rows \(\mathbf{h}_{j}=\operatorname{clr}\left(\mathbf{e}_{j}\right), j=1, \ldots, D-1\) satisfying \(\mathbf{H}_{\mathrm{D}} \mathbf{H}_{\mathrm{D}}^{\top}=\mathbf{I}_{\mathrm{D}-1}\) and \(\mathbf{H}_{\mathrm{D}}^{\top} \mathbf{H}_{\mathrm{D}}= \mathbf{G}_{\mathrm{D}}\) and \(\mathbf{G}_{\mathrm{D}}\) is the \(D\)-dimensional centering matrix \(\mathbf{G}_{\mathrm{D}}=\mathbf{I}_{\mathrm{D}}-D^{-1} \mathbb{1}_{\mathrm{D}} \mathbb{1}_{\mathrm{D}}^{\top}=\mathbf{G}_{\mathrm{D}}^{2}, \mathbf{I}_{\mathrm{D}}\) is the identity matrix of dimension ( \(D \times D\) ), and \(\mathbb{1}_{\mathrm{D}}\) a ( \(D \times 1\) ) vector of ones. Due to the isometric isomorphism established between the Aitchison and the Euclidean geometry by the clr and the ilr transformations, the Aitchison inner product, distances
and metrics coincide with their Euclidean counterparts on the transformed quantities, such that for any \(\mathbf{c}, \mathbf{c}^{\prime} \in \mathbb{S}^{D}\)
\[
d_{A}\left(\mathbf{c}, \mathbf{c}^{\prime}\right)=d_{E}\left(\operatorname{clr}(\mathbf{c}), \operatorname{clr}\left(\mathbf{c}^{\prime}\right)\right)=d_{E}\left(\operatorname{ilr}(\mathbf{c}), \operatorname{ilr}\left(\mathbf{c}^{\prime}\right)\right)
\]
with \(d_{E}\) the Euclidean distance. The ilr transformation is equivalent to the logitfunction used in logistic regression if \(D=2\). When \(D>2\), infinitely many orthonormal basis systems exist and the concrete choice has a crucial impact on the interpretation of the projected data. Particular choices of an orthonormal basis are coordinate representations using (i) balances (Egozcue and Pawlowsky-Glahn, 2005) in which each balancing element can be interpreted as the normalised log-ratio of the geometric means (centers) of two groups, and (ii) pivot coordinates (Fišerová and Hron, 2011; Hron et al., 2017). Balances originate from a sequential binary partition method, which involves dividing the composition into two parts. The \(j\)-th ilr coefficient using pivot coordinates can be expressed as
\[
\operatorname{ilr}_{j}(\mathbf{c})=\sqrt{\left(\frac{D-j}{D-j+1}\right)} \log \left\{\frac{c_{j}}{\sqrt[D-j]{\prod_{k=j+1}^{D} c_{k}}}\right\}, j=1, \ldots, D-1
\]
(Fišerová and Hron, 2011). The initial ilr coefficient is akin to the first clr coefficient, scaled by \(\sqrt{ } D /(D-1)\), and is easily interpreted as the log-ratio of the respective component to the geometric mean. In contrast, interpreting subsequent coefficients is more complex. To address this, the literature proposes generalised pivot coordinates using permuted compositions and symmetric pivot coordinates (see e.g. Kynčlová et al., 2017; Hron et al., 2021).

While the above transformations allow for a representation of the compositionvalued marks in coordinates in a Euclidean space, the underlying log-operations are undefined for zero values. In what follows, we assume that the proportions for all \(D\) parts are non-zero. However, as zeros might be present in some marked point process scenarios when potentially not all components are observed at each location, we also provide a treatment of different transformations in the presence of structural zeros, see Section 1 in the supplementary material.

\subsection*{3.3 First-order tools for composition-valued marks}

We first review first-order mark characteristics for the composition-valued marks. General characteristics for a strictly-positive sample \(\mathbf{c}_{1}, \ldots, \mathbf{c}_{n}\) of compositions commonly applied in the literature include the geometric center, i.e. the closed geometric mean,
\[
\operatorname{cen}(\mathbf{c})=\frac{1}{n} \odot \bigoplus_{j=1}^{n} \mathbf{c}_{j}=\operatorname{clr}^{-1}\left(\frac{1}{n} \sum_{j=1}^{n} \operatorname{clr}\left(\mathbf{c}_{j}\right)\right)
\]
where \(\operatorname{clr}^{-1}(\tilde{\mathbf{c}})=\operatorname{cls}(\exp (\tilde{\mathbf{c}}))\). It serves as the mean composition-valued mark. Further, the variation of the composition-valued marks can be described through the variation matrix \(\mathbf{T}\) with elements \(t_{j l}=\operatorname{Var}\left[\log \left(c_{j} / c_{l}\right)\right], j, l=1, \ldots, D\), and the total, i.e. metric, variance
\[
\mathrm{m} \operatorname{Var}[\mathbf{c}]=\frac{1}{2 D} \sum_{j, l=1}^{D} t_{j l}=\frac{1}{n-1} \sum_{j=1}^{n} d_{A}^{2}\left(\mathbf{c}_{j}, \operatorname{cen}(\mathbf{c})\right),
\]
serving as a global measure of dispersion Aitchison, 1986; Pawlowsky-Glahn and Egozcue, 2001). Instead of the variation matrix, a normalised variation matrix \(\overline{\mathbf{T}}=0.5 \cdot \mathbf{T}\) can be used.

In the presence of zero components, alternative measures of centrality include the spatial median (Brown, 1983), the graph median (Sharp, 2006) and the Fréchet mean of the \(\alpha\)-transformed composition (Tsagris et al., 2011).

\subsection*{3.4 Componentwise summary characteristics for compositionvalued marks}

We now define novel componentwise characteristics for composition-valued marks analogously to the mark characteristics of Section 2.2. These are useful to look at the spatial behaviour of all components individually, while we propose characteristics for the whole composition in the next subsection. We here explicitly focus on the case where each point is augmented by exactly one composition. Extensions to multivariate settings are outlined in Section 6 of the Supplementary material. Recall that \(\mathbf{c}(\circ)\) and \(\mathbf{c}(\mathbf{r})\) denote the composition-valued marks for a pair of points at locations \(\circ\) and \(\mathbf{r}\) at distance \(\|\mathbf{r}\|=r\), and \(\psi(\mathbf{c})\) is the transformed composition-valued mark with \(j\) th element or component \(\psi_{j}(\mathbf{c})\), where \(\psi: \mathbb{S}^{D} \mapsto \mathbb{R}^{\tilde{D}}\). We then define
\[
\nabla_{\mathfrak{t}_{f}}^{\psi, j l}(r)=\mathbb{E}_{\circ, r}\left[\mathfrak{t}_{f}^{\psi, j l}\left(\psi_{j}(\mathbf{c}(\circ)), \psi_{l}(\mathbf{c}(\mathbf{r}))\right)\right]
\]
where \(\mathfrak{t}_{f}^{\psi, j l}\) denotes a test function specific to the transformed composition-valued marks. The test functions of Table 1 can be employed as \(\mathfrak{t}_{f}^{\psi, j l}\); for clarity, they are re-expressed for the transformed marks \(\psi(\mathbf{c})\) in Table 2. Further, we let \(\nabla_{\mathfrak{t}_{f}}^{\psi, j l}\) stand for the limiting case of \(\nabla_{\mathfrak{t}_{f}}^{\psi, j l}(r)\) when \(r \rightarrow \infty\), i.e.
\[
\nabla_{\mathfrak{t}_{f}}^{\psi, j l}=\int_{\mathbb{R}^{\tilde{D}}} \int_{\mathbb{R}^{\tilde{D}}} \mathfrak{t}_{f}^{\psi, j l}\left(\psi_{j}\left(c_{1}\right), \psi_{l}\left(c_{2}\right)\right) \varpi\left(\mathrm{d} \psi_{j}\left(c_{1}\right)\right) \varpi\left(\mathrm{d} \psi_{l}\left(c_{2}\right)\right)
\]
and define the \(\mathfrak{t}_{f}^{\psi, j l}\)-correlation function \(\kappa_{\mathfrak{t}_{f}}^{\psi, j l}(r)\) as
\[
\kappa_{\mathfrak{t}_{f}}^{\psi, j l}(r)=\nabla_{\mathfrak{t}_{f}}^{\psi, j l}(r) / \nabla_{\mathfrak{t}_{f}}^{\psi, j l}
\]

\begin{table}
\begin{tabular}{|l|l|l|l|l|}
\hline Name for \(\mathfrak{t}_{f}^{\psi, j l}\) & Test function \(\mathfrak{t}_{f}^{\psi, j l}\) & Normalising factor \(\nabla_{\mathrm{t}_{f}}^{\psi, j l}\) & Notation for \(\nabla_{\mathfrak{t}_{f}}^{\psi, j l}(r)\) & Notation for \(\kappa_{\mathrm{t}_{f}}^{\psi, j l}(r)\) \\
\hline \(\mathfrak{t}_{1}^{\psi, j l}\) & \(\psi_{j}(\mathbf{c}(\circ)) \psi_{l}(\mathbf{c}(\mathbf{r}))\) & \(\mu_{j}^{\psi} \mu_{l}^{\psi}\) & \(\tau_{j l}^{\psi}(r)\) & \(\kappa_{j l}^{\psi}(r)\) \\
\hline \(\mathfrak{t}_{2}^{\psi, j l}\) & \(\psi_{j}(\mathbf{c}(\circ))\) & \(\mu_{j}^{\psi}\) & \(\tau_{j \bullet}^{\psi}(r)\) & \(\kappa_{j_{\bullet}}^{\psi}(r)\) \\
\hline \(\mathfrak{t}_{3}^{\psi, j l}\) & \(\psi_{l}(\mathbf{c}(\mathbf{r}))\) & \(\mu_{l}^{\psi}\) & \(\tau_{\bullet l}^{\psi}(r)\) & \(\kappa_{\bullet l}^{\psi}(r)\) \\
\hline \(\mathfrak{t}_{4}^{\psi, j l}\) & \(0.5\left(\psi_{j}(\mathbf{c}(\circ))-\psi_{l}(\mathbf{c}(\mathbf{r}))\right)^{2}\) & \(\zeta_{j l}^{\psi}\) & \(\gamma_{j l}^{\psi}(r)\) & \(\gamma_{j l}^{\psi, \mathrm{n}}(r)\) \\
\hline \(\mathfrak{t}_{5}^{\psi, j l}\) & \(\left(\psi_{j}(\mathbf{c}(\circ))-\mu_{j}^{\psi}\right)\left(\psi_{l}(\mathbf{c}(\mathbf{r}))-\mu_{l}^{\psi}\right)\) & \(\sigma_{j l}^{\psi}\) & \(\iota_{j l}^{\psi, \text { Shi }}(r)\) & \(I_{j l}^{\psi, \text { Shi }}(r)\) \\
\hline \(\mathfrak{t}_{6}^{\psi, j l}\) & \(\left(\psi_{j}(\mathbf{c}(\circ))-\mu_{j}^{\psi}(\mathbf{r})\right)\left(\psi_{l}(\mathbf{c}(\mathbf{r}))-\mu_{l}^{\psi}(\mathbf{r})\right)\) & \(\sigma_{j l}^{\psi}\) & \(\iota_{j l}^{\psi, \mathrm{Sch}}(r)\) & \(I_{j l}^{\psi, \text { Sch }}(r)\) \\
\hline
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 2: Specifications for the componentwise test function \(\mathfrak{t}_{f}^{\psi, j l}\left(\psi_{j}(\mathbf{c}(\circ)), \psi_{l}(\mathbf{c}(\mathbf{r}))\right)\) with \(\mu_{j}^{\psi}\) denoting the unconditional mark mean and \(\mu_{j}^{\psi}(\mathbf{r})\) the conditional mark mean of the \(j\)-th element of the \(\psi\)-transformed composition-valued marks, where \(\zeta_{j l}^{\psi}=0.5\left[\sigma_{j j}^{\psi}+\sigma_{l l}^{\psi}+\left(\mu_{j}^{\psi}-\mu_{l}^{\psi}\right)^{2}\right]\) with \(\sigma_{j l}^{\psi}\) denoting the covariance of the \(j\)-th and the \(l\)-th element of the \(\psi\)-transformed composition-valued mark.}
\end{table}
where (11) is again replaced by the corresponding variance in the case of the last two test functions. The choice of the test function \(\mathfrak{t}_{f}^{\psi, j l}\) determines the focus of the analysis, as the different characteristics highlight different properties of marks (see Section 2.2).

Although in the following we restrict our discussion to only some characteristics, the same principle applies to all test functions in Table 2. Consider as an example test function \(\mathfrak{t}_{1}^{\psi, j l}\). It yields a componentwise conditional mean product of marks \(\tau_{j l}^{\psi}(r)\),
\[
\tau_{j l}^{\psi}(r)=\mathbb{E}_{\circ, r}\left[\psi_{j}(\mathbf{c}(\circ)) \psi_{l}(\mathbf{c}(\mathbf{r}))\right],
\]
which describes the average product of two components of the \(\psi\)-transformed marks for any pair of points as a function of the distance \(r\).

Restricting to pairs of the \(j\)-th element of \(\psi(\mathbf{c})\) at two locations, \(j=l\), and specifying \(\psi\) in (13) as log-ratio between parts \(j_{1}\) and \(j_{2}\), to directly focus on their relative contribution to the total, yields the mean pairwise product of mark \(\log\)-ratios \(\tau_{j j}^{\mathrm{lr}}(r)\),
\[
\tau_{j j}^{\operatorname{lr}}(r)=\mathbb{E}_{\circ, r}\left[\log \left(\frac{c_{j_{1}}(\circ)}{c_{j_{2}}(\circ)}\right) \log \left(\frac{c_{j_{1}}(\mathbf{r})}{c_{j_{2}}(\mathbf{r})}\right)\right]
\]
for \(j\) indexing the \(\tilde{D}=D^{2}\) ordered pairs \(\left(j_{1}, j_{2}\right)=(1,1), \ldots,(1, D), \ldots,(D, D)\). This characteristic can highlight if the log-ratios at nearby points are correlated and their products thus tend to be larger (smaller) than expected under independence, \(\nabla_{\mathfrak{t}_{1}}^{\operatorname{lr}, j j}\). It is easy to see that \(\tau_{j j}^{\mathrm{lr}}\) equals the squared mean \(\left(\mu_{j}^{\mathrm{lr}}\right)^{2}\) of the log-ratios of the \(j_{1}\)-th and \(j_{2}\)-th parts under independent marks.

Similarly, the componentwise log-ratio \(\mathbf{r}\)-mark functions \(\tau_{j_{\bullet}}^{\mathrm{lr}}\) and \(\tau_{\bullet j}^{\mathrm{lr}}\) result from taking the conditional expectation of either \(\mathfrak{t}_{2}^{\mathrm{lr}, j j}\) or \(\mathfrak{t}_{3}^{\mathrm{lr}, j j}\), which both coincide with the
mean of the log-ratio of the \(j\)-th part, \(\mu_{j}^{\mathrm{lr}}\), in the case of independent marks. Reflecting the mean behaviour of the transformed mark at either the first or second point, deviations of the empirical curves from \(\mu_{j}^{\mathrm{lr}}\) at some distances indicate the presence of spatial structure in the mark parts, i.e. changes in the average mark ratios if a second point is present at distance \(r\). As such, both quantities could be helpful tools to identify potential interrelations of the points and marks.

Substituting \(\mathfrak{t}_{4}^{\psi}\) for \(\mathfrak{t}_{1}^{\psi}\) in (13) yields a generic componentwise mark variogram
\[
\gamma_{j l}^{\psi}(r)=\mathbb{E}_{\circ, r}\left[0.5\left(\psi_{j}(\mathbf{c}(\circ))-\psi_{l}(\mathbf{c}(\mathbf{r}))\right)^{2}\right]
\]

The precise form again depends on the specific choice of \(\psi\). Using a transformation into log-ratios and \(j=l\) leads to a componentwise log-ratio mark variogram \(\gamma_{j j}^{\mathrm{lr}}\) for composition-valued marks defined by
\[
\gamma_{j j}^{\operatorname{lr}}(r)=\mathbb{E}_{\circ, r}\left[\frac{1}{2}\left(\log \left(\frac{c_{j_{1}}(\circ)}{c_{j_{2}}(\circ)}\right)-\log \left(\frac{c_{j_{1}}(\mathbf{r})}{c_{j_{2}}(\mathbf{r})}\right)\right)^{2}\right]
\]
for \(j\) indexing \(\left(j_{1}, j_{2}\right)=(1,1), \ldots,(1, D), \ldots,(D, D)\). Similar to the classic mark variogram, \(\gamma_{j j}^{\mathrm{lr}}(r)\) concerns the average pairwise variation of the \(j\)-th log-ratio of \(\psi(\mathbf{c})\) at two distinct points at distance \(r\) and tends to the variance \(\sigma_{j j}^{\mathrm{lr}}\) of the log-ratios of the \(j_{1}\)-st and \(j_{2}\)-nd parts for \(r \rightarrow \infty\). The mark variogram can be used to investigate the heterogeneity among the marks, i.e. if the proportions of the specific parts for any pair of points are on average more similar in value for small distances. We note that for each \(r\) all of the above log-ratio characteristics could be stored in local ( \(D \times D\) ) matrices with \(D(D-1)\) log-ratios of different parts in the off-diagonal entries and zeros for the log-ratios of all parts with themselves on its diagonal. In particular, collecting all \(\gamma_{j j}^{\mathrm{lr}}(r)\) into a local mark variogram matrix \(\boldsymbol{\Gamma}^{\mathrm{lr}}(r)\) is similar in spirit to a local variation matrix \(\mathbf{T}(r)\) which captures the spatial dispersion of the composition-valued mark at the distance \(r\).

While the above log-ratio characteristics are most useful for autocovariances and autocorrelations, the clr and ilr transformations allow for both auto- and cross-characteristic formulations. Choosing \(\psi\) to denote the clr transformation, evaluation of the conditional expectation of \(\mathfrak{t}_{1}^{\psi, j l}\) and \(\mathfrak{t}_{4}^{\psi, j l}\) yields the conditional mean product of clr marks, \(\tau_{j l}^{\text {clr }}(r)\), and the clr mark variogram, \(\gamma_{j l}^{\text {clr }}\), defined by
\[
\tau_{j l}^{\mathrm{clr}}(r)=\mathbb{E}_{\circ, r}\left[\log \left(\frac{c_{j}(\circ)}{g(\mathbf{c})(\circ)}\right) \cdot \log \left(\frac{c_{l}(\mathbf{r})}{g(\mathbf{c})(\mathbf{r})}\right)\right]
\]
and
\[
\gamma_{j l}^{\mathrm{clr}}(r)=\mathbb{E}_{\circ, r}\left[\frac{1}{2}\left(\log \left(\frac{c_{j}(\circ)}{g(\mathbf{c})(\circ)}\right)-\log \left(\frac{c_{l}(\mathbf{r})}{g(\mathbf{c})(\mathbf{r})}\right)\right)^{2}\right]
\]
respectively. Due to the construction of (17) and (18), both functions describe the average spatial association/variation of the \(j\)-th and \(l\)-th parts relative to the geometric mean, including both auto- (for \(j=l\) ) and cross-relations (for \(j \neq l\) ). In particular, cross-relations might be useful to analyse the association between different parts, e.g. the proportions of two distinct parasite species for neighbouring trees. Both quantities allow for similar interpretations as their log-ratio counterparts and could provide useful information on the distributional behaviour of the clr transformed parts. Inserting the corresponding terms into (11) immediately implies that \(\tau_{j l}^{\text {clr }}(r)\) converges to the product of means \(\mu_{j}^{\mathrm{clr}} \cdot \mu_{l}^{\mathrm{clr}}\) of the clr transformed parts \(j, l\) for \(r \rightarrow \infty\) as
\[
\begin{aligned}
\nabla_{\mathrm{t}_{1}}^{\mathrm{clr}, j l} & =\int_{\mathbb{R}} \int_{\mathbb{R}} \operatorname{clr}_{j}(\mathbf{c}(\circ)) \operatorname{clr}_{l}(\mathbf{c}(\mathbf{r})) \varpi\left(\mathrm{d} \mathrm{clr}_{j}(\mathbf{c}(\circ)) \varpi\left(\mathrm{d} \mathrm{clr}_{j}(\mathbf{c}(\mathbf{r}))\right)\right. \\
& =\int_{\mathbb{R}} \operatorname{clr}_{j}(\mathbf{c}) \varpi\left(\mathrm{d} \mathrm{clr}_{j}(\mathbf{c})\right) \int_{\mathbb{R}} \operatorname{clr}_{l}(\mathbf{c}) \varpi\left(\mathrm{d} \mathrm{clr}_{l}(\mathbf{c})\right) \\
& =\mu_{j}^{\mathrm{clr}} \mu_{l}^{\mathrm{clr}}
\end{aligned}
\]

Likewise, for \(\gamma_{j l}^{\text {clr }}(r)\), inserting \(\mathfrak{t}_{4}^{\text {clr, } j l}\) into (11) yields
\[
\begin{aligned}
\nabla_{\mathrm{t}_{4}}^{\mathrm{clr}, j l} & =\int_{\mathbb{R}} \int_{\mathbb{R}} 0.5\left(\operatorname{clr}_{j}(\mathbf{c}(\circ))-\operatorname{clr}_{l}(\mathbf{c}(\mathbf{r}))\right)^{2} \varpi\left(\mathrm{~d} \mathrm{clr}_{j}(\mathbf{c}(\circ))\right) \varpi\left(\mathrm{d} \mathrm{clr}_{l}(\mathbf{c}(\mathbf{r}))\right) \\
& =0.5\left[\int_{\mathbb{R}}\left(\operatorname{clr}_{j}(\mathbf{c})\right)^{2} \varpi\left(\mathrm{~d} \mathrm{clr}_{j}(\mathbf{c})\right)+\int_{\mathbb{R}}\left(\operatorname{clr}_{l}(\mathbf{c})\right)^{2} \varpi\left(\mathrm{~d} \mathrm{clr}_{l}(\mathbf{c})\right)\right. \\
& \left.-2 \int_{\mathbb{R}} \int_{\mathbb{R}} \operatorname{clr}_{j}(\mathbf{c}(\circ)) \operatorname{clr}_{l}(\mathbf{c}(\mathbf{r})) \varpi\left(\mathrm{d} \mathrm{clr}_{j}(\mathbf{c}(\circ))\right) \varpi\left(\mathrm{d} \mathrm{clr}_{l}(\mathbf{c}(\mathbf{r}))\right)\right] \\
& =0.5\left[\sigma_{j j}^{\mathrm{clr}}+\left(\mu_{j}^{\mathrm{clr}}\right)^{2}+\sigma_{l l}^{\mathrm{clr}}+\left(\mu_{l}^{\mathrm{clr}}\right)^{2}-2 \mu_{j}^{\mathrm{clr}} \mu_{l}^{\mathrm{clr}}\right] \\
& =0.5\left[\sigma_{j j}^{\mathrm{clr}}+\sigma_{l l}^{\mathrm{clr}}+\left(\mu_{j}^{\mathrm{clr}}-\mu_{l}^{\mathrm{clr}}\right)^{2}\right] \\
& =: \zeta_{j l}^{\mathrm{clr}}
\end{aligned}
\]

The clr characteristics for all \(D\) parts can efficiently be stored in local ( \(D \times D\) ) matrices including the local clr mark variogram matrix \(\boldsymbol{\Gamma}^{\mathrm{clr}}(r)=\left[\gamma_{j l}^{\mathrm{clr}}(r)\right]_{j, l=1, \ldots, D}\). Similarly, using ilr coordinates yields
\[
\tau_{j l}^{\mathrm{ilr}}(r)=\mathbb{E}_{\circ, r}\left[\operatorname{ilr}_{j}(\mathbf{c}(\circ)) \operatorname{ilr}_{l}(\mathbf{c}(\mathbf{r}))\right]
\]
and
\[
\gamma_{j l}^{\mathrm{ilr}}(r)=\mathbb{E}_{\circ, r}\left[\frac{1}{2}\left(\operatorname{ilr}_{j}(\mathbf{c}(\circ))-\operatorname{ilr}_{l}(\mathbf{c}(\mathbf{r}))\right)^{2}\right]
\]
where \(\mathrm{ilr}_{j}\) is the \(j\)-th ilr coordinate of the composition-valued marks and \(\operatorname{ilr}_{l}\) analogous.

\subsection*{3.5 Compositional summary characteristics for compositionvalued marks}

Instead of describing the spatial behaviour of one component of \(\psi\)-transformed compositions using a componentwise test function specification as covered in Table 2, we now discuss an extension which allows to evaluate the spatial properties of whole \(D\) part compositions. This allows to assess the variation and correlation of the complete composition of e.g. nearby trees or business sectors. Such compositional summary characteristics, which summarise the variation and interrelation between the \(D\)-part compositions for a pair of points as a function of the interpoint distance \(r\), can be derived using central concepts from the Aitchison geometry. In particular, defining \(\nabla_{\mathfrak{t}_{f}}^{\mathbb{S}}(r)\) to denote the conditional mean of the test function \(\mathfrak{t}_{f}^{S}\) of a \(D\)-part composition-valued mark \(\mathbf{c}\) at locations \(\circ\) and \(\mathbf{r}\) where \(\|\circ-\mathbf{r}\|=r\), and \(\kappa_{\mathfrak{t}_{f}}^{S}(r)\) as
\[
\kappa_{\mathfrak{t}_{f}}^{\mathbb{S}}(r)=\frac{\nabla_{\mathfrak{t}_{f}}^{\mathbb{S}}(r)}{\nabla_{\mathfrak{t}_{f}}^{\mathbb{S}}},
\]
where \(\nabla_{\mathfrak{t}_{f}}^{\mathrm{S}}\) extends (11) to
\[
\nabla_{\mathfrak{t}_{f}}^{\mathbb{S}}=\int_{\mathbb{S}^{D}} \int_{\mathbb{S}^{D}} \mathfrak{t}_{f}^{\mathbb{S}}(\mathbf{c}(\circ), \mathbf{c}(\mathbf{r})) \varpi(\mathrm{d} \mathbf{c}(\circ)) \varpi(\mathrm{d} \mathbf{c}(\mathbf{r}))
\]
and denotes the conditional expectation of the test function for \(r \rightarrow \infty\), allows to define different compositional mark summary characteristics. An overview of different test functions \(\mathfrak{t}_{f}^{S}\) and the corresponding characteristics is provided in Table 3 where here only those test functions are considered which allow a one-to-one relation to the corresponding componentwise test functions of Table 2. We note that instead of applying (11) to compute \(\nabla_{\mathfrak{t}_{f}}^{\mathrm{S}}\) for the compositional versions of Schlather's and Shimantani's \(I\) functions, \(\nabla_{\mathfrak{t}_{f}}^{S}\) is set to \(\sigma_{\mathbf{c}}^{2}\),
\[
\sigma_{\mathbf{c}}^{2}=\omega \sum_{j=1}^{\tilde{D}} \zeta_{j j}^{\psi}
\]
where \(\omega=1 / 2 D\) if a transformation into logratios is applied and \(\omega=1\) under ilr and clr transformations to allow for a close analogy to Moran's \(I\). Like the componentwise mark characteristics, all of these test functions can be used to highlight particular aspects of the distributional properties. However, providing different insights into the underlying structure of the composition-valued marks, the componentwise and the compositional mark characteristics both provide useful tools to investigate the individual contribution of the components to the results and to assess the overall variation and correlation. Noting the isometry of the Aitchison norm and distance to their Euclidean counterpart versions of the clr or ilr transformed compositions, the test functions of

\begin{table}
\begin{tabular}{lllll}
\hline Name & Test & Normalising & Notation & Notation \\
for \(\mathfrak{t}_{f}^{S}\) & function \(\mathfrak{t}_{f}^{S}\) & factor \(\nabla_{\mathfrak{t}_{f}}^{S}\) & for \(\nabla_{\mathfrak{t}_{f}}^{S}(r)\) & for \(\kappa_{\mathfrak{t}_{f}}^{S}(r)\) \\
\hline \(\mathfrak{t}_{1}^{S}\) & \(\langle\mathbf{c}(\circ), \mathbf{c}(\mathbf{r})\rangle_{A}\) & \(\mu_{\mathbf{c}}^{2}\) & \(\tau_{\mathbf{c c}}(r)\) & \(\kappa_{\mathbf{c c}}(r)\) \\
\(\mathfrak{t}_{4}^{S}\) & \(0.5\|\mathbf{c}(\circ) \ominus \mathbf{c}(\mathbf{r})\|_{A}^{2}\) & \(\sigma_{\mathbf{c}}^{2}\) & \(\gamma_{\mathbf{c c}}(r)\) & \(\gamma_{\mathbf{c c}}^{n}(r)\) \\
\(\mathfrak{t}_{5}^{S}\) & \(\langle\mathbf{c}(\circ) \ominus \operatorname{cen}(\mathbf{c}), \mathbf{c}(\mathbf{r}) \ominus \operatorname{cen}(\mathbf{c})\rangle_{A}\) & \(\sigma_{\mathbf{c}}^{2}\) & \(L_{\mathbf{c c}}^{S \mathrm{Si}}(r)\) & \(I_{\mathbf{c c}}^{S h \mathrm{i}}(r)\) \\
\(\mathfrak{t}_{6}^{\mathbb{S}}\) & \(\langle\mathbf{c}(\circ) \ominus \operatorname{cen}(\mathbf{c})(\mathbf{r}), \mathbf{c}(\mathbf{r}) \ominus \operatorname{cen}(\mathbf{c})(\mathbf{r})\rangle_{A}\) & \(\sigma_{\mathbf{c}}^{2}\) & \(L_{\mathbf{c c}}^{S \mathrm{Sc}}(r)\) & \(I_{\mathbf{c c}}^{S \mathrm{ch}}(r)\) \\
\hline
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 3: Compositional test function specifications with \(\mu_{\mathbf{c}}^{2}=\omega \sum_{j=1}^{\tilde{D}} \mu_{j}^{\psi} \cdot \mu_{j}^{\psi}\) and \(\sigma_{\mathbf{c}}^{2}= \omega \sum_{j=1}^{\tilde{D}} \sigma_{j j}^{\psi}\) with \(\omega=1 / 2 D\) if a transformation into logratios is applied and and \(\omega=1\) under ilr and clr transformations, cen(c) is the center of (10), and cen(c)(r) is the conditional center computed over all pairs of points at a distance \(\|\mathbf{r}\|=r\).}
\end{table}

Table 2 and Table 3 are related as follows. The compositional conditional expectation of the inner product of marks \(\tau_{\mathbf{c c}}(r)\) can be constructed by using \(\mathfrak{t}_{1}^{\mathbb{S}}\), which is related to the componentwise test function \(\mathfrak{t}_{1}^{\psi, j l}\) through
\[
\begin{array}{ccl}
\mathfrak{t}_{1}^{\mathrm{s}} \quad\langle\mathbf{c}(\circ), \mathbf{c}(\mathbf{r})\rangle_{A}= & \langle\operatorname{clr}(\mathbf{c}(\circ)), \operatorname{clr}(\mathbf{c}(\mathbf{r}))\rangle_{E} \\
=\sum_{j=1}^{D} \operatorname{clr}_{j}(\mathbf{c}(\circ)) \operatorname{clr}_{j}(\mathbf{c}(\mathbf{r}))= & \sum_{j=1}^{D} \mathfrak{t}_{1}^{\mathrm{clr}, j j} \\
=\sum_{j=1}^{D-1} \operatorname{ilr}_{j}(\mathbf{c}(\circ)) \operatorname{ilr}_{j}(\mathbf{c}(\mathbf{r}))= & \sum_{j=1}^{D-1} \mathfrak{t}_{1}^{\mathrm{ilr}, j j} \\
\stackrel{B}{2 D} \sum_{j_{1}} \sum_{j_{2}} \log \left(\frac{c_{j_{1}}(\circ)}{c_{j_{2}}(\circ)}\right) \log \left(\frac{c_{j_{1}}(\mathbf{r})}{c_{j_{2}}(\mathbf{r})}\right)= & \frac{1}{2 D} \sum_{j=1}^{D^{2}} \mathfrak{t}_{1}^{\mathrm{lr}, j j}
\end{array}
\]
with \(\langle\cdot, \cdot\rangle_{E}\) denoting the Euclidean inner product. Thus
\[
\tau_{\mathbf{c c}}(r)=\sum_{j=1}^{D} \tau_{j j}^{\mathrm{clr}}(r)=\sum_{j=1}^{D-1} \tau_{j j}^{i l r}(r)=\frac{1}{2 D} \sum_{j=1}^{D^{2}} \tau_{j j}^{\mathrm{lr}}(r)
\]
which allows to evaluate which individual components contribute to the overall mark characteristic. The corresponding compositional mark correlation function \(\kappa_{\mathbf{c c}}\) follows
by normalising \(\tau_{\mathbf{c c}}\) by \(\mu_{\mathbf{c}}^{2}\),
\[
\begin{aligned}
\mu_{\mathbf{c}}^{2} & =\int_{\mathbb{S}^{D}} \int_{\mathbb{S}^{D}}\langle\mathbf{c}(\circ), \mathbf{c}(\mathbf{r})\rangle_{A} \varpi(\mathrm{~d} \mathbf{c}(\circ)) \varpi(\mathrm{d} \mathbf{c}(\mathbf{r})) \\
& =\omega \int_{\mathbb{R}^{\tilde{D}}} \int_{\mathbb{R}^{\tilde{D}}}\langle\psi(\mathbf{c}(\circ)), \psi(\mathbf{c}(\mathbf{r}))\rangle_{E} \varpi(\mathrm{~d} \psi(\mathbf{c}(\circ))) \varpi(\mathrm{d} \psi(\mathbf{c}(\mathbf{r}))) \\
& =\omega \sum_{j=1}^{\tilde{D}} \int_{\mathbb{R}^{\tilde{D}}} \int_{\mathbb{R}^{\tilde{D}}} \psi_{j}(\mathbf{c}(\circ)), \psi_{j}(\mathbf{c}(\mathbf{r})) \varpi\left(\mathrm{d} \psi_{j}(\mathbf{c}(\circ))\right) \varpi\left(\mathrm{d} \psi_{j}(\mathbf{c}(\mathbf{r}))\right) \\
& =\omega \sum_{j=1}^{\tilde{D}} \mu_{j}^{\psi} \cdot \mu_{j}^{\psi}
\end{aligned}
\]
the limiting case which follows by substituting \(\mathfrak{t}_{1}^{\mathbb{S}}\) for \(\mathfrak{t}_{f}^{\mathbb{S}}\) in (22) with \(\omega=1 / 2 D\) under the lr-transformation and \(\omega=1\) if a transformation into ilr or clr coordinates is applied. Hence, for independent composition-valued marks, \(\kappa_{\mathbf{c c}}\) becomes constant one.

Selection of \(\mathfrak{t}_{4}^{\mathbb{S}}\) results in a compositional mark variogram, \(\gamma_{\mathbf{c c}}(r)\),
\[
\gamma_{\mathbf{c c}}(r)=\mathbb{E}_{\circ, r}\left[\frac{1}{2}\|\mathbf{c}(\circ) \ominus \mathbf{c}(\mathbf{r})\|_{A}^{2}\right]=\mathbb{E}_{\circ, r}\left[\frac{1}{2} d_{A}(\mathbf{c}(\circ), \mathbf{c}(\mathbf{r}))^{2}\right]
\]

By the decomposition of the Aitchison squared distance into the sum of the squared distances between the clr or ilr transformed marks as stated in (9), the compositional mark variogram can be decomposed into the sum of the componentwise mark variogram terms. In particular, we have that
\[
\gamma_{\mathbf{c c}}(r)=\sum_{j=1}^{D} \gamma_{j j}^{\mathrm{clr}}(r)=\sum_{j=1}^{D-1} \gamma_{j j}^{\mathrm{ilr}}(r)
\]

This functional quantity describes the average dispersion between the \(D\)-part compositionvalued marks for a focal and a second point at distance \(r\). If composition-valued marks are correlated for small distances, the compositional mark variogram \(\gamma_{\mathbf{c c}}(r)\) will show a corresponding decrease for small spatial distances. The decomposition allows to evaluate which mark components contribute how to the overall compositional mark variogram
\(\gamma_{\mathbf{c c}}(r)\). Again using (22), the mark variogram under independent marks is equal to \(\sigma_{\mathbf{c}}^{\mathbb{S}}\),
\[
\begin{aligned}
\sigma_{\mathbf{c}}^{\mathbb{S}} & =0.5 \int_{\mathbb{S}^{D}} \int_{\mathbb{S}^{D}}\|\mathbf{c}(\circ) \ominus \mathbf{c}(\mathbf{r})\|_{A}^{2} \varpi(\mathrm{~d} \mathbf{c}(\circ)) \varpi(\mathrm{d} \mathbf{c}(\mathbf{r})) \\
& =0.5 \omega \int_{\mathbb{R}^{\tilde{D}}} \int_{\mathbb{R}^{\tilde{D}}}\|\psi(\mathbf{c}(\circ)) \ominus \psi(\mathbf{c}(\mathbf{r}))\|_{E}^{2} \varpi(\mathrm{~d} \psi(\mathbf{c}(\circ))) \varpi(\mathrm{d} \psi(\mathbf{c}(\mathbf{r}))) \\
& =0.5 \omega \sum_{j=1}^{\tilde{D}} \int_{\mathbb{R}^{\tilde{D}}} \int_{\mathbb{R}^{\tilde{D}}}\left(\psi_{j}(\mathbf{c}(\circ))-\psi_{j}(\mathbf{c}(\mathbf{r}))\right)^{2} \varpi\left(\mathrm{~d} \psi_{j}(\mathbf{c}(\circ))\right) \varpi\left(\mathrm{d} \psi_{j}(\mathbf{c}(\mathbf{r}))\right) \\
& =\omega \sum_{j=1}^{\tilde{D}} \zeta_{j \dot{j}}^{\psi}
\end{aligned}
\]

Test functions \(\mathfrak{t}_{5}^{\mathbb{S}}\) and \(\mathfrak{t}_{6}^{\mathbb{S}}\) yield adaptations of Schlather's and Shimatani's \(I\) functions (denoted by \(\iota_{\mathbf{c c}}^{\mathrm{Sch}}\) and \(I_{\mathbf{c c}}^{\mathrm{Sch}}\) and \(\iota_{\mathbf{c c}}^{\mathrm{Shi}}\) and \(I_{\mathbf{c c}}^{\mathrm{Shi}}\) ), respectively, which can be useful tools to investigate the spatial autocorrelation of the composition-valued marks. Writing \(\overline{\mathbf{c}}=\mathbf{c}-\operatorname{cen}(\mathbf{c})\) to denote the centered composition, \(\iota_{\mathbf{c c}}^{\text {Shi }}\) can be decomposed into the weighted sum over the componentwise Shimantani's functions \(\iota_{\mathbf{c c}}^{\psi, j l}\) analogously to 24 , using the equivalences of the Aitchinson inner product of (23). The unnormalised \(\iota_{\mathbf{c c}}^{\mathrm{Sch}}\) function of Schlather can be obtained from its componentwise versions in a similar manner by applying a centering of the composition by the conditional center cen(c)(r). The decomposition of \(\iota_{\mathbf{c c}}^{\mathrm{Sch}}(r)\) into individual contributions is also analogous to that for \(\gamma_{\mathbf{c c}}(r)\).

\subsection*{3.6 Extensions to composition-valued marks with total information}

Finally, extensions of the proposed mark characteristics to combinations of relative, i.e. composition-valued, and absolute point-specific information are outlined. Adapting the results of Pawlowsky-Glahn et al. (2015) to the present context, consider \(\tilde{\mathbf{c}} \in \mathbb{R}^{D}\) such that \(\mathbf{c}=\operatorname{cls}(\tilde{\mathbf{c}}) \in \mathbb{S}^{D}\) and denote by \(y=\sum_{j=1}^{D} \tilde{c}_{j}\) the total. We then call \(\boldsymbol{\eta}=(y, \mathbf{c})\) a mixed mark with real-valued and composition-valued components \(y\) and \(\mathbf{c}\) living on \(\mathbb{T}=\mathbb{R}_{+} \times \mathbb{S}^{D}, \mathbb{R}_{+}=[0, \infty)\). Focusing on the process \(\left\{x_{i}, \boldsymbol{\eta}\left(x_{i}\right)\right\}\) instead of \(\left\{x_{i}, \mathbf{c}\left(x_{i}\right)\right\}\), all the above mark characteristics can be extended to mixed marks on \(\mathbb{T}\). To this end, consider first the scalars \(y, y^{\prime} \in \mathbb{R}_{+}\)and let \(\oplus_{+}\)and \(\oplus_{+}\)denote the plus-perturbation and plus-powering operations, respectively, where \(y \oplus_{+} y^{\prime}=y y^{\prime}\) and \(\xi \odot_{+} y=y^{\xi}\) for \(\xi \in \mathbb{R}\). Further, denote by \(\left\langle y, y^{\prime}\right\rangle_{+}=\left\langle\log (y), \log \left(y^{\prime}\right)\right\rangle_{E}\) the plus-inner product and by \(d_{+}\left(y, y^{\prime}\right)=\operatorname{abs}\left(\log (y)-\log \left(y^{\prime}\right)\right)\) the plus-distance on \(\mathbb{R}_{+}\). The above results can then be used to establish a vector space structure on \(\mathbb{T}\) with \(\mathbb{T}\)-perturbation \(\boldsymbol{\eta} \oplus_{\mathbb{T}} \boldsymbol{\eta}= \left(y \oplus_{+} y^{\prime}, \mathbf{c} \oplus \mathbf{c}^{\prime}\right)\) and \(\mathbb{T}\)-powering operations \(\xi \oplus_{\mathbb{T}} \boldsymbol{\eta}^{\prime}=\left(y^{\xi}, \xi \oplus \mathbf{c}\right)\) where \(\boldsymbol{\eta}, \boldsymbol{\eta}^{\prime} \in \mathbb{T}\) and \(\xi \in \mathbb{R}\) as before. Both variation and correlation related mark characteristics from Section 3.5
can be redefined through the \(\mathbb{T}\)-inner product and \(\mathbb{T}\)-squared distance
\[
\left\langle\boldsymbol{\eta}, \boldsymbol{\eta}^{\prime}\right\rangle_{\mathbb{T}}=\left\langle\mathbf{c}, \mathbf{c}^{\prime}\right\rangle_{A}+\beta\left\langle y, y^{\prime}\right\rangle_{+}
\]
and
\[
d_{\mathbb{T}}\left(\boldsymbol{\eta}, \boldsymbol{\eta}^{\prime}\right)=d_{A}\left(\mathbf{c}, \mathbf{c}^{\prime}\right)+\beta\left(\operatorname{abs}\left(\log (y)-\log \left(y^{\prime}\right)\right)\right),
\]
respectively, where \(\beta\) is a weight which can e.g. be chosen as one or as the ratio of the variances of the composition \(\mathbf{c}\) and of the total \(y\) (Happ and Greven, 2018). Using the above extensions and reformulating the compositional mark variogram for mixed composition- and real-valued marks, we obtain
\[
\gamma_{\mathbf{c c}, y}(r)=\mathbb{E}_{\circ, r}\left[\frac{1}{2} d_{A}(\mathbf{c}(\circ), \mathbf{c}(\mathbf{r}))^{2}+\beta\left(\frac{1}{2} d_{+}(y(\circ), y(\mathbf{r}))\right)^{2} .\right]
\]

The other summary characteristics of Table 3 can be extended analogously.

\subsection*{3.7 Estimation}

Recalling the representation of \(\nabla_{\mathrm{t}_{f}}^{\psi, j l}(r)\) as the ratio of the two second-order product density functions \(\varrho_{t_{f}}^{\psi,(2)}(r)\) and \(\varrho^{(2)}(r), \nabla_{t_{f}}^{\psi, j l}(r)\) can be estimated in close analogy to classic spatial point processes by
\[
\widehat{\nabla_{\mathrm{t}_{f}}^{\psi, j l}}(r)=\widehat{\varrho_{\mathrm{t}_{f}}^{\psi, j l,(2)}}(r) / \widehat{\varrho_{\mathrm{t}_{f}}^{(2)}(r)}
\]
where
\[
\widehat{\varrho_{t_{f}}^{\psi, j l,(2)}}(r)=\frac{1}{2 \pi r \nu(W)} \sum_{x_{1}, x_{2} \in W}^{\neq} \mathfrak{t}_{f}^{\psi, j l}\left(\psi_{j}\left(\mathbf{c}\left(x_{1}\right)\right), \psi_{l}\left(\mathbf{c}\left(x_{2}\right)\right)\right) \mathfrak{K}_{b}\left(\left\|x_{1}-x_{2}\right\|-r\right)
\]
and
\[
\widehat{\varrho^{(2)}}(r)=\frac{1}{2 \pi r \nu(W)} \sum_{x_{1}, x_{2} \in W}^{\neq} \mathfrak{K}_{b}\left(\left\|x_{1}-x_{2}\right\|-r\right)
\]

Here, \(\mathfrak{K}_{b}\) denotes a kernel function of bandwidth \(b\) and \(\nu(W)\) the area of the observation window \(W\). We note that as (30) and (31) are estimated using the same estimation principle, an edge correction factor can be ignored in both expressions (Illian et al., 2008). Similarly, an estimator of \(\kappa_{\mathfrak{t}_{f}}^{\psi, j l}\) of (12) can be obtained from normalizing \(\widehat{\nabla_{\mathfrak{t}_{f}}^{\psi, j l}}(r)\) by \(\widehat{\nabla_{\mathfrak{t}_{f}}^{\psi, j l}}\),
\[
\widehat{\kappa^{\psi, j l}}{\boldsymbol{t}_{f}}(r)=\widehat{\nabla_{\boldsymbol{t}_{f}}^{\psi, j l}}(r) / \widehat{\nabla_{\mathfrak{t}_{f}}^{\psi, j l}}
\]
where \(\widehat{\nabla_{\mathfrak{t}_{f}}^{\psi, j l}}\) in (32) can be estimated analogously to the scalar case from the transformed marks of the \(n\) points \(x_{1}, \ldots, x_{n}\) by
\[
\widehat{\nabla_{\mathfrak{t}_{f}}^{\psi, j l}}=\frac{1}{n^{2}} \sum_{i=1}^{n} \sum_{h=1}^{n} \mathfrak{t}_{f}^{\psi, j l}\left(\psi_{j}\left(\mathbf{c}\left(x_{i}\right)\right), \psi_{l}\left(\mathbf{c}\left(x_{h}\right)\right)\right)
\]
(Illian et al., 2008). For example, specifying \(\psi_{j}\) as the \(j\)-th component of the clr transformation of \(\mathbf{c}\), i.e. \(\log \left(c_{j}(\circ) / g(\mathbf{c})\right)\), an estimator of the clr mark variogram \(\gamma_{j l}^{\mathrm{clr}}\) of (18) can be obtained from the ratio of second-order density functions \(\widehat{\varrho_{\mathfrak{t}_{4}} \widehat{\mathrm{cl}, j l,(2)}}(r) / \widehat{\varrho_{\mathfrak{t}_{4}}^{(2)}(r)}\),
\[
\widehat{\gamma_{j l}^{\mathrm{clr}}}=\frac{\sum_{x_{1}, x_{2} \in W}^{\neq} 0.5\left(\log \left(\frac{c_{j}}{g(\mathbf{c})}\right)\left(x_{1}\right)-\log \left(\frac{c_{l}}{g(\mathbf{c})}\right)\left(x_{2}\right)\right)^{2} \mathfrak{K}_{b}\left(\left\|x_{1}-x_{2}\right\|-r\right)}{\sum_{x_{1}, x_{2} \in W}^{\neq} \mathfrak{K}_{b}\left(\left\|x_{1}-x_{2}\right\|-r\right)}
\]
as constant terms cancel. Likewise, considering a transformation into logratios, an estimator of the mean pairwise product of mark \(\log\)-ratios \(\tau_{j j}^{\operatorname{lr}}(r)\) of (14) is obtained as
\[
\widehat{\tau_{j j}^{\operatorname{lr}}}=\frac{\sum_{x_{1}, x_{2} \in W}^{\neq}\left(\log \left(\frac{c_{j_{1}}}{c_{j_{2}}}\right)\left(x_{1}\right) \cdot \log \left(\frac{c_{j_{1}}}{c_{j_{2}}}\right)\left(x_{2}\right)\right) \mathfrak{K}_{b}\left(\left\|x_{1}-x_{2}\right\|-r\right)}{\sum_{x_{1}, x_{2} \in W}^{\neq} \mathfrak{K}_{b}\left(\left\|x_{1}-x_{2}\right\|-r\right)}
\]

Considering the unmarked case first and applying the Campbell theorem Chiu et al., 2013) we have
\[
\mathbb{E}\left[\widehat{\varrho^{(2)}}(r)\right]=\int \mathfrak{K}_{b}(s) \varrho^{(2)}(r+b s) \mathrm{d} s
\]

Noting that \(\mathbb{E}\left[\widehat{\varrho^{(2)}}(r)\right] \rightarrow \varrho^{(2)}\) as \(b \rightarrow 0\) it follows that (31) is an unbiased estimator for \(b \rightarrow 0\) of the second-order product density function. Analogously, \(\widehat{\varrho_{\mathrm{t}_{f}}^{\psi, j l,(2)}}\) of (30) can be shown to be an unbiased estimator of \(\varrho_{\mathrm{t}_{f}}^{\psi, j l,(2)}\) for \(b \rightarrow 0\) by applying the Campbell theorem to the marked case (Daley and Vere-Jones, 2003) such that \(\mathbb{E}\left[\widehat{\varrho_{\mathrm{t}_{f}}^{\psi, j l,(2)}}(r)\right] \rightarrow \varrho_{\mathrm{t}_{f}}^{\psi, j l,(2)}\). As both (30) and (31) yield unbiased estimators, (29) yields a ratio-unbiased estimator for \(b \rightarrow 0\).

\section*{4 Test of random labeling hypothesis for compositionvalued marked point processes}

To test for deviations from the null hypothesis of random labels, i.e. marks that are i.i.d. and thus independent of each other and the points, we adopt global envelope tests.

These are non-parametric tests based on \(s\) simulations of the test statistic under the null model, originally introduced by Myllymäki et al. (2017) to solve multiple testing problems in spatial statistics. In the case of the random labeling hypothesis, the simulations can be obtained simply by permuting the marks of the points (e.g. Myllymäki et al., 2015). In the case of composition-valued marks, we take the same approach, i.e. permute the composition-valued marks.

Thus, first, the marks are permuted \(s\) times. The next step then is to compute the test statistic from the observed marked point pattern \(\left\{\left(x_{i}, \mathbf{c}\left(x_{i}\right)\right)\right\}_{i=1}^{n}\) and the \(s\) simulated patterns with permuted composition-valued marks. Let \(\vartheta_{1}(r)\) stand for the empirical functional test statistic computed from the observed pattern and let \(\vartheta_{2}(r), \ldots, \vartheta_{s+1}(r)\) be the \(s\) functional test statistics computed from the \(s\) simulated patterns. Then a Monte Carlo test is done based on \(\vartheta_{1}(r), \ldots, \vartheta_{s+1}(r)\). If the functional test statistics can be ordered from the least extreme to the most extreme, a Monte Carlo \(p\)-value can be computed for the test in a similar manner as in the classical Monte Carlo test (Barnard, 1963). Here, to order the statistics, we use the extreme rank length (ERL) measure (Myllymäki et al., 2017; Mrkvička et al., 2020) as a particular instance of rankbased measures which also allow for the graphical interpretation of the test in terms of a global envelope. Please refer to the publications citet above and Myllymäki and Mrkvička (2023, Appendix A) for the definition of ERL and a discussion of alternative rank measures.

More precisely, the idea of the global envelope test is the following: Let us denote by \(E_{i}, i=1, \ldots, s+1\), the measure associated with the \(i\)-th functional test statistic. Let \(\prec\) be an ordering for the measures \(E_{i}\) such that \(E_{i} \prec E_{j}\) whenever \(\vartheta_{i}\) is more extreme than \(\vartheta_{j}\) with respect to the measure \(E\). The critical value \(E_{(\alpha)}\) under a given significance level \(\alpha\) can then be found as the largest \(E_{i}\) which satisfies
\[
\sum_{i=1}^{s+1} \mathbb{1}\left(E_{i}<E_{(\alpha)}\right) \leqslant \alpha(s+1)
\]

Let us then denote by \(I_{(\alpha)}\) the index set of the test statistics \(\vartheta_{i}\) that are less or as extreme as \(E_{(\alpha)}\) as measured by their associated \(E_{i}\). Then the \(100(1-\alpha) \%\) global envelope is the band given by the two functions
\[
\vartheta_{(\alpha)}^{l}(r)=\min _{i \in I_{(\alpha)}} \vartheta_{i}(r)
\]
and
\[
\vartheta_{(\alpha)}^{u}(r)=\max _{i \in I_{(\alpha)}} \vartheta_{i}(r)
\]

If \(\vartheta_{1}(r)\) goes outside of the envelope \(\left(\vartheta_{(\alpha)}^{l}(r), \vartheta_{(\alpha)}^{u}(r)\right)\) for any of its argument values \(r\), there is evidence to reject the null hypothesis at the given significance level \(\alpha\). Further, the values of \(r\) where \(\vartheta_{1}(r)\) leaves the envelope show the reasons of the rejection of the
test. There is a one-to-one correspondence between the graphical interpretation by the global envelope and the Monte Carlo \(p\)-value of the test, given by
\[
p=\frac{1}{s+1}\left\{1+\sum_{i=2}^{s+1} \mathbf{1}\left(E_{i} \prec E_{1}\right)\right\}
\]

That is, assuming that there are no pointwise ties in \(\vartheta_{i}(r), i=1, \ldots, s+1\), with probability 1 , the empirical statistic \(\vartheta_{1}(r)\) leaves the envelope if and only if \(p \leqslant \alpha\), and \(\vartheta_{(\alpha)}^{l}(r) \leqslant \vartheta_{1}(r) \leqslant \vartheta_{(\alpha)}^{u}(r)\) if \(p>\alpha\) (Mrkvička et al., 2022, Theorem 1). The size of the test based on the \(p\)-value (33), and thus the global envelope, is \(\alpha\) when the test statistics \(\vartheta_{i}(r)\) can be strictly ordered and \(\alpha(s+1)\) is an integer (Myllymäki et al., 2017, Lemma 1).

In practice, the functional test statistics are estimators of functional summary characteristics computed on a chosen finite but dense set of argument values \(r\). If the random labeling hypothesis concerns all \(D\) parts of the composition-valued marks, the functional test statistics could be specified by \(\vartheta(r)=\widehat{\kappa_{\mathbf{t}_{f}}^{\mathrm{S}}}(r)\) for any of the functional test statistics in Table 3. Such a test can point out distances \(r\) which are responsible for the potential rejection of the test, but does not focus on which components in the composition-valued marks show the strongest dependence.

Alternatively, the functional test statistic can be constructed from the componentwise mark summary characteristics using the combining procedure of Myllymäki and Mrkvička (2023, Appendix B). For example, for the componentwise mark variograms, i.e. \(\gamma_{j j}^{\psi}(r), j=1, \ldots, \tilde{D}\), and a set of \(d r\)-values, the test vector is constructed from the estimated mark variograms \(\widehat{\gamma_{j j}^{\psi}}(r), j=1, \ldots, \tilde{D}\) :
\[
\boldsymbol{\vartheta}=\left(\left({\widehat{\gamma^{\psi}}}_{11}\left(r_{1}\right), \ldots,{\widehat{\gamma^{\psi}}}_{11}\left(r_{d}\right)\right), \ldots,\left({\widehat{\gamma^{\psi}}}_{\tilde{D} \tilde{D}}\left(r_{1}\right), \ldots,{\widehat{\gamma^{\psi}}}_{\tilde{D} \tilde{D}}\left(r_{d}\right)\right)\right)
\]

This test summarizes the information from all the components of the compositionvalued marks using the chosen componentwise characteristics. This test holds the global significance level for the complete test vector (34) and it can point out both the components \(j\) and distances \(r\) which are responsible for the potential rejection of the test.

\section*{5 Applications}

We use our new mark characteristics to analyse two marked point patterns from forestry and urban economics. In particular, focusing on our extensions of three prominent characteristics from the literature with each of these addressing a different aspect of the marks, we utilized the mark variogram, the conditional mean product of the marks and Shimantani's I function to investigate the variation, association and autocorrelation of

\begin{table}
\begin{tabular}{lrrrrrr}
\hline Mark & Min & 1st quantile & Mean & Median & 3rd quantile & Max \\
\hline\(h_{b} / h\) & 0.11 & 0.30 & 0.39 & 0.40 & 0.48 & 0.81 \\
\(h_{c r} / h\) & 0.19 & 0.52 & 0.61 & 0.60 & 0.70 & 0.89 \\
\hline\(h\) & 1.36 & 3.19 & 6.30 & 8.10 & 12.50 & 26.2 \\
\hline
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 4: Summary statistics for the relative height from the ground to the first living branch of the tree \(\left(h_{b} / h\right)\), the relative height of the crown \(\left(h_{c r} / h\right)\), and the total height \((h)\) in meters for the Finnish tree data.}
\end{table}
the composition-valued marks. The forestry data is an example of a spatial point pattern with 2-part composition-valued marks and additional absolute information (totals), while the urban economics data includes compositional marks with 4 parts.

\subsection*{5.1 Application to tree data with crown-to-base proportion}

The data on forest tree stands under study originates from a forest development study of managed, uneven-aged Norway spruce forests conducted under the ERIKA research project at the Natural Resources Institute Finland (Luke) (Eerikäinen et al., 2007, 2014; Saksa and Valkonen, 2011). Several plots of size \(40 \mathrm{~m} \times 40 \mathrm{~m}\) were recorded in southern Finland. Here we analyse a plot with 349 trees located in Vesijako (see Figure 1). The tree locations and associated tree characteristics, including the total height of the tree, \(h\), and the height from the ground to the first living branch (of the crown), \(h_{b}\), were recorded for all trees with \(h \geqslant 1.3 \mathrm{~m}\), which corresponds to the height at which the diameter at breast height (DBH) of a tree is measured. The height of the crown was obtained by \(h_{c r}=h-h_{b}\). In what follows, we call \(h_{b}\) 'base' for short. Instead of the absolute values of \(h_{c r}\) and \(h_{b}\) (which clearly depend on the individual age of the trees), we considered the relative base height \(h_{b} / h\) (with geometric mean 0.37 ), relative crown height \(h_{c r} / h\) (with geometric mean 0.58 ) and the total height \(h\). A summary of these three marks (two relational, one absolute) is provided in Table 4. The total tree heights varied up to 26.2 meters with an average crown proportion of \(61 \%\). The spatial distribution of the marks and the corresponding ilr coordinates are depicted in Figure 1. Crown proportions tended to be rather large, with only some small values (see top right panel of Figure 1). Compared to the crown-to-base composition, there is greater variation in the tree heights (bottom left). The ilr coordinates (bottom right) corresponding to the log-ratio transformation of the crown-to-base compositions support the above impressions: Most of the ilr coordinates are positive (indicated in red), which occurs when crown proportions are larger than base proportions.

Initially we considered only the crown-to-base ratios through their ilr coordinates and conducted a separate analysis where the total height was examined at its original scale. Next, to extend the composition-valued mark by the absolute height information, we additionally included the log transformation of the total heights in a vector-valued

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-24.jpg?height=875&width=875&top_left_y=391&top_left_x=627}
\captionsetup{labelformat=empty}
\caption{Figure 1: Spatial distribution of the marks and ilr-transformed composition-valued marks for the Finnish forest stand data. Base (top left) and crown (top right) proportions and absolute height in metres (bottom left) per tree with the diameter of the discs proportional to the mark values. The ilr-transformed crown-to-base proportions (bottom right) with negative values shown in blue and positive values highlighted in red.}
\end{figure}
mark in our computations according to Section 3.6. For both steps of our analysis and each mark characteristic, we computed the \(95 \%\) global envelope under the random labeling hypothesis (see Section 4 for details), based on 3000 permutations. While the interplay of tree crowns and total heights, and competition between neighbouring trees have been of interest in different studies (Hegyi, 1974; Gavrikov et al., 1993; Hui et al., 2018, Pitkänen et al., 2022), we are not aware of studies on the interdependencies of crown-to-base proportions and, additionally, the total heights over space. In particular, different from the existing approaches, the proposed rescaling of the absolute information into relative proportions allows for the characterisation of the structural properties of the marks without being affected by any heterogeneity in the age or types of the trees under study.

The top row of Figure 2 shows the compositional mark variogram \(\gamma_{\mathbf{c c}}\) (left), the conditional mean product of marks \(\tau_{\mathbf{c c}}\) (central) and Shimantani's \(\iota_{\mathbf{c c}}\) (right) together with \(95 \%\) global envelopes. We note that for \(D=2\) as in this application, all three compositional characteristics coincide with their componentwise analogues using the
ilr transformation. While the first two empirical characteristics for the crown-to-base composition leave the global envelope for some distances \(r\), the third one is completely inside the global envelope. The mark variogram suggests that the average dispersion of the crown-to-base log-ratios is smaller than expected under the random labeling hypothesis for any pair of points at interpoint distances of about \(r=6 \mathrm{~m}\). This would correspond to above random similarity in the crown-to-base log-ratios for neighbouring trees at intermediate distances. The results for \(\tau_{\mathbf{c c}}\) suggest that the product of the

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-25.jpg?height=951&width=1503&top_left_y=818&top_left_x=311}
\captionsetup{labelformat=empty}
\caption{Figure 2: Selected compositional mark summary characteristics and \(95 \%\) global envelopes (shaded) based on 3000 permutations of marks computed from the ilr transformed crown-to-base proportions (top row) and the total height (bottom row): mark variogram \(\gamma_{\mathbf{c c}}\) (top left), the conditional mean product of marks \(\tau_{\mathbf{c c}}\) (top central), Shimantani's \(\iota_{\mathbf{c c}}^{\text {Shi }}\) (top right), mark variogram \(\gamma_{y}\) (bottom left), the conditional mean scalar product of marks \(\tau_{y}\) (bottom central) and Shimantani's \(\iota_{y}^{\text {Shi }}\) (bottom right). The dashed line corresponds to the mean function under the random labeling hypothesis and the solid curve to the test function on the observed data. Distances are given in meters.}
\end{figure}
crown-to-base log-ratios of two points at distance \(r \approx 1 \mathrm{~m}\) apart from each other tend to be smaller than expected under the random labeling hypothesis. Recall that small values of \(\tau_{\mathbf{c c}}(r)\) occur for distance \(r\) if the transformed tree compositions for any two points at a distance \(r\) are more different, e.g. trees with large crown proportions
(i.e. positive coordinates) are surrounded by trees with small crown proportions (i.e. negative coordinates). This finding might be explained by crown competition and growth restrictions due to space limitations for closely neighbouring trees. Reinspecting the ilr scores of Figure 1, positive values (red, i.e. large crown to small base ratio) occur in close distance to negative ilr coordinates (blue, i.e. small crown to large base ratio) which could explain the results.

Comparing the findings with the results for the total information depicted in the bottom row of Figure 2, both the conditional mean product of marks \(\tau_{y}\) and Shimantani's \(\iota_{y}^{\text {Shi }}\) show clear negative deviations from the global envelopes for distances \(r<1.25\) m . These findings indicate that on average large trees are surrounded by smaller trees at shorter distances and vice versa implying negative autocorrelation potentially due to competition. On the other hand, the mark variogram \(\gamma_{y}\) is completely within the \(95 \%\) global envelope, even though it is rather close to the lower boundary for small \(r\). Thus, no significance dependence was detected using the squared difference of the tree heights as the test function.

Additionally, specifying the weight \(\beta\) introduced in Section 3.6 as the ratio of the variances for the ilr-transformed composition and the log-transformed totals, we computed all three mark characteristics following the concepts outlined in Section 3.6. Accounting for the total information in the analysis of the crown-to-base composition, all three characteristics are completely covered within the global envelopes. This result means that when taking both mark components \(\mathbf{c}\) and \(y\) jointly into account, the marks do not show any significant spatial dependence or autocorrelation.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-26.jpg?height=448&width=1481&top_left_y=1572&top_left_x=320}
\captionsetup{labelformat=empty}
\caption{Figure 3: Selected compositional mark summary characteristics and envelopes (shaded) based on 3000 simulations computed from the ilr transformed crown-to-base proportions and log transformed total height as mixed mark with \(\beta=0.57\). Mark variogram \(\gamma_{\mathbf{c c}, y}\) (left), the conditional mean inner product of marks \(\tau_{\mathbf{c c}, y}\) (central) and Shimantani's \(\iota_{\mathbf{c c}, y}^{\mathrm{Shi}}\) (right). The dashed line corresponds to the mean function under the random labeling hypothesis and the solid curve to the test function on the observed data. Distances are given in meters.}
\end{figure}

\subsection*{5.2 Application to Spanish municipalities data with local business sector compositions}

As a second application, we considered data from the National Statistics Institute of Spain (INE) on a de-composition of the local economy into four different sectors at municipality level. The data at hand was generated using a data query at the official webpage (www.ine.es) and derives directly from information collected in the Spanish central business register. The decomposition into the four distinct sectors was available for Spanish municipalities with at least 1,000 inhabitants and refers to the total number of all economic actors, e.g. companies, with location in a given municipality. To protect against inconsistencies in the data assignment and potential problems with multiple spatially wide spreading business locations, each company was matched by INE in a pre-processing step to exactly one municipality using the address information of the corporate headquarter. In subsequent steps, each company was categorised into one of the four business sectors
(a) industry (including extractive and manufacturing industries, energy and water supply, sanitation activities, waste management and decontamination),
(b) construction,
(c) commerce (including wholesale, retail trade, repair of motor vehicles and motorcycles, transport and storage, hostelry) and
(d) services (including communication, financial and insurance services, administrative activities and support services, education, health and social services, artistic, recreational and entertainment activities)
according to its main economic activity with reference date 1st January 2022.
Next, restricting the data to the Spanish provinces of Albacete, Cuenca, Cuidad Real and Toledo yielded a sample of 66 municipalities with complete information on the business sector decomposition. The four selected provinces are located southeast of Madrid and belong to the geographic region of La Mancha, the Spanish Plateau which is characterised by a homogeneous climate and strong similarity in its local population densities. Due to these characteristics, La Mancha has been used as a particular example of a homogeneous spatial point process in the literature (see e.g. Glass and Tobler, 1971; Ripley, 1977; Chiu et al., 2013). All collected information was then considered as a marked spatial point process by treating the attached coordinates as points and the closed four-part composition as point attribute. A visualisation of the point pattern at hand with the corresponding composition-valued marks shown as pie charts is depicted in Figure 4. The observed configuration of the marked points reflects a clear tendency of clustering for the point locations in combination with some variation over the individual pie charts, which indicate a clear predominance of the commerce

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-28.jpg?height=639&width=873&top_left_y=403&top_left_x=625}
\captionsetup{labelformat=empty}
\caption{Figure 4: Distribution of 4-part business sector composition on the Spanish Plateau for municipalities with at least 1000 inhabitants. Pie charts show the local decomposition of the economy into the business sector proportions of industry (blue), construction (green), commerce (orange) and service (red).}
\end{figure}
and service sectors within the four-part compositions. While the heterogeneity of the business sector decomposition seems to be maximal at larger interpoint distances, the composition seems to become more homogeneous for closely neighbouring points.

\begin{table}
\begin{tabular}{lrrrrrr}
\hline Sector (in \%) & Min & 1st quantile & Mean & Median & 3rd quantile & Max \\
\hline Industry & 2.81 & 7.35 & 8.83 & 9.56 & 11.72 & 19.78 \\
Construction & 6.69 & 11.96 & 15.61 & 15.90 & 18.49 & 29.30 \\
Commerce & 31.15 & 38.28 & 41.73 & 42.02 & 44.96 & 56.07 \\
Service & 18.18 & 27.11 & 31.47 & 32.52 & 37.08 & 56.72 \\
\hline
\end{tabular}
\captionsetup{labelformat=empty}
\caption{Table 5: Summary statistics for the closed 4 -part business sector composition computed from the Spanish business sector data for 66 municipalities with at least 1000 inhabitants on the Spanish Plateau.}
\end{table}

This observed variation of the marks is also supported by the numerical summary statistics of the business sector proportions reported in Table 5, which again reflect a clear predominance of the commerce and service sectors contrasted with only smaller proportions of the industry and construction sectors. The geometric means of the four parts highlight clear differences between the sectors industry ( 0.09 ), construction ( 0.16 ), commerce (0.43) and services (0.32).

For the composition-valued marks, we computed the same three compositional mark summary characteristics with global envelopes as before to investigate the spatial joint

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-29.jpg?height=440&width=1468&top_left_y=433&top_left_x=326}
\captionsetup{labelformat=empty}
\caption{Figure 5: Compositional mark summary characteristics computed from the 4-part economic sector compositions: Mark variogram \(\gamma_{\mathbf{c c}}\) (left), the conditional mean of the product of marks \(\tau_{\mathbf{c c}}\) (central) and Shimantani's \(\iota_{\mathbf{c c}}^{\text {Shi }}\) (right). The grey areas are \(95 \%\) global envelopes constructed from 3000 simulations under the random labeling hypothesis. The dashed line corresponds to the mean function under the random labeling hypothesis and the solid curve to the test function on the observed data.}
\end{figure}
variation, association and autocorrelation of the complete 4 -part composition (see Figure 5). While global envelopes for the mark variogram \(\gamma_{\mathbf{c c}}\) (left) and Shimantani's \(\iota_{\mathbf{c c}}^{\mathrm{Shi}}\) (right) suggest deviations from the independent mark assumption for some distances \(r\), the conditional mean of the inner product of marks \(\tau_{\mathbf{c c}}\) (central) is completely covered by the global envelopes. The mark variogram \(\gamma_{\mathbf{c c}}\) indicates that the average dispersion between the transformed 4 -part compositions of any pairs of points is greater than expected under the independent mark assumption at distances of around 0.3 units. For distances \(r \leqslant 0.2\) units, the mark variogram is smaller than expected under the random labeling hypothesis, although the empirical function stays within the envelope. This suggests similarity of the marks for pairs of nearby municipalities (although not significant), with increasing variability as the distances between point locations become larger, at least until about 0.3 units, where there is the most data. This might be explained by a strong variation in and clustering of the contribution of sectors such as e.g. tourism or banking to the local economy which, in turn, would affect the relative size of the service and commerce sectors. For Shimantani's \(\iota_{\mathbf{c c}}^{\text {Shi }}\) the findings suggest positive conditional spatial autocorrelation of the compositions at any nearby points with \(r<0.15\) units, and negative autocorrelation for \(r \approx 0.3\) units, consistent with the results of the variogram.

To investigate the effect of each of the four parts on the above results, we additionally computed the componentwise clr mark variograms \(\gamma_{j j}^{\mathrm{clr}}\) and componentwise Shimantani's \(c_{j j}^{\text {clr,Shi }}\) functions and plotted the results. Of all four \(\gamma_{j j}^{\text {clr }}\) functions depicted in Figure 6, only the mark variogram of clr(services) highlights deviations from the independent mark setting. This would suggest that the service sector proportions are

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-30.jpg?height=923&width=970&top_left_y=433&top_left_x=575}
\captionsetup{labelformat=empty}
\caption{Figure 6: Componentwise mark variograms \(\gamma_{j j}^{\text {clr }}\) computed from the clr transformed 4 -part economic sector composition. The grey bands represent \(95 \%\) global envelopes constructed from 3000 simulations under the random labeling hypothesis. The dashed line corresponds to the mean function under the random labeling hypothesis and the solid curve to the test function on the observed data.}
\end{figure}
more heterogeneous at larger distances which is consistent with the visual impression from Figure 4. By contrast, except for the clr-transformed construction sector, parts of all empirical \(\iota_{j j}^{\text {clr }}\),Shi functions are outside the \(95 \%\) envelopes for some distances \(r\) (see Figure 7). While we found positive autocorrelations for clr(commerce) and clr(services) at smaller distances \(r \approx 0.1\) units, both clr(industry) and clr(services) are below the \(95 \%\) envelopes at distances \(r \approx 0.3\) units corresponding to negative autocorrelation. This again would relate to similarity among the business sectors for closely neighbouring municipalities and an increasing heterogeneity with increasing distances, with some differences in (the strength of) this pattern between sectors.

\section*{6 Conclusion}

Combining methodological concepts for compositional data and spatial point processes, this paper introduces a novel class of composition-valued marked spatial point processes.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/06eca5fb-2b3b-4c6a-9645-429e929b34f7-31.jpg?height=920&width=970&top_left_y=436&top_left_x=575}
\captionsetup{labelformat=empty}
\caption{Figure 7: Componentwise Shimantani's \(\iota_{j j}^{\mathrm{clr}, \text { Shi }}\) and envelopes based on 3000 simulations computed from the clr transformed 4-part economic sector compositions. The dashed line corresponds to the mean function under the random labeling hypothesis and the solid curve to the test function on the observed data.}
\end{figure}

Our proposed set of different (functional) mark summary characteristics allows to decide on the mark independence assumption and investigate the pairwise dependencies of a new type of marked spatial point process. All developments are formalized through extended test functions, which generalise well-known interpretations to the present context. Transforming the composition-valued marks to the Euclidean space, the proposed tools can build on established methods for real-valued marks and can borrow strength from existing computational implementations.

Allowing for the characterisation of the spatial variation and association between both the complete composition-valued marks as well as their distinct compositional parts, the proposed extensions are helpful tools to highlight different aspects of the mark pattern. While the overall measures characterise the global interdependencies of the marks as a function of the interpoint distance, their componentwise counterparts provide useful insights into the contribution of the distinct parts to the overall results.

Apart from methods for purely composition-valued marks, we covered extensions to mixed marks including both composition-valued and absolute information. These
methods also allow to decompose vector-valued marks into the absolute and the relative information contained therein, highlighting patterns also in the relative information and avoiding that results for the vector-valued marks are mainly driven by the absolute information. For a combined analysis of both pieces of information, we introduced weights into the overall scalar product, which control for the differences in variation between both types of marks (Happ and Greven, 2018). We have here covered the case of spatial point processes with composition-valued marks, which can be seen as a special case of spatial point processes with more general object-valued marks.

Further extensions in this direction might also include alternative non-scalar marks such as density-valued or shape-valued marks. Note that as a by-product of our developed methods, we also showed how to handle vector-valued marks and derived both componentwise and (full-vector) compositional summary characteristics as well as their relationship, a result which is of independent interest in its own right.

\section*{Acknowledgements}

The authors gratefully acknowledge financial support through the German Research Foundation and Research Council of Finland. Matthias Eckardt and Sonja Greven were funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - project numbers 467634837 (Walter Benjamin grant for ME) and 459422098 (SG), respectively. Mari Myllymäki was financially supported by the Research Council of Finland (Grant numbers 295100, 327211) and the European Union - NextGenerationEU in the Research Council of Finland project (Grant number 348154) under flagship ecosystem for Forest-Human-Machine Interplay - Building Resilience, Redefining Value Networks and Enabling Meaningful Experiences (UNITE) (Grant numbers 337655 and 357909).

\title{
Supplement of On spatial point processes with compositionvalued marks
}

\section*{Transformations in the presence of structural zeros}

Where zeros are not structural (genuine) but e.g. due to small samples or rounding, they are commonly replaced or imputed from the data (see e.g. Martín-Fernández et al., 2003; Lubbe et al., 2021). By contrast, general approaches for structural zeros include amalgamations of the zero components into larger nonzero groups, separate treatment of zero and nonzero components, or restriction of the data to completely observed compositions only. Alternative approaches for structural zeros include representations of the data in linear or nonlinear geometric spaces through either combinations of logand power transformations or projections onto the sphere using square root transformations. In the linear geometry framework, early contributions include the folded power (Atkinson, 1985) and Box-Cox (Aitchison, 1986, Rayens and Srinivasan, 1991) transformations, which both tend to the alr-transformation if the power parameter \(p \rightarrow 0\). However, these transformations are problematic for some samples and certain properties of the composition are not well preserved (Barceló et al., 1996). Instead, the \(\alpha\)-transformation introduced by Tsagris et al. (2011) and further investigated by Tsagris (2015) and Tsagris et al. (2016) is defined through the map \(\alpha: \mathbb{S}^{D} \rightarrow \mathbb{R}^{D-1}\) with \(\alpha(\mathbf{c})=\mathbf{H}_{\mathrm{D}} \mathbf{u}_{\alpha}\) where \(\mathbf{u}_{\alpha}=\left(D \cdot\left(\operatorname{cls}(\alpha \odot \mathbf{c})-\mathbb{1}_{D}\right)\right) / \alpha\). The \(\alpha\)-transformation tends to ilr-coordinates when \(\alpha \rightarrow 0\) and a linear transformation if \(\alpha \rightarrow 1\). We note that similar ideas were also developed by Greenacre (2009a b) within the context of correspondence analysis. While the \(\alpha\)-transformation is well-defined for \(\alpha>0\), it maps the data into a codomain that is a subspace \(\mathbb{A}_{\alpha}^{D-1}\) of \(\mathbb{R}^{D-1}\) with \(\lim _{\alpha \rightarrow 0} \mathbb{A}_{\alpha}^{D-1} \rightarrow \mathbb{R}^{D-1}\) given by
\[
\mathbb{A}^{D-1}=\left\{\mathbf{H}_{\mathrm{D}} \mathbf{u}_{\alpha} \left\lvert\,-\frac{1}{\alpha} \leqslant u_{j, \alpha} \leqslant \frac{D-1}{\alpha}\right., \sum_{j=1}^{D} u_{j, \alpha}=0\right\}
\]
(Tsagris and Stewart, 2020). To overcome this limitation, Clarotto et al. (2022) proposed to use a centered \(\alpha\)-transformation (c \(\alpha \mathrm{t}\) ),
\[
\operatorname{c\alpha t}(\mathbf{c})=\left(\alpha^{-1}\left(c_{1}^{\alpha}-\frac{1}{D} \sum_{j=1}^{D} c_{j}^{\alpha}\right), \cdots, \alpha^{-1}\left(c_{D}^{\alpha}-\frac{1}{D} \sum_{j=1}^{D} c_{j}^{\alpha}\right)\right)
\]
with sum-to-zero constraint or an isometric \(\alpha\)-transformation (i \(\alpha \mathrm{t}\) ) where \(\mathrm{i} \alpha \mathrm{t}(\mathbf{c})= \mathrm{c} \alpha \mathrm{t}(\mathbf{c}) \mathrm{H}_{\mathrm{D}}^{\top}\) and \(\operatorname{c} \alpha \mathrm{t}(\mathbf{c})=\alpha^{-1}\left(\mathbf{G}_{\mathrm{D}} \mathbf{c}^{\alpha}\right)\). Both, i \(\alpha \mathrm{t}\) and \(\mathrm{c} \alpha \mathrm{t}\) tend to ilr and clr, respectively, when \(\alpha \rightarrow 0\).

Different from the linear space formulation, some parts of the literature considered square root transformations, which project the data onto a ( \(D-1\) )-dimensional (hyper-)sphere. While any such transformation allows for structural zeros, it imposes
a nonlinear geometry and knowledge in directional data analysis techniques is required to interpret the results (Scealy and Welsh, 2011; Scealy et al., 2015; Wang et al., 2007).

\section*{Extensions to multitype point processes with two distinct compositionvalued marks}

While the above methods considered the analysis of univariate point processes with one composition-valued mark, we next discuss extensions to \(k\)-variate point processes \(\left(\mathbf{x}_{1}, \ldots, \mathbf{x}_{k}\right)\) with two distinct associated composition-valued marks \(\left(\mathbf{c}^{a}, \mathbf{c}^{b}\right)\) on \(\mathbb{R}^{2} \times \mathbb{S}^{D_{1}} \times \mathbb{S}^{D_{2}}\) such that each component \(\mathbf{x}_{k}=\left\{x_{i},\left(\mathbf{c}^{a}\left(x_{i}\right), \mathbf{c}^{b}\left(x_{i}\right)\right)\right\}_{i=1}^{n_{k}}\) is a set of \(n_{k}\) points with two compositional marks. Such extended characteristics might be useful tools to investigate e.g. the association of one composition-valued mark, say a base-to-trunk composition, for different tree species and also to explore the associations and distributional characteristics of different compositions over space, say a tree and a soil composition. Extending the above methods to this setting not only allows to highlight particular aspects of the spatial behaviour of one composition, but also the crossassociation or cross-variation between different components of a multitype point process and/or different compositions.

Extending (13) to two distinct compositions yields a generic cross-composition conditional expectation of the product of \(\psi\)-transformed marks \(\tau_{j l}^{\psi, a b}\),
\[
\tau_{j l}^{\psi, a b}(r)=\mathbb{E}_{\circ, r}\left[\psi_{j}\left(\mathbf{c}^{a}(\circ)\right) \psi_{l}\left(\mathbf{c}^{b}(\mathbf{r})\right)\right]
\]
where \(\psi_{j}\left(\mathbf{c}^{a}(\circ)\right)\) and \(\psi_{l}\left(\mathbf{c}^{b}(\mathbf{r})\right)\) are the \(j\)-th and \(l\)-th elements of the \(\psi\)-transformed compositions \(\mathbf{c}^{a}\) and \(\mathbf{c}^{b}\). Under independence of the two compositions at distance \(r\), (35) tends to the product of means \(\mu_{j}^{\psi, a} \cdot \mu_{l}^{\psi, b}\). Similarly, reformulation of (15) leads to a generic cross-composition mark variogram \(\gamma_{j l}^{\psi, a b}(r)\) as
\[
\gamma_{j l}^{\psi, a b}(r)=\mathbb{E}_{\circ, r}\left[0.5 \cdot\left(\psi_{j}\left(\mathbf{c}^{a}(\circ)\right)-\psi_{l}\left(\mathbf{c}^{b}(\mathbf{r})\right)\right)^{2}\right]
\]

In settings where there are additional integer-valued marks, e.g. the kind of tree (spruce, beech etc.), i.e. the process can be considered as multivariate, we can extend (35) above to integer-valued marks, say points of type \(p\) and \(q\), yielding a crosscomposition cross-type mean product of marks \(\tau_{j l, h w}^{a b, p q}\),
\[
\tau_{j l, p q}^{\psi, a b}(r)=\mathbb{E}_{\circ, r}\left[\psi_{j}\left(\mathbf{c}_{p}^{a}(\circ)\right) \psi_{l}\left(\mathbf{c}_{q}^{b}(\mathbf{r})\right]\right.
\]
where \(\psi_{j}\left(\mathbf{c}_{p}^{a}\right)\) and \(\psi_{j}\left(\mathbf{c}_{q}^{b}\right)\) are the \(j\)-th part, respectively \(l\)-th part, of the \(\psi\)-transformed compositions \(\mathbf{c}^{a}\) and \(\mathbf{c}^{b}\) for points of type \(p\) and \(q\), respectively.

\section*{References}

Aitchinson, J. (1983): "Principal component analysis of compositional data," Biometrika, 70, 57-65.

Aitchinson, J. and S. Shen (1980): "Logistic-normal distributions: Some properties and uses," Biometrika, 67, 261-272.

Aitchison, J. (1986): The Statistical Analysis of Compositional Data, GBR: Chapman \& Hall.
- (2001): "Simplicial inference," in Algebraic methods in statistics and probability (Notre Dame, IN, 2000), Amer. Math. Soc., Providence, RI, vol. 287 of Contemp. Math., 1-22.

Atkinson, A. C. (1985): Plots, Transformations, and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis, Oxford: Clarendon Press.

Baddeley, A. (2010): Handbook of Spatial Statistics, CRC Press, chap. Multivariate and Marked Point Processes, 371-402, Chapman \& Hall/CRC Handbooks of Modern Statistical Methods.

Barceló, C., V. Pawlowsky, and E. Grunsky (1996): "Some aspects of transformations of compositional data and the identification of outliers," Mathematical Geology, 28, 501-518.

Barnard, G. A. (1963): "Discussion on 'The spectral analysis of point processes' (by M. S. Bartlett)," Journal of the Royal Statistical Society: Series B, 25, 264-296.

Billheimer, D., P. Guttorp, and W. F. Fagan (2001): "Statistical Interpretation of Species Composition," Journal of the American Statistical Association, 96, 12051214.

Brown, B. M. (1983): "Statistical Uses of the Spatial Median," Journal of the Royal Statistical Society. Series B (Methodological), 45, 25-30.

Chayes, F. (1960): "On correlation between variables of constant sum," Journal of Geophysical Research (1896-1977), 65, 4185-4193.

Chiu, S. N., D. Stoyan, W. S. Kendall, and J. Mecke (2013): Stochastic Geometry and Its Applications, John Wiley \& Sons, third ed.

Clarotto, L., D. Allard, and A. Menafoglio (2022): "A new class of \(\alpha\) transformations for the spatial analysis of Compositional Data," Spatial Statistics, 47, 100570.

Comas, C., P. Delicado, and J. Mateu (2008): "Analysing spatial point pat- terns with associated functional data," in Statistics for Spatio-temporal Modelling. Proceedings of the 4th International Workshop on Spatio-temporal Modelling (METMA-4), ed. by 157-163.
- (2011): "A second order approach to analyse spatial point patterns with functional marks," TEST, 20, 503-523.

Comas, C., L. Mehtätalo, and J. Miina (2013): "Analysing space-time tree interdependencies based on individual tree growth functions," Stoch. Environ. Res. Risk. A., 27, 1673-1681.

Cressie, N. (1993): Statistics for Spatial Data, Wiley.
Daley, D. and D. Vere-Jones (2003): An Introduction to the theory of point processes. Volume I, Springer, Berlin-Heidelberg.

Diggle, P. J. (1986): "Displaced amacrine cells in the retina of a rabbit: analysis of a bivariate spatial point pattern," Journal of Neuroscience Methods, 18, 115-125.

Eckardt, M., C. Comas, and J. Mateu (2024): "Summary characteristics for multivariate function-valued spatial point process attributes," International Statistical Review, n/a, https://doi.org/10.1111/insr.12582.

Eckardt, M. and J. Mateu (2019a): "Analysing Multivariate Spatial Point Processes with Continuous Marks: A Graphical Modelling Approach," International Statistical Review, 87, 44-67.
- (2019b): "Partial characteristics for marked spatial point processes," Environmetrics, 30, e2565.

Eckardt, M. and M. Moradi (2024a): "Marked Spatial Point Processes: Current State and Extensions to Point Processes on Linear Networks," Journal of Agricultural, Biological and Environmental Statistics, 29, 346-378.
- (2024b): "Rejoinder on 'Marked Spatial Point Processes: Current State and Extensions to Point Processes on Linear Networks'," Journal of Agricultural, Biological and Environmental Statistics, 29, 405-416.

Eerikäinen, K., J. Miina, and S. Valkonen (2007): "Models for the regeneration establishment and the development of established seedlings in uneven-aged, Norway spruce dominated forest stands of southern Finland," Forest Ecology and Management, 242, 444-461.

Eerikäinen, K., S. Valkonen, and T. Saksa (2014): "Ingrowth, survival and height growth of small trees in uneven-aged Picea abies stands in southern Finland," Forest Ecosystems, 1, 5.

Egozcue, J. J. and V. Pawlowsky-Glahn (2005): "Groups of Parts and Their Balances in Compositional Data Analysis," Mathematical Geology, 37, 795-828.

Fišerová, E. and K. Hron (2011): "On the Interpretation of Orthonormal Coordinates for Compositional Data," Math. Geosci., 43, 455.

Gavrikov, V. L., P. Y. Grabarnik, and D. Stoyan (1993): "Trunk-Top Relations in a Siberian Pine Forest," Biometrical Journal, 35, 487-498.

Ghorbani, M., O. Cronie, J. Mateu, and J. Yu (2021): "Functional marked point processes: a natural structure to unify spatio-temporal frameworks and to analyse dependent functional data," TEST, 30, 529-568.

Glass, L. and W. R. Tobler (1971): "General: Uniform Distribution of Objects in a Homogeneous Field: Cities on a Plain," Nature, 233, 67-68.

Greenacre, M. (2009a): "Log-Ratio Analysis Is a Limiting Case of Correspondence Analysis," Mathematical Geosciences, 42, 129.
- (2009b): "Power transformations in correspondence analysis," Computational Statistics \& Data Analysis, 53, 3107-3116.

Guan, Y. (2006): "Tests for Independence between Marks and Points of a Marked Point Process," Biometrics, 62, 126-134.

Guan, Y. and D. R. Afshartous (2007): "Test for independence between marks and points of marked point processes: a subsampling approach," Environmental and Ecological Statistics, 14, 101-111.

Happ, C. and S. Greven (2018): "Multivariate Functional Principal Component Analysis for Data Observed on Different (Dimensional) Domains," Journal of the American Statistical Association, 113, 649-659.

Harkness, R. D. and V. Isham (1983): "A Bivariate Spatial Point Pattern of Ants' Nests," Journal of the Royal Statistical Society. Series C (Applied Statistics), 32, 293-303.

Hegyi, F. (1974): "A simulation model for managing jack-pine stands," in Growth Models for Tree and Stand Simulation, ed. by J. Fries, Royal College of Forestry, 74-90.

Ho, L. P. and D. Stoyan (2008): "Modelling marked point patterns by intensitymarked Cox processes," Statistics \& Probability Letters, 78, 1194-1199.

Hron, K., M. Engle, P. Filzmoser, and E. Fišerová (2021): "Weighted Symmetric Pivot Coordinates for Compositional Data with Geochemical Applications," Math. Geosci., 53, 655-674.

Hron, K., P. Filzmoser, P. de Caritat, E. Fišerová, and A. Gardlo (2017): "Weighted Pivot Coordinates for Compositional Data and Their Application to Geochemical Mapping," Math. Geosci., 49, 797-814.

Huang, T., G. Saporta, and H. Wang (2021): A Spatial Durbin Model for Compositional Data, Cham: Springer, 471-488.

Hui, G., Y. Wang, G. Zhang, Z. Zhao, C. Bai, and W. Liu (2018): "A novel approach for assessing the neighborhood competition in two different aged forests," Forest Ecology and Management, 422, 49-58.

Illian, J., A. Penttinen, H. Stoyan, and D. Stoyan (2008): Statistical Analysis and Modelling of Spatial Point Patterns, John Wiley \& Sons, New York.

Isham, V. (1985): "Marked point processes and their correlations," in Spatial processes and spatial time series analysis, Publications des Facultes Universitaires Saint-Louis Brussels, 63-75.

Kynčlová, P., K. Hron, and P. Filzmoser (2017): "Correlation Between Compositional Parts Based on Symmetric Balances," Mathematical Geosciences, 49, 777796.

Leininger, T. J., A. E. Gelfand, J. M. Allen, and J. A. Silander (2013): "Spatial Regression Modeling for Compositional Data With Many Zeros," Journal of Agricultural, Biological, and Environmental Statistics, 18, 314-334.

Lotwick, H. W. and B. W. Silverman (1982): "Methods for Analysing Spatial Processes of Several Types of Points," Journal of the Royal Statistical Society. Series B (Methodological), 44, 406-413.

Lubbe, S., P. Filzmoser, and M. Templ (2021): "Comparison of zero replacement strategies for compositional data with large numbers of zeros," Chemometrics and Intelligent Laboratory Systems, 210, 104248.

Martín-Fernández, J. A., C. Barceló-Vidal, and V. Pawlowsky-Glahn (2003): "Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation," Mathematical Geology, 35, 253-278.

Mateu-Figueras, G., V. Pawlowsky-Glahn, and J. Egozcue (2011): The Principle of Working on Coordinates, John Wiley \& Sons, Ltd, chap. 3, 29-42.

Moran, P. A. P. (1950): "Notes on continuous stochastic phenomena," Biometrika, 37, 17-23.

Mrkvička, T., M. Myllymäki, M. Jílek, and U. Hahn (2020): "A One-Way ANOVA Test for Functional Data with Graphical Interpretation," Kybernetika, 56, 432-458.

Mrkvička, T., M. Myllymäki, M. Kuronen, and N. N. Narisetty (2022): "New methods for multiple testing in permutation inference for the general linear model," Statistics in Medicine, 41, 276-297.

Myllymäki, M. and T. Mrkvička (2023): "GEt: Global envelopes in R," arXiv:1911.06583 [stat.ME].

Myllymäki, M., T. Mrkvička, P. Grabarnik, H. Seijo, and U. Hahn (2017): "Global envelope tests for spatial processes," Journal of the Royal Statistical Society: Series B, 79, 381-404.

Myllymäki, M., P. Grabarnik, H. Seijo, and D. Stoyan (2015): "Deviation test construction and power comparison for marked spatial point patterns," Spatial Statistics, 11, 19-34.

Pawlowsky, V. (1984): "On spurious spatial covariance between variables of constant sum," Sciences de la terre. Informatique géologique, 107-113.
- (1986): "Räumliche Strukturanalyse und Schätzung ortsabhängiger Kompositionen mit Anwendungsbeispielen aus der Geologie," Ph.D. thesis, Freie Universität Berlin, Fachbereich Geowissenschaft.

Pawlowsky-Glahn, V. and J. Egozcue (2001): "Geometric approach to statistical analysis on the simplex," Stochastic Environmental Research and Risk Assessment, 15, 384-398.

Pawlowsky-Glahn, V., J. J. Egozcue, and D. Lovell (2015): "Tools for compositional data with a total," Statistical Modelling, 15, 175-190.

Pawlowsky-Glahn, V. and A. O. Ricardo (2004): Geostatistical Analysis of Compositional Data, vol. 7 of IAMG studies in mathematical geology, United Kingdom: Oxford University Press.

Pawlowsky-Glahn, V. and A. Buccianti (2011): Compositional Data Analysis, John Wiley \& Sons, Ltd.

Penttinen, A. and D. Stoyan (1989): "Statistical Analysis for a Class of Line Segment Processes," Scandinavian Journal of Statistics, 16, 153-168.

Penttinen, A., D. Stoyan, and H. M. Henttonen (1992): "Marked point processes in forest statistics," Forest Science, 38, 806-824.

Pitkänen, T., S. Bianchi, and A. Kangas (2022): "Quantifying the effects of competition on the dimensions of Scots pine and Norway spruce crowns," International Journal of Applied Earth Observation and Geoinformation, 112, 102941.

Rayens, W. S. and C. Srinivasan (1991): "Box-Cox transformations in the analysis of compositional data," Journal of Chemometrics, 5, 227-239.

Ripley, B. D. (1976): "The second-order analysis of stationary point processes," Journal of Applied Probability, 13, 255-266.
- (1977): "Modelling Spatial Patterns," Journal of the Royal Statistical Society. Series B, 39, 172-212.

Saksa, T. and S. Valkonen (2011): "Dynamics of seedling establishment and survival in uneven-aged boreal forests," Forest Ecology and Management, 261, 14091414.

Scealy, J. L., P. de Caritat, E. C. Grunsky, M. T. Tsagris, and A. H. Welsh (2015): "Robust Principal Component Analysis for Power Transformed Compositional Data," Journal of the American Statistical Association, 110, 136-148.

Scealy, J. L. and A. H. Welsh (2011): "Regression for compositional data by using distributions defined on the hypersphere," Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73, 351-375.

Schlather, M. (2001): "On the second-order characteristics of marked point processes," Bernoulli, 7, 99-117.

Schlather, M., P. Riberio, and P. Diggle (2004): "Detecting Dependence between Marks and Locations of Marked Point Processes," Journal of the Royal Statistical Society, Series B (Methodological), 66, 79-93.

Sharp, W. E. (2006): "The graph median-A stable alternative measure of central tendency for compositional data sets," Mathematical Geology, 38, 221-229.

Shimatani, K. (2002): "Point Processes for Fine-Scale Spatial Genetics and Molecular Ecology," Biometrical Journal, 44, 325-352.

Stoyan, D. (1984a): "Correlations of the Marks of Marked Point Processes - Statistical Inference and Simple Models," Elektronische Informationsverarbeitung und Kybernetik, 20, 285-294.
- (1984b): "On Correlations of Marked Point Processes," Mathematische Nachrichten, 116, 197-207.
- (1987): "Statistical Analysis of Spatial Point Processes: A Soft-Core Model and Cross-Correlations of Marks," Biometrical Journal, 29, 971-980.

Stoyan, D. and H. Stoyan (1994): Fractals, Random Shapes, and Point Fields : Methods of Geometrical Statistics, Chichester, New York: Wiley.

Stoyan, D. and O. Wälder (2000): "On variograms in point process statistics, II: Models for markings and ecological interpretation," Biometrical Journal, 42, 171-187.

Tolosana Delgado, R. (2006): "Geostatistics for constrained variables: positive data, compositions and probabilities. Applications to environmental hazard monitoring," Ph.D. thesis, University of Girona.

Tsagris, M. (2015): "Regression analysis with compositional data containing zero values," Chilean Journal of Statistics, 6, 47-57.

Tsagris, M., S. Preston, and A. T. Wood (2011): A data-based power transformation for compositional data.

Tsagris, M., S. Preston, and A. T. A. Wood (2016): "Improved Classification for Compositional Data Using the \(\alpha\)-transformation," Journal of Classification, 33, 243-261.

Tsagris, M. and C. Stewart (2020): "A folded model for compositional data analysis," Australian \& New Zealand Journal of Statistics, 62, 249-277.

Van Lieshout, M. N. M. and A. J. Baddeley (1999): "Indices of Dependence Between Types in Multivariate Point Patterns," Scandinavian Journal of Statistics, 26, 511-532.
van Lieshout, M. N. M. v. (2006): "A J-Function for Marked Point Patterns," Annals of the Institute of Statistical Mathematics, 58, 235.

Wang, H., Q. Liu, H. M. Mok, L. Fu, and W. M. Tse (2007): "A hyperspherical transformation forecasting model for compositional data," European Journal of Operational Research, 179, 459-468.

Wiegand, T. and K. A. Moloney (2013): Handbook of Spatial Point-Pattern Analysis in Ecology, Chapman and Hall/CRC.

Zhang, T. and Q. Zhuang (2014): "On the local odds ratio between points and marks in marked point processes," Spatial Statistics, 9, 20-37.