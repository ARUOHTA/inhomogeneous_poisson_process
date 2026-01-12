\title{
The Statistical Analysis of Compositional Data
}

\author{
By J. Aitchison \\ University of Hong Kong \\ [Read before the Royal Statistical Society at a meeting organized by the Research Section on Wednesday, 13th January, 1982, Professor R. N. Curnow in the Chair]
}

\begin{abstract}
Summary
The simplex plays an important role as sample space in many practical situations where compositional data, in the form of proportions of some whole, require interpretation. It is argued that the statistical analysis of such data has proved difficult because of a lack both of concepts of independence and of rich enough parametric classes of distributions in the simplex. A variety of independence hypotheses are introduced and interrelated, and new classes of transformed-normal distributions in the simplex are provided as models within which the independence hypotheses can be tested through standard theory of parametric hypothesis testing. The new concepts and statistical methodology are illustrated by a number of applications.
\end{abstract}

\section*{1. Introduction}

There are many practical problems for which the positive simplex
\[
\mathbb{S}^{d}=\left\{\left(x_{1}, \ldots, x_{d}\right): x_{i}>0(i=1, \ldots, d), x_{1}+\ldots+x_{d}<1\right\},
\]
forms the whole, or a major component, of the sample space. For such problems, concepts of independence must often play an important role in any form of statistical analysis. The simplex, however, has proved to be an awkward space to handle statistically; the difficulties appear to lie in the scarcity of meaningful definitions of independence and of measures of dependence and in the absence of satisfactory parametric classes of distributions on \(\mathbb{S}^{d}\). It is the aim of this paper to introduce a number of concepts of independence in the simplex, to relate these to some existing concepts, and to develop within the framework of rich new parametric classes of distributions appropriate statistical methods of analysis.

To motivate all the concepts introduced and to provide illustrations of the statistical methodology developed we shall use data sets in two very different areas of application, geology and consumer demand analysis. We hope that the expert reader will see these examples for what they are, attempts at providing potential statistical insights into these and similar disciplines rather than presumptuous criticism by a novice of interpretations already placed on the particular data sets.

Geology. The geological literature abounds with problems of the interpretation of chemical, mineral and fossil compositions of rock and sediment specimens. Each composition of each specimen is a set of some three to twenty proportions summing to unity and so can be represented by a point in an appropriately dimensioned simplex. We concentrate on three published geological data sets chosen to illustrate, as simply as possible, various aspects of our analysis.

Example 1: Skye lavas. Thompson, Esson and Duncan (1972), in their Table 2, give the chemical compositions of 32 basalt specimens from the Isle of Skye in the form of percentages of 10 major oxides. A typical percentage vector in \(\mathbb{S}^{9}\) is thus

\begin{tabular}{llllllllll}
\(\mathrm{SiO}_{2}\) & \(\mathrm{Al}_{2} \mathrm{O}_{3}\) & \(\mathrm{Fe}_{2} \mathrm{O}_{3}\) & MgO & CaO & \(\mathrm{Na}_{2} \mathrm{O}\) & \(\mathrm{K}_{2} \mathrm{O}\) & \(\mathrm{TiO}_{2}\) & \(\mathrm{P}_{2} \mathrm{O}_{5}\) & MnO \\
46.31 & 14.18 & 12.32 & 12.74 & 9.62 & 2.51 & 0.34 & 1.53 & 0.16 & 0.18
\end{tabular}

For this example we shall discuss classes of parametric models for describing the experienced
pattern of variability, investigate the adequacy of such models and test a number of independence hypotheses for such sets of proportions.

Example 2: Glacial tills in North-Central New York. As part of a study of the composition of glacial till samples Kaiser (1962) presents, within his Table 1, the percentage compositions in terms of four pebble types, together with the total pebble count, of 93 till samples. Typical sample information thus takes the form

\begin{tabular}{ccccc} 
& \multicolumn{2}{c}{ Percentage composition } & & Total \\
Red sandstone & Grey sandstone & Crystalline & Miscellaneous & pebbles \\
67.2 & 31.5 & 0.3 & 1.0 & 387
\end{tabular}

In addition to the composition in \(\mathbb{S}^{3}\) we have here an abundance or size associated with each sample. Interest may then be in the extent, if any, to which composition depends on size.

Example 3: Arctic lake sediments. Coakley and Rust (1968) give, in their Table 1, the compositions in terms of sand, silt and clay percentages of 39 sediment samples at different water depths in an Arctic lake, with typical entry

\begin{tabular}{cccc}
\multicolumn{4}{r}{ Sediment composition in percentages } \\
Sand & Silt & Water \\
\(10 \cdot 5\) & \(55 \cdot 4\) & Clay & depth \((\mathrm{m})\) \\
& & \(34 \cdot 1\) & \(49 \cdot 4\)
\end{tabular}

Of interest here is the question of quantifying the extent to which water depth is explanatory of compositional pattern.

An appreciation of the difficulty imposed by this confinement of data points, such as the compositions in the above examples, to a simplex is inherent in the comments of Pearson (1897) on spurious correlations, and in geological circles the difficulty has since become known as the constant or bounded sum problem and the problem of closed arrays. As our analysis unfolds we shall cite various attempts to overcome this difficulty, and, in identifying reasons for limited success, we shall discover a means of overcoming most of the difficulties.

Consumer demand analysis. An important aspect of the study of consumer demand is the analysis of household budget surveys, in which attention focuses on expenditures on a number of mutually exclusive and exhaustive commodity groups and their relations to total expenditure, income, type of housing, household composition, and so on.

Example 4: Hong Kong household expenditure budgets. The set of household expenditure data available to us is from a pilot selection of 199 Hong Kong households, used as a preparatory study for a large-scale household expenditure survey by the Hong Kong Census and Statistics Department. From this set we have for simplicity selected subgroups of 41 and 42 households in two low-cost housing categories A and B. For each household information is available on number of persons, household composition, total household income, and monthly expenditures in nine commodity/service groups. The contents of these commodity groups are fully defined in the monthly Consumer Price Index Report of the Hong Kong Census and Statistics Department. To keep our illustrative analysis simple we have avoided the problem of zero components by combining two pairs of commodity groups to obtain the following seven: (1) housing, (2) fuel and light, (3) foodstuffs, (4) transport and vehicles, (5) tobacco, alcohol and miscellaneous goods, (6) services, (7) clothing, footwear and durable goods, and by replacing the few remaining zero expenditures in these groups by HK \(\$ 0.05\), half the lowest recordable expenditure.

In the investigation of such data the pattern or composition of expenditures, the proportions of total expenditure allocated to the commodity groups, can be shown to play a central role, and indeed some economists (Working, 1943; Leser, 1976; Deaton, 1978; Deaton and Muellbauer, 1980) have investigated such a budget-share approach. Since each pattern of expenditures is again represented by a point in the simplex, questions such as "To what extent
does the pattern of expenditure depend on the total amount spent?" and "Are there some commodity groups which are given priority in the allocation of expenditure?" obviously require adequate models to describe patterns of variability in the simplex and careful definitions of independence structure in the simplex for their satisfactory resolution.

\section*{2. Parametric Classes of Distributions on \(\mathbb{S}\) d}

\subsection*{2.1. Fundamental Operations on Compositions}

As a first step towards the introduction of new classes of distributions and independence concepts we establish a suitable terminology and notation for certain mathematical operations in the simplex which help in the study and manipulation of compositional data.

Spaces and vectors. Let \(\mathbb{R}^{d}\) denote \(d\)-dimensional real space, \(\mathbb{P}^{d}\) its positive orthant and \(\mathbb{S}^{d}\) its positive simplex (1.1). The symbols, \(\mathbf{w}, \mathbf{x}\) and \(\mathbf{y}\) are reserved for vectors in \(\mathbb{P}^{d}, \mathbb{S}^{d}\) and \(\mathbb{R}^{d}\), respectively, although we shall occasionally have to use other symbols for such vectors. Any vector or point \(\mathbf{x}\) in \(\mathbb{S}^{d}\) is termed a composition and any collection of such vectors, compositional data. We use the symbol \(x_{d+1}\) always in the sense
\[
x_{d+1}=1-x_{1}-\ldots-x_{d}
\]
to denote the fill-up value. The notation \(\mathbf{x}^{(c)}=\left(x_{1}, \ldots, x_{c}\right)\) allows focusing on leading subvectors with the dimension of the subvector indicated by the superscript. Thus \(\mathbf{x}^{(c)}\) with \(c<d\) is a subvector of \(\mathbf{x}\) or equivalently \(\mathbf{x}^{(d)}\), and \(\mathbf{x}^{(d+1)}\) is the augmented \(\mathbf{x}\) vector ( \(x_{1}, \ldots, x_{d}, x_{d+1}\) ). The subvector ( \(x_{c+1}, \ldots, x_{d+1}\) ) obtained by deletion of \(\mathbf{x}^{(c)}\) from \(\mathbf{x}^{(d+1)}\) is denoted by \(\mathbf{x}_{(c)}\). We use \(T\left(\mathbf{x}^{(c)}\right)\) to denote the sum \(x_{1}+\ldots+x_{c}\) of the elements of any vector or subvector, such as \(\mathbf{x}^{(c)}\).

Basis of a composition. In our household expenditure example the \(d\)-dimensional budgetshare composition \(\mathbf{x}^{(d+1)}\) is derived from the actual amounts spent \(\mathbf{w}^{(d+1)}\) on the \(d+1\) commodity groups through an operation \(C: \mathbb{P}^{d+1} \rightarrow \mathbb{S}^{d}\) defined by \(\mathbf{x}^{(d+1)}=C\left(\mathbf{w}^{(d+1)}\right)\) where \(x_{i}=w_{i} / T\left(\mathbf{w}^{(d+1)}\right)(i=1, \ldots, d+1)\). For convenience we term such a vector \(\mathbf{w}^{(d+1)} \in \mathbb{P}^{d+1}\), when it exists, the basis of the composition \(\mathbf{x}^{(d+1)}\).

Subcomposition. Often in the study of geochemical compositions attention is directed towards the relative proportions of a few oxides. For example, a popular diagrammatic representation treats the relative proportions
\[
\left(\mathrm{CaO}, \mathrm{Na}_{2} \mathrm{O}, \mathrm{~K}_{2} \mathrm{O}\right) /\left(\mathrm{CaO}+\mathrm{Na}_{2} \mathrm{O}+\mathrm{K}_{2} \mathrm{O}\right)
\]
in \(\mathbb{S}^{2}\) as triangular coordinates in a CNK ternary diagram. We can formalize this process of focusing on a subset of components as follows. Any subvector, such as \(\mathbf{x}^{(c)}\), of a composition \(\mathbf{x}^{(d+1)}\) can play the role of a basis in \(\mathbb{P}^{c}\) for a composition \(C\left(\mathbf{x}^{(c)}\right)\) in \(\mathbb{S}^{c-1}\). Such a composition is termed a subcomposition \(C\left(\mathbf{x}^{(c)}\right)\) of \(\mathbf{x}^{(d+1)}\).

Amalgamation. In a household expenditure enquiry there may be reasons for combining some commodity groups, to form new amalgamated groups. If we suppose that the composition has been ordered in such a way that combinations are between neighbouring components, the formal general process can be set out as follows. Let the integers \(c_{0}, \ldots . c_{k+1}\) satisfy
\[
0=c_{0}<c_{1}<\ldots<c_{k}<c_{k+1}=d+1
\]
and define
\[
t_{j}=x_{c_{j-1}+1}+\ldots+x_{c_{j}} \quad(j=1, \ldots, k+1) .
\]

Then \(\mathbf{t}^{(k+1)} \in \mathbb{S}^{k}\) and so is a \(k\)-dimensional composition which we term an amalgamation of \(\mathbf{x}^{(d+1)}\). It is obvious that the transformation from \(\mathbf{x}^{(d+1)}\) to \(\mathbf{t}^{(k+1)}\) can be represented by a matrix operation \(\mathbf{t}^{(k+1)}=\mathbf{A} \mathbf{x}^{(d+1)}\) from \(\mathbb{S}^{d}\) to \(\mathbb{S}^{k}\), where \(\mathbf{A}\) consists of 0 s and 1 s , with a single 1 in each column.

Partition. The amalgamation just discussed involves a separation of the vector \(\mathbf{x}^{(d+1)}\) into \(k+1\) subvectors. When considering such an amalgamation we may often be interested also in
the \(k+1\) subcompositions associated with these subvectors. The \(j\) th such subcomposition, \(\mathbf{s}_{j} \in \mathbb{S}^{d_{j}}\) where \(d_{j}=c_{j}-c_{j-1}-1\), has components
\[
s_{j r}=x_{c_{j-1}+r} / t_{j} \quad\left(r=1, \ldots, d_{j}+1\right)
\]
where the \(\left(d_{j}+1\right)\) th component is the fill-up value. An extremely useful feature is that the transformation from \(\mathbb{S}^{d}\) to
\[
\mathbb{S}^{k} \times \prod_{j=1}^{k+1} \mathbb{S}^{d_{j}}
\]
specified by
\[
P\left(\mathbf{x}^{(d+1)}\right)=\left(\mathbf{t} ; \mathbf{s}_{1}, \ldots, \mathbf{s}_{k+1}\right)
\]
is one-to-one, with Jacobian \(D \mathbf{x}^{(d)} / D\left(\mathbf{t} ; \mathbf{s}_{1}, \ldots, \mathbf{s}_{k+1}\right)=t_{1}^{d_{1}} \ldots t_{k+1}^{d_{k+1}}\) and with inverse \(P^{-1}\) given by \(x_{c_{j-1}+r}=t_{j} s_{j r}\left(r=1, \ldots, d_{j} ; j=1, \ldots, k+1\right)\). We shall refer to \(P\left(\mathbf{x}^{(d+1)}\right)\) as a partition of order \(k\) of the composition \(\mathbf{x}^{(d+1)}\). Thus a partition directs attention to an amalgamation together with its associated subcompositions.

Independence notation. In discussing statistical independence we use the notation of Dawid (1979). Thus \(C\left(\mathbf{x}^{(c)}\right) \| C\left(\mathbf{x}_{(c)}\right)\) denotes independence of the two subcompositions, and \(C\left(\mathbf{x}^{(c)}\right) \Perp C\left(\mathbf{x}_{(c)}\right) \mid T\left(\mathbf{x}^{(c)}\right)\) denotes their conditional independence, given the sum, \(x_{1}+\ldots+x_{c}\). We use \(\frac{1}{} \mathbf{w}^{(d+1)}\) to indicate that \(\mathbf{w}^{(d+1)}\) consists of independent components.

\subsection*{2.2. The Dirichlet Class}

Undoubtedly the only familiar class of distributions on \(\mathbb{S}^{d}\) is the Dirichlet class with typical member \(D^{d}(\alpha)\) having density function
\[
\prod_{i=1}^{d+1} x_{i}^{\alpha_{i}-1} / \Delta(\alpha) \quad\left(\mathbf{x}^{(d)} \in \mathbb{S}^{d}\right)
\]
where \(\alpha\) or \(\alpha^{(d+1)} \in \mathbb{P}^{d+1}\) is a \((d+1)\)-vector parameter and
\[
\Delta(\alpha)=\Gamma\left(\alpha_{1}\right) \ldots \Gamma\left(\alpha_{d+1}\right) / \Gamma\left(\alpha_{1}+\ldots+\alpha_{d+1}\right)
\]
is the Dirichlet function. A major obstacle to its use in the statistical analysis of compositional data is that it seldom, if ever, provides an adequate description of actual patterns of variability of compositions. The reasons for this are not difficult to find. First, the isoprobability contours of every Dirichlet distribution with \(\alpha_{i}>1(i=1, \ldots, d+1)\) are convex, and so the Dirichlet class must fail to describe obviously concave data patterns such as in Fig. 1. More importantly, the Dirichlet class has so much independence structure built into its definition that it represents, not a convenient modelling class for compositional data but the ultimate in independence hypotheses. This strong independence structure stems from a well-known relationship between the Dirichlet and gamma classes, which can be expressed in the terminology of compositional data as follows.
D1. Any Dirichlet composition in \(\mathbb{S}^{d}\) can be expressed as the composition of a basis of \(d+1\) independent gamma-distributed quantities, each with the same scale parameter.
There are many ways of expressing the strong internal independence structure of \(D^{d}(\alpha)\) without reference to a conceptual external basis. For our purposes here we can collect most of these into a single general result concerning any partition of a Dirichlet composition.
D2. If \(\mathbf{x}^{(d+1)}\) is \(D^{d}(\alpha)\) then, for partition (2.6), \(\mathbf{t}\left\|\mathbf{s}_{1}\right\| \ldots \leq \mathbf{s}_{k+1}\), with \(\mathbf{t}\) of \(D^{k}(\gamma)\) form and \(\mathbf{s}_{j}\) of \(D^{d_{j}}\left(\boldsymbol{\beta}_{j} \gamma_{j}\right)\) form \((j=1, \ldots, k+1)\) where \(P\left(\boldsymbol{\alpha}^{(d+1)}\right)=\left(\boldsymbol{\gamma} ; \boldsymbol{\beta}_{1}, \ldots, \boldsymbol{\beta}_{k+1}\right)\).
We shall show later the relevance of these two properties to various concepts of independence in the simplex.

The realization that the Dirichlet class leans so heavily towards independence has prompted a number of authors (Connor and Mosimann, 1969; Darroch and James, 1974; Mosimann, 1975b; James and Mosimann, 1980; James, 1981) to search for generalizations of

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0201cf8f-8652-425d-8516-e7a333b99c5c-05.jpg?height=700&width=814&top_left_y=174&top_left_x=417}
\captionsetup{labelformat=empty}
\caption{Fig. 1. A concave data set and the 95 per cent prediction region of a fitted additive logistic normal distribution. The points are the subcompositions \(C\left(\mathrm{Na}_{2} \mathrm{O}+\mathrm{K}_{2} \mathrm{O}, \mathrm{Fe}_{2} \mathrm{O}_{3}, \mathrm{MgO}\right)\) of 23 aphyric Si-poor basalt-benmoreites from the AFM diagram of Fig. 7 of Thomson, Essen and Duncan (1972).}
\end{figure}
the Dirichlet class with less structure. Their efforts have met with only limited success and it remains an open problem to find a useful parametric class of distributions on \(\mathbb{S}^{d}\) which contains the Dirichlet class but also contains distributions which do not satisfy any of the simplex independence properties already appearing in the literature or to be introduced in this paper.

In our view the way out of the impasse is simply to travel by a different route, escaping from the awkward constrictions of \(\mathbb{S}^{d}\) into the wide open spaces of \(\mathbb{R}^{d}\) through suitably selected transformations between \(\mathbb{S}^{d}\) and \(\mathbb{R}^{d}\).

\subsection*{2.3. Transformed Normal Classes}

The idea of inducing a tractable class of distributions over some awkward sample space from a proven and well-established class over some simpler space is at least a century old. McAlister (1879), faced with the "awkward" sample space \(\mathbb{P}^{1}\), saw that if he considered \(y\) in \(\mathbb{R}^{1}\) to be \(N\left(\mu, \sigma^{2}\right)\) then the transformation \(x=\exp (y)\) would induce a useful "expnormal" distribution \(\Lambda\left(\mu, \sigma^{2}\right)\) on \(\mathbb{P}^{1}\) : he, of course, expressed the idea in terms of the inverse, logarithmic, transformation and we are stuck with the name lognormal. Over the century there has been a continuing interest in transformations to normality, intensified in recent years following the work of Box and Cox (1964) and the increasing availability of tests of multinormality, as in Andrews, Gnanadesikan and Warner (1973). It seems surprising therefore that the idea of moving from multinormal distributions \(N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) on \(\mathbb{R}^{d}\) to a class \(f N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) of distributions on \(\mathbb{S}^{d}\) by a suitable transformation \(f: \mathbb{R}^{d} \rightarrow \mathbb{S}^{d}\) has been so slow to emerge. Our surprise must be even greater when one such transformation, the additive logistic transformation \(a_{d}: \mathbb{R}^{d} \rightarrow \mathbb{S}^{d}\) defined in Table 1, is already heavily exploited in other areas of statistical activity, such as logistic discriminant analysis (Cox, 1966; Day and Kerridge, 1967; Anderson, 1972) and in the analysis of binary data (Cox, 1970).

Aitchison and Shen (1980) have identified as the logistic-normal class those distributions induced on \(\mathbb{S}^{d}\) from the class of \(N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) distributions on \(\mathbb{R}^{d}\) by the transformation \(a_{d}\). The earliest explicit mention of this class we have traced is in a personal communication to Johnson and Kotz (1972, p. 20) by Obenchain, who does not seem to have developed the idea

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1
Elementary logistic transformations from \(\mathbb{R}^{d}\) to \(\mathbb{S}^{d}\)}
\begin{tabular}{|l|l|l|}
\hline Name and notation & Specification & Inverse \\
\hline Additive \(a_{d}\) & \(x_{i}\left\{1+\sum_{j=1}^{d} \exp \left(y_{i}\right)\right\}= \begin{cases}\exp \left(y_{i}\right) & (i=1, \ldots, d) \\ 1 & (i=d+1)\end{cases}\) & \(y_{i}=\log \frac{x_{i}}{x_{d+1}}\) \\
\hline Multiplicative \(m_{d}\) & \(x_{i} \prod_{j=1}^{i}\left\{1+\exp \left(y_{i}\right)\right\}= \begin{cases}\exp \left(y_{i}\right) & (i=1, \ldots, d) \\ 1 & (i=d+1)\end{cases}\) & \(y_{i}=\log \frac{x_{i}}{1-\sum_{j=1}^{i} x_{j}}\) \\
\hline Hybrid \(h_{d}\) & \[
\begin{gathered}
x_{1}=\exp \left(y_{1}\right) /\left\{1+\exp \left(y_{1}\right)\right\} \\
x_{i}\left\{1+\sum_{j=1}^{i-1} \exp \left(y_{i}\right)\right\}\left\{1+\sum_{j=1}^{i} \exp \left(y_{j}\right)\right\} \\
=\exp \left(y_{i}\right),(i=2, \ldots, d) \\
x_{d+1}=1 /\left\{1+\sum_{j=1}^{d} \exp \left(y_{j}\right)\right\}
\end{gathered}
\] & \[
y_{i}=\log \frac{x_{1}}{1-x_{1}}
\]
\[
\begin{aligned}
& y_{i}=\log \frac{x_{i}}{\left(1-\sum_{j=1}^{i-1} x_{j}\right)\left(1-\sum_{j=1}^{i} x_{j}\right)} \\
& (i=2, \ldots, d)
\end{aligned}
\] \\
\hline
\end{tabular}
\end{table}
further. Aitchison and Shen (1980) cite a number of earlier implicit uses, particularly as a vehicle for the description of prior and posterior distributions of vectors of multinomial probabilities which are naturally confined to a suitably dimensioned simplex. Leonard (1973) started a thorough investigation of this use of the class over simplex parameter spaces. The first use of the class for describing patterns of variability of data appears to be for probabilistic data in a medical diagnostic problem by Aitchison and Begg (1976) and for compositional data by Aitchison and Shen (1980), who discuss a number of useful properties and demonstrate the simplicity of its application in a variety of problems. Our interest here in logistic-normal distributions is in their membership of a wider class of transformed normal distributions on the simplex and their use in relation to the independence concepts of subsequent sections.

The additive logistic transformation \(a_{d}\) is by no means the only transformation from \(\mathbb{R}^{d}\) to \(\mathbb{S}^{d}\), and may be quite unsuited to particular investigations. Table 1 gives two other elementary transformations, the multiplicative logistic \(m_{d}\) and the hybrid logistic \(h_{d}\). All three transformations \(a_{d}, m_{d}, h_{d}\) have Jacobian \(D \mathbf{x} / D \mathbf{y}\) given by \(x_{1} x_{2} \ldots x_{d+1}\). We shall see that such elementary transformations can act as the building blocks of much more complicated transformations. An obvious comment is that the exponential function used in the definitions is not an essential feature; it could be replaced by any one-to-one transformation from \(\mathbb{R}^{1}\) to \(\mathbb{P}^{1}\), though there are few transformations as tractable.

There are two main ways of building further useful transformations.
Linear transformation method. The fact that the \(N^{d}\) class on \(\mathbb{R}^{d}\) is closed under the group of non-singular linear transformations implies that any one of the elementary transformednormal classes on \(\mathbb{S}^{d}\) will have a related closure property (Aitchison and Shen, 1980). In practical terms this means that we could replace \(\mathbf{y}^{\text {(d) }}\) by \(\mathbf{Q y}{ }^{\text {(d) }}\), with \(\mathbf{Q}\) non-singular, in any one of the elementary transformations and formally obtain a new transformation but with the assurance that we are remaining within the same class of distributions on \(\mathbb{S}^{d}\). For example, with \(a_{d}\) and
\[
q_{i i}=1(i=1, \ldots, d), \quad q_{i, i+1}=-1(i=1, \ldots, d-1), \quad q_{i j}=0
\]
otherwise, we obtain a new transformation
\[
x_{i}\left\{1+\sum_{k=1}^{d} \exp \left(\sum_{j=k}^{d} y_{j}\right)\right\}=\exp \left(\sum_{j=i}^{d} y_{j}\right), \quad y_{i}=\log \left(x_{i} / x_{i+1}\right) \quad(i=1, \ldots, d) .
\]
involving ratios of adjacent components of the composition.
Partition transformation method. When a partition \(P\left(\mathbf{x}^{(d+1)}\right)\), as defined in (2.6), is under consideration a relevant transformation from \(\mathbb{R}^{d}\) to \(\mathbb{S}^{d}\) may be constructed as follows. Let
\[
f_{0}: \mathbb{R}^{d} \rightarrow \mathbb{S}^{k}, \quad f_{j}: \mathbb{R}^{d_{j}} \rightarrow \mathbb{S}^{d_{j}} \quad(j=1, \ldots, k+1)
\]
be any \(k+2\) suitably dimensioned elementary transformations from Table 1. The compound \(\mathbf{f}=\left(f_{0} ; f_{1}, \ldots, f_{k+1}\right)\) is a one-to-one transformation
\[
\mathbf{f}: \mathbb{R}^{d}=\mathbb{R}^{k} \times \prod_{j=1}^{k+1} \mathbb{R}^{d_{j}} \rightarrow \mathbb{S}^{k} \times \prod_{j=1}^{k+1} \mathbb{S}^{d_{j}} .
\]
and the inverse transformation \(P^{-1}\) then takes us further on to \(\mathbb{S}^{d}\) to complete a transformation \(P^{-1} \mathbf{f}\) from \(\mathbb{R}^{d}\) to \(\mathbb{S}^{d}\). We denote this resultant transformation shortly by ( \(f_{0} ; f_{1}, \ldots, f_{k+1}\) ). The choice of \(f_{j}(j=0, \ldots, k+1)\) from among the appropriately dimensioned elementary transformations obviously offers a multitude of transformations from \(\mathbb{R}^{d}\) to \(\mathbb{S}^{d}\). The choice in any particular application should clearly depend on the situation under investigation.

\section*{3. Validity of Transformed-Normal Models}

Any statistical weapon designed to overcome such a resistant fortress as the simplex is unlikely to gain acceptance before undergoing proving tests as to its suitability to the terrain.

Goodness-of-fit tests. If \(\mathbf{x}^{(d+1)}\) follows a \(f N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) distribution in \(\mathbb{S}^{d}\) then \(\mathbf{y}^{(d)}=f^{-1}\left(\mathbf{x}^{(d+1)}\right)\) follows a \(N^{d}(\mu, \Sigma)\) distribution in \(\mathbb{R}^{d}\). We can thus test the goodness of fit of any transformednormal class to a compositional data set by applying the now extensive battery of multivariate normal tests, as for example in Andrews, Gnanadesikan and Warner (1973), to the transformed data set.

For \(d\)-dimensional compositional data sets we have applied Kolmogorov-Smirnov and Cramér-von Mises tests in their Stephens (1974) versions to all \(d\) marginal distributions, to all \(\frac{1}{2} d(d-1)\) bivariate angle distributions, and to the distribution of \(d\)-dimensional radii. For the Skye lava compositions we have tested in this way both the additive \(a N^{9}\) and the multiplicative \(m N^{9}\) logistic-normal models. For the additive version not a single one of the battery of 92 tests gives a significant indication of non-normality at the 5 per cent significance level; for the multiplicative version only one of the marginal tests gives evidence of any departure from normality, at the 1 per cent significance level. Application of the battery of tests to another 20 data sets of different geological types similarly encourages the view that transformed-normal distributions may have an important practical role to play in the analysis of compositional data.

The ability of transformed-normal distributions to cope with concave data sets in \(\mathbb{S}^{d}\) is illustrated in Fig. 1 where the 95 per cent prediction region of a fitted additive logistic-normal distribution, constructed by transformation of the corresponding elliptical region in \(\mathbb{R}^{2}\), neatly contains the data points.

Two caveats are worth recording. First, testing for multivariate normality and trying to detect outliers are two highly interrelated activities (Gnanadesikan and Kettenring, 1972); delicate judgements may occasionally have to be made between rejection of an apparent outlier to justify multivariate normal modelling and retention of suspect data with consequently more complex modelling. Secondly, in multivariate normal regression modelling, multivariate normality of the vector residuals, not of the regressand vectors, is the hypothesis under scrutiny. Thus in Example 3 the sediment compositions show significant departure from
logistic-normality, whereas in the appropriate regression analysis on the explanatory water depth, reported later in Section 7.3, the residuals survive such scrutiny.

Genesis models. Many of the natural and sampling processes by which compositions are determined are extremely complex; see, for example, the description by Chayes (1971, p. 44) for some geological sampling. Just as some support for normal and lognormal modelling can be provided by additive and multiplicative central limit theorems so we can postulate a process of random modifications to compositions which lead, through central limit theory arguments, to transformed-normal distributions for compositions. The underlying concept is that of a perturbation \(\mathbf{w}^{(d+1)} \in \mathbb{P}^{d+1}\), whose effect on a composition \(\mathbf{x}^{(d+1)} \in \mathbb{S}^{d}\) is to produce a perturbed composition
\[
\mathbf{w} \circ \mathbf{x}=C\left(w_{1} x_{1}, \ldots, w_{d+1} x_{d+1}\right) .
\]

Successive perturbations \(\mathbf{w}_{[1]}, \mathbf{w}_{[2]}, \ldots\) on an initial composition \(\mathbf{x}_{[0]}\) produce a sequence of compositions \(\mathbf{x}_{[1]}, \mathbf{x}_{[2]}, \ldots\), related by \(\mathbf{x}_{[r]}=\mathbf{w}_{[r]}{ }^{\circ} \mathbf{x}_{[r-1]}(r=1,2, \ldots)\) and satisfying
\[
\log \left(x_{r j} / x_{r, d+1}\right)=\log \left(x_{0 j} / x_{0, d+1}\right)+\sum_{i=1}^{r} \log \left(w_{i j} / w_{i, d+1}\right) .
\]

It is then clear that suitable conditions on the perturbations could lead, for large \(r\), to approximately additive logistic-normal or \(a N^{d}\) distributions for \(\mathbf{x}_{[r]}\).

\section*{4. Extrinsic Analysis of Independence}

\subsection*{4.1. Introduction}

We distinguish between two forms of structural analysis of compositional data:
(1) extrinsic analysis, where compositions in \(\mathbb{S}^{d}\) have been derived, or are conceptualized as arising, from bases in \(\mathbb{P}^{d+1}\) and interest is in the relation of composition to basis;
(2) intrinsic analysis, where there is no basis and so interest is not directed outside the simplex but in the composition per se.
In this section we consider two independence concepts of extrinsic analysis.
One general point should first be made. It will be obvious that most of the independence concepts introduced and their properties could be presented in a weaker moment form involving correlations. Since a main aim is to develop tests of hypotheses within transformed normal models, where independence and zero correlation coincide, we have not considered it worthwhile to interrupt the narrative to draw such fine distinctions when they exist.

\subsection*{4.2. Compositional Invariance}

In Examples 2 and 4 the compositions arise from actual bases in the form of quantities of different types of pebbles and expenditures in different commodity groups. Questions such as "Is pebble-type composition independent of the abundance of the pebbles?" and "To what extent is the pattern of household expenditure dependent on total expenditure?" direct us towards investigation of the relationship between the composition \(\mathbf{x}=C(\mathbf{w})\) and the total size \(t=T(\mathbf{w})\) of a basis \(\mathbf{w} \in \mathbb{P}^{d+1}\). This leads naturally to the following independence concept.

Definition: compositional invariance of a basis. A basis \(\mathbf{w} \in \mathbb{P}^{d+1}\) is compositionally invariant if \(C(\mathbf{w}) \| T(\mathbf{w})\).

This concept has appeared under a variety of guises: as the Lukacs condition in a characterization of the Dirichlet distribution (Mosimann, 1962), as additive isometry in the analysis of biological shape and size (Mosimann, 1970, 1975a, b), as proportion invariance in the study of \(F\)-independence (Darroch and James, 1974).

The development of a satisfactory parametric test of compositional invariance seems to have been delayed by two model-building deficiencies of the multivariate lognormal class \(\Lambda^{d+1}(\boldsymbol{\mu}, \boldsymbol{\Omega})\), a natural first-thought contender for the role of modelling the variability of bases in \(\mathbb{P}^{d+1}\).
(1) If \(\mathbf{w}\) is \(\Lambda^{d+1}(\boldsymbol{\mu}, \boldsymbol{\Omega})\) there is no simple, tractable form for the distribution of \(T(\mathbf{w})\) and so investigation of \(C(\mathbf{w}) \Perp T(\mathbf{w})\) is difficult.
(2) A multivariate lognormal basis \(\mathbf{w}\) can be compositionally invariant only if \(\mathbf{w}\) has a degenerate, one-dimensional distribution with covariance matrix \(\boldsymbol{\Omega}=\operatorname{cov}(\log \mathbf{w})\) a scalar multiple of the matrix \(\mathbf{U}_{d+1}\) consisting of unit elements and so of rank 1 (Mosimann, 1975b).
Thus not only from a point of view of tractability but also on logical grounds, study of compositional invariance within multivariate lognormal modelling of the basis is doomed to failure. Since non-degenerate compositional invariance is obviously a logical possibility the message to the practical statistician is clear: he must do better in his modelling. With transformed normal classes the answer is easy. Since interest is in \(\mathbf{x} \Perp \boldsymbol{t}\) we need not insist on finding an elegant model for the joint distribution of ( \(t, \mathbf{x}\) ) but concentrate on the conditional distribution \(p(\mathbf{x} \mid t)\) using a transformed normal regression form such as \(f N^{d}(\boldsymbol{\alpha}+\boldsymbol{\beta} t, \boldsymbol{\Sigma})\) or \(f N^{d}(\boldsymbol{\alpha}+\boldsymbol{\beta} \log t, \Sigma)\). Then compositional invariance is simply the parametric hypothesis \(\boldsymbol{\beta}=\mathbf{0}\). Moreover, testing this hypothesis on a data set consisting of \(n\) bases, and hence of \(n\) pairs of corresponding compositions and sizes, is standard methodology in multivariate analysis of dispersion (Morrison, 1976, Chapter 5). This regression approach seems appropriate since we would surely want, in the event of rejecting the hypothesis of compositional invariance, to study the basis further by trying to describe the nature of the dependence of composition on size.

Glacial tills. We have tested compositional invariance for the 93 pebble samples of Example 2 in both the \(a N^{3}(\boldsymbol{\alpha}+\boldsymbol{\beta} t, \Sigma)\) and \(a N^{3}(\boldsymbol{\alpha}+\boldsymbol{\beta} \log t, \Sigma)\) models with very similar results. Using the generalized likelihood ratio criterion as in Morrison (1976, p. 222) we obtain values 2.74 and 3.05 for the test statistics, each to be compared against 7.81 , the upper 5 per cent \(\chi^{2}\) (3) point. Thus there is no evidence against compositional invariance in these glacial tills. Two comments should be made. First, while two of the marginal tests indicate evidence of departure from additive logistic normality the other tests show no such evidence. Secondly, zero components in 14 of the samples were replaced by proportions 0.0005 , half the lowest recorded value, before analysis. We shall return to this problem of zeros in Section 7.4.

Household expenditure budgets. Incorporating compositional analysis directly into the analysis of household budgets has many advantages and provides opportunities for new forms of investigation. Modelling as above with \(p(\mathbf{x} \mid t)\) of \(a N^{d}(\boldsymbol{\alpha}+\boldsymbol{\beta} \log t, \boldsymbol{\Sigma})\) form has interesting consequences. First, the sometimes troublesome Engel aggregation condition (Brown and Deaton, 1972, p. 1163) that, for each household, total expenditure should equal the sum of all commodity expenditures, is automatically satisfied. Secondly, the hypothesis of compositional invariance, \(\boldsymbol{\beta}=\mathbf{0}\), has a direct interpretation in terms of the income elasticities \(e_{i}=\partial \log w_{i} / \partial \log t\) of demand \((i=1, \ldots, d+1)\), if for simplicity we identify household total expenditure with household income. In expectation terms \(\beta_{i}=e_{i}-e_{d+1}(i=1, \ldots, d)\), so that compositional invariance corresponds to equality of all \(d+1\) income elasticities. Thirdly, whether or not there is compositional invariance, the modelling can clearly be extended to a full consumer demand analysis by the incorporation of commodity prices and other explanatory variables such as household type and household composition into the mean parameter of the \(a N^{d}\) distribution. Indeed such an extension can be shown to be identical with the Houthakker (1960) indirect addilog model of consumer demand (Brown and Deaton, 1972, equation 115).

There is, however, an important extra flexibility in the present compositional approach, for we are not restricted to the additive logistic transformation but could equally use other forms, for example, directed towards the investigation of whether households place priorities in allocation of expenditures on some commodity groups.

In the above discussion we have identified household total expenditure \(t\) with household income \(s\). This is not an essential feature of the modelling since we could approach it through the conditioning
\[
p(s, t, \mathbf{x})=p(s) p(t \mid s) p(\mathbf{x} \mid s, t)
\]
with perhaps the reasonable assumption that \(\mathbf{x} \Perp s \mid t\) leading to the above focus on \(p(\mathbf{x} \mid t)\).

Application of the test of compositional invariance gives observed values of \(36 \cdot 4\) and \(39 \cdot 0\) for the test statistics for household types A and B respectively, each to be compared against upper \(\chi^{2}(6)\) values, and hence highly significant. Thus for both types A and B the hypothesis of compositional invariance is firmly rejected, not surprisingly when we recall that the hypothesis is equivalent to the equality of the income elasticities for all commodity groups. More interestingly, from the estimated values of \(\beta_{i}\) the relationship \(\beta_{i}=e_{i}-e_{d+1}\) provides us with an ordering of the commodity groups in terms of increasing magnitude of income elasticity, that is in conventional economic jargon from necessity to increasing luxury groups. For household type A this ordering is as follows: housing; fuel and light; foodstuffs; transport and vehicles; alcoholic drinks, tobacco and miscellaneous goods; services; clothing, footwear and durable goods. For household type B the ordering is identical except that the groups 4 and 5 are interchanged. While these orderings seem reasonable for Hong Kong it should be clear that any satisfactory analysis must involve the introduction of concomitant explanatory variables such as household size and the use of data from the eventual household expenditure survey rather than from specially selected pilot households. We hope to report on a more detailed analysis elsewhere.

\subsection*{4.3. Basis Independence}

Even when no basis actually exists a number of authors, conscious of the difficulties of defining independence concepts for compositions, have seen a method of escape through the relating of the compositional property to that of independence of an imaginary basis. Their various forms of this idea can be simply expressed as follows.

Definition: basis independence. A composition \(\mathbf{x}^{(d)} \in \mathbb{S}^{d}\) is said to have basis independence if there exists a basis \(\mathbf{w}^{(d+1)} \in \mathbb{P}^{d+1}\) with \(\Perp \mathbf{w}^{(d+1)}\) and such that \(\mathbf{x}^{(d)}=C\left(\mathbf{w}^{(d+1)}\right)\).

Since every Dirichlet-distributed composition has basis independence, by property D1 of Section 2.2, the Dirichlet class has obviously no fruitful role to play in the investigation of this independence property.

Attention has concentrated on assessing null correlations, the spurious correlations that would arise in the raw proportions solely from the process of forming proportions from conceptual, independent basis measurements, and subsequently on comparing sample correlations against these null values (Chayes, 1960, 1962, 1971; Mosimann, 1962; Chayes and Kruskal, 1966; Darroch, 1969). Many awkward features and pitfalls of this direct correlational approach have been pointed out: see, for example, Aitchison (1981a) who, after emphasizing the limitations of inferences about bases from compositions imposed by the fact that a composition \(\mathbf{x}^{(d+1)}\) determines a basis \(\mathbf{w}^{(d+1)}=t \mathbf{x}^{(d+1)}\) only up to a multiplicative factor \(t\), provides an overall test by showing that basis independence is associated with a particularly simple covariance structure of logratios of the raw proportions:
\[
\operatorname{cov}\left\{\log \left(\mathbf{x}^{(d)} / x_{d+1}\right)\right\}=\operatorname{diag}\left\{\lambda_{1}, \ldots, \lambda_{d}\right\}+\lambda_{d+1} \mathbf{U}_{d}, \quad\left(\lambda_{i}>0, i=1, \ldots, d+1\right),
\]
where \(\mathbf{U}_{d}\) is the \(d \times d\) matrix of units. Even a simplified approach, however, has merit only so long as it proves impossible to provide an equivalent intrinsic concept. Since we have now discovered a simple way of defining the illusive concept of almost-independence within the composition itself we proceed immediately to this new concept.

\section*{5. Intrinsic Analysis: Complete Subcompositional Independence}

It has long been appreciated that there must be at least one pair of correlated components in any composition \(\mathbf{x}^{(d+1)}\). An obvious first problem in studying independence in \(\mathbb{S}^{d}\) is therefore to find a structure which most closely approaches the unattainable goal of \(\Perp \mathbf{x}^{(d+1)}\). The following definition embodies such a concept.

Definition: complete subcompositional independence. A composition \(\mathbf{x}^{(d+1)}\) has complete
subcompositional independence if, for each possible partition of \(\mathbf{x}^{(\mathbf{d}+\mathbf{1})}\), the set of all its subcompositions is independent.

Every Dirichlet composition has complete subcompositional independence, by D2 of Section 2.2. Note also that complete subcompositional independence is automatically satisfied by any composition of dimension \(d=1\) or 2 , since partitions involve one-component subvectors such as \(x_{1}\) which have trivial subcompositions such as \(C\left(x_{1}\right)=1\).

For a composition \(\mathbf{x}^{(d+1)}\) with complete subcompositional independence, \(C\left(\mathbf{x}^{(b)}\right) \Perp C\left(\mathbf{x}_{(c)}\right)\) for \(b \leqq c\). Moreover, since every subcomposition based on a two-dimensional subvector such as ( \(x_{1}, x_{2}\) ) is a function only of the ratio \(x_{1} / x_{2}\), complete subcompositional independence implies independence of every pair of ratios \(x_{i} / x_{j}\) and \(x_{k} / x_{i}\) with \(i, j, k, l\) all different and, a fortiori, of the logratios \(\log \left(x_{i} / x_{j}\right)\) and \(\log \left(x_{k} / x_{i}\right)\). This implication can be fully expressed in terms of the special form for the covariance structure
\[
\boldsymbol{\Sigma}_{H}=\operatorname{cov}\left\{\log \left(\mathbf{x}^{(d)} / x_{d+1}\right)\right\}=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{d}\right)+\lambda_{d+1} \mathbf{U}_{d},
\]
where \(\lambda^{(d+1)}\) has the following interpretations:
\[
\lambda_{i}=\operatorname{cov}\left\{\log \left(x_{j} / x_{i}\right), \log \left(x_{k} / x_{i}\right)\right\}, \quad \lambda_{i}+\lambda_{j}=\operatorname{var}\left\{\log \left(x_{i} / x_{j}\right)\right\}
\]
where \(i, j, k\) are unequal. This attractive form for the covariance structure suggests that additive logistic-normal modelling may be useful. This approach is further encouraged by the easily proved equivalence result, that, for an additive logistic-normal composition, complete subcompositional independence and covariance structure (5.1) are equivalent.

The similarity of (5.1) to (4.1) confirms that we have found an intrinsic counterpart of the doubtful extrinsic concept of basis independence. The difference lies only in the restrictions placed on \(\lambda^{(d+1)}\), the positivity in form (4.1) being relaxed to the extent that \(\lambda^{(d+1)}\) need only ensure positive-definiteness of form (5.1).

Within this framework of the \(a N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) class for the composition \(\mathbf{x}^{(d+1)}\), testing for complete subcompositional independence becomes testing the parametric hypothesis that the covariance structure is of form (5.1). Note that this hypothesis places \(\frac{1}{2} d(d-1)-1\) constraints on the parameters. No exact test of the hypothesis has been found but the familiar Wilks (1938) asymptotic generalized likelihood ratio test gives a reasonable substitute. This compares
\[
n\left\{\log \left(\left|\hat{\mathbf{N}}_{H}\right| /\left|\hat{\mathbf{\Sigma}}_{M}\right|\right)+\operatorname{trace}\left(\hat{\mathbf{\Sigma}}_{H}^{-1} \mathbf{V}\right)-d\right\},
\]
where \(\mathbf{V}\) is the sample covariance matrix of the transformed vector \(\log \left(\mathbf{x}^{(d)} / x_{d+1}\right)\) and \(\hat{\mathbf{\Sigma}}_{H}\) and \(\hat{\boldsymbol{\Sigma}}_{M}\) are the maximum likelihood estimates of \(\boldsymbol{\Sigma}\) under the hypothesis and model, against the appropriate upper percentile of \(\chi^{2}\left\{\frac{1}{2} d(d-1)-1\right\}\). The estimate \(\hat{\mathbf{\Sigma}}_{M}\) is simply \(\mathbf{V}\) but the computation of \(\hat{\mathbf{\Sigma}}_{H}\) requires a suitable numerical maximization procedure. We have used a modification of the Marquardt (1963) mixture of Newton-Raphson and steepest ascent methods, exploiting the special forms taken by \(\left|\boldsymbol{\Sigma}_{\boldsymbol{H}}\right|, \boldsymbol{\Sigma}_{\boldsymbol{H}}^{-1}\) and the positive-definiteness constraint. The details are tedious and unimportant to our context: any reader interested may obtain a program in BASIC from the author.

Skye lavas. For the Skye lava data of Example 1 with \(n=32\) and \(d=9\) we obtain the value 325 for the test quantity (5.3) to be compared against upper \(\chi^{2}(35)\) values, with consequent sound rejection of the hypothesis of complete subcompositional independence.

\section*{6. Intrinsic Analysis: Partition of Order One}

\subsection*{6.1. Introduction}

In their considerations of geochemical compositions geologists almost invariably concentrate on a few low-dimensional subcompositions, often with some amalgamation and represented in ternary diagrams such as AFM for \(\mathrm{C}\left(\mathrm{Na}_{2} \mathrm{O}+\mathrm{K}_{2} \mathrm{O}, \mathrm{Fe}_{2} \mathrm{O}_{3}, \mathrm{MgO}\right)\). Such partial analyses inevitably raise questions about possible loss of information and one relevant form of analysis is to ask the extent of the dependence of the subcomposition on other aspects of the
complete composition. We suspect that an underlying reason for some of the subcompositional approaches has been the absence of suitable and readily available methodology for their undoubtedly special multivariate problems with a consequential need to project down into dimensions which can be inspected by eye. We hope that transformed multinormal modelling on the simplex will encourage full multivariate analyses of geochemical data. It should also throw some light on the validity of past choices, and the optimization of future choices, of subcompositions. More positively, with this methodology and with the concepts of intrinsic independence about to be introduced, it may be possible for the geologist to formulate his questions about subcompositions more precisely. For example, if he wishes to ask what factors affect the relative proportions of iron and manganese oxides in specimens, part of his investigation must concern the relationship of the subcomposition \(C\left(\mathrm{FeO}+\mathrm{Fe}_{2} \mathrm{O}_{3}, \mathrm{MnO}\right)\) to the other aspects of the whole composition. There may, of course, be other contributory factors external to the composition such as water content. We shall see later that these could be investigated within a multivariate regression model for compositional data. Here we concentrate only on compositional factors.

As nothing more than an illustration of the analytical possibilities we consider for Example 1 the popular AFM subcomposition, actually used by Thompson, Esson and Duncan (1972); it is then natural to reorder the components, make a division of the complete vector as follows
\[
\left(\mathrm{A}=\mathrm{Na}_{2} \mathrm{O}+\mathrm{K}_{2} \mathrm{O}, \mathrm{~F}=\mathrm{Fe}_{2} \mathrm{O}_{3}, \mathrm{M}=\mathrm{MgO} \mid \mathrm{MnO}, \mathrm{P}_{2} \mathrm{O}_{5}, \mathrm{TiO}_{2}, \mathrm{CaO}, \mathrm{Al}_{2} \mathrm{O}_{3}, \mathrm{SiO}_{2}\right)
\]
and thus direct interest to this partition of order one of the composition now in \(\mathbb{S}^{8}\).
More generally then our interest is in a partition ( \(\mathbf{x}^{(c)}, \mathbf{x}_{(c)}\) ) of \(\mathbf{x}^{(d+1)}\) and in the extent of interdependence of the amalgamation \(\mathbf{t}=\left\{T\left(\mathbf{x}^{(c)}\right), T\left(\mathbf{x}_{(c)}\right)\right\}=(t, 1-t)\) and the associated left and right subcompositions \(\mathbf{s}_{1}=C\left(\mathbf{x}^{(c)}\right)\) and \(\mathbf{s}_{2}=C\left(\mathbf{x}_{(c)}\right)\). We can form altogether ten independence hypotheses, falling into four types (i) \(\mathbf{s}_{1} \Perp \mathbf{s}_{2} \mid t\); (ii) \(\mathbf{s}_{1} \| t\); (iii) \(\mathbf{s}_{1} \|\) ( \(\mathbf{s}_{2}, t\) ); (iv) \(\mathbf{s}_{1} \Perp \mathbf{s}_{2} \Perp t\); types (i)-(iii) each have two other obvious versions. Note that, by D2, the Dirichlet class satisfies all these ten independence properties. Only type (iii), in its versions \(\mathbf{s}_{1} . \Perp\) ( \(\mathbf{s}_{2}\), t) and \(\mathbf{s}_{2} \Perp\left(\mathbf{s}_{1}, t\right)\), has been previously studied, following its introduction by Connor and Mosimann (1969) under the name of neutrality. In any particular application only some subset of the ten independence hypotheses is likely to be relevant and it is clearly not practicable to consider here all possible selections of such independence hypotheses. We have therefore chosen to concentrate on six hypotheses; these, we believe, are appropriate to a large number of applications, can be fully illustrated by the application specified above, and display interesting relationships which throw light on the concept of neutrality.

\subsection*{6.2. Related Concepts of Independence}

For convenience of reference the definitions of the six forms of independence are set out formally in Table 2, their implication relationships are completely summarized in the Venn diagram of Fig. 2, and a lattice of interest in our illustrative application is shown in Fig. 3. Our main purpose in the text is then to motivate the concepts, to describe modelling within which tests can be devised and to provide a rationale for the multiple-hypothesis testing situation of the lattice.

Subcompositional invariance. In the relation of a composition to its basis the concept of compositional invariance, independence of the composition \(C(\mathbf{w})\) and the total size \(T(\mathbf{w})\) of the basis \(\mathbf{w}\) as defined in Section 4.2, plays an important role. There is a simple and useful intrinsic counterpart of this concept for subcompositions, namely subcompositional invariance, defined as independence of a subcomposition from the share of the available unit which is taken up by its components. Thus \(\mathbf{s}_{1}\) has subcompositional invariance, denoted by \(\mathscr{I}_{1}\), when \(\mathbf{s}_{1} \| t\). There is, of course, another possible subcompositional invariance associated with the partition, namely \(\mathbf{s}_{2} \Perp 1-t\) or equivalently \(\mathbf{s}_{2} \Perp t\), and denoted by \(\mathscr{I}_{2}\).

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2
Some forms of independence for the partition \(\left(\mathbf{x}^{(c)}, \mathbf{x}_{(c)}\right)\) of \(\mathbf{x}^{(d+1)}\)}
\begin{tabular}{|l|l|l|}
\hline Notation & Definition & Parametric hypothesis \\
\hline \multicolumn{3}{|l|}{Subcompositional invariance} \\
\hline \(\mathscr{I}_{1}\) & \(C\left(\mathbf{x}^{(c)}\right) \| T\left(\mathbf{x}^{(c)}\right)\) & \(\boldsymbol{\beta}_{1}=\mathbf{0}\) \\
\hline \(\mathscr{I}_{2}\) & \(C\left(\mathbf{x}_{(c)}\right) T\left(\mathbf{x}_{(c)}\right)\) & \(\beta_{2}=0\) \\
\hline \multicolumn{3}{|l|}{Conditional subcompositional independence} \\
\hline \(\mathscr{6}\) & \(C\left(\mathbf{x}^{\left({ }^{c}\right)}\right) C\left(\mathbf{x}_{(c)}\right) \mid T\left(\mathbf{x}^{(c)}\right)\) & \(\boldsymbol{\Sigma}_{12}=\mathbf{0}\) \\
\hline \multicolumn{3}{|l|}{Neutrality} \\
\hline \(\sigma_{1}\) (left) & \(C\left(\mathbf{x}^{(c)}\right) \mathbf{x}_{(c)}\) & \(\beta_{1}=0, \Sigma_{12}=0\) \\
\hline \(\mathcal{A}_{2}\) (right) & \(C\left(\mathbf{x}_{(c)}\right) \mathbf{x}^{(c)}\) & \(\boldsymbol{\beta}_{\mathbf{2}}=\mathbf{0}, \boldsymbol{\Sigma}_{\mathbf{1 2}}=\mathbf{0}\) \\
\hline \multicolumn{3}{|l|}{Partition independence} \\
\hline \(\mathscr{P}^{\mathscr{P}}\) & || \(\left\{C\left(\mathbf{x}^{(c)}\right), C\left(\mathbf{x}_{(c)}\right), T\left(\mathbf{x}^{(c)}\right)\right\}\) & \(\beta_{1}=0, \beta_{2}=0, \Sigma_{12}=0\) \\
\hline
\end{tabular}
\end{table}

Conditional subcompositional independence. The subcompositional invariances \(\mathscr{I}_{1}\) and \(\mathscr{I}_{2}\) are not concerned with the relationship of \(\mathbf{s}_{1}\) and \(\mathbf{s}_{2}\). A question of some interest concerning the two subcompositions \(\mathbf{s}_{1}\) and \(\mathbf{s}_{2}\), if, for example, \(\mathscr{I}_{1}\) and \(\mathscr{I}_{2}\) do not hold, is whether their dependence on each other may be only through the total amounts \(t\) and \(1-t\) being assigned to each. This leads naturally to the concept of conditional subcompositional independence defined as \(\mathbf{s}_{1} \| \mathbf{s}_{2} \mid t\) and denoted by \(\mathscr{C}\). We note that this hypothesis is symmetric in \(\mathbf{s}_{1}\) and \(\mathbf{s}_{2}\) so that \(\mathscr{C}\) requires no distinguishing suffices in contrast to \(\mathscr{I}_{1}\) and \(\mathscr{I}_{2}\).

Neutrality. Connor and Mosimann (1969) introduced the concept of neutrality which in our notation may be expressed as \(C\left(\mathbf{x}_{(c)}\right) \Perp \mathbf{x}^{(c)}\). This question of whether the subcomposition on the right is independent of the entire subvector on the left was motivated by a biological problem of whether turtle scutes compete for space along the plastron during their development. The concept has been the source of a number of developments by Darroch and James (1974), Darroch and Ratcliff (1970, 1971, 1978), James (1975), James and Mosimann (1980), Mosimann (1975a, b), but much of the statistical analysis of neutrality has been hampered because until recently no parametric class of distributions on the simplex had been found rich enough to accommodate both neutrality and non-neutrality.

Since there is a one-to-one transformation between \(\mathbf{x}^{(c)}\) and \(\left(\mathbf{s}_{1}, t\right)\), neutrality as defined above can be expressed as \(\mathbf{s}_{2} \|\left(\mathbf{s}_{1}, t\right)\). We term this neutrality on the right and denote it by \(\mathcal{N}_{2}\), to distinguish it from \(\mathscr{N}_{1}\), neutrality on the left where the independence property \(\mathbf{s}_{1} \Perp\left(\mathbf{s}_{2}, t\right)\) involves the relationship of the subcomposition on the left to the entire subvector on the right. Since \(\mathbf{s}_{1} \| t\) and \(\mathbf{s}_{1}\left\|s_{2} \mid t \Leftrightarrow \mathbf{s}_{1}\right\|\left(\mathbf{s}_{2}, t\right)\) we obtain the very simple relationships \(\mathscr{I}_{1} \cap \mathscr{C}=\mathscr{N}_{1}\), \(\mathscr{I}_{2} \cap \mathscr{C}=\mathscr{N}_{2}\). These, together with other similar relationships, are recorded in Fig. 2. Subcompositional invariance and conditional subcompositional independence are weaker forms of independence than neutrality and may thus be appropriate forms for investigation in situations where neutrality is rejected.

Partition independence. We have been discussing above various forms of independence involving \(\mathbf{s}_{1}, \mathbf{s}_{2}\) and \(t\), and it is natural to go to the ultimate form \(\mathbf{s}_{1} \Perp \mathbf{s}_{2} \Perp t\). We term this partition independence, denote it by \(\mathscr{P}\), and note the relation \(\mathscr{I}_{1} \cap \mathcal{N}_{2}=\mathscr{P}\) depicted in Fig. 2

Note that for \(d=1\) all the independence properties introduced are trivially satisfied. For \(d=2\) and partition ( \(x_{1}, x_{2} \mid x_{3}\) ), satisfaction of \(\mathscr{C}, \mathscr{I}_{2}\) and \(\mathscr{N}_{2}\) is again automatic; for the partition \(\left(x_{1} \mid x_{2}, x_{3}\right)\) the concepts are identical with \(\mathscr{C}=\mathscr{A}_{2}=\mathscr{N}_{2}\). It is only for \(d \geqq 3\) that we have a real distinction between the various concepts.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0201cf8f-8652-425d-8516-e7a333b99c5c-14.jpg?height=540&width=637&top_left_y=179&top_left_x=525}
\captionsetup{labelformat=empty}
\caption{Fig. 2. Diagrammatic representation of the relationships between independence properties for a partition of order one.}
\end{figure}

\subsection*{6.3. Modelling and Testing}

The problem we now face is how to model the partition ( \(t\); \(\mathbf{s}_{1}, \mathbf{s}_{2}\) ), and hence the original composition, in such a way that the independence hypotheses just discussed become appropriate parametric hypotheses. Since \(\mathscr{C}\) involves conditioning on \(t\) it is natural to try to accommodate all the hypotheses within a conditional model for \(\left(\mathbf{s}_{1}, \mathbf{s}_{2} \mid t\right)\). For example, we can adopt additive logistic modelling for \(\mathbf{s}_{1}\) and \(\mathbf{s}_{2}\) with mean vector parameters dependent on \(t\) or some transform of \(t\). With
\[
\mathbf{y}_{1}=a_{c-1}^{-1}\left(\mathbf{s}_{1}\right)=\log \left\{\mathbf{s}_{1}^{(c-1)} / s_{1 c}\right\}, \quad \mathbf{y}_{2}=a_{d-c}^{-1}\left(\mathbf{s}_{2}\right)=\log \left\{\mathbf{s}_{2}^{(d-c)} / s_{2, d-c+1}\right\}
\]
and \(z=\log \{t /(1-t)\}\) we can take our model \(M\) with conditional model for \(\left(\mathbf{y}_{1}, \mathbf{y}_{2} \mid z\right)\) of the following form:
\[
N^{d-1}\left\{\left[\begin{array}{l}
\alpha_{1}+\beta_{1} z \\
\alpha_{2}+\beta_{2} z
\end{array}\right], \quad\left[\begin{array}{ll}
\Sigma_{11} & \Sigma_{12} \\
\Sigma_{21} & \Sigma_{22}
\end{array}\right]\right. \text {. }
\]

All the independence hypotheses considered are then easily identified with constraints on the parameters \(\boldsymbol{\beta}_{1}, \boldsymbol{\beta}_{2}, \boldsymbol{\Sigma}_{12}\). For example, \(\mathscr{I}_{2}\) requires \(\mathbf{y}_{2} \Perp z\) and so has parametric counterpart \(\boldsymbol{\beta}_{2}=\mathbf{0}\); and \(\mathscr{N}_{2}\) requires the further condition \(\mathbf{y}_{1} \Perp \mathbf{y}_{2} \mid z\) or \(\boldsymbol{\Sigma}_{12}=\mathbf{0}\) and so is identical to the parametric hypotheses \(\boldsymbol{\beta}_{\mathbf{2}}=\mathbf{0}, \boldsymbol{\Sigma}_{12}=\mathbf{0}\). All these parametric counterparts are listed for convenience beside the definitions in Table 2.

Since the hypotheses under test impose linear constraints on mean vector and simple restrictions on covariance matrices the generalized likelihood ratio test statistic again takes the form (5.3) with approximate critical values given through asymptotic theory as upper \(\chi^{2}\) percentiles with appropriate degrees of freedom \(q_{H}\) for hypotheses \(H\). The derivation of \(\hat{\Sigma}_{H}\) and \(q_{H}\) for the various hypotheses and of \(\hat{\Sigma}_{M}\) is routine; for easy reference we provide the computational forms in Table 3.

\subsection*{6.4. Testing a Lattice of Hypotheses: An Application}

If only one of the independence hypotheses already discussed is under scrutiny then the appropriate test procedure set out in Section 6.3 applies. If, however, we have under investigation a number of the hypotheses then we must consider more carefully our strategy, such as order of testing. In the lattice of hypotheses set out in Fig. 3 for the partition (6.1) of the Skye lava compositions, the model is at the highest level with hypotheses at deeper levels corresponding to more and more constraints on the parameters. Viewed from the bottom of

Table 3
Maximum likelihood estimates of \(\mathbf{\Sigma}\) associated with independence hypotheses

\begin{tabular}{|l|l|l|}
\hline Hypothesis H or model \(M\) & Maximum likelihood estimate of \(\mathbf{\Sigma}\) with submatrices in the order \(\boldsymbol{\Sigma}_{11}, \boldsymbol{\Sigma}_{12}, \boldsymbol{\Sigma}_{22}\) & Degrees of freedom \(q_{H}\) \\
\hline M & \(\hat{\boldsymbol{\Sigma}}_{11}, \hat{\boldsymbol{\Sigma}}_{12}, \hat{\boldsymbol{\Sigma}}_{22}\) & \\
\hline \(\mathscr{I}_{1}\) & \(\boldsymbol{S}_{11}, \hat{\boldsymbol{\Sigma}}_{12}, \hat{\boldsymbol{\Sigma}}_{22}\) & \(c-1\) \\
\hline \(\mathscr{I}_{2}\) & \(\hat{\boldsymbol{\Sigma}}_{11}, \hat{\boldsymbol{\Sigma}}_{12}, \mathbf{S}_{22}\) & \(d-c\) \\
\hline 4 & \(\hat{\Sigma}_{11}, 0, \hat{\Sigma}_{22}\) & \((c-1)(d-c)\) \\
\hline \(F_{1}\) & \(\mathbf{S}_{11}, \mathbf{0}, \boldsymbol{\Sigma}_{22}\) & \((c-1)(d-c+1)\) \\
\hline . \(1_{2}\) & \(\hat{\boldsymbol{\Sigma}}_{11}, \mathbf{0}, \boldsymbol{S}_{22}\) & \(c(d-c)\) \\
\hline \(\mathscr{P}\) & \(\mathbf{S}_{11}, \mathbf{0}, \mathbf{S}_{22}\) & \(c(d-c)+c-1\) \\
\hline \multicolumn{3}{|l|}{Required matrix computations} \\
\hline \multicolumn{3}{|l|}{\(n S_{i j}=\sum_{r=1}^{n}\left(y_{i r}-\bar{y}_{i}\right)\left(y_{j r}-\bar{y}_{j}\right) \quad(i, j=1,2) ; \quad n S_{z z}=\sum_{r=1}^{n}\left(z_{r}-\bar{z}\right)^{2} ;\)} \\
\hline \multicolumn{3}{|l|}{\(n S_{i z}=\sum_{r=1}^{n}\left(y_{i r}-\bar{y}_{i}\right)(z,-\bar{z}) ; \quad \hat{\beta}_{i}=S_{i z} / S_{z z},(i=1,2) ;\)} \\
\hline \multicolumn{3}{|l|}{\(\hat{\Sigma}_{i j}=S_{i j}-\hat{\beta}_{i} \hat{\beta}_{j} S_{z z},(i, j=1,2)\).} \\
\hline
\end{tabular}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0201cf8f-8652-425d-8516-e7a333b99c5c-15.jpg?height=783&width=511&top_left_y=1001&top_left_x=538}
\captionsetup{labelformat=empty}
\caption{Fig. 3. Lattice for Skye lava analysis showing values of test statistics with associated degrees of freedom in brackets.}
\end{figure}
the lattice the hypothesis \(\mathscr{P}\) is the simplest explanation of the relationship of \(C\left(\mathbf{x}^{(c)}\right), C\left(\mathbf{x}_{(c)}\right)\) and \(T\left(\mathbf{x}^{(c)}\right)\), namely mutual independence. As we move up the lattice, for example to \(\mathcal{N}_{1}\), we have to introduce more parameters, namely \(\boldsymbol{\beta}_{2}\), to provide an explanation of the pattern of variability, and to \(\mathscr{C}\) further parameters, namely \(\boldsymbol{\beta}_{1}\).

For such multiple-hypothesis testing, a sensible approach is to adopt the simplicity postulate of Jeffreys (1961, p. 47): in order to move from a simple explanation, such as \(\mathscr{P}\), to a more complex explanation, such as \(\mathcal{N}_{1}\), we require to reject the simpler explanation through
an appropriate significance test. In other words, to justify the introduction of more parameters, we require a mandate, provided by significant rejection, to allow us to move to a higher level in the lattice. Thus our procedure would involve the following steps. First test \(\mathscr{P}\) within \(M\). If we cannot reject \(\mathscr{P}\) then there is nothing to justify moving from the simple explanation \(\mathscr{P}\). If we reject \(\mathscr{P}\) then we move up to the next level, testing each of \(\mathscr{N}_{1}\) and \(\mathscr{N}_{2}\) within \(M\). If we cannot reject both then we have a feasible explanation at this level. If we reject both then we move to a test of \(\mathscr{C}\) within \(M\), and so on. Note that the tests are all of a hypothesis \(H\) within \(M\) and the mechanism of these tests has already been described in Section 6.3.

For our geochemical partition the values of the test statistics with their bracketed degrees of freedom are shown at the appropriate nodes of the lattice. All the hypotheses of the lattice are rejected at significance levels well below 0.1 per cent. However we care to interpret the lattice, the \(C(\mathrm{~A}, \mathrm{~F}, \mathrm{M})\) subcomposition has clearly neither subcompositional invariance nor is it conditionally independent of the complementary subcomposition. Further analysis, not reported here, shows that it is also not (absolutely) independent, defined as \(\mathbf{s}_{1} \| \mathbf{s}_{2}\), of the complementary subcomposition. Thus any analysis of AFM which subsumes that this subcomposition is independent of other apsects of the composition is surely suspect.

\section*{7. Further Aspects of Intrinsic Analysis}

\subsection*{7.1. Partial Subcompositional Independence}

In the extrinsic approach to compositional structure some geologists, for example Sarmanov and Vistelius (1959), consider forms of partial basis independence under such terms as concretionary and metasomatic. These have satisfactory intrinsic counterparts whose form we can now indicate briefly in terms of a partition ( \(\mathbf{x}^{(c)}, \mathbf{x}_{(c)}\) ) or ( \(t ; \mathbf{s}_{1}, \mathbf{s}_{2}\) ) of order one.

Definition: partial subcompositional independence restricted by \(\mathbf{x}^{(c)}\). A composition \(\mathbf{x}^{(d+1)}\) has partial subcompositional independence restricted by \(\mathbf{x}^{(c)}\) if \(\mathbf{s}_{1} \| \mathbf{s}_{2}\) and \(\mathbf{s}_{2}\) has complete subcompositional independence within \(\mathbb{S}^{d-c}\).

Since the amalgamation \(\mathbf{t}=(t, 1-t)\) is not involved in the definition we can investigate such partial subcompositional independence within a model for the joint distribution of \(\left(\mathbf{s}_{1}, \mathbf{s}_{2}\right)\). Taking this to be of transformed normal form \(\left\{a_{c-1}^{-1}\left(\mathbf{s}_{1}\right), a_{d-c}^{-1}\left(\mathbf{s}_{2}\right)\right\}\) and hence with covariance matrix
\[
\Sigma=\operatorname{cov}\left\{\log \left(x_{i} / x_{c}\right)(i=1, \ldots, c-1) ; \log \left(x_{c+i} / x_{d+1}\right)(i=1, \ldots, d-c)\right\}
\]
we can specify partial subcompositional independence as the parametric hypothesis
\[
\boldsymbol{\Sigma}_{12}=\mathbf{0}, \quad \boldsymbol{\Sigma}_{22}=\operatorname{diag}\left(\lambda_{c+1}, \ldots, \lambda_{d}\right)+\lambda_{d+1} \mathbf{U}_{d-c}
\]
in term of the obvious partitioning of \(\mathbf{\Sigma}\). Such a formulation brings this form of independence within the scope of the test procedures developed in Section 6. Moreover, the fact that partial subcompositional independence is seen as the conjunction of two less stringent hypotheses, unconditional subcompositional independence \(\mathbf{s}_{1} \| \mathbf{s}_{2}\) and complete subcompositional independence of \(\mathbf{s}_{2}\) within \(\mathbb{S}^{d-c}\), open up another means of probing compositional structure through a lattice approach.

\subsection*{7.2. Independence up to Level \(c\)}

There are a number of situations, where a specific ordering of the \(d+1\) components has been made and already embodied in \(\mathbf{x}^{(d+1)}\), and where interest is in considering independence properties for partitions of order one at a sequence of levels \(c\). Since we consider here only independence in the form \(\mathscr{C}, \mathscr{I}_{2}\) and \(\mathscr{N}_{2}\) we drop the suffix 2 to allow us to emphasize the level \(c\) at which division has been made. Thus \(\mathscr{C}_{c}, \mathscr{I}_{c}, \mathscr{N}_{c}\) denote \(\mathscr{C}_{,} \mathscr{I}_{2}, \mathscr{N}_{2}\) at level \(c\). We recall the basic relation \(\mathscr{C}_{c} \cap \mathscr{I}_{c}=\mathscr{N}_{c}\), and, for any one of these hypotheses, say \(H_{c}\), define the corresponding concept \(H^{c}\) up to level \(c\) as follows.

Definition: independence property up to level c. A composition \(\mathbf{x}^{(\boldsymbol{d}+\mathbf{i})}\) has independence property \(H\) up to level \(c\) if \(H_{k}\) holds for \(k=1, \ldots, c\).

It follows from the relationship that \(\mathscr{C}^{c} \cap \mathscr{I}^{c}=\mathscr{N}^{c}\). For the special case when \(c=d-1\) (or equivalently \(d\) ) we use the term complete.

Definition: complete independence property. A composition \(\mathbf{x}^{(d+1)}\) possesses the complete independence property \(H\) if \(H^{d-1}\) holds.

Thus, for example, complete neutrality (Connor and Mosimann, 1969) requires \(C\left(\mathbf{x}_{(c)}\right) \| \mathbf{x}^{(c)}\) for \(c=1, \ldots, d-1\). The investigation of neutrality at different levels is best pursued in terms of the multiplicative logistic transformation. This approach has been adopted by Aitchison (1981b) to provide a suitable parametric statistical framework within which to test \(\mathscr{N}_{c} \mathscr{N}^{c}\) and lattices of hypotheses involving these. Adopting a \(m N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) model we see that the hypotheses \(\mathscr{N}_{c} \mathscr{N}^{c}\) and \(\mathscr{N}^{d-1}\) correspond to the following covariance matrix structures
\[
\left[\begin{array}{cc}
\boldsymbol{\Sigma}_{11} & \mathbf{0} \\
\mathbf{0} & \boldsymbol{\Sigma}_{22}
\end{array}\right], \quad\left[\begin{array}{cc}
\operatorname{diag}\left(\sigma_{11}, \ldots, \sigma_{c c}\right) & \mathbf{0} \\
\mathbf{0} & \boldsymbol{\Sigma}_{22}
\end{array}\right], \quad \operatorname{diag}\left(\sigma_{11}, \ldots, \sigma_{d d}\right),
\]
where \(\Sigma_{11}\) is of order \(c \times c\). The numbers of constraints imposed by the three hypotheses are \(c(d-c), c\left\{d-\frac{1}{2}(c+1)\right\}\) and \(\frac{1}{2} d(d-1)\). If \(\mathbf{V}\) is the estimated covariance matrix associated with \(d\) dimensional vectors \(\mathbf{y}_{r}(r=1, \ldots, n)\) defined in the \(m_{d}\) entry of Table 1 then the test statistics (Aitchison, 1981b) associated with \(\mathscr{N}_{c} \mathscr{N}^{c}\) and \(\mathscr{N}^{d-1}\) are again of form (5.3) with
\[
\begin{aligned}
& \left|\widehat{\mathbf{\Sigma}}_{M}\right|=|\mathbf{V}|, \quad\left|\widehat{\mathbf{\Sigma}}_{c}\right|=\left|\mathbf{V}_{11}\right| .\left|\mathbf{V}_{22}\right| \\
& \left|\widehat{\mathbf{\Sigma}}_{c}\right|=v_{11} \ldots v_{c c}\left|\mathbf{V}_{22}\right|, \quad\left|\widehat{\mathbf{\Sigma}}_{c}\right|=v_{11} \ldots v_{d d}
\end{aligned}
\]
where \(\mathbf{V}_{11}, \mathbf{V}_{22}\) are obvious submatrices of \(\mathbf{V}\) in a ( \(c, d-c\) ) partitioning and \(v_{i j}\) is the ( \(i, j\) )th element of \(\mathbf{V}\). Note that in all these tests the term trace \(\left(\hat{\boldsymbol{\Sigma}}_{H}^{-1} \hat{\boldsymbol{\Sigma}}_{M}\right)=d\) so that the test statistic reduces to \(n \log \left(\left|\hat{\boldsymbol{\Sigma}}_{H}\right| /\left|\hat{\boldsymbol{\Sigma}}_{M}\right|\right)\).

Since \(\mathscr{H}_{c}\) and \(\mathscr{N}_{c}(c>1)\) are quite distinct hypotheses we might expect \(\mathscr{I}^{c}\) and \(\mathscr{N}^{c}\) to be distinct and, since \(\mathcal{N}^{c} \subset \mathscr{I}^{c}\), to be able to devise a model for which \(\mathscr{I}^{c}\) holds but \(\mathscr{N}^{c}\) does not. We have failed to produce such a model and are beginning to conjecture that, within the framework of transformed normal modelling, \(\mathscr{I}^{c} \equiv \mathscr{N}^{c}\), though so far we have failed to prove the conjecture.

That there is a distinction between \(\mathscr{C}^{c}\) and \(\mathscr{S}^{c}\) can be readily seen for the case \(d=3\). Since \(\mathscr{C}_{1}\) and \(\mathscr{C}_{3}\) are trivially satisfied for any compositional distribution, model (6.2) with \(c=2\) and \(\sigma_{12}=0\) supports \(\mathscr{C}_{2}\) and hence complete conditional subcompositional independence, whereas \(\mathscr{N}_{2}\) does not hold unless \(\beta_{2}=0\). We have not so far found any practical problem to which the idea of \(\mathscr{C}^{c}\) seems relevant and have not therefore pursued the modelling problem further.

Skye lavas. From the strong rejection of complete subcompositional independence and the neutrality hypotheses \(\mathscr{N}_{2}\) there can be little surprise in discovering that tests of neutrality associated with an ordering such as (6.1) of the entire compositional vector lead to rejections. Simply as an illustrative example for numerical comparison therefore we show in Table 4

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4
Test results for neutrality hypotheses for Skye lavas}
\begin{tabular}{ccc}
\hline Level & Test statistic & Degrees of freedom \\
\hline 1 & 58.0 & 7 \\
2 & 135.6 & 13 \\
3 & 182.7 & 18 \\
4 & 227.0 & 22 \\
5 & 275.0 & 25 \\
6 & 283.6 & 27 \\
7 & 290.0 & 28 \\
\hline
\end{tabular}
\end{table}
values of the test statistics, as described above for testing \(\mathcal{N}^{c}\) up to all possible levels \(c\) for this ordering together with the corresponding degrees of freedom at each of the seven levels. Since the comparison is against upper chi-squared values at the degrees of freedom shown, the neutrality hypotheses up to all levels for this ordering are strongly rejected. To those who regard hypothesis-testing as a means towards arriving at a model for subsequent analyses we reiterate the important fact that rejection of all these hypotheses still leaves transformed normal models on the simplex as possible describers of patterns of variability of non-neutral compositional data.

\subsection*{7.3. Compositional Regression Models}

On finding a subcomposition \(\mathbf{s}_{1}\), such as AFM, dependent on complementary aspects \(t\) and \(\mathbf{s}_{2}\) of the complete composition we may wish to assess the conditional distribution \(p\left(\mathbf{s}_{1} \mid t, \mathbf{s}_{2}\right)\). This aspect of estimation, essentially regression analysis in transformed normal modelling, has already been illustrated for sediment compositions by Aitchison and Shen (1980) and need not be detailed here. Rather we present briefly examples where we may wish to explore the dependence of a composition \(\mathbf{x}^{(d+1)} \in \mathbb{S}^{d}\) on concomitant information \(z\).

Arctic lake sediments. In Example 3 for each Arctic lake composition the associated depth \(z\) is provided. There is, through the transformed normality approach, an obvious way of modelling to allow investigation of the dependence of composition on depth, namely to take \(p\left(\mathbf{x}^{(d+1)} \mid z\right)=f N^{d}(\mathbf{g}(z), \mathbf{\Sigma})\), where the regression function \(\mathbf{g}(z)\) can be investigated in the usual multivariate regression form. In our model \(M\) we have taken \(f\) to be \(a_{2}\) and \(\mathbf{g}(z)\) to include terms in \(z, z^{2}, \log z\) and \((\log z)^{2}\). We have then worked through a lattice of increasingly complex hypotheses along the lines of Section 6.4, and found that linear regression is certainly rejected, but that hypotheses of the form \(\mathbf{g}(z)=\boldsymbol{\alpha}+\boldsymbol{\beta} \log z\), or quadratic regression \(\mathbf{g}(z)=\boldsymbol{\alpha}+\boldsymbol{\beta} z+\boldsymbol{\gamma} z^{2}\) are equally good fits and cannot be rejected. Moreover the residuals based on either of these fitted regressions pass the complete battery of multivariate normal tests.

Household budgets. As another illustration of the simplicity of regression techniques we might extend the model of Section 4.2 to include the possibility of compositional dependence on household size, for example with the regression function of the form
\[
\alpha+\beta \log \text { (total expenditure) }+\gamma \log \text { (household size). }
\]

If we then investigate the lattice with nodes at \(\boldsymbol{\beta}=\mathbf{0}, \boldsymbol{\gamma}=\mathbf{0}\), at \(\boldsymbol{\beta}=\mathbf{0}\) and at \(\boldsymbol{\gamma}=\mathbf{0}\) we find that the hypothesis \(\boldsymbol{\gamma}=\mathbf{0}\) is the only one that cannot be rejected. Moreover, fitting of this accepted regression function leaves residuals which survive the battery of goodness-of-fit tests.

\subsection*{7.4. The Problem of Zero Components}

Throughout the paper attention has been confined to the strictly positive simplex. The reason is the obvious one that we cannot take logarithms of zero. And yet zero components do occur in a number of applications, for example, when a household spends nothing on the commodity group "tobacco and alcohol" or a rock specimen contains "no trace" of a particular mineral. In the absence of a one-to-one monotonic transformation between the real line and its non-negative subset the problem of zeros is unlikely ever to be satisfactorily resolved. A similar problem occurs in lognormal modelling and, as there, ad hoc solutions naturally depend on the frequency and nature of the zeros.

If there are only a few zeros of the no-trace type then replacement by positive values smaller than the smallest traceable amounts will allow an analysis. In such circumstances it will always be wise to perform a sensitivity analysis to determine the effect that different zero replacement values have on the conclusions of the analysis. For example, in the investigation of compositional invariance in glacial tills in Section 4.2 we replaced 14 zero proportions by 0.0005 obtaining the value 3.05 for the test statistic in the \(a N^{3}(\boldsymbol{\alpha}+\boldsymbol{\beta} \log t, \Sigma)\) modelling. For other replacement values \(0 \cdot 001,0 \cdot 00025,0 \cdot 00001\) and \(0 \cdot 000001\) the values of the test statistic are \(3 \cdot 93\),
\(2 \cdot 46,2 \cdot 01\) and \(1 \cdot 54\) all leading to the same conclusion of no evidence against compositional invariance at the 5 per cent significance level.

If there is a moderate number of real zeros it may be worth considering the device of threeparameter lognormal modelling (Aitchison and Brown, 1957, p. 14), whereby a constant, either known or to be estimated, is added to every observation. One compositional counterpart would be to apply the transformations, not to \(\mathbf{x}^{(d+1)}\), but to \(C\left(\mathbf{x}^{(d+1)}+\boldsymbol{\tau}^{(d+1)}\right)\) where \(\tau^{(d+1)}\) is either chosen or estimated. For example, for the case \(d=1\) and an additive logistic model we are considering the model with \(\log \left\{\left(x+\tau_{1}\right) /\left(1+\tau_{2}-x\right)\right.\) of \(N^{1}\left(\mu, \sigma^{2}\right)\) form, which is a fourparameter lognormal model of Johnson (1949). Clearly if \(\tau^{(d+1)}\) has to be estimated there are substantial estimation and interpretation problems even for small \(d\).

If there is a substantial number of zeros mostly in a few components and if amalgamations of components are ruled out, then some form of conditional modelling separating out the zero may be possible. For example, if the zeros are confined to the last component then the conditional distribution of \(C\left(\mathbf{x}^{(d)}\right)\) on \(x_{d+1}\) might be modelled by taking \(\log \left(\mathbf{x}^{(d-1)} / x_{d}\right)\) to be \(N^{d-1}\left(\boldsymbol{\alpha}+\boldsymbol{\beta} x_{d+1}, \Sigma\right)\) with the marginal distribution of \(x_{d+1}\) having a mass probability at zero and \(\log \left\{x_{d+1} /\left(1-x_{d+1}\right)\right\}\) following \(N^{1}\left(\mu, \sigma^{2}\right)\) for \(x_{d+1}>0\).

\subsection*{7.5. Partitions of Higher Order}

For a partition of order one we saw there are ten different independence hypotheses and that careful selection of hypotheses relevant to the practical problem is of primary importance to a successful analysis. The choice of relevant hypotheses for a higher order partition \(\left(\mathbf{t} ; \mathbf{s}_{1}, \ldots, \mathbf{s}_{k+1}\right)\) is even more crucial and we have no space to discuss it at length here. A brief look at a partition of order 2 should, however, indicate the potentialities of transformed normal modelling.

Suppose that for a partition ( \(\mathbf{t} ; \mathbf{s}_{1}, \mathbf{s}_{2}, \mathbf{s}_{3}\) ) of order 2 we wish to investigate the extent of subcompositional invariance with respect to the sums \(t_{1}, t_{2}, t_{3}\) and also whether the amalgamation \(\mathbf{t}^{(3)}=\left(t_{1}, t_{2}, t_{3}\right)\) displays complete neutrality. If we model in terms of the transformed partition ( \(\mathbf{z} ; \mathbf{y}_{1}, \mathbf{y}_{2}, \mathbf{y}_{3}\) ) of ( \(m_{k} ; a_{d_{1}}, a_{d_{2}}, a_{d_{3}}\) ) type we might use conditional modelling \(p(\mathbf{z}) p\left(\mathbf{y}_{1}, \mathbf{y}_{2}, \mathbf{y}_{3} \mid \mathbf{z}\right)\) with \(p\left(\mathbf{y}_{1}, \mathbf{y}_{2}, \mathbf{y}_{3} \mid \mathbf{z}\right)\) of multinormal form
\[
N^{d-1}\left\{\left[\begin{array}{llll}
\alpha_{1}+\beta_{1} z & \Sigma_{11} & \Sigma_{12} & \Sigma_{13} \\
\alpha_{2}+\beta_{2} z & \Sigma_{21} & \Sigma_{22} & \Sigma_{23} \\
\alpha_{3}+\beta_{3} z & \Sigma_{31} & \Sigma_{32} & \Sigma_{33}
\end{array}\right]\right\}
\]
and \(p(\mathbf{z})\) of \(N^{2}(\gamma, \boldsymbol{\Omega})\) form. It must now be clear that there could be a large number of hypotheses of interest.

As an example of a simple lattice approach we refer to Fig. 4 where forms of hypotheses of total subcompositional invariance, ( \(\mathbf{s}_{1}, \mathbf{s}_{2}, \mathbf{s}_{3}\) ) t or parametrically \(\boldsymbol{\beta}_{h}=\mathbf{0}(h=1,2,3)\), and of complete neutrality of \(\mathbf{t}^{(3)}\), namely \(\omega_{12}=0\), are brought together. The testing of such a lattice is straightforward following the lines of Section 6.4. It is also clear that the total subcompositional invariance hypothesis could be broken into interesting hypotheses such as \(\boldsymbol{\beta}_{1}=\mathbf{0}\) at a higher level of the lattice. Note that in the selection of the transformation we used \(m\) for the amalgamation since interest was in complete neutrality. Had an objective been to study neutrality within the subcompositions then \(m\) transformations could have replaced the \(a\) transformations actually used in the modelling.

\section*{8. Discussion}

There remain many loose ends to our transformed normal package. We hope that discussion in the Society tonight will reveal many statistical fingers anxious to tie up, to add to, even to repack, the package and to address it for delivery to new areas of application. The

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/0201cf8f-8652-425d-8516-e7a333b99c5c-20.jpg?height=885&width=739&top_left_y=181&top_left_x=497}
\captionsetup{labelformat=empty}
\caption{Fig. 4. Lattice for testing subcompositional invariance and neutrality for subcompositional shares.}
\end{figure}
following collection of random thoughts on the current state of the package is little more than an attempt to draw attention to topics of personal interest.
(i) We have dealt only with one-way compositions. There are problems where the components fall naturally into a two-way classification. It would be of interest to discuss problems of this type and the means of analysing them.
(ii) Of our three elementary transformations in Table 1 we have used only \(a_{d}\) and \(m_{d}\). Are there any applications where \(h_{d}\) is essential? What other transformations between \(\mathbb{R}^{d}\) and \(\mathbb{S}^{d}\) might find applications? To what extent will it be necessary to widen the class of transformations, as suggested by Aitchison and Shen (1980), through the Box and Cox (1964) approach, with \(y_{i}=\left\{\left(x_{i} / x_{d+1}\right)^{\lambda}-1\right\} / \lambda(i=1, \ldots, d)\) and \(\lambda\) being estimated from the compositional data?
(iii) Although the immediate relationship to multivariate normality usually ensures the carry-over of existing techniques, such as discriminant analysis, to compositional data some care is needed to check the validity of this transfer. For example, reduction of the compositional dimension through the use of principal components based on \(\boldsymbol{\Sigma}=\operatorname{cov}\left\{\log \left(\mathbf{x}^{(d)} / x_{d+1}\right)\right\}\) might seem a hopeful technique until it is realized that trace \((\boldsymbol{\Sigma})\) is not invariant under a permutation of the components \(x_{1}, \ldots, x_{d+1}\). A substantial modification to standard principal component analysis is required to restore the desirable invariance property.
(iv) One embarrassment of the transformed normal approach is the galaxy of possible models it offers. For example, in our discussion of right neutrality \(\mathcal{F}_{c}\) in Section 6.3 either the model (6.2) with \(\boldsymbol{\beta}_{2}=\mathbf{0}, \boldsymbol{\Sigma}_{12}=\mathbf{0}\) or the model \(m N^{d}(\boldsymbol{\mu}, \boldsymbol{\Sigma})\) of Section 7.2 with \(\boldsymbol{\Sigma}_{12}=\mathbf{0}\) could be used. Although the problem of choice between models here is no different from similar problems in other areas of statistics, tests of these separate classes could prove troublesome because of the dimension of the parameter space. One possible line of investigation might be the examination of how close the models are in the same way as Aitchison and Shen (1980) considered the closeness of logistic-normal and Dirichlet classes.
(v) Conjecture about the potential of the transformed normal approach to the analysis of the structure of geological compositions is a fascinating subject. Geologists, for example Chayes (1971), assure us that the study of correlations in compositions is essential to their understanding and yet it appears difficult to pinpoint their precise hypotheses of interest. There seems little doubt that the package can play a useful rôle in descriptive geostatistics, such as in classification, but can the fundamental hypotheses of compositional structure now be specified within the concepts of this paper?
(vi) Compositional data obviously occur in areas other than the geological and economic applications cited here; for example, in developmental biology if we wish to explore how the shape (composition) of a linear organism relates to size, the model used for the study of compositional invariance will obviously play a rôle. There are also problems with simplex sample spaces where the data are not compositions; for example, probabilistic data in \(\mathbb{S}^{d}\) occur in the analysis of subjective performance of inferential tasks (Aitchison, 1981c). More complex product sample spaces, such as \(\mathbb{S}^{d} \times \mathbb{R}^{c}\) or \(\mathbb{S}^{d} \times \mathbb{P}^{c}\), also arise, as in medical diagnosis (Aitchison and Begg, 1976), and succumb to the transformed normal technique.
(vii) There are still many distributional problems to be resolved. For example, although we have in \(m N^{d}\) a model for the investigation of complete right neutrality and in a separate \(m N^{d}\) model applied to the reversal of the vector \(\mathbf{x}^{(d+1)}\) a means of investigating left neutrality, we have been unable to find a class of models which will accommodate both forms of neutrality as parametric hypotheses and will also have non-neutral members. Thus the battle of the statistical knights who search for the holy grail of a parametric class which will include the highly structured Dirichlet distributions and all forms of dependent distributions, is obviously not over. We hope, however, that transformed normal distributions may sharpen their lances and encourage the search.

\section*{Acknowledgements}

I had the good fortune of meeting at conferences in Trieste and Sydney during 1980 three seasoned campaigners in the simplex, J. N. Darroch, W. Kruskal and J. E. Mosimann. The first draft of this paper owed much to discussions with them, with other participants at these conferences and with my Hong Kong colleagues during the formative period of some of the ideas. As usual the final version has been greatly improved by the diverse, but always penetrating and constructive, comments of four referees. I am also grateful to Colin Greenfield, Commissioner of Census and Statistics in Hong Kong, for making available the pilot household budget data.

\section*{References}

Aitchison, J. (1981a). A new approach to null correlations of proportions. J. Math. Geol., 13, 175-189.
- (1981b). Distributions on the simplex for the analysis of neutrality. In Statistical Distributions in Scientific Work (C. Taillie, G. P. Patil and B. Baldessari, eds), vol. 4, pp. 147-156. Dordrecht, Holland: D. Reidel Publishing Company.
- (1981c). Some distribution theory related to the analysis of subjective performance in inferential tasks. In Statistical Distributions in Scientific Work (C. Taillie, G. P. Patil and B. Baldessari, eds), vol. 5, pp. 363-386. Dordrecht, Holland: D. Reidel Publishing Company.
Aitchison, J. and Begg, C. B. (1976). Statistical diagnosis when the cases are not classified with certainty. Biometrika, 63, 1-12.
Aitchison, J. and Brown, J. A. C. (1957). The Lognormal Distribution. Cambridge University Press.
Aitchison, J. and Shen, S. M. (1980). Logistic-normal distributions: some properties and uses. Biometrika, 67, 261-272.
Anderson, J. A. (1972). Separate sample logistic discrimination. Biometrika, 59, 19-35.
Andrews, D. F., Gnanadesikan, R. and Warner, J. L. (1973). Methods for assessing multivariate normality. In Multivariate Analysis III (P. R. Krishnaiah, ed) pp. 95-116. New York: Academic Press.
Box, G. E. P and Cox, D. R. (1964). The analysis of transformations (with discussion). J. R. Statist. Soc. B, 26, 211-252.
Brown, A. and Deaton, A. S. (1972). Models of consumer behaviour. Econ. J., 82, 177-268.
Chayes, F. (1960). On correlations between variables of constant sum. J. Geophys. Res., 65, 4185-4193.
- (1962). Numerical correlation and petrographic variation. J. Geol., 70, 440-452.
- (1971). Ratio Correlation. University of Chicago Press.

Chayes, F. and Kruskal, W. (1966). An approximate statistical test for correlations between proportions. J. Geol., 74, 692-702.
Coakley, J. P. and Rust, B. R. (1968). Sedimentation in an Arctic lake. J. Sedimentary Petrology, 38, 1290-1300.
Connor, R. J. and Mosimann, J. E. (1969). Concepts of independence for proportions with a generalization of the Dirichlet distribution. J. Amer. Statist. Assoc., 64, 194-206.
Cox, D. R. (1966). Some procedures associated with the logistic qualitative response curve. In Research Papers in Statistics: Festschrift for J. Neyman (F. N. David, ed), pp. 57-71. New York: Wiley.
- (1970). The Analysis of Binary Data. London: Methuen.

Darroch, J. N. (1969). Null correlations for proportions. J. Math. Geol., 3, 467-483.
Darroch, J. N. and James, I. R. (1974). F-independence and null correlations of continuous, bounded-sum, positive variables. J. R. Statist. Soc. B, 36, 467-483.
Darroch, J. N. and Ratcliff, D. (1970). Null correlations for proportions II. J. Math. Geol., 2, 307-312.
- (1971). A characterization of the Dirichlet distribution. J. Amer. Statist. Assoc., 66, 641-643.
- (1978). No-association of proportions. J. Math. Geol., 10, 361-368.

Dawid, A. P. (1979). Conditional independence in statistical theory (with discussion). J. R. Statist. Soc. B, 41, 1-31.
Day, N. E. and Kerridge, D. F. (1967). A general maximum likelihood discriminant. Biometrics, 23, 313-323.
Deaton, A. S. (1978). Specification and testing in applied demand analysis. Econ. J., 88, 524-536.
Deaton, A. S. and Muellbauer, J. (1980). Economics and Consumer Behavior. New York: Cambridge University Press.
Gnanadesikan, R. and Kettenring, J. R. (1972). Robust estimates, residuals and outlier detection with multiresponse data. Biometrics, 28, 81-124.
Houthakker, H. S. (1960). Additive preferences. Econometrica, 28, 244-254.
James, I. R. (1975). Multivariate distributions which have beta conditional distributions. J. Amer. Statist. Assoc., 70, 681-684.
- (1981). Distributions associated with neutrality properties for random proportions. In Statistical Distributions in Scientific Work (C. Taillie, G. P. Patil and B. Baldessari, eds), vol. 4, pp. 125-136. Dordrecht, Holland, D. Reidel Publishing Company.
James, I. R. and Mosimann, J. E. (1980). A new characterization of the Dirichlet distribution through neutrality. Ann. Statist., 8, 183-189.
Jeffreys, H. (1961). Theory of Probability. Oxford University Press.
Johnson, N. L. (1949). Systems of frequency curves generated by methods of translation. Biometrika, 36, 149-176.
Johnson, N. L. and Kotz, S. (1972). Distributions in Statistics. Continuous Multivariate Distributions. Boston: Houghton Mifflin.
Kaiser, R. F. (1962). Composition and origin of glacial till Mexico and Kasoog quadrangles, New York. J. Sedimentary Petrology, 32, 502-513.
Leonard, T. (1973). A Bayesian method for histograms. Biometrika, 59, 581-589.
Leser, C. E. V. (1976). Income, household size and price changes 1953-1973. Oxford Bull. Econ. Statist., 38, 1-10.
McAlister, D. (1879). The law of the geometric mean. Proc. R. Soc., 29, 367.
Marquardt, D. W. (1963). An algorithm for least-squares estimation of non-linear parameters. J. SIAM, 11, 431-441.
Morrison, D. F. (1976). Multivariate Statistical Methods. New York: McGraw-Hill.
Mosimann, J. E. (1962). On the compound multinomial distribution, the multivariate \(\beta\)-distribution and correlations among proportions. Biometrika, 49, 65-82.
- (1970). Size allometry: size and shape variables with characterizations of the lognormal and generalized gamma distributions. J. Amer. Statist. Assoc., 65, 630-645.
- (1975a). Statistical problems of size and shape. I. Biological applications and basic theorems. In Statistical Distributions in Scientific Work (G. P. Patil, S. Kotz and J. K. Ord, eds), pp. 187-217, Dordrecht, Holland: D. Reidel Publishing Company.
- (1975b). Statistical problems of size and shape. II. Characterizations of the lognormal, gamma and Dirichlet distributions. In Statistical Distributions in Scientific Work (G. P. Patil, S. Kotz and J. K. Ord, eds), pp. 219-239. Dordrecht, Holland: D. Reidel Publishing Company.
Pearson, K. (1897). Mathematical contributions to the theory of evolution. On a form of spurious correlations which may arise when indices are used in the measurement of organs. Proc. R. Soc., \(\mathbf{6 0 ,} \mathbf{4 8 9 - 4 9 8 .}\)
Sarmanov, O. V. and Vistelius, A. B. (1959). On the correlation of percentage values. Doklady Akad. Nauk. SSSR, 126, 22-25.
Stephens, M. A. (1974). EDF statistics for goodness of fit and some comparisons. J. Amer. Statist. Assoc., 69, 730-737.
Thompson, R. N., Esson, J. and Duncan, A. C. (1972). Major element chemical variation in the Eocene lavas of the Isle of Skye, Scotland. J. Petrology, 13, 219-253.
Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. Ann. Math. Statist., 9, 60-62.
Working, H. (1943). Statistical laws of family expenditure. J. Amer. Statist. Assoc., 38, 43-56.