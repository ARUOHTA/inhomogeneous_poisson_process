\section*{MACROECOLOGICAL METHODS}
![](https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-01.jpg?height=260&width=274&top_left_y=231&top_left_x=227)

\title{
Accounting for imperfect detection and survey bias in statistical analysis of presence-only data
}

\author{
Robert M. Dorazio*
}

\begin{abstract}
Southeast Ecological Science Center, US Geological Survey, 7920 NW 71st Street, Gainesville, FL 32653, USA
\end{abstract}
*Correspondence: Robert Dorazio, US Geological Survey, Southeast Ecological Science Center, 7920 NW 71 Street, Gainesville, FL 32653, USA.
E-mail: bdorazio@usgs.gov

\begin{abstract}
Aim During the past decade ecologists have attempted to estimate the parameters of species distribution models by combining locations of species presence observed in opportunistic surveys with spatially referenced covariates of occurrence. Several statistical models have been proposed for the analysis of presence-only data, but these models have largely ignored the effects of imperfect detection and survey bias. In this paper I describe a model-based approach for the analysis of presence-only data that accounts for errors in the detection of individuals and for biased selection of survey locations.

Innovation I develop a hierarchical, statistical model that allows presence-only data to be analysed in conjunction with data acquired independently in planned surveys. One component of the model specifies the spatial distribution of individuals within a bounded, geographic region as a realization of a spatial point process. A second component of the model specifies two kinds of observations, the detection of individuals encountered during opportunistic surveys and the detection of individuals encountered during planned surveys.

Main conclusions Using mathematical proof and simulation-based comparisons, I demonstrate that biases induced by errors in detection or biased selection of survey locations can be reduced or eliminated by using the hierarchical model to analyse presence-only data in conjunction with counts observed in planned surveys. I show that a relatively small number of high-quality data (from planned surveys) can be used to leverage the information in presence-only observations, which usually have broad spatial coverage but may not be informative of both occurrence and detectability of individuals. Because a variety of sampling protocols can be used in planned surveys, this approach to the analysis of presence-only data is widely applicable. In addition, since the point-process model is formulated at the level of an individual, it can be extended to account for biological interactions between individuals and temporal changes in their spatial distributions.
\end{abstract}

\section*{Keywords}

Ecological niche model, N-mixture model, predictive biogeography, site occupancy model, spatial point process, species distribution model.

\section*{INTRODUCTION}

Species distribution models (SDMs) are used to predict the spatial distribution of a species, often as a function of spatially varying covariates. SDMs have a variety of uses (Elith \& Leathwick, 2009), but predicting the occurrence of a species within its potential geographic range is perhaps the most common (Scott et al., 2002). These predictions are particularly
relevant in conservation or management problems that require assessments of the effects of environmental modifications on a species' geographic distribution (Cabeza et al., 2004).

During the past decade ecologists have attempted to estimate the parameters of a SDM using only locations where a species is observed. These so-called presence-only observations are often made during opportunistic surveys and are recorded in museum collections or online databases. To fit a SDM, the presence-only
observations are supplemented with measures of potential covariates of species occurrence stored in geographic information systems or other geographic databases. These covariate measurements are generally available at a grid of locations that spans the study area.

Analyses of presence-only data, which include the presence locations and spatially referenced covariate measurements, are motivated, in part, by the difficulty and expense of conducting planned surveys of natural populations. As the region of interest grows in size and spatial complexity, so does the number of surveys required to predict species occurrence accurately. Because presence-only data usually have broad spatial coverage, they provide attractive sources of information for SDMs.

A variety of statistical models have been used to analyse presence-only data (Elith et al., 2006). Many of these models are strictly appropriate for the analysis of data observed in planned surveys, where presence or absence of individuals is observable at each survey location. The predictive accuracy of these models generally suffers when individuals are present at locations that have not been surveyed, as with presence-only data. In this case alternative models that allow for the presence of individuals at unsurveyed locations are generally more accurate in predicting the spatial distribution of species occurrences. These models include Maxent (Phillips et al., 2006; Elith et al., 2010), models of case-augmented binary outcomes (Lee et al., 2006; Lele \& Keim, 2006) and models of spatial point patterns (Warton \& Shepherd, 2010).

In the latter class of models, presence-only locations are modelled as a realization of a spatial point process - specifically, a Poisson process wherein the first-order intensity function specifies a SDM. This approach offers several advantages over analyses of presence-only data. First, the parameters of a pointprocess model are invariant to spatial scale, so the model can be used to predict the abundance or occurrence of individuals for any subregion located within the region of interest. In contrast, the parameters of other models (Maxent and case-augmented binary regression) vary with the spatial resolution of the data. Technical equivalences have recently been established between the parameters of Poisson process models and those of Maxent models (Fithian \& Hastie, 2013; Renner \& Warton, 2013) and case-augmented binary regression models (Dorazio, 2012). Specifically, as the spatial resolution of the data is increased, the scale-dependent parameters of the latter two models converge to those of Poisson process models; therefore, it can be argued that spatial point-process models provide a conceptual unification for classes of models whose parameters are not invariant to spatial scale.

One limitation of the point-process model of Warton \& Shepherd (2010) is that it fails to account for the effects of errors in the detection of individuals. This limitation is important because nearly all surveys of natural populations, including opportunistic surveys that produce presence-only observations, are prone to detection errors (Yoccoz et al., 2001; Chen et al., 2013). In addition, Dorazio (2012) proved that failure to account for imperfect detectability in models of presence-only data induces bias in estimates of SDMs when the covariates of
abundance are not distinct and stochastically independent of the covariates of detectability (see also Lahoz-Monfort et al., 2014). Bias in estimates of SDMs can also occur if the presence-only locations are unrepresentative of the region of interest (Phillips et al., 2009; Yackulic et al., 2013). For example, this source of bias might be present if survey locations are selected based on their accessibility or convenience. To alleviate these biases, Chakraborty et al. (2011) and Fithian \& Hastie (2013) proposed models based on a location-dependent thinning of a Poisson process; however, the parameters of these models are not identifiable unless the covariates of abundance are distinct and linearly independent of the covariates of detectability.

In this paper I propose a hierarchical, statistical model that allows presence-only data to be analysed in conjunction with data acquired independently in planned surveys. Previously, such data have been used only to validate the predictions of SDMs fitted to presence-only data (Newbold et al., 2010; Peterman et al., 2011; Gormley et al., 2013). Here I show that jointly modelling both kinds of data allows the parameters of a SDM to be estimated while accounting for the effects of imperfect detection and survey bias. In the hierarchy of the model one component specifies the SDM as a spatial point process and a second component, which is conditional on the first, specifies two kinds of observations: presence-only locations in opportunistic surveys and location-specific counts in planned surveys. I establish a set of restrictive conditions under which the parameters of a SDM can be estimated from presence-only data alone; however, I also show that these conditions may safely be ignored when the counts observed in planned surveys are analysed in conjunction with the presence-only data. I show that this approach is widely applicable owing to the variety of sampling protocols - including site-occupancy sampling - that can be used to conduct planned surveys.

\section*{MODELLING THE SPATIAL DISTRIBUTION OF A SPECIES AND ITS DETECTION IN OPPORTUNISTIC AND PLANNED SURVEYS}

In this section I describe a hierarchical model composed of two components. One component specifies the spatial distribution of individuals within a bounded geographic region that is relevant in the context of some scientific or management-related problem. This component is the true, but unknown, SDM and specifies how the limiting, expected density of individuals which will be given a mathematically precise definition - varies geographically as a function of one or more spatially varying covariates (such as measures of habitat quality).

The second component specifies models for two kinds of observations: (1) detections of individuals encountered during opportunistic surveys (i.e. presence-only observations) and (2) detections of individuals encountered during planned surveys of locations selected using a prescribed sampling design. These observations depend on several factors, such as the area of the surveyed locations, the number of individuals present at these locations, the effects of spatially varying covariates on an observer's detection ability and heterogeneity among observers in
detection skills. In the following subsections I describe both components of the hierarchical model.

\section*{Spatial point-process model of individual locations}

I assume that the spatial distribution of individuals - or, more correctly, of the activity centres of mobile individuals - may be modelled using a Poisson point process. This model has been used in the analysis of spatial capture-recapture data (Efford, 2004; Borchers \& Efford, 2008; Dorazio, 2013) and in the analysis of presence-only locations (Warton \& Shepherd, 2010). The latter analysis was intended as a SDM, but it did not consider the ramifications of detection errors or sources of survey bias. The hierarchical model proposed here is intended to correct this deficiency.

To formulate the SDM, I consider individuals that reside within a bounded geographic region \(B \subset \mathbb{R}^{2}\), where \(\mathbb{R}^{2}\) denotes the real plane. The choice of \(B\) is often driven by issues related to science, management or conservation of a species. I assume that the activity centres of these individuals are a realization of a Poisson point process parameterized by a first-order intensity function \(\lambda(\boldsymbol{s})\), where \(\boldsymbol{s}\) denotes a location (point) in \(B\). In the context of SDMs, \(\lambda(\boldsymbol{s})\) denotes the limiting, expected (E) density of individuals (number of individuals per unit area) at location \(\boldsymbol{s}\). That is, for a small region \(\mathrm{d} \boldsymbol{s}\) of area \(A(\mathrm{~d} \boldsymbol{s})\) centred at \(\boldsymbol{s}\) :
\(\lambda(\boldsymbol{s})=\lim _{A(\mathrm{~d} \Omega) \rightarrow 0} \mathrm{E}\{N(\mathrm{~d} \boldsymbol{s})\} / A(\mathrm{~d} \boldsymbol{s})\)
for the Poisson point process \(N\) (Cressie \& Wikle, 2011). The total number \(N(B)\) of individuals in the region is a Poisson random variable that depends on the mean intensity of the process over \(B, \mu(B)=\int_{B} \lambda(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\). In other words, \(N(B) \sim \operatorname{Poisson}(\mu(B))\). Furthermore, if \(B\) is partitioned into a set of disjoint (non-overlapping) subregions, say \(C_{1} \mathrm{U} \ldots \mathrm{U} C_{K}=B\), the number of individuals in each subregion also is a Poisson random variable. In other words, \(N\left(C_{k}\right) \sim \operatorname{Poisson}\left(\mu\left(C_{k}\right)\right)\), where \(\mu\left(C_{k}\right)=\int_{C_{k}} \lambda(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\) (Møller \& Waagepetersen, 2004; Illian et al., 2008). The number of individuals present in one subregion is also independent of the number in another subregion, a property that will be exploited in other components of the hierarchical model.

To specify a SDM, \(\boldsymbol{\lambda}(\boldsymbol{s})\) is formulated as a \(\log\)-linear function of unknown parameters \(\boldsymbol{\beta}\) and location-specific regressors \(\boldsymbol{x}(\boldsymbol{s})\) as follows: \(\log (\boldsymbol{\lambda}(\boldsymbol{s}))=\boldsymbol{\beta}^{\prime} \boldsymbol{x}(\boldsymbol{s})\) (the prime indicates the transpose of a matrix or vector.) The regressors \(\boldsymbol{x}(\boldsymbol{s})\) may be computed from covariates measured at a grid of locations that span \(B\) (e.g. measures of habitat quality). Therefore, accurate estimates of \(\boldsymbol{\beta}\) may be used to predict the abundance or occurrence of individuals for any location or any subregion within \(B\).

Let \(n\) denote the unknown total number of individuals in region \(B\). I have established that
\[
\operatorname{Pr}(N(B)=n)=\exp (-\mu(B))(\mu(B))^{n} / n!
\]

Let \(S_{i} \in B\) denote a random variable for the activity centre of the \(i\) th individual residing in region \(B(i=1, \ldots, n)\). Conditional on \(n\), the joint probability density of the \(n\) activity centres is
\[
f\left(\boldsymbol{s}_{1}, \boldsymbol{s}_{2}, \ldots, \boldsymbol{s}_{n} \mid n\right)=\prod_{i=1}^{n} \lambda\left(\boldsymbol{s}_{i}\right) / \mu(B)
\]
(Møller \& Waagepetersen, 2004; Illian et al., 2008). If \(n\) and \(\left\{\boldsymbol{s}_{i}\right\}\) were observable, their joint density
\(g\left(\boldsymbol{s}_{1}, \boldsymbol{s}_{2}, \ldots, \boldsymbol{s}_{n}, n\right)=\frac{\exp (-\mu(B))}{n!} \prod_{i=1}^{n} \lambda\left(\boldsymbol{s}_{i}\right)\),
(obtained by combining equations 1 and 2 ) could be used as a likelihood function for estimating \(\boldsymbol{\beta}\) (Warton \& Shepherd, 2010). However, only some of the \(n\) individuals are actually observed owing to errors in detection and survey biases. The second component of the hierarchical model (described in the following two subsections) allows \(\boldsymbol{\beta}\) to be estimated while accounting for unobserved individuals.

\section*{Detections of individuals in opportunistic surveys}

I make several assumptions in modelling the detections of individuals encountered during opportunistic surveys. First, I assume that each individual is detected independently by an observer with probability \(p(\boldsymbol{s})\) that depends only on the location \(\boldsymbol{s}\) of an individual. Therefore, \(p(\boldsymbol{s})\) includes both the observer's detection ability and choice of survey location. For example, \(p(\boldsymbol{s})\) will be zero if an observer does not survey individuals at location \(\boldsymbol{s}\) because the location is inaccessible. I also assume that the presence of other individuals encountered during an opportunistic survey has no effect on an observer's detection rate. Instead, the probability of detecting an individual located at \(\boldsymbol{s}\) is assumed to be a logit-linear function of unknown parameters \(\boldsymbol{\alpha}\) and location-specific regressors \(\boldsymbol{w}(\boldsymbol{s})\) as follows: \(\operatorname{logit}(p(\boldsymbol{s}))=\boldsymbol{\alpha}^{\prime} \boldsymbol{w}(\boldsymbol{s})\). The regressors \(\boldsymbol{w}(\boldsymbol{s})\) are assumed to be computable at all locations in \(B\) and may include features thought to influence an observer's detection ability or choice of survey location. For example, distance to the nearest road could be included in \(\boldsymbol{w}(\boldsymbol{s})\) if the observer is known to have chosen the survey locations based on their proximity to roads.

Although the identity of each observer potentially could be used as a regressor of \(p(\boldsymbol{s})\), especially in cases with few observers, most presence-only data include many observers and their identities are not always available. Therefore, I make the simplifying assumption that \(\boldsymbol{\alpha}\) contains a single intercept parameter. Although differences in \(p(\boldsymbol{s})\) may exist due to the differing abilities of observers, I assume here that these differences are small in comparison with the effects of spatially varying covariates on an observer's detection ability and choice of survey location.

A final assumption concerns the effects of movements of individuals around their activity centres. Unlike spatial capturerecapture data, presence-only data contain only a single observation (and location) for each individual. This leaves no information for estimating the magnitude of individual movements; therefore, I assume that individuals are detected at their activity centres, recognizing that this assumption is likely to be violated if individuals are highly mobile or possess large territories (that is, large relative to the size of region \(B\) ).

Given these assumptions, the detections of individuals in opportunistic surveys may be modelled as a location-dependent thinning of the Poisson point process described earlier. Let \(Y\) denote a binary random variable whose observed value indicates whether an individual residing in region \(B\) is detected \((Y=1)\) or not ( \(Y=0\) ) during incidental surveys of the region. I assume \(Y_{i} \mid \boldsymbol{s}_{i} \sim \operatorname{Bernoulli}\left(p\left(\boldsymbol{s}_{i}\right)\right)\), which specifies a marked point-process model of \(\left(y_{i}, \boldsymbol{s}_{i}\right)(i=1, \ldots, n)\). However, the individuals detected in opportunistic surveys all have the same mark (that is, \(Y=1\) for these individuals). Therefore, let \(M(B)\) denote a random variable for the number of individuals that are present and detected during opportunistic surveys of region \(B\). It is easily shown that
\[
\operatorname{Pr}(M(B)=m)=\exp (-v(B))(v(B))^{m} / m!
\]
where \(\boldsymbol{v}(B)=\int_{B} \lambda(\boldsymbol{s}) p(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\) and that \(M\) is a Poisson process resulting from an independent thinning of the original process \(N\) with probability \(1-p(s)\) (Møller \& Waagepetersen, 2004; Illian et al., 2008). It follows that \(v(B)=\mathrm{E}(M(B))\) and that the joint density of \(m\) and the locations of individuals detected in the opportunistic surveys is
\[
h\left(\boldsymbol{s}_{1}, \boldsymbol{s}_{2}, \ldots, \boldsymbol{s}_{m}, m\right)=\frac{\exp (-v(B))}{m!} \prod_{i=1}^{m} \lambda\left(\boldsymbol{s}_{i}\right) p\left(\boldsymbol{s}_{i}\right)
\]
where the first \(m\) of \(n\) locations are assumed to correspond to those of detected individuals. \({ }^{1}\) Note that the integral required to compute \(\nu(B)\) cannot be evaluated in closed form. In practice \(\nu(B)\) is approximated as a Riemann sum by partitioning \(B\) into a sufficiently fine grid.

The joint density of the observed data (i.e. equation 4 ) may be used as a likelihood function \(L(\boldsymbol{\beta}, \boldsymbol{\alpha})\) for estimating the parameters \(\boldsymbol{\beta}\) and \(\boldsymbol{\alpha}\). In later subsections I establish the conditions required for the identification of parameters based on this likelihood function.

\section*{Detections of individuals in planned surveys}

Several survey protocols allow the abundance of imperfectly detected individuals to be estimated from data observed in spatial sample units. Examples of these protocols include double-observer surveys, removal surveys and replicated pointcount surveys (Royle \& Dorazio, 2008, chapter 8). In practice the choice of protocol often depends on the species and on the methods used to detect or capture individuals.

Suppose one of these protocols is used to supplement the information obtained through opportunistic surveys. The idea is to select a representative sample of region \(B\) using a planned design and to survey individuals at each sample location using a protocol that is informative of both abundance and detectability. The number of locations in this sample will usually be small in comparison with the number of presence-only loca-

\footnotetext{
\({ }^{1}\) The order of the \(n\) observations can be changed without loss of generality.
}
tions observed in opportunistic surveys. To illustrate this approach, I describe here a model of replicated point counts at each location, though any of the survey protocols mentioned earlier can be applied without loss of generality.

For the purposes of sampling, assume that region \(B\) is partitioned into a finite number of disjoint sample units (e.g. quadrats, strip transects, etc.). Let \(C_{1}, C_{2}, \ldots, C_{K}\) denote a representative (but not necessarily random) sample of these units, and let \(A\left(C_{1}\right), A\left(C_{2}\right), \ldots, A\left(C_{K}\right)\) denote their areas. Under the assumptions of the SDM, the number \(N\left(C_{k}\right)\) of individual activity centres in sample unit \(C_{k}\) is a Poisson random variable with mean \(\mu\left(C_{k}\right)=\int_{C_{k}} \lambda(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\). If the values of spatially varying regressors \(x(\boldsymbol{s})\) are constant within a sample unit, which is often true in practice, the Poisson mean may be expressed as a function of the area of the unit and the fixed value of the regressors \(\boldsymbol{x}\left(C_{k}\right)\) as follows: \(\mu\left(C_{k}\right)=A\left(C_{k}\right) \exp \left(\boldsymbol{\beta}^{\prime} \boldsymbol{x}\left(C_{k}\right)\right)\). Thus, the expected number of individual activity centres in each sample unit increases with its area.

Suppose \(J_{k}\) independent point-count surveys of individuals in unit \(C_{k}\) are completed during a period where the number of individuals present in \(C_{k}\) remains constant. Independence can be achieved by using different observers or methods of detection, or by allowing sufficient time to elapse between successive surveys. These surveys yield a vector \(y_{k}=\left(y_{k 1}, y_{k 2}, \ldots, y_{k j}\right)^{\prime}\) of counts that contains the number of individuals detected during each survey of unit \(C_{k}\). Generally speaking, the number of individuals available to be detected in a sample unit depends on the locations of individual activity centres and on the spatial extent of individual movements about their activity centres. If the survey protocol is repeated over time, the number of individuals that are present and available to be detected in a sample unit can be estimated (Chandler et al., 2011). However, for many species temporal replication of the survey protocol is often physically or logistically unfeasible owing to sampling constraints. In the absence of such replication I assume that only individuals whose activity centres lie within a sample unit are available to be detected; thus \(N\left(C_{k}\right)\) denotes a random variable for the number of individuals present and available to be detected in unit \(C_{k}\). This assumption may not be valid for highly mobile species (large birds or mammals), but it will be satisfied for species whose movements are more limited (plants, clams, some amphibians, small nesting birds, etc.).

To specify the effects of imperfect detection on the observed point counts, I use the product-binomial model proposed by Royle (2004):
\[
\operatorname{Pr}\left(\boldsymbol{Y}_{k}=\boldsymbol{y}_{k} \mid N\left(C_{k}\right)=n_{k}\right)=\prod_{j=1}^{J_{k}}\binom{n_{k}}{y_{k j}} p_{k j}^{y_{k j}}\left(1-p_{k j}\right)^{n_{k}-y_{k j}}
\]
where \(p_{k j}\) is the conditional probability of detecting an individual given that it is present and available to be observed during the \(j\) th survey of sample unit \(C_{k}\). The detection probability \(p_{k j}\) may be specified as a function of unknown parameters and covariate measurements that differ among sample units or surveys. These covariates may include those used to model detections in opportunistic surveys, with the exception of the
covariates used to specify the effects of survey bias. The latter covariates are obviously unnecessary, since the locations of sample units in planned surveys are prescribed by design. I assume here that \(p_{k j}\) is a logit-linear function of unknown parameters \(\boldsymbol{\gamma}\) and regressors \(\boldsymbol{v}\left(C_{k}\right)\) whose values are constant within unit \(C_{k}\) : \(\operatorname{logit}\left(p_{k j}\right)=\boldsymbol{\gamma}^{\boldsymbol{\prime}} \boldsymbol{v}\left(C_{k}\right)\). This function implies that the probability of detecting an individual during planned surveys can differ among sample units (i.e., spatially) but not among surveys within a unit.

The unconditional probability of the point counts observed in sample unit \(C_{k}\) is obtained by marginalizing the joint density of \(n_{k}\) and \(\boldsymbol{y}_{k}\) over the admissible values of \(n_{k}\), as follows:
\[
\begin{aligned}
\operatorname{Pr}\left(\boldsymbol{Y}_{k}=\boldsymbol{y}_{k}\right)= & \sum_{n_{k}=\max \left(y_{k}\right)}^{\infty} \frac{\exp \left(-\mu\left(C_{k}\right)\right)\left(\mu\left(C_{k}\right)\right)^{n_{k}}}{n_{k}!} \\
& \prod_{j=1}^{J_{k}}\binom{n_{k}}{y_{k j}} p_{k j}^{y_{k j}}\left(1-p_{k j}\right)^{n_{k}-y_{k j}} .
\end{aligned}
\]

In practice, the infinite upper limit of summation is replaced with an integer that is sufficiently large to ignore the contributions of additional terms. A likelihood function for estimating the parameters \(\boldsymbol{\beta}\) and \(\boldsymbol{\gamma}\) from the point counts is
\[
L(\boldsymbol{\beta}, \boldsymbol{\gamma})=\prod_{k=1}^{K} \operatorname{Pr}\left(\boldsymbol{Y}_{k}=\boldsymbol{y}_{k}\right)
\]
which stems from the independence of abundances and counts among sample units. If two or more point counts are observed in each sample unit, the parameters of this model ( \(\boldsymbol{\beta}\) and \(\boldsymbol{\gamma}\) ) are identifiable (Royle, 2004). If only one point count is observed in each sample unit, the parameters of the model may be identifiable if some of the covariate measurements used as abundance regressors are distinct (or at least linearly independent) from the covariate measurements used as regressors of detection probability (Sólymos et al., 2012).

\section*{ESTIMATING SDMS FROM DETECTIONS OF INDIVIDUALS IN OPPORTUNISTIC SURVEYS}

In this section I identify requirements for estimating SDMs using only presence-only observations. First, if detection probability \(p\) is assumed to be constant at all locations (that is, if \(\boldsymbol{\alpha}=\alpha_{0}=\operatorname{logit}(p)\) ), then \(\beta_{0}\) and \(\alpha_{0}\) are not identifiable and the SDM is not estimable. This result is easily deduced from the log-likelihood function of this model's parameters:
\[
\begin{aligned}
\log \{L(\boldsymbol{\beta}, p)\}= & -\int_{B} \lambda(\boldsymbol{s}) p \mathrm{~d} \boldsymbol{s}+\sum_{i=1}^{m} \log \left(\lambda\left(\boldsymbol{s}_{i}\right) p\right) \\
= & -p \int_{B} \exp \left(\boldsymbol{\beta}_{0}+\tilde{\boldsymbol{\beta}}^{\prime} \tilde{\boldsymbol{x}}(\boldsymbol{s})\right) \mathrm{d} \boldsymbol{s} \\
& +m\left(\boldsymbol{\beta}_{0}+\log (p)\right)+\sum_{i=1}^{m} \tilde{\boldsymbol{\beta}^{\prime}} \tilde{\boldsymbol{x}}\left(\boldsymbol{s}_{i}\right) \\
= & -\exp \left(\boldsymbol{\beta}_{0}+\log (p)\right) \int_{B} \exp \left(\tilde{\boldsymbol{\beta}}^{\prime} \tilde{\boldsymbol{x}}(\boldsymbol{s})\right) \mathrm{d} \boldsymbol{s} \\
& +m\left(\boldsymbol{\beta}_{0}+\log (p)\right)+\sum_{i=1}^{m} \tilde{\boldsymbol{\beta}} \tilde{\boldsymbol{x}}\left(\boldsymbol{s}_{i}\right)
\end{aligned}
\]
where \(\tilde{\boldsymbol{\beta}} \tilde{\boldsymbol{x}}(\boldsymbol{s})\) corresponds to terms in the log-linear predictor of \(\lambda(s)\) that do not include \(\beta_{0}\). In this log-likelihood function the parameters \(\beta_{0}\) and \(p\) appear only as a sum \(\left(\beta_{0}+\log (p)\right)\); therefore, only the sum of these parameters is identified, and unique estimates of \(\beta_{0}\) and \(p\) cannot be obtained.

More generally, consider models of presence-only data that assume spatial variation in both \(\lambda\) and \(p\). For these models the log-likelihood function is
\[
\log (L(\boldsymbol{\beta}, \boldsymbol{\alpha}))=-\int_{B} \lambda(\boldsymbol{s}) p(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}+\sum_{i=1}^{m} \log \left(\lambda\left(\boldsymbol{s}_{i}\right) p\left(\boldsymbol{s}_{i}\right)\right)
\]
(equal to the logarithm of the right-hand side of equation 4, ignoring terms that don't include parameters). Mathematical proofs of parameter identifiability based on equation 6 are challenging because the integral cannot be evaluated in closed form and because neither the product of \(\lambda(\boldsymbol{s})\) and \(p(\boldsymbol{s})\), which equals
\(\lambda(\boldsymbol{s}) p(\boldsymbol{s})=\frac{\exp \left(\boldsymbol{\beta}^{\prime} \boldsymbol{x}(\boldsymbol{s})+\boldsymbol{\alpha}^{\prime} \boldsymbol{w}(\boldsymbol{s})\right)}{1+\exp \left(\boldsymbol{\alpha}^{\prime} \boldsymbol{w}(\boldsymbol{s})\right)}\),
nor its logarithm can be expressed as a linear combination of \(\boldsymbol{\beta}\) and \(\boldsymbol{\alpha}\).

One exception occurs when the probabilities of detecting individuals are relatively low (say, \(p(s)<0.2\) ) at all locations in \(B\). In this case, \(\lambda(\boldsymbol{s}) p(\boldsymbol{s}) \doteq \exp \left(\boldsymbol{\beta}^{\prime} \boldsymbol{x}(\boldsymbol{s})+\boldsymbol{\alpha}^{\prime} \boldsymbol{w}(\boldsymbol{s})\right)\), so the parameters \(\beta_{0}\) and \(\alpha_{0}\) only occur as a sum ( \(\beta_{0}+\alpha_{0}\) ) and again are not identifiable. Similarly, other elements of \(\boldsymbol{\beta}\) and \(\boldsymbol{\alpha}\) are not identified when the regressors \(\boldsymbol{x}\) and \(\boldsymbol{w}\) are linearly dependent (Fithian \& Hastie, 2013). This implies that in circumstances where \(p\) is relatively low at all locations in \(B\), the regressors of \(\lambda\) and \(p\) should include covariates that are distinct and not strongly correlated (either positively or negatively).

To assess the identifiability of parameters more generally for models in which \(\lambda\) and \(p\) are assumed to vary among locations, I derived the Fisher information matrix \(\boldsymbol{I}(\boldsymbol{\theta})\) for the \(\log\) likelihood function of equation 6, where \(\boldsymbol{\theta}=\left(\boldsymbol{\beta}^{\prime}, \boldsymbol{\alpha}^{\prime}\right)^{\prime}\). Recall that \(\boldsymbol{I}(\boldsymbol{\theta})\) equals the negative of the expected value of the Hessian matrix obtained by taking second partial derivatives of the loglikelihood function. If the Fisher information matrix has full rank (i.e. if \(\boldsymbol{I}(\boldsymbol{\theta})\) is positive definite and invertible), the parameters in \(\boldsymbol{\theta}\) are identifiable (Bowden, 1973). Therefore, parameter identifiability may be assessed by examining the condition number of \(\boldsymbol{I}(\boldsymbol{\theta})\), which equals the ratio of the largest to the smallest eigenvalues of this matrix.

The elements of \(\boldsymbol{I}(\boldsymbol{\theta})\) cannot be expressed in closed form, but they can be evaluated for fixed values of the parameters ( \(\boldsymbol{\beta}\) and \(\boldsymbol{\alpha}\) ) and the regressors ( \(\boldsymbol{x}\) and \(\boldsymbol{w}\) ) using numerical integration (see Appendix S1 in Supporting Information). The eigenvalues of \(\boldsymbol{I}(\boldsymbol{\theta})\) can then be calculated from these integral approximations. If \(\boldsymbol{I}(\boldsymbol{\theta})\) is positive definite, all of its eigenvalues are positive; however, if \(\boldsymbol{I}(\boldsymbol{\theta})\) is positive semi-definite (not full rank), the ratio of smallest to largest eigenvalues (reciprocal of the condition number) will be zero, though numerical evaluations of this ratio may not equal zero exactly. Therefore, this ratio provides a useful diagnostic for assessing whether the parameters of a SDM are identifiable.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-06.jpg?height=1402&width=834&top_left_y=174&top_left_x=171}
\captionsetup{labelformat=empty}
\caption{Figure 1 Spatial distributions of a covariate of individual density (upper panel) and a covariate of detection probability (lower panel).}
\end{figure}

To illustrate the utility of the Fisher information matrix, I compared two different presence-only models. The SDM of both models was identical - a single spatially varying covariate (Fig. 1, upper panel) was used to predict \(\lambda(\boldsymbol{s})\) as follows:
\[
\log (\lambda(s))=\log (8000)+0.5 x(s) .
\]

The covariate measurements were centred and scaled so that \(x\) had zero mean and unit variance. The two models of presenceonly data differed in their specification of spatially varying detection probabilities. In one model a single covariate \(w\) (Fig. 1, lower panel) whose values were computed independently of \(x\) was used to predict \(p(s)\) as follows:
\[
\operatorname{logit}(p(s))=\alpha_{0}-1.0 w(s) .
\]

The covariate measurements were centred and scaled so that \(w\) had zero mean and unit variance. The parameter \(\alpha_{0}\) was

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-06.jpg?height=788&width=819&top_left_y=174&top_left_x=1052}
\captionsetup{labelformat=empty}
\caption{Figure 2 Assessment of parameter identifiability for two point-process models of presence-only data. One model is based on independent regressors of abundance and detection (solid line); the other model is based on identical regressors (dashed line). Upper panel: reciprocal of the condition number of the Fisher information matrix is plotted as a function of \(\alpha_{0}\) (the logit-scale detection parameter). Lower panel: expected number of individuals detected in region \(B\) is plotted as a function of \(\alpha_{0}\).}
\end{figure}
assigned values ranging from -5 to 5 so that detection probabilities at the average value of the covariate ranged from very low values (near zero) to very high values (near one). In the second model the covariate \(x\), which was used to predict \(\lambda(\boldsymbol{s})\), was also used to predict \(p(\boldsymbol{s})\), that is:
\[
\operatorname{logit}(p(s))=\alpha_{0}-1.0 x(s) .
\]

The values of \(\boldsymbol{\alpha}_{0}\) and \(\boldsymbol{\alpha}_{1}\) used in both models of the detection probabilities were identical. The second model was intended to be representative of a worst-case scenario wherein the regressors of \(\lambda(s)\) and \(p(s)\) were perfectly correlated.

Whereas the expected number of individuals detected in region \(B(\nu(B))\) increased from 242 to 35,470 over the assumed range of \(\alpha_{0}\) values, the difference in \(\nu(B)\) between models for any single value of \(\alpha_{0}\) was minor in comparison (Fig. 2). In contrast, the ratio of smallest to largest eigenvalues of \(\boldsymbol{I}(\boldsymbol{\theta})\) signalled striking differences in the identifiability of the parameters of these two models. For example, regardless of the value of \(\alpha_{0}\), the ratio of smallest to largest eigenvalues was nearly constant (approximately zero) when identical regressors were assumed for \(\lambda(s)\) and \(p(s)\) (Fig. 2). For the other model, which assumed independent regressors, the ratio of the smallest to the largest eigenvalues was highest at \(\alpha_{0}=0\) (corresponding to a detection probability of 0.5 at the average value of the covariate) and declined to approximately zero at extreme values of \(\alpha_{0}(<-4\) or \(>4)\). The reason for this pattern is clear. At extreme values of \(\alpha_{0}\), detection probabilities are nearly constant over \(B\) [i.e. \(p(s)\) is
nearly zero (if \(\alpha_{0}<-4\) ) or one (if \(\alpha_{0}>4\) ) at every location in \(B\) ]. At these extremes the parameters \(\beta_{0}\) and \(\alpha_{0}\) are not identifiable and cannot be estimated uniquely, as established earlier.

This example illustrates that presence-only data may contain only limited information about the parameters of a SDM. If at least some of the regressors used to specify spatial differences in \(\lambda\) and \(p\) are not independent, SDMs will be unidentified. The parameters of SDMs can also be difficult to estimate if detection probabilities are either very small or very large at all locations. However, in these circumstances the ratio of smallest to largest eigenvalues of \(\boldsymbol{I}(\boldsymbol{\theta})\) provides a useful diagnostic. This ratio is also invariant to the magnitude of \(\beta_{0}\) (see Appendix S1) and is therefore relevant to the analysis of presence-only observations from populations of all sizes.

\section*{ESTIMATING SDMS FROM DETECTIONS OF INDIVIDUALS IN OPPORTUNISTIC AND PLANNED SURVEYS}

One way of overcoming the limited information in presenceonly observations is to analyse these data in conjunction with the detections of individuals in planned surveys. A joint analysis of these two kinds of data leverages the spatial coverage of presence-only observations, which is usually large, with the strength of information about abundance and detection in planned surveys.

As an example of this approach, consider a joint analysis of presence-only data and replicated point-count survey data. The likelihood functions for these two data sets (i.e. equations 4 and 5) both include the parameter \(\boldsymbol{\beta}\), which determines the SDM. Observations in opportunistic and planned surveys are obtained independently, so these two data sets can be analysed together using the following likelihood function:
\[
L(\boldsymbol{\beta}, \boldsymbol{\alpha}, \boldsymbol{\gamma})=L(\boldsymbol{\beta}, \boldsymbol{\alpha}) \times L(\boldsymbol{\beta}, \boldsymbol{\gamma}) .
\]

Estimates of \(\boldsymbol{\beta}\) obtained by maximizing equation 10 should have less bias and greater precision than estimates obtained by maximizing either equation 4 or equation 5 separately. In this section I describe simulation studies that illustrate the benefits of a joint analysis of data observed in opportunistic and planned surveys. I designed the simulation studies to build upon the models of presence-only data described in the previous section and to demonstrate that addition of planned survey data (in this case, point counts) can reduce or eliminate the estimation problems induced by parameter unidentifiability.

\section*{Design of simulation studies}

In the simulation studies presence-only data were simulated using the two models described earlier. That is, in both models \(\lambda(\boldsymbol{s})\) was specified as a function of \(x(\boldsymbol{s})\) using equation 7 . In one model \(p(\boldsymbol{s})\) was specified as a function of \(w(\boldsymbol{s})\) (using equation 8 ); in the other model \(p(\boldsymbol{s})\) was specified as a function of \(x(\boldsymbol{s})\) (using equation 9). The parameter vector \(\boldsymbol{\alpha}\) was assigned the same value in both models: \(\boldsymbol{\alpha}=(-1,-1)^{\prime}\).

The locations of individuals were simulated within a square region \(B\) using the covariate measurements shown in Fig. 1. The spatial distribution of \(\lambda\) and of the expected densities of detections (one for each model) are shown in Fig. 3. Based on these distributions, the expected number of individuals present in region \(B\) was 35,857 , and the expected number of individuals detected in region \(B\) equalled either 11,251 (for the model where detection depended on \(w\) ) or 7922 (for the model where detection depended on \(x\) ). Each simulated presence-only data set (i.e. the number and locations of detected individuals) was computed as a realization of a thinned Poisson process, as described earlier.

To simulate observations from point-count surveys, region \(B\) was partitioned into a grid of 10,000 square quadrats of equal size. Samples of \(K=50,100,200,400\) or 800 quadrats were selected randomly (without replacement), and \(J=4\) independent point-count surveys were conducted in each of these quadrats. The sample sizes were deliberately chosen to be small ( \(0.5-8 \%\) of the total number of quadrats) relative to the expected number of presence-only locations to simulate the usual situation where presence-only observations greatly outnumber the locations visited in planned surveys. The covariates of detectability were identical to those used in the models of presence-only data, i.e. \(v(C)=\int_{C} w(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\) for one model, and \(v(C)=\int_{C} x(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\) for the other model. The logit-scale detection parameters were identical in both point-count models: \(\gamma=(0\), \(-1)^{\prime}\). Note that the parameters \(\alpha_{1}\) and \(\gamma_{1}\) were assigned the same value to ensure that the effects of covariates on detection probability were the same regardless of whether an individual was encountered during an opportunistic survey or a planned survey.

The spatial distributions of the expected number of individuals per quadrat and of the expected number of individuals detected per quadrat during a single survey are shown in Fig. 4. The similarity between these distributions and those shown in Fig. 3 is not surprising given the similarity in the detection probability models assumed for opportunistic and planned surveys. Each simulated set of point counts was computed by aggregating the realized locations of individuals in \(B\) into quadrats, by selecting a random sample of these quadrats, and by taking \(J\) independent binomial draws from the individuals present in each sampled quadrat. This sequence of steps produced \(J\) point counts for each sampled quadrat.

A total of 1000 data sets containing both presence-only observations and point counts were simulated for each of the two models. Maximum-likelihood estimates of \(\boldsymbol{\beta}\) (parameters of the SDM) were computed for each simulated data set by fitting the model of presence-only data (using equation 6 ), the model of point counts (using equation 5) and the model of combined data (presence-only observations and point counts) (using equation \(10)\). Operating characteristics of the estimators (such as bias and variance) were estimated from these simulation results and compared to illustrate the inferential benefits of including both presence-only observations and point counts in the analysis.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-08.jpg?height=2028&width=805&top_left_y=169&top_left_x=176}
\captionsetup{labelformat=empty}
\caption{Figure 3 Spatial distributions of the expected density of individuals (upper panel) and the expected densities of detections in opportunistic surveys as a function of covariate \(w\) (middle panel) or covariate \(x\) (lower panel). Superimposed points indicate a single realization from each distribution.}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-08.jpg?height=2141&width=820&top_left_y=171&top_left_x=1046}
\captionsetup{labelformat=empty}
\caption{Figure 4 Spatial distributions of the expected number of individuals per sample unit (upper panel) and the expected number of individuals detected per survey as a function of covariate \(w\) (middle panel) or covariate \(x\) (lower panel).}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-09.jpg?height=848&width=1113&top_left_y=186&top_left_x=210}
\captionsetup{labelformat=empty}
\caption{Figure 5 Operating characteristics of maximum-likelihood estimators of the species distribution model parameters \(\beta_{0}\) and \(\beta_{1}\). Detections of individuals were assumed to depend on covariate \(w\), which was independent of the covariate used in the species distribution model. Symbols indicate that estimates were obtained by fitting the model of presence-only data (triangles), the model of point counts (circles) or the model of combined data (squares). Error bars indicate 95\% confidence intervals.}
\end{figure}

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/c0478a30-1dbf-4efd-ae4d-58e69105daaa-09.jpg?height=850&width=1117&top_left_y=1116&top_left_x=210}
\captionsetup{labelformat=empty}
\caption{Figure 6 Operating characteristics of maximum-likelihood estimators of the species distribution model parameters \(\beta_{0}\) and \(\beta_{1}\). Detections of individuals were assumed to depend on covariate \(x\), which was identical to the covariate used in the species distribution model. Symbols indicate that estimates were obtained by fitting the model of presence-only data (triangles), the model of point counts (circles) or the model of combined data (squares). Error bars indicate 95\% confidence intervals.}
\end{figure}

Appendix S2 contains the R code (R Core Team, 2014) that was used to generate and analyse the simulated data sets.

\section*{Results of simulation studies}

For the model in which detections of individuals depended only on the covariate \(w\), presence-only based estimators of \(\boldsymbol{\beta}\) appeared to have negligibly small bias (Fig. 5). When presenceonly observations and point counts were analysed together, the uncertainty in estimates of \(\beta_{0}\) was reduced considerably, while uncertainty in the estimates of \(\beta_{1}\) was unchanged.

Conspicuously different results were obtained for the model in which detections of individuals depended on the same covariate used in the SDM (i.e. covariate \(x\) ). For this model presence-only based estimators of \(\beta_{0}\) and \(\beta_{1}\) were strongly biased and highly variable (Fig. 6). This result was expected because the model parameters are not well identified, as established earlier. More surprising were the inferential benefits obtained by including point counts in the analysis. Adding point counts from as few as 50 quadrats ( \(0.5 \%\) of the sample frame) to the analysis dramatically reduced the bias and the uncertainty in the estimates of \(\beta_{0}\) and \(\beta_{1}\). The estimators appeared to exhibit
virtually no bias in point-count samples of 200 or more quadrats.

Additional simulation-based comparisons (not shown) using different values of \(\beta_{0}\) and \(J\) produced results that were qualitatively similar to those illustrated in Figs 5 and 6, i.e. bias and uncertainty of presence-only based estimators of SDMs can be reduced considerably by analysing presence-only observations and point counts together in a single model. An important caveat, however, is that \(J=4\) replicate point-count surveys in each quadrat were required to achieve these inferential improvements. Bias and uncertainty of the presence-only based estimators of SDMs were not reduced if only one point-count survey \((J=1)\) was conducted per quadrat. The reason, of course, is that the parameters of the point-count model cannot be identified when \(J=1\) unless some of the covariate measurements used as abundance regressors are distinct (or at least linearly independent) from the covariate measurements used as regressors of detection probability (Sólymos et al., 2012).

\section*{DISCUSSION}

Several statistical models have been proposed for the analysis of presence-only data, but they have largely ignored the effects of imperfect detectability and survey bias (Dorazio, 2012; Lahoz-Monfort et al., 2014). In this paper I have shown that the bias in estimates of SDMs induced by detection errors or survey bias can be reduced or eliminated by modelling presence-only data in conjunction with counts observed in planned surveys. In this modelling approach the SDM is specified using a spatial point process and the observable data (presence-only locations and counts) are specified conditional on a realization of this process.

If this approach is adopted but only the presence-only data are analysed, the parameters of a SDM may not be identifiable. For example, if the probability of detecting an individual is identical at all locations, the parameters of a SDM are not identifiable. On the other hand, if the detection probability differs among locations, the parameters of a SDM can be estimated if some of the regressors used to specify spatial differences in individual density are linearly independent of the regressors of detection probability. This restriction is similar to that identified by Fithian \& Hastie (2013), who proposed a different thinned Poisson process model for the analysis of presence-only data. In practice, I have shown that the condition number of the Fisher information matrix (Appendix S1) can be used to assess the identifiability of a model's parameters.

The simulation study I conducted illustrates how including counts from planned surveys in the analysis of presence-only data generally improves inferences about the parameters of a SDM. In this way a relatively small number of high-quality data - which are informative about both detection and abundance of individuals - can be used to leverage the information in presence-only observations. This approach follows a recommendation of Phillips \& Elith (2013) that additional data should be collected if estimates of absolute species occurrences are desired.

Several benefits, both conceptual and practical, stem from using a point-process model as a SDM. As shown earlier, the parameters of a point-process model are defined on an areal basis and are invariant to spatial scale; thus, predictions may be computed for any summary of the spatial distribution of individuals within the region of interest. For example, under the assumptions of the Poisson process model that I described, the expected abundance \(\mathrm{E}(N(C))\) and occurrence \(\operatorname{Pr}(N(C)>0)\) of individuals in any subregion \(C \subset B\) are defined implicitly as follows:
\(\mathrm{E}(N(C))=\mu(C)\)
\[
\operatorname{Pr}(N(C)>0)=1-\exp (-\mu(C))
\]
where \(\mu(C)=\int_{C} \lambda(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\) is a function of the model's parameters \(\boldsymbol{\beta}\) and the values of the spatially varying regressors \(\boldsymbol{x}\) in \(C\). Phillips \& Elith, (2013, p. 1411) assert that the Poisson process model cannot be used to estimate the probability of presence, which equals \(\operatorname{Pr}(N(C)>0)\). That claim is obviously incorrect. Furthermore, by defining occurrence probability as a function of \(\mu(C)\), any prediction of occurrence probability increases automatically with an increase in the area of subregion \(C\). This property is not shared by models of occurrence, such as Maxent, where occurrence probability is defined per grid cell and depends on the spatial resolution used in the analysis (Renner \& Warton, 2013). For the same reason conventional siteoccupancy models (MacKenzie et al., 2002; Tyre et al., 2003) are limited for use as SDMs, i.e. the interpretation of occurrence probability in site-occupancy models depends on the spatial resolution used in the analysis. Because of this limitation, Efford \& Dawson (2012) suggested that site-occupancy models be modified to account for spatial resolution and for the movements of individuals within their home ranges.

One practical benefit of using a point-process model as a SDM is that it clarifies how the location-specific measurements and interpolations of the covariates should be used in the analysis. In the SDM literature the number and spatial resolution of covariate measurements and interpolations needed in statistical analyses of presence-only data have been subjects of considerable debate (Pearce \& Boyce, 2006; Guisan et al., 2007; Elith \& Leathwick, 2009; Barbet-Massin et al., 2012). However, as Warton \& Shepherd (2010) and I have shown, the role of the covariate values is clear when using a point-process model as a SDM. Both measured and interpolated values of covariates over the entire region of interest \(B\) are used when \(B\) is partitioned to approximate the integral in the likelihood function as a Riemann sum. The spatial resolution required for accurate approximation is easily determined by increasing the resolution until the maximized value of the likelihood function stabilizes.

In this paper I used a model of point counts to illustrate the benefits of adding information from planned surveys to the analysis of presence-only data. Similarly, information from alternative sampling protocols, such as double-observer sampling, removal sampling or even capture-recapture sampling, can be used to improve the analysis of presence-only data. These
counts are then modelled conditional on the latent (unobserved) abundance of individuals within a unit (say, \(C_{k}\) ). This idea is an extension of Royle's (2004) \(N\)-mixture model of point counts and has been used to analyse many kinds of spatially referenced counts (Royle \& Dorazio, 2008, Chapter 8). The key is to exploit the variation in abundance among sample units and to specify the effects of the size, location, and habitat characteristics of each unit in the model of abundance. If abundance within sample units is not too high, this approach also can be applied using a site-occupancy survey protocol. For example, suppose \(J_{k}(>1)\) independent presence-absence surveys are conducted in unit \(C_{k}\). Let \(Z_{k j}\) denote a binary random variable whose observed value indicates whether one or more individuals were detected ( \(Z_{k j}=1\) ) or not ( \(Z_{k j}=0\) ) during the \(j\) th survey of sample unit \(C_{k}\). As with conventional site-occupancy models, \(Z_{k j}\) is modelled conditional on the latent presence of individuals in unit \(C_{k}\) as follows:
\[
Z_{k j} \mid N\left(C_{k}\right)=n_{k} \sim \operatorname{Bernoulli}\left(q_{j k} I\left(n_{k}>0\right)\right)
\]
where \(q_{k j}\) is the conditional probability of detection during the \(j\) th survey given that one or more individuals are present in unit \(k\) ( \(I(e)\) denotes the indicator function, which equals one if Boolean argument \(e\) is true and zero if \(e\) is false). Unlike conventional site-occupancy models, the probability of occurrence in unit \(C_{k}\) (say, \(\psi\left(C_{k}\right)\) ) is automatically scaled for the size of the sample unit because \(\psi\left(C_{k}\right)=\operatorname{Pr}\left(N\left(C_{k}\right)>0\right)=1-\exp \left(-\mu\left(C_{k}\right)\right)\) and \(\mu\left(C_{k}\right)=\int_{C_{k}} \lambda(\boldsymbol{s}) \mathrm{d} \boldsymbol{s}\) increases with the area of \(C_{k}\). The parameters of this site-occupancy model are therefore invariant to spatial scale. An alternative approach, proposed by Royle \& Nichols (2003), is to specify a model wherein detection probability \(q_{j k}\) is assumed to increase with abundance \(n_{k}\) as follows: \(q_{k j}=1-\left(1-p_{k j}\right)^{n_{k}}\). The estimable parameter \(p_{j k}\) is the probability of detection per individual during the \(j\) th survey. In this approach \(Z_{k j}\) is modelled conditional on the latent abundance of individuals in unit \(C_{k}\) as follows:
\[
Z_{k j} \mid N\left(C_{k}\right)=n_{k} \sim \operatorname{Bernoulli}\left[1-\left(1-p_{k j}\right)^{n_{k}}\right] .
\]

This site-occupancy model is closely related to Royle's (2004) model of point counts. In fact, the two models provide mathematically equivalent estimators of site occupancy (Dorazio, 2007).

\section*{CONCLUSIONS AND RECOMMENDATIONS}

Estimates of SDMs obtained in analyses of presence-only data are vulnerable to biases induced by detection errors and survey bias. These biases can be reduced or eliminated by: (1) including covariates of detection and survey bias in the presence-only model and (2) analysing counts observed in planned surveys in conjunction with the presence-only data. These planned surveys must include a sampling protocol that is informative of both abundance and detectability. For some protocols (doubleobserver, removal or capture-recapture sampling) this requirement is implied; in others (independent point counts or
presence-absence samples) two or more replicates are needed at each survey location to obtain multiple observations of the individuals present at each location.

The benefits of planned surveys in species distribution modelling may depend on the species. For example, planned surveys of highly mobile species (large birds and mammals) can be problematic if movements of individuals make them available for detection in more than one sample unit. Combining presence-only data with observations in planned surveys is more likely to benefit species whose movements are more limited relative to the size of a sample unit.

\section*{Extensions of SDMs}

The hierarchical model developed in this paper can be extended to address a variety of inference problems in species distribution modelling. For example, several authors (Guisan \& Thuiller, 2005; Goodsoe \& Harmon, 2012; Higgins et al., 2012; Kissling et al., 2012; Wisz et al., 2013) have noted that biological interactions between individuals (e.g. competitive or predator-prey interactions) should be included in SDMs along with the effects of spatially varying habitat characteristics. Spatial point-process models have been used to infer the effects of competitive interactions within and among species of plants or ants (Högmander \& Särkkä, 1999; Grabarnik \& Särkkä, 2004; Wiegand et al., 2007a, b), but to my knowledge point-process models have not been formulated to specify the effects of these interactions mechanistically. Doing this, while accounting for species- and location-specific differences in detectability and survey bias, will no doubt present challenges for the analyst. Another important class of inference problems requires the construction of dynamic (i.e. space-time) point-process models. Much of the SDM literature is limited to the construction of static (purely spatial) models, but the need for dynamic, process-based models is clearly evident and growing (Guisan \& Thuiller, 2005; Elith \& Leathwick, 2009; Dormann et al., 2012). Models are needed to predict changes in the spatial distribution of species due to changes in climate, habitat, levels of disturbance and abundance of non-indigenous species. Point-process models are also needed to predict changes in the spatial distribution of exotic species following their initial introduction to a region of interest. Considerable data may be required to fit these dynamic models, but recent progress with invasive plant species (Balderama et al., 2012) illustrates the feasibility of this approach. I anticipate that the construction of point-process models driven by scientific needs and biological mechanisms will greatly enhance the inferences and predictions of future SDMs.

\section*{ACKNOWLEDGEMENTS}

Editorial suggestions from the editor, C. Yackulic, and an anonymous referee greatly improved an earlier draft of this article. Any use of trade, product or firm names is for descriptive purposes only and does not imply endorsement by the US government.

\section*{REFERENCES}

Balderama, E., Schoenberg, F.P., Murray, E. \& Rundel, P.W. (2012) Application of branching models in the study of invasive species. Journal of the American Statistical Association, 107, 467-476.
Barbet-Massin, M., Jiguet, F., Albert, C.H. \& Thuiller, W. (2012) Selecting pseudo-absences for species distribution models: how, where and how many? Methods in Ecology and Evolution, 3, 327-338.
Borchers, D.L. \& Efford, M.G. (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics, 64, 377-385.
Bowden, R. (1973) The theory of parametric identification. Econometrica, 41, 1069-1074.
Cabeza, M., Araújo, M.B., Wilson, R.J., Thomas, C.D., Cowley, M.J.R. \& Moilanen, A. (2004) Combining probabilities of occurrence with spatial reserve design. Journal of Applied Ecology, 41, 252-262.
Chakraborty, A., Gelfand, A.E., Wilson, A.M., Latimer, A.M. \& Silander, J.A. (2011) Point pattern modelling for degraded presence-only data over large regions. Applied Statistics, 60, 757-776.
Chandler, R.B., Royle, J.A. \& King, D.I. (2011) Inference about density and temporary emigration in unmarked populations. Ecology, 92, 1429-1435.
Chen, G., Kéry, M., Plattner, M., Ma, K. \& Gardner, B. (2013) Imperfect detection is the rule rather than the exception in plant distribution studies. Journal of Ecology, 101, 183191.

Cressie, N. \& Wikle, C.K. (2011) Statistics for spatio-temporal data. John Wiley and Sons, Hoboken, NJ.
Dorazio, R.M. (2012) Predicting the geographic distribution of a species from presence-only data subject to detection errors. Biometrics, 68, 1303-1312.
Dorazio, R.M. (2007) On the choice of statistical models for estimating occurrence and extinction from animal surveys. Ecology, 88, 2773-2782.
Dorazio, R.M. (2013) Bayes and empirical Bayes estimators of abundance and density from spatial capture-recapture data. PLoS ONE, 8, e84017.
Dormann, C.F., Schymanski, S.J., Cabral, J., Chuine, I., Graham, C., Hartig, F., Kearney, M., Morin, X., Römermann, C., Schröder, B. \& Singer, A. (2012) Correlation and process in species distribution models: bridging a dichotomy. Journal of Biogeograpy, 39, 2119-2131.
Efford, M. (2004) Density estimation in live-trapping studies. Oikos, 106, 598-610.
Efford, M.G. \& Dawson, D.K. (2012) Occupancy in continuous habitat. Ecosphere, 3, art. 32, http://dx.doi.org/10.1890/ES1100308.1.

Elith, J. \& Leathwick, J.R. (2009) Species distribution models: ecological explanation and prediction across space and time. Annual Review of Ecology, Evolution, and Systematics, 40, 677697.

Elith, J., Graham, C.H., Anderson, R.P. et al. (2006) Novel methods improve prediction of species' distributions from occurrence data. Ecography, 29, 129-151.
Elith, J., Phillips, S.J., Hastie, T., Dudik, M., Chee, Y.E. \& Yates, C.J. (2010) A statistical explanation of MaxEnt for ecologists. Diversity and Distributions, 17, 43-57.
Fithian, W. \& Hastie, T. (2013) Finite-sample equivalence in statistical models for presence-only data. Annals of Applied Statistics, 7, 1917-1939.
Goodsoe, W. \& Harmon, L.J. (2012) How do species interactions affect species distribution models? Ecography, 35, 811-820.
Gormley, A.M., Forsyth, D.M., Griffioen, P., Lindeman, M., Ramsey, D.S.L., Scroggie, M.P. \& Woodford, L. (2013) Using presence-only and presence-absence data to estimate the current and potential distributions of established invasive species. Journal of Applied Ecology, 48, 25-34.
Grabarnik, P. \& Särkkä, A. (2004) Modelling the spatial structure of forest stands by multivariate point processes with hierarchical interactions. Ecological Modelling, 220, 1232-1240.
Guisan, A. \& Thuiller, W. (2005) Predicting species distribution: offering more than simple habitat models. Ecology Letters, 8, 993-1009.
Guisan, A., Graham, C.H., Elith, J. et al. (2007) Sensitivity of predictive species distribution models to change in grain size. Diversity and Distributions, 13, 332-340.
Higgins, S.I., O'Hara, R.B. \& Römermann, C. (2012) A niche for biology in species distribution models. Journal of Biogeography, 39, 2091-2095.
Högmander, H. \& Särkkä, A. (1999) Multitype spatial point patterns with hierarchical interactions. Biometrics, 55, 10511058.

Illian, J., Penttinen, A., Stoyan, H. \& Stoyan, D. (2008) Statistical analysis and modelling of spatial point patterns. John Wiley and Sons, Chichester.
Kissling, W.D., Dormann, C.F., Groeneveld, J., Hickler, T., Kühn, I., McInerny, G.J., Montoya, J.M., Römermann, C., Schiffers, K., Schurr, F.M., Singer, A., Svenning, J.C., Zimmermann, N.E. \& O'Hara, R.B. (2012) Towards novel approaches to modelling biotic interactions in multispecies assemblages at large spatial extents. Journal of Biogeography, 39, 21632178.

Lahoz-Monfort, J.J., Guillera-Arroita, G. \& Wintle, B.A. (2014) Imperfect detection impacts the performance of species distribution models. Global Ecology and Biogeography, 23, 504515.

Lee, A.J., Scott, A.J. \& Wild, C.J. (2006) Fitting binary regression models with case-augmented samples. Biometrika, 93, 385397.

Lele, S.R. \& Keim, J.L. (2006) Weighted distributions and estimation of resource selection probability functions. Ecology, 87, 3021-3028.
MacKenzie, D.I., Nichols, J.D., Lachman, G.B., Droege, S., Royle, J.A. \& Langtimm, C.A. (2002) Estimating site occupancy rates when detection probabilities are less than one. Ecology, 83, 2248-2255.

Møller, J. \& Waagepetersen, R.P. (2004) Statistical inference and simulation for spatial point processes. Chapman and Hall, Boca Raton, FL.
Newbold, T., Reader, T., El-Gabbas, A., Berg, W., Shohdi, W.M., Zalat, S., El Din, S.B. \& Gilbert, F. (2010) Testing the accuracy of species distribution models using species records from a new field survey. Oikos, 119, 1326-1334.
Pearce, J.L. \& Boyce, M.S. (2006) Modelling distribution and abundance with presence-only data. Journal of Applied Ecology, 43, 405-412.
Peterman, W.E., Crawford, J.A. \& Kuhns, A.R. (2011) Using species distribution and occupancy modeling to guide survey efforts and assess species status. Journal for Nature Conservation, 21, 114-121.
Phillips, S.J. \& Elith, J. (2013) On estimating probability of presence from use-availability or presence-background data. Ecology, 94, 1409-1419.
Phillips, S.J., Anderson, R.P. \& Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling, 190, 231-259.
Phillips, S.J., Dudik, M., Elith, J., Graham, C.H., Lehmann, A., Leathwick, J. \& Ferrier, S. (2009) Sample selection bias and presence-only distribution models: implications for background and pseudo-absence data. Ecological Applications, 19, 181-197.
R Core Team (2014) R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna.
Renner, I.W. \& Warton, D.I. (2013) Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics, 69, 274-281.
Royle, J.A. (2004) \(N\)-mixture models for estimating population size from spatially replicated counts. Biometrics, 60, 108-115.
Royle, J.A. \& Dorazio, R.M. (2008) Hierarchical modeling and inference in ecology. Academic Press, Amsterdam.
Royle, J.A. \& Nichols, J.D. (2003) Estimating abundance from repeated presence-absence data or point counts. Ecology, 84, 777-790.
Scott, J.M., Heglund, P.J., Morrison, M.L., Haufler, J.B., Raphael, M.G., Wall, W.A. \& Samson, F.B. (2002) Predicting species occurrences: issues of accuracy and scale. Island Press, Washington, DC.
Sólymos, P., Lele, S. \& Bayne, E. (2012) Conditional likelihood approach for analyzing single visit abundance survey data in the presence of zero inflation and detection error. Environmetrics, 23, 197-205.
Tyre, A.J., Tenhumberg, B., Field, S.A., Niejalke, D., Parris, K. \& Possingham, H.P. (2003) Improving precision and reducing bias in biological surveys: estimating false-negative error rates. Ecological Applications, 13, 1790-1801.

Warton, D.I. \& Shepherd, L.C. (2010) Poisson point process models solve the 'pseudo-absence problem' for presence-only data in ecology. Annals of Applied Statistics, 4, 1383-1402.
Wiegand, T., Gunatilleke, S. \& Gunatilleke, N. (2007a) Species associations in a heterogenous Sri Lankan dipterocarp forest. The American Naturalist, 170, E77-E95.
Wiegand, T., Gunatilleke, S., Gunatilleke, N. \& Okuda, T. (2007b) Analyzing the spatial structure of a Sri Lankan tree species with multiple scales of clustering. Ecology, 88, 30883102.

Wisz, M.S., Pottier, J., Kissling, W.D. et al. (2013) The role of biotic interactions in shaping distributions and realised assemblages of species: implications for species distribution modelling. Biological Reviews, 88, 15-30.
Yackulic, C.B., Chandler, R., Zipkin, E.F., Royle, J.A., Nichols, J.D., Grant, E.H.C. \& Veran, S. (2013) Presence-only modelling using MAXENT: when can we trust the inferences? Methods in Ecology and Evolution, 4, 236-243.
Yoccoz, N.G., Nichols, J.D. \& Boulinier, T. (2001) Monitoring of biological diversity in space and time. Trends in Ecology and Evolution, 16, 446-453.

\section*{SUPPORTING INFORMATION}

Additional supporting information may be found in the online version of this article at the publisher's web-site.

Appendix S1 Fisher information matrix.
Appendix S2 R code (R Core Team, 2014) used to simulate and analyse data in opportunistic and planned surveys.

\section*{BIOSKETCH}

Robert M. Dorazio is a Research Statistician at the US Geological Survey's Southeast Ecological Science Center. He also holds a Courtesy Associate Professorship in the Department of Statistics at the University of Florida. His research is motivated primarily by statistical inference problems that arise in the general areas of population dynamics, community ecology and conservation biology. He is also interested in developing the theory and practice of adaptive decision-making in problems of natural resource management.

Editor: Niklaus Zimmermann