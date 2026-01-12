\title{
\title{
Point process models for presence-only analysis
}
}

\author{
lan W. Renner \({ }^{\mathbf{1}}\) *, Jane Elith \({ }^{\mathbf{2}}\), Adrian Baddeley \({ }^{\mathbf{3}}\), William Fithian \({ }^{\mathbf{4}}\), Trevor Hastie \({ }^{\mathbf{4}}\), Steven J. Phillips \({ }^{\mathbf{5}}\), Gordana Popovic \({ }^{\mathbf{6}}\) and David I. Warton \({ }^{\mathbf{6}}\) \\ \({ }^{1}\) School of Mathematical and Physical Sciences, The University of Newcastle, University Drive, Callaghan, NSW 2308, Australia; \({ }^{2}\) School of BioSciences, The University of Melbourne, Parkville, Vic. 3010, Australia; \({ }^{3}\) Department of Mathematics \& Statistics, Curtin University, GPO Box U1987, Perth, WA 6845, Australia; \({ }^{4}\) Department of Statistics, Stanford University, 390 Serra Mall, Stanford, CA 94303, USA; \({ }^{5} 2201\) 4th Street, Boulder, CO 80304, USA; and \({ }^{6}\) School of Mathematics and Statistics and Evolution \& Ecology Research Centre, The University of New South Wales, Sydney, NSW 2052, Australia
}

\begin{abstract}
Summary
1. Presence-only data are widely used for species distribution modelling, and point process regression models are a flexible tool that has considerable potential for this problem, when data arise as point events.
2. In this paper, we review point process models, some of their advantages and some common methods of fitting them to presence-only data.
3. Advantages include (and are not limited to) clarification of what the response variable is that is modelled; a framework for choosing the number and location of quadrature points (commonly referred to as pseudoabsences or 'background points') objectively; clarity of model assumptions and tools for checking them; models to handle spatial dependence between points when it is present; and ways forward regarding difficult issues such as accounting for sampling bias.
4. Point process models are related to some common approaches to presence-only species distribution modelling, which means that a variety of different software tools can be used to fit these models, including maxent or generalised linear modelling software.
\end{abstract}

Key-words: Cox processes, Gibbs processes, maxent, pseudo-absences, species distribution modelling

\section*{Introduction}

Species distribution modelling (SDM) provides a framework for determining the distribution of a species' habitat as a function of environmental variables and is a highly researched topic of interest to ecologists, biologists and climate change scientists. Often, the best available species data come in the form of a list of reported presence locations of a species without any corresponding information about where a species is absent. This type of data is known as 'presence-only' data (Pearce \& Boyce 2006) and can be found in museums, atlases and herbaria. A researcher interested in exploring the relationship between a species and the environment is faced with the question of which methods to choose, and there have been calls for unification of SDM concepts (Elith \& Leathwick 2009; Aarts, Fieberg \& Matthiopoulos 2012). Here, using language common to the SDM literature, we aim to progress understanding of methods by focussing on emerging knowledge of the links between point process models (PPMs), regression and MAXENT.

Presence-only data typically arise as point events - a set of point locations where a species has been observed. In the statis-

\footnotetext{
*Correspondence author. E-mail: ian.renner@newcastle.edu.au
}
tical literature, a set of point events (in which the location and number of points is random) is known as a point process. Spatial statistics literature provides a suite of tools for modelling point processes (Cressie 1993; Diggle 2003), but only recently have point process models been proposed as a natural way for analysing species presence-only data in a regression framework (Warton \& Shepherd 2010; Chakraborty et al. 2011). PPMs are closely connected to methods already in widespread use in ecology such as MAXENT (Aarts, Fieberg \& Matthiopoulos 2012; Fithian \& Hastie 2013; Renner \& Warton 2013), some implementations of logistic regression (Baddeley et al. 2010; Warton \& Shepherd 2010) and estimation of resource selection functions (Aarts, Fieberg \& Matthiopoulos 2012; McDonald et al. 2013). PPMs enjoy particular benefits in interpretation and implementation. Consequently, we believe PPMs are a natural choice of analysis method for presence-only SDM, when the data arise as point events.

In this paper, we review PPMs for species distribution modellers, and different methods for fitting them. We show that the point process viewpoint resolves a number of important questions regarding presence-only data, including exactly what target quantity is being modelled, how to select background or pseudo-absence points (named quadrature points in the PPM literature), what assumptions are made and how they can be
checked. It also suggests a natural way to deal with biases. We provide a worked example to demonstrate how to fit PPMs with different methods.

\section*{Example - distribution of Eucalyptus sparsifolia}

We will proceed by briefly describing an example data set to be used throughout the paper to illustrate key ideas. The data set comprises 230 presence-only locations of Eucalyptus sparsifolia within the Greater Blue Mountains World Heritage Area (GBMWHA) and a surrounding \(100-\mathrm{km}\) buffer zone (Fig. 1), a \(86227-\mathrm{km}^{2}\) area near Sydney, Australia (NSW Office of Environment and Heritage 2012). This species is known to be abundant and broadly distributed across this region (Hager \& Benson 2010). Maps of environmental variables were available over the study region, and the goals of analysis were to map the distribution of Eucalyptus sparsifolia and identify its key environmental correlates.

The presence records were entirely from incidental sightings of the species collated by the responsible state department since 1972. These records were cleaned prior to analysis to remove data from systematically sampled transects, and records with high location errors, leaving only records of opportunistic sightings whose point location was known to a reasonable \((1 \mathrm{~km})\) degree of accuracy. As these observations are reported as locations, not counts in transects or grid cells, they are best described as point locations in continuous space, which motivates the use of point process models for analysis.

Environmental predictors selected as likely relevant to the distribution of the species and its records are minimum and maximum annual temperature, annual rainfall, number of fires

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/32261777-e61f-4cb5-87f9-d15254cb1f65-02.jpg?height=845&width=606&top_left_y=1630&top_left_x=300}
\captionsetup{labelformat=empty}
\caption{Fig. 1. Locations of 230 presence-only Eucalyptus sparsifolia observations within 100 km of the Greater Blue Mountains World Heritage Area.}
\end{figure}
since 1943 and a categorical soil variable. A description of the soil categories is presented in Section 1 of the Appendix S1. All variables were available at \(100-\mathrm{m}\) grid cell resolution (Renner et al. 2015).

\section*{Point process models}

Presence-only data consist of a set of locations \(\mathbf{s}_{P}=\left\{s_{1}, s_{2}, \ldots, s_{m}\right\}\) at which a species has been observed in some region \(\mathcal{A}\). While methods for fitting PPMs are closely related to common regression models, particularly generalised linear models (GLMs, McCullagh \& Nelder 1989), PPMs are posed differently. A regression model is typically used when the object of interest is a random variable \(y_{i}\), for which we model the mean \(\mu_{i}\) as a function of covariates \(\mathbf{x}_{i}\). By contrast, the objects of primary interest in a PPM are the spatial locations of presence points \(\mathbf{s}_{P}\) - that is the focus is on where the points were observed. We model the locations in \(\mathbf{s}_{P}\) jointly with the number of presences \(m\) and characterise them via the intensity or limiting expected number of presence points per unit area \(\lambda(s)\). The link to regression comes because we typically model \(\lambda(s)\) as a function of covariates \(\mathbf{x}(s)\) measured throughout the study region \(\mathcal{A}\).

The first advantage of PPMs, before looking any further, is greater clarity about what exactly is being modelled (Aarts, Fieberg \& Matthiopoulos 2012; Dorazio 2012). The target of interest, intensity, is not a probability and is instead a measure of abundance - the number of presence records per unit area. Thus, it need not have an upper bound of one (Aarts, Fieberg \& Matthiopoulos 2012). Further, intensity as defined is a function of only two quantities - spatial patterning in the presenceonly data, and the spatial measurement units. Changing the spatial units from kilometres to metres should change intensity proportionally (decreasing by a factor of \(1000^{2}\) ).

It should be emphasised that in most instances, the intensity \(\lambda(s)\) does not reflect the expected abundance per unit area of a species; rather, it reflects the expected abundance of species reportings. It can typically only be used to make inferences about relative patterns in species abundance (Fithian \& Hastie 2013).

\section*{THE POISSON CASE - NO SPATIAL DEPENDENCE}

The simplest type of PPM of use in presence-only analysis is an inhomogeneous Poisson point process (hereafter referred to as a Poisson PPM), in which we assume (a) point events are independent of each other, which can be shown to imply that the total number of points in the study region is a Poisson random variable, and (b) that the intensity \(\lambda(s)\) varies spatially (and so is indexed by location \(s\) ). We will further assume it varies according to environmental conditions \(\mathbf{x}(s)\).

Assumption (a), in which point locations are independent, is a restrictive assumption which often is not satisfied by pres-ence-only data. Methods to check the independence assumption are described in Section 'Software for fitting point process models' and methods to fit models that account for dependence are described in Section 'Spatial dependence in point processes'.

Assumption (b) is often refined to a loglinearity assumption, where we model intensity as a loglinear function of environmental covariates:
\(\ln \lambda(s)=\mathbf{x}(s)^{\prime} \beta\),
(eqn 1)
where \(\beta=\left\{\beta_{1}, \ldots, \beta_{p}\right\}\) is a vector that contains the parameters corresponding to the \(p\) environmental covariates \(\mathbf{x}(s)\). Loglinearity is a natural assumption because it ensures that intensity is a non-negative quantity, and it is the canonical link for Poisson data. While the form of this equation is loglinear, it can readily capture nonlinear relationships between intensity and the environment, for example, using quadratic and interaction terms, smoothed functions in generalised additive models (GAMs, Hastie \& Tibshirani 1990) or via kernel regression (Guan 2008; Baddeley et al. 2012).

Loglinear models of the form of (1) are commonly fitted to count data, and one way to fit a PPM is in fact to break the data into grid cells and fit a Poisson loglinear model to the counts of presence points in each grid cell. However, such a model has the potential to lose information from the data during aggregation from point location to the grid cell level (Renner \& Warton 2013). It may, however, be helpful to readers to think of a Poisson PPM as like a model for count data the key distinction being that the data for analysis are the set of point locations of presences, rather than counts in grid cells.

\section*{REGRESSION MODELS OF PRESENCE-BACKGROUND DATA}

Here we pause to briefly discuss connections with a common practice in the SDM literature, presence-background (PB) regression, also referred to as pseudo-absence regression (e.g. Chefaoui \& Lobo 2008; Phillips et al. 2009; Barbet-Massin et al. 2012) and more recently as 'naïve regression' (Fithian \& Hastie 2013). This models presence ( \(y=1\) ) and background (treated as \(y=0\) ) with logistic regression methods usually used for presence-absence data. This approach can be understood as being 'naïve' essentially because of a mismatch between the model being fitted and the data that were collected - the presences ( \(y=1\) ) are the raw data for which we wish to specify a model, and the background points \((y=0)\) are a fabrication. PB regression was motivated by the need to model distributions of species for which survey (PA) data were unavailable and, early on, by a lack of suitable alternatives. It has remained popular, perhaps because of the examples in which this approach seems to work reasonably well compared with other methods, and the fact that many ecologists are already familiar with regression methods. The approach is somewhat ad hoc. Some users apply arbitrary weights to the background samples, for pragmatic rather than statistically based reasons (Elith et al. 2006). The fitted quantity is interpreted as a relative likelihood of presence, with an unknown scaling linking it to the true probability of presence.

One of the most challenging steps in fitting a PB regression model is selection of the background points. A common choice has been to select a large number (thousands) of points across
the landscape of interest (Elith et al. 2006). Other more complex schemes include identifying points more likely to represent a true absence, or at least avoiding presence points (Engler, Guisan \& Rechsteiner 2004), or trying to specify an optimal number of background points (or presence-background ratio) for different methods (Barbet-Massin et al. 2012). These are usually based on an idea about data structure or on evidence that it performs better than an alternative in particular case studies or simulations, but without stronger statistical justification. These have created some confusion among users regarding which approach is the best to use. We think that efforts to clarify background sampling schemes under the naïve model are misdirected and that much is to be gained by changing to the point process viewpoint. As will be seen later, this viewpoint provides a solid statistical framework for understanding the role of background points and for deciding how many, placed where, are sufficient.

Warton \& Shepherd (2010) and Fithian \& Hastie (2013) discuss problems with using naïve logistic regression and its various extensions for presence-only data. One problem is scale dependence- the scale of PB regression predictions is meaningless since the predictions change as more background points are added to the sample. But Warton \& Shepherd (2010) and Baddeley et al. (2010) showed that PB regression can be understood as an approximation to fitting a Poisson point process model and that the latter resolves the scale dependence issue, and many issues with background choice.

\section*{SPATIAL DEPENDENCE IN POINT PROCESSES}

An underlying assumption of the Poisson PPMs (and most PB regression methods) is that data are conditionally independent given the covariates; that is, the similarities in intensity in nearby regions are fully explained by the environmental and sampling covariates in the model. This is, however, often not the case, and failing to account for spatial dependence can significantly alter conclusions (Dormann 2007). Common examples of spatial dependence are clustering through dispersal or social aggregation. Spatial dependence may also be induced by failing to measure environmental variables which act to make presence patterns in regions close together seem more similar than those further apart.

The two most common classes of point process models for species dependence are Gibbs and Cox processes.

Gibbs or 'interaction' processes relax the independence assumption by assuming interactions between sets of points. One useful example is area-interaction processes (Widom \& Rowlinson 1970; Baddeley \& van Lieshout 1995), which assume interactions among all points within a distance of \(2 r\). These interactions can be understood as introducing an additional term to the intensity function (conditional on the observed presence points):
\(\ln \lambda(s)=\mathbf{x}(s)^{\prime} \beta+t_{s}\left(\mathbf{s}_{P}\right) \theta\)
(eqn 2)
where \(\theta\) is an interaction parameter (positive values implying clustering of points) and \(t_{s}(\mathbf{y})\) is the area of a disc of radius \(r\)
centred at the location \(s\) that does not intersect with similar discs centred around each of the presence points \(\mathbf{s}_{P}\).

Area-interaction models are potentially useful for SDMs because they are capable of modelling either clustering or inhibition. They also have some biological justification; for example, interaction radius \(r\) could be considered as a maximum dispersal distance, with intensity increasing at locations within a distance \(r\) of a known presence because of the chance of establishment from that presence point. Beyond area-interaction processes, there are a number of other types of Gibbs process, in particular processes that involve pairwise interaction between points (for a list, see Baddeley \& Turner 2005, Section 9 of the Appendix S1).

An alternative way to deal with clustering and the effects of unmeasured covariates is by fitting a Cox process, the most common example of which is the spatial log-Gaussian Cox process (LGCP) (Møller, Syversveen \& Waagepetersen 1998). This can be understood as a point process analogue of a generalised linear mixed model with a random intercept that is normally distributed.

The intensity \(\lambda(s)\) in a LGCP is a function not just of environmental variables, but also of a stochastic Gaussian process \(\xi(s)\) :
\(\ln \lambda(s)=\mathbf{x}(s)^{\prime} \beta+\xi(s)\).

Here, \(\xi(s)\) is a spatial Gaussian process with zero mean, and a covariance function that depends on the distance between observations, such that observations closer together in space are assumed to be more positively correlated than those further apart. The \(\xi(s)\) can be understood as an unmeasured covariate which is associated with the distribution of the species. Conditional on this latent process, the point events are assumed to be inhomogeneous Poisson. In other words, it is assumed that any spatial dependence in the data is entirely captured by \(\xi(s)\).

\section*{Fitting a point process model}

This section provides an overview of the process of fitting a PPM, with more detailed information about software and example code provided in the Appendix S1. We use the example Eucalyptus sparsifolia data introduced previously in illustratory analyses.

\section*{FITTING POISSON POINT PROCESSES}

The most common approach to fitting a Poisson PPM is to maximise the log-likelihood function (Cressie 1993), which can be written as:
\(l\left(\beta ; \mathbf{s}_{P}\right)=\sum_{i=1}^{m} \ln \lambda\left(s_{i}\right)-\int_{\mathcal{A}} \lambda(s) d s\).

For a derivation of this log-likelihood, and how it differs from the homogeneous Poisson point process case, see Appendix S1 (Section 2). The integral in this expression can be interpreted as the expected number of presence points in the whole
study region \(\mathcal{A}\), and it is the approximation of this quantity that is the main challenge in model fitting.

Estimation of parameters in the inhomogeneous Poisson PPM is not straightforward, because the integral in (4) does not have a closed form and must be approximated in some way. A standard way to approximate this integral is through the use of numerical integration, otherwise known as 'quadrature' (Davis \& Rabinowitz 1984). The general idea is to choose a set of 'quadrature points' at which the intensity function is evaluated, and these evaluations are then combined as a weighted sum to estimate the integral. Common examples of quadrature methods are Riemann sums and the trapezoidal rule. Irrespective of the quadrature method used, the likelihood can then be written as:
\[
\begin{aligned}
l\left(\beta ; \mathbf{s}_{P}\right) & \approx \sum_{i=1}^{m} \ln \lambda\left(s_{i}\right)-\sum_{j=1}^{m+n} w_{j} \lambda\left(s_{j}\right) \\
& =\sum_{j=1}^{m+n} w_{j}\left(y_{j} \ln \lambda\left(s_{j}\right)-\lambda\left(s_{j}\right)\right)
\end{aligned}
\]
where \(\quad \mathbf{w}=\left\{w_{1}, \ldots, w_{m+n}\right\} \quad\) are quadrature weights, \(\mathbf{s}_{0}=\left\{s_{m+1}, \ldots, s_{m+n}\right\}\) are quadrature points and \(y_{j}=\frac{1}{w_{j}}\) for presence points \((j=1, \ldots, n)\) and \(y_{j}=0\) otherwise. The quadrature weights \(w_{j}\) can be understood as applying a spatial scaling, so that the response being modelled (intensity, \(\boldsymbol{\lambda}\) ) has spatial units not observational units. For example, in analysing the Eucalyptus sparsifolia data, we modelled the expected number of presences per square kilometre. Broadly speaking, \(w_{j}\) represents the area of the neighbourhood around the point \(s_{j}\), found after partitioning the study region \(\mathcal{A}\) into neighbourhoods around each point (including presences and quadrature points).

Equation 6 is due to Berman \& Turner (1992), and it reexpresses the likelihood as a Poisson likelihood with observation weights \(w_{j}\). The significance of this result is that it makes Poisson PPMs relatively straightforward to fit - they can be fitted using any standard glm software, such as the glm function in R (R Development Core Team 2010), using the Poisson family and appropriate observation weightings. This can be done in just a few lines of code, as illustrated in the Appendix S1, although specialised packages are available that offer enhanced options - the spatstat package on R (Baddeley \& Turner 2005) has a suite of model-checking tools, and the ppmlasso package adds a LASSO penalty for improved predictive performance.

The connection to Poisson GLMs in equation 6 has enabled equivalence results between Poisson PPMs and PB regression in large samples (Baddeley et al. 2010; Warton \& Shepherd 2010) and more recently MAXENT (Aarts, Fieberg \& Matthiopoulos 2012; Fithian \& Hastie 2013; Renner \& Warton 2013), offering some insight into these methods. For example, the scale dependence of PB regression can be understood as arising due to the omission of appropriate observation weights \(w_{j}\). The equivalence with MAXENT implies that mAXENT software can be used to fit a Poisson PPM. We will discuss the capabilities of these alternate methods for fitting Poisson PPMs later.

\section*{HOW TO CHOOSE QUADRATURE POINTS}

A first step in fitting a PPM is selection of quadrature points, to allow estimation of the PPM likelihood. This is an equivalent issue to that of background selection (Section 'Regression models of presence-background data'). However, our choice of the term 'quadrature points' in this Section reflects a desire to pose the question of their choice as a quadrature problem, which clarifies their role in analysis and provides a framework for their selection.

From the point process viewpoint, quadrature points are merely a device to estimate an integral. This is true across the various methods that can be used to fit a Poisson PPM (Section Software for fitting point process models, see also Fithian \& Hastie 2013). Being such a device, the important question becomes: How many points, placed where, will be sufficient to accurately estimate the likelihood? This is primarily a question for which the number and location of presence records are irrelevant; hence, ideas of matching number of presence points and sampling far away from presence points are not relevant, except when computational efficiency is an issue, which then requires specifically designed schemes (see later).
Assuming that the extent of the study area is pre-determined by the analyst, two simple strategies to select the location of quadrature points within that extent are (i) to choose them on a regular mesh (e.g. at regularly spaced intervals in each cardinal direction) or (ii) to choose them randomly. Alternative strategies whereby the density of quadrature points increases with environmental variability may be helpful in scenarios where computation time is slow as these would lead to a smaller data set with negligible loss in accuracy. This is related to the ideas of importance sampling and adaptive quadrature (Davis \& Rabinowitz 1984), which seem so far under-utilised in the SDM literature.

The number of quadrature points should be sufficient for an accurate estimate of the likelihood, which will lead to a stable model that is approximately invariant across repeat samples of the points. To determine an appropriate number of quadrature points, we advise that analysts check that sufficient accuracy is
achieved by increasing the number of quadrature points until there is little appreciable change in model fit or in predictive performance (Phillips \& Dudík 2008). For instance, for Fig. 2a, we fitted Poisson PPMs to Eucalyptus sparsifolia data repeatedly halving the spacing between quadrature points, selected across a regular mesh. This was done using the ppmlasso package, which has a function specifically designed to perform this operation. The log-likelihood converged at a \(1-\mathrm{km}\) spacing, which required more than 86000 quadrature points for our example (Fig. 2a), noticeably more than common software defaults (e.g. 10000 in maxent).

Undersampling of quadrature points will lead to error in our coefficient estimates, although this error may still be small compared to the coefficients' overall standard errors. If quadrature points have been randomly sampled (rather than using a regular mesh), we can easily quantify this error by considering what might happen under repeated sampling of quadrature points. For example, in Fig. 2b, models were fitted with an increasing number of randomly chosen quadrature points, and results replicated 30 times to study how much results varied when using different sets of random quadrature points (code available in Section 5 of the Appendix S1). Figure 2b suggests that 100000 or more randomly chosen quadrature points would be needed to reliably estimate the maximised log-likelihood - with significant variation from one set of quadrature points to another before this point. In particular, if using 10000 random quadrature points, the maximised log-likelihood varied over a range of more than 50 for different samples of quadrature points.

Alternatively, the error introduced by quadrature can be estimated analytically quite easily. The integral was estimated as an average of estimated intensities at \(n\) random points (multiplied by \(|\mathcal{A}|\) ), so its uncertainty can be estimated using the formula for the standard error of a sample mean, \(\frac{|\mathcal{A}| \sigma}{\sqrt{n}}\). In our example data set, given an initial fit using 10000 random quadrature points, the standard deviation of estimated intensities at the quadrature points was \(s=0 \cdot 0103\), yielding an estimated standard error of \(\frac{86227 \times 0.0103}{\sqrt{10000}}=8.89\). If we desire an

\begin{figure}
\captionsetup{labelformat=empty}
\caption{Fig. 2. Checking for likelihood convergence as the number of quadrature points changes: (a) Using a rectangular mesh of quadrature points at different spatial resolutions, as available in the ppmlasso package; (b) using random sets of quadrature points and progressively increasing the sample size, as estimated using downweighted Poisson regression models. It appears that there is little benefit in analysing the data at a spatial resolution finer than 1 km (a), or with more than 100000 quadrature points (b).}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/32261777-e61f-4cb5-87f9-d15254cb1f65-05.jpg?height=514&width=1090&top_left_y=1926&top_left_x=495}
\end{figure}
estimate of the log-likelihood to be within a standard error of two of its true value, we can estimate the required number of quadrature points as \(|\mathcal{A}|^{2} s^{2} / 2^{2}=197745\). This corresponds well with the results of Fig. 2b.

The precise number of points needed for a sufficiently accurate estimate should vary with the roughness of the intensity surface (hence the difficulty of the integration problem). Thus, we may expect to need more quadrature points when environmental data are measured at a finer resolution (equivalently, smaller grid cell sizes) or broader spatial extent.
'Quadrature thinking' also leads to other emphases. First, if computational efficiency is of key importance, it could be prudent to sample fewer quadrature points (with higher corresponding quadrature weights) in areas where the species is unlikely to occur, which should have negligible impact on the intensity and the likelihood of the fitted model. However, this needs to be done with care because it relies on correct identification and inclusion in the model of the covariates causing species absence from certain parts of the landscape. An example of where this may be appropriate is in telemetry studies, where an individual's location is typically strongly associated with distance from the last observed location; hence, quadrature points far from that location make negligible likelihood contribution (Warton \& Aarts 2013).

Secondly, the quadrature viewpoint tends to place more emphasis on specifying the model in a way that deals with biases, rather than fiddling with quadrature points. Hence with quadrature thinking, one is more likely to accommodate bias by explicitly specifying covariates in the model (e.g. Warton, Renner \& Ramp 2013; Fithian et al. 2015) than through the selection of quadrature points (e.g. target-group background as in Phillips et al. 2009). For example, in our Eucalyptus sparsifolia model, we added two predictors related to site accessibility in order to model observer bias (distance from main roads and distance from urban areas). These two observer bias variables were included to try to account for spatial patterning in presence locations due to behaviour of observers rather than behaviour of the study species. By then making predictions at a common level of observer bias (e.g. distance equals zero), we can map E. sparsifolia distribution controlling for observer bias (Warton, Renner \& Ramp 2013; Fithian et al. 2015).

Finally, most methods for fitting Poisson PPMs attach weights to the quadrature points (the \(w_{j}\) in equation 5) for a scale-invariant estimate of the log-likelihood that is comparable across sets of quadrature samples of different size. This can have advantages for model fitting and interpretation (Warton \& Shepherd 2010).

\section*{CHECKING ASSUMPTIONS}

A suite of diagnostic tools are available in the point process literature that can be used to ground-truth model assumptions. Just as with ordinary regression models, residual analysis (Baddeley et al. 2005) can be used to assess adequacy of the model for intensity, in particular by checking for a spatial trend in residuals. The assumption of independence among point locations can be checked using Ripley's
\(K\)-function (Ripley 1977) and its generalisations (Baddeley, Møller \& Waagepetersen 2000).

Consider, for example, a Poisson PPM fitted to the E. sparsifolia data. A check of the independence assumption suggests significant clustering of points at radii \(<10 \mathrm{~km}\) (Fig. 3a), so we fit an area-interaction model (see Section 3 of the Appendix S1 for details). The resulting model fit exhibits some pattern (Fig. 3b), but the magnitude of the residuals is not exceedingly large. Cumulative residuals for increasing longitude and latitude ( \(x\) and \(y\) ) do not significantly deviate from Monte Carlo simulation envelopes (Fig. 3c-d), so the model fit may be deemed sufficiently appropriate.

Other useful diagnostic features demonstrated in Section 3 of the Appendix S1 include influence, leverage and partial residual plots, all derived in direct analogy to how they are used in generalised linear models (as in Baddeley et al. 2013). All of these diagnostic plots are relatively easy to produce using the spatstat package.

When checking the assumption of independence of presence points, an alternative approach is to treat models that incorporate spatial dependence as the default, to fit such models, and then study the level of dependence in the subsequent model fit, as in Fig. 4a. This approach makes particular sense as an approach if one is expecting spatial dependence a priori, and data of sufficient quality to see the spatial dependence signal. The fitting of such models is discussed below.

\section*{FITTING POINT PROCESSES WITH SPATIAL DEPENDENCE}

Both Gibbs and Cox processes are more difficult to fit than Poisson PPMs, although for different reasons.

Gibbs processes are difficult to fit by maximum likelihood because their specification involves a proportionality constant that is difficult to estimate. One workaround is to use the Poisson process likelihood in place of the true likelihood (Besag 1977), that is, use (5) as a pseudo-likelihood. This is often done because the Gibbs likelihood has a complex form, but by using the Poisson pseudo-likelihood instead, estimates are readily available via GLM. This approach is implemented in spatstat and ppmlasso (see Sections 3 and 4 of the Appendix S1 for detailed code), with a number of types of interaction processes to choose from in spatstat, as specified using the interaction argument. One issue to be aware of when using this approach is that traditional methods of likelihood-based inference (such as likelihood-based standard errors, likelihood ratio tests, AIC) may no longer apply, because parameters have not been estimated by maximum likelihood.

Cox process models are difficult to fit because they involve an unobserved Gaussian random process (the \(\xi(s)\) in equation 3), and so maximum likelihood estimation would involve complex integrals. These models are hierarchical and are therefore naturally suited to a Bayesian hierarchical approach to estimation as in Illian, Møller \& Waagepetersen (2009) and Chakraborty et al. (2011). Other estimation techniques include composite likelihood (Guan 2006) and weighted esti-

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/32261777-e61f-4cb5-87f9-d15254cb1f65-07.jpg?height=1790&width=1663&top_left_y=247&top_left_x=203}
\captionsetup{labelformat=empty}
\caption{Fig. 4. Estimation of spatial dependence using a log-Gaussian Cox process model: (a) Mean and (b) standard deviation of the random Gaussian field, and (c) a posterior mean and \(95 \%\) credible interval for the interaction coefficient, computed using the dist function on r-inla. The mean is significantly larger than the standard deviation in some regions, and the interaction coefficient in (c) has a prediction interval greater than zero at small distances, both of which suggest clustering in the data beyond that explained by covariates.}
\end{figure}
mating equations (Guan \& Shen 2010), as implemented in spatstat.
Two methods of Cox process estimation that are currently popular, and which both have been implemented under a Bayesian framework, are the integrated nested Laplace approximation (INLA, implemented in the r-inla package, Rue, Martino \& Chopin 2009) and a Markov chain Monte Carlo (MCMC) method (implemented in the lGCp package,

Taylor et al. 2013). The INLA method seeks to calculate the integrals by a set of carefully chosen approximations. It is generally fast compared to MCMC methods, which can be quite time-consuming but have potentially greater accuracy. For a comparison of these, see Taylor \& Diggle (2014).

Sections 7 and 8 of the Appendix S1 have example code for fitting a Cox process using both the r-inla and lgCP packages. We found these packages much more difficult for
the practitioner to use than Poisson PPMs, with specialist guidance often required. The main gain from this additional effort is the possibility of making valid inferences from the model, taking into account uncertainty, under the assumption of spatial dependence. This is otherwise harder to achieve without resorting to resampling, as below.

\section*{BLOCK RESAMPLING FOR INFERENCE IN THE PRESENCE OF SPATIAL DEPENDENCE}

A complication that can arise when modelling spatially dependent point processes is that likelihood-based standard errors may be estimated to be too small. This problem arises when using the pseudo-likelihood approach to fit Gibbs processes or when failing to account for interpoint dependence correctly. This issue does not arise in fitting Cox processes, if correctly specified, a significant advantage of that approach.

Similarly, but irrespective of what model is fitted to data, standard cross-validation measures of out-of-sample prediction error can be over-optimistic in the presence of dependence, potentially leading us to erroneously prefer models that overfit to local structure in the data (Wenger \& Olden 2012). Block resampling techniques can deal with short-range interpoint dependence nonparametrically.

Suppose that interpoint dependence is strong for nearby locations, but weak at radius \(r\). Then if we tile our geographic domain into \(c\) rectangular blocks of size \(r \times r\), the data falling in one block are approximately independent of the data in other blocks. If we believe this, then we can obtain accurate standard errors using a bootstrap algorithm that resamples whole blocks with replacement (Efron \& Tibshirani 1993), as in Slavich et al. (2014). Likewise, we can obtain estimates of out-of-sample prediction error by cross-validation where blocks are assigned whole to each fold (Burman, Chow \& Nolan 1994), as in Wenger \& Olden (2012), Pearson et al. (2013) and Warton, Renner \& Ramp (2013).

For example, the model for Fig. 5c involved a LASSO penalty which was estimated by 5 -fold cross-validation using blocks of \(32 \mathrm{~km} \times 32 \mathrm{~km}\). Section 4 of the Appendix S1 has example code for block cross-validation with the ppmlasso package.

\section*{SOFTWARE FOR FITTING POINT PROCESS MODELS}

The main software packages currently available for fitting point process models, and their key differences in properties, are summarised in Table 1. All are available in r . In the Appendix S1, we have developed short tutorials stepping the user through analysis of the Eucalyptus sparsifolia data using each of these packages, and we encourage new users to work through these resources when deciding on an approach to analysis of point process data.

The most established package for point process modelling is spatstat, whose main advantages are the extensive suite of diagnostic tools, and its ability to simulate data from a given point process model. The ppmlasso package was written to be spatstat-compatible, so it inherits many useful diagnostic
tools, while adding a couple of functions of particular interest for SDMs - the ability to regularise parameter estimates (i.e. shrink them towards zero to reduce variance, Hastie, Tibshirani \& Friedman 2009) using the LASSO or elastic net, and functions to guide the user regarding quadrature point choice, as in Fig. 2a. Standard errors are not returned in ppmlasso output, since these become very approximate when using a LASSO or other regularisation approach in parameter estimation.

Because point process models can be fitted using standard glm software, one can entirely avoid using specialised point process software, but the onus is on the user to ensure that quadrature points have been chosen in sufficient numbers and locations and that the appropriate quadrature weights have been assigned. The main advantage of this approach is that the user has greater control over what type of model is fitted - for example, as well as GLM, one could use any of its extensions (GAM, elastic net, CART, . . ). Along these lines, Fithian \& Hastie (2013) proposed a simple algorithm based on weighted logistic regression that can be used to estimate slope coefficients in a Poisson PPM, referred to as infinitely weighted logistic regression (IWLR). But IWLR only gives a solution proportional to a point process, because the intercept term and hence the scale of the log-likelihood is arbitrary. A modification of this idea, which we call 'downweighted Poisson regression' (DWPR), is proposed in Section 5 of the Appendix S1. Given a random set of pseudo-absences, this reduces the steps of assigning quadrature weights and fitting the model to just a few lines. When using DWPR, the intercept and hence the likelihood are estimable, so we can use this technique to look at questions like how many quadrature points to select (as done in Fig. 2b).

Given that MAXENT has recently been shown to be proportional to a Poisson PPM (Aarts, Fieberg \& Matthiopoulos 2012; Fithian \& Hastie 2013; Renner \& Warton 2013), maxent software can also be thought of as a means of fitting a Poisson PPM. This does, however, require some departures from the default maxent software settings, as described in Section 6 of the Appendix S1. Using maxent to fit a Poisson PPM has the advantages that it is familiar to many users, has a lot of nice mapping features and is easy to use. Example features include the interactive 'explain' tool which visualises the link between the mapped prediction and the model at any selected point, and maps that show whether environments in new regions or times of interest are within the training range of the modelled data. maxent can be run from r using the package dismo, which allows streamlined data preparation and modelling, and the potential to use block resampling along the lines of Wenger \& Olden (2012). A short tutorial for fitting PPMs using dismo is presented in Section 6 of the Appendix S1.

A key issue, however, with implementations of Poisson PPM using GLM and maxent software is the lack of assumption checking tools. One should not take 'on faith' the assumption that there is no spatial dependence in the data beyond that explained by environmental variables included in the model. A related issue is the lack of a capacity to fit models that account for point interactions using GLM or maxent.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/32261777-e61f-4cb5-87f9-d15254cb1f65-09.jpg?height=1722&width=1100&top_left_y=241&top_left_x=203}
\captionsetup{labelformat=empty}
\caption{Fig. 5. Predicted intensity of Eucalyptus sparsifolia observations for a Poisson point process model fitted using (a) spatstat or (b) maxent; (c) an area-interaction model fitted using ppmlasso; and (d) a log-Gaussian Cox process fitted using r-inla. Note that spatstat and maxent results look quite similar, as they fit similar models, and that the ppmlasso and rinla results, which account for spatial dependence, highlight additional areas of relatively high intensity further west.}
\end{figure}

Fitting Cox processes in a Bayesian framework, as in R-INLA and lGCP, has the advantage that spatial dependence can be estimated and accounted for in a flexible and statistically efficient way. However, the additional complexity of estimating the latent field results in much slower computation times as compared to competitors. Careful selection of the respective prior distribution for the Gaussian random field is important to avoid overfitting (Illian et al. 2013). But there is currently little guidance about prior selection - such guidelines are currently being developed for R -INLA.

\section*{Analysis of Eucalyptus sparsifolia}

When fitting point process models to Eucalyptus sparsifolia data, results were broadly similar across methods of fitting

Poisson PPMs (Fig. 5a-b), with some variation due to different decisions being made in default implementations (e.g. spatstat does not apply a LASSO penalty, whereas maxent does). But as identified previously, there is evidence of spatial clustering between presence points in close proximity - at distances less than 10 km , the inhomogeneous \(K\)-function in Fig. 3a strays above its simulation envelope, and the point interaction coefficient from a fitted Cox process is significantly above zero (Fig. 4c). Further evidence of dependence can be seen in the mean of latent Gaussian field from the Cox process fit, which was sometimes large relative to its standard deviation (Fig. 4a-b). Maps produced by models which account for this spatial dependence (Fig. 5c and d) highlight areas of relatively higher intensity further west.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1. Summary table of software properties}
\begin{tabular}{|l|l|l|l|l|l|l|l|}
\hline Property & SPATSTAT & PPMLASSO & IWLR & DWPR & MAXENT & R-INLA & LGCP \\
\hline Regularisation & × & v & ✓ & ✓ & \(\boldsymbol{v}^{1}\) & × & × \\
\hline Standard errors & \(\boldsymbol{v}^{2}\) & × & \(v^{2}\) & \(v^{2}\) & × & v & \(\nu\) \\
\hline Variable importance plots & × & × & × & × & v & × & × \\
\hline Diagnostic plots & \(\nu\) & \(\nu\) & × & × & × & × & × \\
\hline Spatial dependence & v & v & × & × & × & v & v \\
\hline Nonlinearity (e.g. smoothers) & v & v & v & v & v & v & \(\nu\) \\
\hline Scale invariant & \(\nu\) & \(\nu\) & × & \(\nu\) & \(\boldsymbol{v}^{\mathbf{3}}\) & ✓ & ✓ \\
\hline
\end{tabular}
\end{table}
\({ }^{1}\) LASSO only.
\({ }^{2}\) For Poisson models only.
\({ }^{3}\) Raw output only.

The analysis is largely consistent with pre-existing knowledge of the distribution of Eucalyptus sparsifolia, thought to prefer 'low nutrient soils, but some on medium and high nutrient soils, over a wide range of rainfall' (Hager \& Benson 2010). The maxent response curve for rainfall (Fig. 6a) indicates suitability over a wide rainfall range, and several models (Poisson PPM and area-interaction models produced by ppmlasso, and maxent) either dropped all rainfall terms or assessed them as relatively unimportant in describing the distribution of E. sparsifolia (results not shown). The species appears most strongly associated with low-nutrient high-quartz sedimentary soils and has low intensity in volcanic soils. See Section 9 of the Appendix S1 for a list of coefficient estimates.

Minimum annual temperature is strongly associated with E. sparsofila distribution, yet not mentioned by Hager \& Benson (2010). The quadratic term is significantly negative in models produced by ppmlasso and r-inla, consistent with the response curve produced by maxent (Fig. 6a). This variable has implications for climate change projections, suggesting a substantial decrease in E. sparsofila intensity at the southern end of its range under warming scenarios (Fig. 6a).

A key distinction in the models produced by the different methods is in the number of significant variables (Section 9 of the Appendix S1). The Poisson PPM and area-interaction model fitted by ppmlasso added 24 and 18 nonzero terms in the model, respectively, while the Cox process model produced by r-inla added only seven, five of which were soil indicators. One possible explanation is that there may have been collinearity between the Gaussian random field and environmental predictors, dampening the environmental signal - such a 'spatial confounding' effect is seen elsewhere in spatial statistics, and adjustments can account for it when fitting spatial generalised linear mixed models (Hodges \& Reich 2010; Hughes \& Haran 2013). Extension of these ideas to address spatial confounding in Cox process models is a potential avenue for future research.

\section*{Extensions}

To this point, the focus has been on spatial point processes, to describe a set of point locations in space. A number of poten-
tial variations on the method may be of interest to species distribution modellers.

A time stamp is often available with presences as well as their point location. It may be of interest to study the patterning of points jointly in space and time, thus fitting a spatio-temporal point process (Cressie \& Wikle 2011). This seems especially relevant in telemetry (Hooten et al. 2013), where one would expect strong temporal dependence in spatial patterning (with individuals tending to be found near their last known location); thus, there is a strong case for modelling such point events jointly in space and time, in order to tease apart habitat preference from habitat availability. Another example where spatiotemporal modelling is of clear interest is in the study of invasive species not at equilibrium (Hooten \& Wikle 2008).

Sometimes presences are observed along a network rather than in a region in space. For example, when modelling roadkill events (Ramp et al. 2005), presence points occur along a road network. Poisson PPMs can be fitted to point events arising along networks relatively easily; however, methods for studying and accounting for dependence in point events along networks are a little explored topic (Baddeley, Jammalamadaka \& Nair 2014).

Simultaneously modelling data from multiple species has potential from a number of standpoints. Fithian et al. (2015) showed how estimating observer bias simultaneously across species can improve outcomes, making use of the idea that the sources of observer bias can often be reasonably assumed constant across species, given that this bias is a property of the observer more so than of the study species. Multispecies models could also be used to study species interaction, by specifying what is known as a marked point process model (Cressie 1993) which explicitly includes terms for species interaction.

An issue with presence-only data is data quality, and there are a number of potential extensions to take into account data of varying quality. For example, typically there is uncertainty in the spatial location of presence points, and locations are often assigned 'accuracy' scores to estimate this which can be accounted for in subsequent modelling (Hefley et al. 2014). Further, environmental data

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/32261777-e61f-4cb5-87f9-d15254cb1f65-11.jpg?height=1522&width=1094&top_left_y=245&top_left_x=206}
\captionsetup{labelformat=empty}
\caption{Fig. 6. Maps produced by maxent software to aid interpretation: (a) The explain tool - the location indicated by the arrow has low intensity, and the plots at the right suggest that the high minimum temperature is largely responsible. (b) Multivariate similarity surface plot (right) and most dissimilar variable plot (left) for climate change predictions - the redder colours in the left panel indicate areas that are most dissimilar to the environmental conditions used to build the model, while the right panel identifies which variables are most responsible for the dissimilarity. It seems minimum temperature may be limiting the distribution near the coast under the climate change model.}
\end{figure}
are also measured with uncertainty, maps of environmental data often being spatially interpolated from (often sparsely distributed) weather stations. This is a type of errors-invariables problem (Carroll et al. 2012), and PPMs for such problems, while rare in the literature, should be a relatively straightforward extension given established methods for errors-in-variables approaches to GLM (Stoklosa et al. 2015). Presence-only data are furthermore subject to imperfect detectability, inducing biased estimates, but this bias can be reduced by building a hierarchical model including both presence-only data and independent pres-ence-absence data (Dorazio 2014).
A key strength of point process models is that they operate at what is usually the most ecologically relevant of sampling levels - the level of the individual. This means they can (in principle) incorporate processes operating at the level of the individual, such as interactions between individual organisms, or covariates that vary across individuals. A limiting factor obvi-
ously is the quality of data at hand, but there is exciting potential in this framework.

\section*{Discussion}

Point process modelling has been introduced as a natural framework for modelling presence-only data that is better understood as point events rather than as data from transects or grid cells.

A particular benefit of the PPM specification that we have highlighted is greater clarity around the issue of how to choose quadrature points, and the possibility of querying the data being analysed to verify that a given choice of quadrature points is appropriate (Fig. 2). We found in Section 'How to choose quadrature points' that for our example data, the number of quadrature points required for sufficient convergence of the log-likelihood (a change of \(<2\) ) depended upon the selection method, but was closer to 100000 than to the 10000 usually
advocated. Perhaps more quadrature points are needed when randomly chosen than when on a regular mesh, suggesting greater efficiency when using a regular mesh design, although at the cost of making it more difficult to quantify uncertainty in the approximation. This is related to classical results from survey sampling, where systematic sampling is often more efficient than random sampling under serial correlation, but for which it is more difficult to quantify uncertainty (Cochran 1946).

Our hope is that approaching the 'pseudo-absence problem' via numerical quadrature will shift attention of analysis away from quadrature point choice and towards where it belongs developing and interpreting a plausible model for intensity as a function of environment and possibly observer bias variables.

PPMs as in this paper can be understood as applying regression methods to point event data, and as such, issues that arise in other areas of regression analysis apply equally well here. For example, there is increasing awareness of a dichotomy between prediction and explanation (Elith \& Leathwick 2009) - many SDM researchers are interested primarily in the prediction problem and hence are mostly concerned with using a method which has good predictive performance. This leads the user down the road towards regularised methods (such as in ppmlasso) or model averaging (Araújo \& New 2007). Others are interested primarily in explanation - identifying key associations between environmental variables and a species. This leads the user to put greater focus on appropriately accounting for different sources of uncertainty, which requires careful consideration of the question of spatial dependence in the data, and potential problems like multicollinearity (Zuur, Ieno \& Elphick 2010).

Care must be taken when interpreting a fitted PPM concerning what is actually being modelled, with particular reference to how the data were collected. For example, in order to interpret intensity as relative abundance of individuals per unit area, the intensity of presence records should be proportional to the intensity of individuals (i.e. abundance) of the species. Contrast this with how observers might collect data and how data are entered into online data bases. First, observers may tend to go to a site and only record one individual even though many are present, so individuals in abundant sites are underrepresented. In that case, really one is modelling relative intensity of occupied sites or, equivalently, relative probability of presence (which opens a can of worms, since the extent of a site may be unknown or vary between observers). Secondly, in data aggregation services such as the Global Biodiversity Information Facility (GBIF, Belbin et al. 2013), records arrive from multiple providers. Duplicate records may represent the same individual, for instance, because the same record has been contributed through multiple channels or because different specimens of the same individual were lodged in different museums or herbaria. These duplicate records do not always have exactly the same coordinates because of variable data handling practices. These two examples illustrate typical problems in dealing with the realities of presence-only data, emphasising that it is important that appropriate data cleaning and model interrogation is considered, and that the interpreta-
tion of the final model not stray far from the data used to construct it.

While this paper has focussed on the merits from a modelfitting perspective of taking a point process approach, there are broader advantages afforded by a coherent modelling framework for presence-only data. From a technical perspective, PPMs can serve as a model for generation of simulated presence-only data, with and without point interactions (the spatstat package makes this quite straightforward). From an analyst's perspective, PPMs offer a way forward regarding the assessment of goodness-of-fit to presence-only data, obviating the need to adapt tools like ROC curves to the presence-only context (Jiménez-Valverde 2012). Apart from a suite of diagnostic plots for goodness-of-fit, likelihoodbased procedures can be used to quantify predictive success, for example Kullback-Leibler distance. From an ecologist's perspective, rather than aggregating to an arbitrary sampling unit, one can specify a model for location and behaviour of individual organisms, for example, incorporate covariate information particular to individuals into analyses, where available. In this respect, the point process framework offers an exciting platform that can be used to study ecological processes.

\section*{Data accessibility}

Locations of Eucalyptus sparsifolia, including the 230 presence-only locations used in this paper, are available from http://www.bionet.nsw.gov.au. Environmental data for Blue Mountains region: DRYAD entry doi:10.5061/ dryad.985s5.

\section*{References}

Aarts, G., Fieberg, J. \& Matthiopoulos, J. (2012) Comparative interpretation of count, presence-absence and point methods for species distribution models. Methods in Ecology and Evolution, 3, 177-187.
Araújo, M.B. \& New, M. (2007) Ensemble forecasting of species distributions. Trends in Ecology \& Evolution, 22, 42-47.
Baddeley, A.J. \& van Lieshout, M.N.M. (1995) Area-interaction point processes. Annals of the Institute of Statistical Mathematics, 47, 601-619.
Baddeley, A. \& Turner, R. (2005) Spatstat: an R package for analyzing spatial point patterns. Journal of Statistical Software, 12, 1-42.
Baddeley, A.J., Møller, J. \& Waagepetersen, R. (2000) Non- and semiparametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54, 329-350.
Baddeley, A.J., Turner, R., Møller, J. \& Hazelton, M. (2005) Residual analysis for spatial point processes. Journal of the Royal Statistical Society, Series B, 67, 617-666.
Baddeley, A., Berman, M., Fisher, N.I., Hardegen, A., Milne, R.K., Schuhmacher, D., Shah, R. \& Turner, R. (2010) Spatial logistic regression and change-ofsupport in Poisson point processes. Electronic Journal of Statistics, 4, 11511201.

Baddeley, A.J., Chang, Y.M., Song, Y. \& Turner, R. (2012) Nonparametric estimation of the dependence of a spatial point process on spatial covariates. Statistics and Its Interface, 5, 221-236.
Baddeley, A., Chang, Y.-M., Song, Y. \& Turner, R. (2013) Residual diagnostics for covariate effects in spatial point process models. Journal of Computational and Graphical Statistics, 22, 886-905.
Baddeley, A., Jammalamadaka, A. \& Nair, G. (2014) Multitype point process analysis of spines on the dendrite network of a neuron. Journal of the Royal Statistical Society: Series C (Applied Statistics), 63, 673694.

Barbet-Massin, M., Jiguet, F., Albert, C.H. \& Thuiller, W. (2012) Selecting pseudo-absences for species distribution models: how, where and how many? Methods in Ecology and Evolution, 3, 327-338.

Belbin, L., Daly, J., Hirsch, T., Hobern, D. \& La Salle, J. (2013) A specialists audit of aggregated occurrence records: an aggregators perspective. ZooKeys, 305, 67-76.
Berman, M. \& Turner, T.R. (1992) Approximating point process likelihoods with GLIM. Journal of the Royal Statistics Society, Series C, 41, 31-38.
Besag, J. (1977) Some methods of statistical analysis for spatial data. Bulletin of the International Statistical Institute, 47, 77-91.
Burman, P., Chow, E. \& Nolan, D. (1994) A cross-validatory method for dependent data. Biometrika, 81, 351-358.
Carroll, R.J., Ruppert, D., Stefanski, L.A. \& Crainiceanu, C.M. (2012) Measurement Error in Nonlinear Models: A Modern Perspective. CRC press, Boca Raton, FL.
Chakraborty, A., Gelfand, A.E., Wilson, A.M., Latimer, A.M. \& Silander, J.A. (2011) Point pattern modelling for degraded presence-only data over large regions. Journal of the Royal Statistical Society, Series C, 60, 757-776.
Chefaoui, R.M. \& Lobo, J.M. (2008) Assessing the effects of pseudo-absences on predictive distribution model performance. Ecological Modelling, 210, 478-486.
Cochran, W.G. (1946) Relative accuracy of systematic and stratified random samples for a certain class of populations. The Annals of Mathematical Statistics, 17, 164-177.
Cressie, N.A.C. (1993) Statistics for Spatial Data. John Wiley \& Sons, New York, NY.
Cressie, N. \& Wikle, C.K. (2011) Statistics for Spatio-temporal Data. John Wiley \& Sons, Hoboken.
Davis, P.J. \& Rabinowitz, P. (1984) Methods of Numerical Integration. Academic Press, Orlando, FL.
Diggle, P. (2003) Statistical Analysis of Spatial Point Patterns. Oxford University Press, New York, NY.
Dorazio, R.M. (2012) Predicting the geographic distribution of a species from presence-only data subject to detection errors. Biometrics, 68, 1303-1312.
Dorazio, R.M. (2014) Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. Global Ecology and Biogeography, 23, 1472-1484.
Dormann, C.F. (2007) Effects of incorporating spatial autocorrelation into the analysis of species distribution data. Global Ecology and Biogeography, 16, 129-138.
Efron, B. \& Tibshirani, R. (1993) An Introduction to the Bootstrap. CRC press, Boca Raton, FL.
Elith, J. \& Leathwick, J.R. (2009) Species distribution models: ecological explanation and prediction across space and time. Annual Review of Ecology, Evolution, and Systematics, 40, 677-697.
Elith, J., Graham, C.H., Anderson, R.P., Dudík, M., Ferrier, S., Guisan, A. et al. (2006) Novel methods improve prediction of species' distributions from occurrence data. Ecography, 29, 129-151.
Engler, R., Guisan, A. \& Rechsteiner, L. (2004) An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. Journal of Applied Ecology, 41, 263-274.
Fithian, W. \& Hastie, T. (2013) Finite-sample equivalence in statistical models for presence-only data. The Annals of Applied Statistics, 7, 1917-1939.
Fithian, W., Elith, J., Hastie, T. \& Keith, D.A. (2015) Bias correction in species distribution models: pooling survey and collection data for multiple species. Methods in Ecology and Evolution, 6, 424-438.
Guan, Y. (2006) A composite likelihood approach in fitting spatial point process models. Journal of the American Statistical Association, 101, 1502-1512.
Guan, Y. (2008) On consistent nonparametric intensity estimation for inhomogeneous spatial point processes. Journal of the American Statistical Association, 103, 1238-1247.
Guan, Y. \& Shen, Y. (2010) A weighted estimating equation approach for inhomogeneous spatial point processes. Biometrika, 97, 867-880.
Hager, T. \& Benson, D. (2010) The Eucalypts of the Greater Blue Mountains World Heritage Area: distribution, classification and habitats of the species of Eucalyptus, Angophora and Corymbia (family Myrtaceae) recorded in its eight conservation reserves. Cunninghamia, 10, 425-444.
Hastie, T. \& Tibshirani, R. (1990) Generalized Additive Models. Chapman \& Hall, Boca Raton.
Hastie, H., Tibshirani, R. \& Friedman, J. (2009) The Elements of Statistical Learning: Data Mining, Inference, and Prediction. Springer, New York, NY.
Hefley, T.J., Baasch, D.M., Tyre, A.J. \& Blankenship, E.E. (2014) Correction of location errors for presence-only species distribution models. Methods in Ecology and Evolution, 5, 207-214.
Hodges, J.S. \& Reich, B.J. (2010) Adding spatially-correlated errors can mess up the fixed effect you love. The American Statistician, 64, 325-334.

Hooten, M.B. \& Wikle, C.K. (2008) A hierarchical Bayesian non-linear spatio-temporal model for the spread of invasive species with application to the Eurasian Collared-Dove. Environmental and Ecological Statistics, 15, 59-70.
Hooten, M.B., Hanks, E.M., Johnson, D.S. \& Alldredge, M.W. (2013) Reconciling resource utilization and resource selection functions. Journal of Animal Ecology, 82, 1146-1154.
Hughes, J. \& Haran, M. (2013) Dimension reduction and alleviation of confounding for spatial generalized linear mixed models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75, 139-159.
Illian, J.B., Møller, J. \& Waagepetersen, R.P. (2009) Hierarchical spatial point process analysis for a plant community with high biodiversity. Environmental and Ecological Statistics, 16, 389-405.
Illian, J.B., Martino, S., Sørbye, S.H., Gallego-Fernández, J.B., Zunzunegui, M., Esquivias, M.P. \& Travis, J.M. (2013) Fitting complex ecological point process models with integrated nested Laplace approximation. Methods in Ecology and Evolution, 4, 305-315.
Jiménez-Valverde, A. (2012) Insights into the area under the receiver operating characteristic curve (AUC) as a discrimination measure in species distribution modelling. Global Ecology and Biogeography, 21, 498-507.
Møller, J., Syversveen, A.R. \& Waagepetersen, R.P. (1998) Log Gaussian Cox processes. Scandinavian Journal of Statistics, 25, 451-482.
McCullagh, P. \& Nelder, J. (1989) Generalized Linear Models. Chapman and Hall, London.
McDonald, L., Manly, B., Huettmann, F. \& Thogmartin, W. (2013) Locationonly and use-availability data: analysis methods converge. Journal of Animal Ecology, 82, 1120-1124.
NSW Office of Environment and Heritage (2012) Atlas of NSW Wildlife database. Data accessed 31/05/2012.
Pearce, J.L. \& Boyce, M.S. (2006) Modelling distribution and abundance with presence-only data. Journal of Applied Ecology, 43, 405-412.
Pearson, R.G., Phillips, S.J., Loranty, M.M., Beck, P.S., Damoulas, T., Knight, S.J. \& Goetz, S.J. (2013) Shifts in Arctic vegetation and associated feedbacks under climate change. Nature Climate Change, 3, 673-677.
Phillips, S.J. \& Dudík, M. (2008) Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography, 31, 161-175.
Phillips, S.J., Dudík, M., Elith, J., Graham, C.H., Lehmann, A., Leathwick, J. \& Ferrier, S. (2009) Sample selection bias and presence-only distribution models: implications for background and pseudo-absence data. Ecological Applications, 19, 181-197.
R Development Core Team (2010) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, AustriaISBN 3-900051-07-0
Ramp, D., Caldwell, J., Edwards, K.A., Warton, D. \& Croft, D.B. (2005) Modelling of wildlife fatality hotspots along the Snowy Mountain Highway in New South Wales, Australia. Biological Conservation, 126, 474-490.
Renner, I.W. \& Warton, D.I. (2013) Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics, 69, 274-281.
Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T., Phillips, S., Popovic, G. \& Warton, D. (2015) Data from: Point process models for presence-only analysis - a review. Dryad Digital Repository, doi: 10.5061/dryad.985s5.
Ripley, B.D. (1977) Modelling spatial patterns (with discussion). Journal of the Royal Statistical Society, Series B, 39, 172-212.
Rue, H., Martino, S. \& Chopin, N. (2009) Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 71, 319-392.
Slavich, E., Warton, D.I., Ashcroft, M.B., Gollan, J.R. \& Ramp, D. (2014) Topoclimate versus macroclimate: how does climate mapping methodology affect species distribution models and climate change projections? Diversity and Distributions, 20, 952-963.
Stoklosa, J., Daly, C., Foster, S.D., Ashcroft, M.B. \& Warton, D.I. (2015) A climate of uncertainty: accounting for error in climate variables for species distribution models. Methods in Ecology and Evolution, 6, 412423.

Taylor, B.M. \& Diggle, P.J. (2014) INLA or MCMC? A tutorial and comparative evaluation for spatial prediction in log-Gaussian Cox processes. Journal of Statistical Computation and Simulation, 84, 2266-2284.
Taylor, B.M., Davies, T.M., Rowlingson, B.S. \& Diggle, P.J. (2013) lgcp: an R package for inference with spatial and spatio-temporal log-Gaussian Cox processes. Journal of Statistical Software, 52, 1-40.
Warton, D.I. \& Aarts, G. (2013) Advancing our thinking in presence-only and used available analysis. Journal of Animal Ecology, 82, 1125-1134.

Warton, D.I. \& Shepherd, L.C. (2010) Poisson point process models solve the "pseudo-absence problem" for presence-only data in ecology. Annals of Applied Statistics, 4, 1383-1402.
Warton, D.I., Renner, I.W. \& Ramp, D. (2013) Model-based control of observer bias for the analysis of presence-only data in ecology. PloS One, 8, e79168.
Wenger, S.J. \& Olden, J.D. (2012) Assessing transferability of ecological models: an underappreciated aspect of statistical validation. Methods in Ecology and Evolution, 3, 260-267.
Widom, B. \& Rowlinson, J.S. (1970) New model for the study of liquid-vapor phase transitions. The Journal of Chemical Physics, 52, 1670-1984.
Zuur, A.F., Ieno, E.N. \& Elphick, C.S. (2010) A protocol for data exploration to avoid common statistical problems. Methods in Ecology and Evolution, 1, 3-14.

Received 20 October 2014; accepted 12 January 2015
Handling Editor: Robert B. O'Hara

\section*{Supporting Information}

Additional Supporting Information may be found in the online version of this article.

Appendix S1. Full description of environmental covariates, derivation of likelihood formulae, and detailed software tutorials.