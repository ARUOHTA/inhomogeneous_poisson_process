\title{
Spatial Modeling with Spatially Varying Coefficient Processes \\ Author(s): Alan E. Gelfand, Hyon-Jung Kim, C. F. Sirmans and Sudipto Banerjee
}

Source: Journal of the American Statistical Association, Jun., 2003, Vol. 98, No. 462 (Jun., 2003), pp. 387-396

Published by: Taylor \& Francis, Ltd. on behalf of the American Statistical Association
Stable URL: https://www.jstor.org/stable/30045248

JSTOR is a not-for-profit service that helps scholars, researchers, and students discover, use, and build upon a wide range of content in a trusted digital archive. We use information technology and tools to increase productivity and facilitate new forms of scholarship. For more information about TSTOR, please contact support@jstor.org.

Your use of the JSTOR archive indicates your acceptance of the Terms \& Conditions of Use, available at https://about.jstor.org/terms

Taylor \& Francis, Ltd. and American Statistical Association are collaborating with JSTOR to digitize, preserve and extend access to Tournal of the American Statistical Association

\title{
Spatial Modeling With Spatially Varying Coefficient Processes
}

\author{
Alan E. Gelfand, Hyon-Jung Kim, C. F. Sirmans, and Sudipto Banerjee
}

\begin{abstract}
In many applications, the objective is to build regression models to explain a response variable over a region of interest under the assumption that the responses are spatially correlated. In nearly all of this work, the regression coefficients are assumed to be constant over the region. However, in some applications, coefficients are expected to vary at the local or subregional level. Here we focus on the local case. Although parametric modeling of the spatial surface for the coefficient is possible, here we argue that it is more natural and flexible to view the surface as a realization from a spatial process. We show how such modeling can be formalized in the context of Gaussian responses providing attractive interpretation in terms of both random effects and explaining residuals. We also offer extensions to generalized linear models and to spatio-temporal setting. We illustrate both static and dynamic modeling with a dataset that attempts to explain (log) selling price of single-family houses.
\end{abstract}

KEY WORDS: Bayesian framework; Multivariate spatial processes; Prediction; Spatio-temporal modeling; Stationary Gaussian process.

\section*{1. INTRODUCTION}

The broad availability of fast, inexpensive computing along with the development of very capable, user-friendly geographic information systems (GIS) software has led to the increased collection of spatial and spatio-temporal data in such diverse fields as real estate/finance, epidemiology, environmetrics/ecology, and communications. In turn, this has fueled increased spatial modeling and data analysis activity within the statistical community.

In many applications, the objective is to build regression models to explain a response variable observed over a region of interest, say \(D\), under the assumption that the responses are spatially associated. That is, whereas some spatial modeling may be accomplished through the mean, it is still anticipated that the responses are dependent and that this dependence becomes stronger as pairs of responses become closer in space. With continuous response that is point referenced, if a normality assumption (perhaps on a transformed scale) seems plausible, this dependence is typically modeled directly using a Gaussian process. The literature here is enormous. The book by Cressie (1993) is perhaps a place to start. In the case of, say, binary or count response, a hierarchical model is often adopted using an exponential family model at the first stage and then introducing normally distributed spatial random effects into the mean structure on a transformed scale related by a link function (see, e.g., Diggle, Tawn, and Moyeed 1998).

In nearly all of this work, the regression coefficients are assumed to be constant across the region. In certain applications, this would not be appropriate. The coefficients may be expected to vary at the local or subregion level. For instance, Assunçao, Gamerman, and Assunçao (1999) introduced a Bayesian space varying parameter model to examine microregion factor productivity and the degree of factor substitution

\footnotetext{
Alan E. Gelfand is a Professor, Institute of Statistics and Decision Sciences, Duke University, Durham, NC 27708-0251. Hyon-Jung Kim is a Research Associate, Department of Mathematical Sciences, University of Oulu, Finland. C. F. Sirmans is Professor of Finance and Director, Center for Real Estate and Urban Economic Studies, University of Connecticut, Storrs, CT 062691041. Sudipto Banerjee is Assistant Professor, Division of Biostatistics, University of Minnesota, Minneapolis, MN 55455. A portion of Gelfand's work was supported by National Science Foundation grant DMS 9971206. Kim was supported in part by National Science Foundation/Environmental Protection Agency grant "Statistical Methods for Environmental Social Sciences" and in part by the Center for Real Estate and Urban Economic Studies.
}
in the Brazilian agriculture. Agarwal, Gelfand, Sirmans, and Thibadeau (2003), in the context of locally stationary spatial modeling, introduced local regression models for factors affecting house price at the school district (and sub-school district) level. A flexible modeling approach for space-varying regression models was developed with a simulation study by Gamerman, Moreira, and Rue (2001). These authors all made the rather restrictive assumption that for a given coefficient, it is constant on specified areal units. The levels of the surface on these units is modeled using independent or conditionally autoregressive specifications. Concerns arise about the arbitrariness of the scale of resolution, the lack of smoothness of the surface, and the inability to interpolate the value of the surface to individual locations. When working with pointreferenced data, it will be more attractive to allow the coefficients to vary by location, to envision a spatial surface for a particular coefficient. For instance, in our application we also model the (log) selling price of single family houses. Customary explanatory variables include the age of the house, the square feet of living area, the square feet of other area, and the number of bathrooms. If the region of interest is a city or greater metropolitan area, then it is evident that the capitalization rate (e.g., for age), will vary across the region. Older houses will have higher value in some parts of the region than in other parts. By allowing the coefficient of age to vary with location, we can remedy the foregoing concerns. With practical interest in mind, say real estate appraisal, we can predict the coefficient for arbitrary properties, not just for those that sold during the period of investigation. Similar issues arise in modeling environmental exposure to a particular pollutant, where covariates might include temperature and precipitation.

One possible approach would be to model the spatial surface for the coefficient parametrically. In the simplest case this would require the rather arbitrary specification of a polynomial surface function; a range of surfaces too limited or inflexible might result. More flexibility could be introduced using a spline surface over two-dimensional space (see, e.g., Luo and Wahba 1998). However, this requires selection of a spline function and determination of the number of and locations of the knots in the
space. Also, with multiple coefficients, a multivariate specification of a spline surface is required. The approach that we adopt here is arguably more natural and at least as flexible. We model the spatially varying coefficient surface as a realization from a spatial process. For multiple coefficients, we use a multivariate spatial process model. In fact, we use a stationary specification in which desired degree of smoothness of process realization can be modeled through the choice of covariance function (e.g., the Matèrn class). Kent (1989) and Stein (1999a) have provided discussions of univariate process; Banerjee and Gelfand (2003), of multivariate process. A nonstationary model results for the data.

We adopt a Bayesian approach for our modeling framework. This is attractive in the proposed setting, because we are specifically interested in inference for the random spatial effects. In particular, we obtain an entire posterior for the spatial coefficient process at both observed and unobserved locations, as well as posteriors for all model parameters. Interpolation for a process that is neither observed nor arising as a residual seems inaccessible in any other framework. For Gaussian responses, although some inference is possible through likelihood methods, it is more limited, and in particular, interval estimation relies on possibly inappropriate asymptotics. For non-Gaussian data, approximations will almost surely be required, possibly in the form that we provide in Section 7, but still inference will be limited.

To clarify interpretation and implementation, we first develop our general approach in the case of a single covariate, and thus we have two spatially varying coefficient processes, one for "intercept" and one for "slope." We then turn to the case of multiple covariates. Because even in the basic multiple regression setting, coefficient estimates typically reveal some strong correlations, the collection of spatially varying coefficient processes is expected to be dependent. Hence we use a multivariate process model. Indeed, we present a further generalization to build a spatial analog of a multivariate regression model (see, e.g., Goldstein 1995). We also consider flexible spatio-temporal possibilities. The previously mentioned real estate setting provides site level covariates whose coefficients are of considerable practical interest and a dataset of single-family home sales from Baton Rouge, LA enables illustration. Except for regions exhibiting special topography, we anticipate that a spatially varying coefficient model will prove more useful than, for instance, a trend surface model. That is, incorporating a polynomial in latitude and longitude into the mean structure would not be expected to serve as a surrogate for allowing the variability across the region of a coefficient for say age or living area of a house.

The article is organized as follows. Section 2 details the Gaussian modeling approach for a single covariate. Section 3 addresses the multiple-covariate case. Section 4 proposes a sequence of varying coefficient specifications in the spatio-temporal setting. Section 5 comments briefly on model comparison. Section 6 presents an example using the aforementioned single-family home sales data. Section 7 concludes with some discussion of the generalized linear model setting.

\section*{2. THE MODELING APPROACH FOR A SINGLE COVARIATE}

Recall the usual Gaussian stationary spatial process model as in, for example, Cressie (1993),
\[
Y(\mathbf{s})=\mu(\mathbf{s})+W(\mathbf{s})+\epsilon(\mathbf{s}),
\]
where \(\mu(\mathbf{s})=\mathbf{x}(\mathbf{s})^{T} \beta\) and \(\epsilon(\mathbf{s})\) is a white noise process, that is, \(\mathrm{E}(\epsilon(\mathbf{s}))=0, \operatorname{var}(\epsilon(\mathbf{s}))=\tau^{2}, \operatorname{cov}\left(\epsilon(\mathbf{s}), \epsilon\left(\mathbf{s}^{\prime}\right)\right)=0\), and \(W(\mathbf{s})\) is a second-order stationary mean 0 process independent of the white noise process; that is, \(\mathrm{E}(W(\mathbf{s}))=0\), \(\operatorname{var}(W(\mathbf{s}))=\sigma^{2}, \operatorname{cov}\left(W(\mathbf{s}), W\left(\mathbf{s}^{\prime}\right)\right)=\sigma^{2} \rho\left(\mathbf{s}, \mathbf{s}^{\prime} ; \phi\right)\), where \(\rho\) is a valid two-dimensional correlation function.

The \(W(\mathbf{s})\) are viewed as spatial random effects, and (1) implicitly defines a hierarchical model. Letting \(\mu(\mathbf{s})=\beta_{0}+ \beta_{1} x(\mathbf{s})\), write \(W(\mathbf{s})=\beta_{0}(\mathbf{s})\) and define \(\tilde{\beta}_{0}(\mathbf{s})=\beta_{0}+\beta_{0}(\mathbf{s})\). Then \(\beta_{0}(\mathbf{s})\) can be interpreted as a random spatial adjustment at location \(\mathbf{s}\) to the overall intercept \(\beta_{0}\). Equivalently, \(\tilde{\beta}_{0}(\mathbf{s})\) can be viewed as a random intercept process. For an observed set of locations \(\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{n}\) given \(\beta_{0}, \beta_{1},\left\{\beta_{0}\left(\mathbf{s}_{i}\right)\right\}\) and \(\tau^{2}\), the \(Y\left(\mathbf{s}_{i}\right)=\beta_{0}+\beta_{1} x\left(\mathbf{s}_{\mathbf{i}}\right)+\beta_{0}\left(\mathbf{s}_{\mathbf{i}}\right)+\epsilon\left(\mathbf{s}_{\mathbf{i}}\right), i=1, \ldots, n\), are conditionally independent. The first-stage likelihood is
\[
\begin{aligned}
L\left(\beta_{0}, \beta_{1},\left\{\beta_{0}\left(\mathbf{s}_{i}\right)\right\}, \tau^{2} ; \mathbf{y}\right)= & \left(\tau^{2}\right)^{-\frac{n}{2}} \exp \left\{-\frac{1}{2 \tau^{2}} \sum\left(Y\left(\mathbf{s}_{i}\right)\right.\right. \\
& \left.\left.-\left(\beta_{0}+\beta_{1} x\left(\mathbf{s}_{i}\right)+\beta_{0}\left(\mathbf{s}_{i}\right)\right)^{2}\right)\right\}
\end{aligned}
\]

In obvious notation, the distribution of \(\boldsymbol{\beta}_{0}=\left(\beta_{0}\left(\mathbf{s}_{1}\right), \ldots\right.\), \(\left.\beta_{0}\left(\mathbf{s}_{n}\right)\right)^{T}\) is
\[
f\left(\beta_{0} \mid \sigma_{0}^{2}, \phi_{0}\right)=N\left(\mathbf{0}, \sigma_{0}^{2} H_{0}\left(\phi_{0}\right)\right),
\]
where \(\left(H_{0}\left(\phi_{0}\right)\right)_{i j}=\rho_{0}\left(\mathbf{s}_{i}-\mathbf{s}_{j} ; \phi_{0}\right)\). For all of the discussion and examples that follow, we adopt the Matern correlation function, \(\rho(\mathbf{h}, \phi) \propto(\gamma(\|\mathbf{h}\|))^{v} K_{v}(\gamma(\|\mathbf{h}\|))\). Here \(K_{v}\) is a modified Bessel function, \(\phi=(\gamma, v)\), where \(\gamma\) is a decay parameter and \(v\) is a smoothness parameter (see Stein 1999a for a more in-depth discussion). With a prior on \(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}\), and \(\phi_{0}\), specification of the Bayesian hierarchical model is completed. Under (2) and (3), we can integrate over \(\boldsymbol{\beta}_{0}\), obtaining the marginal likelihood
\[
\begin{aligned}
& L\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \phi_{0} ; \mathbf{y}\right) \\
& \quad=\left|\sigma_{0}^{2} H_{0}\left(\phi_{0}\right)+\tau^{2} I\right|^{-\frac{1}{2}} \\
& \quad \quad \times \exp \left\{-\frac{1}{2}\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)^{T}\left(\sigma_{0}^{2} H_{0}\left(\phi_{0}\right)+\tau^{2} I\right)^{-1}\right. \\
& \left.\quad \quad \times\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)\right\}
\end{aligned}
\]
where \(\mathbf{x}=\left(x\left(\mathbf{s}_{1}\right), \ldots, x\left(\mathbf{s}_{n}\right)\right)^{T}\).
We note analogies with usual Gaussian random-effects models where \(Y_{i j}=\beta_{0}+\beta_{1} x_{i j}+\alpha_{i}+e_{i j}\), with \(\alpha_{i}\) iid \(N\left(0, \sigma_{\alpha}{ }^{2}\right)\) and \(\epsilon_{i j}\) iid \(N\left(0, \sigma_{\epsilon}{ }^{2}\right)\). In this case, replications are needed to identify (separate) the variance components. Because of the dependence between the \(\beta_{0}\left(s_{i}\right)\), replications are not needed in the spatial case, as (4) reveals. Also, if \(U(\mathbf{s})\) denotes the total error in the regression model, then \(U(\mathbf{s})\) is partitioned into "intercept process" error and "pure" error.

If the Bayesian model is fitted using the marginal likelihood in (4) with simulation-based model fitting, then samples essentially from the posterior \(f\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \phi \mid y\right)\) are obtained. But then samples from \(\boldsymbol{\beta}_{0} \mid \mathbf{y}\) can be obtained one-for-one because
\[
\begin{aligned}
f\left(\boldsymbol{\beta}_{0} \mid \mathbf{y}\right)=\int f\left(\boldsymbol{\beta}_{0} \mid \beta_{0}, \beta_{1}, \tau^{2}\right. & \left., \sigma_{0}^{2}, \phi_{0}, \mathbf{y}\right) \\
& \times f\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \phi_{0} \mid \mathbf{y}\right)
\end{aligned}
\]
where
\[
\begin{aligned}
& f\left(\beta_{0} \mid \beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \phi_{0}, \mathbf{y}\right) \\
& \quad=N\left(\left(\frac{1}{\tau^{2}} I+\frac{1}{\sigma_{0}^{2}} H_{0}^{-1}\left(\phi_{0}\right)\right)^{-1} \frac{1}{\tau^{2}}\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)\right. \\
& \left.\quad\left(\frac{1}{\tau^{2}} I+\frac{1}{\sigma_{0}^{2}} H_{0}^{-1}\left(\phi_{0}\right)\right)^{-1}\right)
\end{aligned}
\]

We can also obtain samples from the posterior of the \(\beta_{0}(\mathbf{s})\) process at a new location, say \(\mathbf{s}_{\text {new }}\), to provide interpolation for the \(\beta_{0}(\mathbf{s})\) surface. Specifically,
\[
f\left(\beta_{0}\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)=\int f\left(\beta_{0}\left(\mathbf{s}_{\text {new }}\right) \mid \boldsymbol{\beta}_{0}, \sigma_{0}^{2}, \phi_{0}\right) f\left(\boldsymbol{\beta}_{0}, \sigma_{0}^{2}, \phi_{0} \mid \mathbf{y}\right)
\]

The first density under the integral is a univariate normal that can be written down directly from the specification of \(\beta_{0}(\mathbf{s})\). For the prediction of \(y\left(\mathbf{s}_{\text {new }}\right)\) given \(\mathbf{y}\), we require
\[
\begin{aligned}
f\left(y\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)=\int & f\left(y\left(\mathbf{s}_{\text {new }}\right) \mid \beta_{0}, \beta_{1}, \beta_{0}\left(\mathbf{s}_{\text {new }}\right), \tau^{2}\right) \\
& \times f\left(\beta_{0}\left(\mathbf{s}_{\text {new }}\right) \mid \beta_{0}, \sigma_{0}^{2}, \phi_{0}\right) \\
& \times f\left(\beta_{0}, \beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \phi_{0} \mid \mathbf{y}\right)
\end{aligned}
\]

The first term under the integral sign is a normal density. Again, it is straightforward to obtain samples from this predictive distribution.

The foregoing development immediately suggests how to formulate a spatially varying coefficient model. Suppose that we write
\[
Y(\mathbf{s})=\beta_{0}+\beta_{1} x(\mathbf{s})+\beta_{1}(\mathbf{s}) x(\mathbf{s})+\epsilon(\mathbf{s}) .
\]

In (8), \(\beta_{1}(\mathbf{s})\) is a second-order stationary mean 0 Gaussian process with variance \(\sigma_{1}^{2}\) and correlation function \(\rho_{1}\left(\cdot ; \phi_{1}\right)\). Also, let \(\tilde{\beta}_{1}(\mathbf{s})=\beta_{1}+\beta_{1}(\mathbf{s})\). Now \(\beta_{1}(\mathbf{s})\) can be interpreted as a random spatial adjustment at location \(\mathbf{s}\) to the overall slope \(\beta_{1}\). Equivalently, \(\tilde{\beta}_{1}(\mathbf{s})\) can be viewed as a random slope process. In effect, we are using an infinite-dimensional function to explain the relationship between \(x(\mathbf{s})\) and \(Y(\mathbf{s})\).

Expression (8) yields an obvious modification of (2) and (3). In particular, the resulting marginalized likelihood becomes
\[
\begin{aligned}
& L\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{1}^{2}, \phi_{1} ; \mathbf{y}\right) \\
& =\left|\sigma_{1}^{2} D_{x} H_{1}\left(\phi_{1}\right) D_{x}+\tau^{2} I\right|^{-\frac{1}{2}} \\
& \quad \quad \times \exp \left\{-\frac{1}{2}\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)^{T}\left(\sigma_{1}^{2} D_{x} H_{1}\left(\phi_{1}\right) D_{x}+\tau^{2} I\right)^{-1}\right. \\
& \left.\quad \quad \times\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)\right\}
\end{aligned}
\]
where \(D_{x}\) is diagonal with \(\left(D_{x}\right)_{i i}=x\left(\mathbf{s}_{i}\right)\). Moreover, with \(\beta_{1}=\left(\beta_{1}\left(\mathbf{s}_{1}\right), \ldots, \beta_{1}\left(\mathbf{s}_{n}\right)\right)^{T}\), we can sample \(f\left(\beta_{1} \mid \mathbf{y}\right)\) and \(f\left(\beta_{1}\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)\) using obvious analogs of (5) and (6).

Note that (8) provides a heterogeneous, nonstationary process for the data regardless of the choice of covariance function for the \(\beta_{1}(\mathbf{s})\) process. Here, \(\operatorname{var}\left(Y(\mathbf{s}) \mid \beta_{0}, \beta_{1}, \tau^{2}, \sigma_{1}^{2}, \phi_{1}\right)= x^{2}(\mathbf{s}) \sigma_{1}^{2}+\tau^{2}\) and \(\operatorname{cov}\left(Y(\mathbf{s}), Y\left(\mathbf{s}^{\prime}\right) \mid \beta_{0}, \beta_{1}, \tau^{2}, \sigma_{1}^{2}, \phi_{1}\right)= \sigma_{1}^{2} x(\mathbf{s}) x\left(\mathbf{s}^{\prime}\right) \rho_{1}\left(\mathbf{s}-\mathbf{s}^{\prime} ; \phi_{1}\right)\). As a result, we observe that in practice, (8) is sensible only if we have \(x(\mathbf{s})>0\). In fact, centering and scaling, usually advocated for better-behaved model fitting, is inappropriate here. With centered \(x(\mathbf{s})\) 's, we would find the likely untenable behavior that \(\operatorname{var}(Y(\mathbf{s}))\) decreases and then increases in \(x(\mathbf{s})\). Worse, for an essentially central \(x(\mathbf{s})\), we would find \(Y(\mathbf{s})\) essentially independent of \(Y\left(\mathbf{s}^{\prime}\right)\) for any \(\mathbf{s}^{\prime}\). Also, scaling the \(x(\mathbf{s})\) 's accomplishes nothing. \(\beta_{1}(\mathbf{s})\) would be inversely rescaled, because the model identifies only \(\beta_{1}(\mathbf{s}) x(\mathbf{s})\).

This leads to concerns regarding possible approximate collinearity of \(\mathbf{x}\), the vector of \(x\left(\mathbf{s}_{i}\right)\) 's, with the vector \(\mathbf{1}\). Expression (9) shows that a badly behaved likelihood will arise if \(\mathbf{x} \approx c \mathbf{1}\). Fortunately, we can reparameterize (8) to \(Y(\mathbf{s})= \beta_{0}^{\prime}+\beta_{1}^{\prime} \tilde{x}(\mathbf{s})+\beta_{1}(\mathbf{s}) x(\mathbf{s})+\epsilon(\mathbf{s})\), where \(\tilde{x}(\mathbf{s})\) is centered and scaled with obvious definitions for \(\beta_{0}^{\prime}\) and \(\beta_{1}^{\prime}\). Now \(\tilde{\beta_{1}}(\mathbf{s})= \beta_{1}^{\prime} / s_{x}+\beta_{1}(\mathbf{s})\), where \(s_{x}\) is the sample standard deviation of the \(x(\mathbf{s})\) 's.

As after (4), we can draw an analogy with usual longitudinal linear growth curve modeling where \(Y_{i j}=\beta_{0}+ \beta_{1} x_{i j}+\beta_{1 i} x_{i j}+\epsilon_{i j}\), that is, a random slope for each individual. Also, \(U(\mathbf{s})\), the total error in the regression model (8), is now partitioned into "slope process" error and "pure" error.

The general specification encompassing (1) and (7) would be
\[
Y(\mathbf{s})=\beta_{0}+\beta_{1} x(\mathbf{s})+\beta_{0}(\mathbf{s})+\beta_{1}(\mathbf{s}) x(\mathbf{s})+\epsilon(\mathbf{s}) .
\]

Expression (10) parallels the usual linear growth curve modeling by introducing both an intercept process and a slope process. The model in (10) requires a bivariate process specification to determine the joint distribution of \(\boldsymbol{\beta}_{0}\) and \(\boldsymbol{\beta}_{1}\). We return to this in Section 3, but, under independence of the processes, (10) can be easily marginalized over \(\beta_{0}\) and \(\beta_{1}\), yielding
\[
\begin{aligned}
& L\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \sigma_{1}^{2}, \phi_{0}, \phi_{1} ; \mathbf{y}\right) \\
& =\left|\sigma_{0}^{2} H_{0}\left(\phi_{0}\right)+\sigma_{1}^{2} D_{x} H_{1}\left(\phi_{1}\right) D_{x}+\tau^{2} I\right|^{-\frac{1}{2}} \\
& \quad \times \exp \left\{-\frac{1}{2}\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)^{T}\right. \\
& \quad \quad \times\left(\sigma_{0}^{2} H_{0}\left(\phi_{0}\right)+\sigma_{1}^{2} D_{x} H_{1}\left(\phi_{1}\right) D_{x}+\tau^{2} I\right)^{-1} \\
& \left.\quad \times\left(\mathbf{y}-\beta_{0} \mathbf{1}-\beta_{1} \mathbf{x}\right)\right\}
\end{aligned}
\]

Again, simulation from \(f\left(\beta_{0} \mid \mathbf{y}\right), f\left(\beta_{1} \mid \mathbf{y}\right), f\left(\beta_{0}\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)\), \(f\left(\beta_{1}\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)\), and \(f\left(y\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)\) is straightforward. The possibility of predicting the spatial surface at arbitrary locations makes a compelling case for a Bayesian inference approach. Classical kriging methods cannot address this problem, because no \(\beta_{0}\left(\mathbf{s}_{i}\right)\) or \(\beta_{1}\left(\mathbf{s}_{i}\right)\) are observed. Also, in (10), the total error has been partitioned into three independent pieces
with obvious interpretation. The relative sizes of the error components at any \(\mathbf{s}_{i}\) can be studied through \(E\left(\beta_{0}\left(\mathbf{s}_{i}\right) \mid \mathbf{y}\right)\), \(E\left(\beta_{1}\left(\mathbf{s}_{i}\right) x\left(\mathbf{s}_{i}\right) \mid \mathbf{y}\right)\), and \(E\left(\epsilon\left(\mathbf{s}_{i}\right) \mid \mathbf{y}\right)\). To compare the relative variability contributions, we could calculate \(E\left(\tau^{2} \mid \mathbf{y}\right)\), \(E\left(\sigma_{0}^{2} \mid \mathbf{y}\right)\), and \(\bar{x}^{2} E\left(\sigma_{1}^{2} \mid \mathbf{y}\right)\), where \(\bar{x}=\sum x\left(\mathbf{s}_{i}\right) / n\). In the case where \(\rho_{0}\) and \(\rho_{1}\) are isotropic, posteriors for the ranges can be obtained to compare spatial range of the intercept process with that of the slope process. The posterior mean surfaces \(E\left(\beta_{0}(\mathbf{s}) \mid \mathbf{y}\right)\) and \(E\left(\beta_{1}(\mathbf{s}) \mid \mathbf{y}\right)\) can be obtained over an arbitrary grid of locations and displayed graphically using standard software.

Bayesian models using the marginal likelihoods in either (4), (9), or (11) are easily fitted using a slice Gibbs sampler as discussed by Agarwal and Gelfand (2001) (see Neal 2002, as well). In particular, under, say (11), we add an auxiliary variable \(U \sim \operatorname{unif}\left(0, L\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \sigma_{1}^{2}, \phi_{0}, \phi_{1} ; \mathbf{y}\right)\right)\). Then the posterior \(f\left(\beta_{0}, \beta_{1}, \tau^{2}, \sigma_{0}^{2}, \sigma_{1}^{2}, \phi_{0}, \phi_{1}, U \mid \mathbf{y}\right) \propto I(U<L) \pi\left(\beta_{0}\right.\), \(\left.\beta_{1}, \tau^{2}, \sigma_{0}^{2}, \sigma_{1}^{2}, \phi_{0}, \phi_{1}\right)\), where \(\pi\) is the prior on the marginalized parameters. The slice Gibbs sampler updates \(U\) with a uniform draw and updates the other parameters with a prior draw subject to the indicator restriction. This algorithm is "off-the-shelf," requiring no tuning. It converges faster than Metropolis alternatives and avoids the autocorrelation problem that often arises with these alternatives.

\section*{3. A MULTIVARIATE SPATIALLY VARYING COEFFICIENT MODEL}

Here we turn to the case of a \(p \times 1\) multivariate covariate vector \(\mathbf{X}(\mathbf{s})\) at location \(\mathbf{s}\), where, for convenience, \(\mathbf{X}(\mathbf{s})\) includes a 1 as its first entry to accommodate an intercept. We generalize (10) to
\[
Y(\mathbf{s})=\mathbf{X}^{T}(\mathbf{s}) \tilde{\beta}(\mathbf{s})+\epsilon(\mathbf{s})
\]
where \(\tilde{\boldsymbol{\beta}}(\mathbf{s})\) is assumed to follow a \(p\)-variate spatial process model. With observed locations \(\mathbf{s}_{1}, \mathbf{s}_{2}, \ldots, \mathbf{s}_{n}\), let \(X^{T}\) be \(n \times n p\) block diagonal having as block for the \(i\) th row \(\mathbf{X}^{T}\left(\mathbf{s}_{i}\right)\). Then we can write \(\mathbf{Y}=X^{T} \tilde{\beta}+\epsilon\), where \(\tilde{\beta}\) is \(n p \times 1\), the concatenated vector of the \(\tilde{\beta}(\mathbf{s})\) and \(\epsilon \sim N\left(0, \tau^{2} I\right)\).

In practice, to assume that the component processes of \(\tilde{\boldsymbol{\beta}}(\mathbf{s})\) are independent is likely inappropriate. That is, in the simpler case of simple linear regression, negative association between slope and intercept is usually seen. (This is intuitive if one envisions overlaying random lines that are likely relative to a fixed scattergram of data points.) The dramatic improvement in model performance when dependence is incorporated is shown in the example of Section 6. To formulate a multivariate Gaussian process for \(\tilde{\boldsymbol{\beta}}(\mathbf{s})\), we require the mean and the cross-covariance function. For the former, following Section 2, we take this to be \(\mu_{\beta}=\left(\beta_{1}, \ldots, \beta_{p}\right)^{T}\). For the latter, we require a valid \(p\)-variate choice. In the sequel we work with a computationally convenient separable choice following Mardia and Goodall (1993) (see also Banerjee and Gelfand 2002 in this regard).

More precisely, let \(\boldsymbol{C}\left(\mathbf{s}, \mathbf{s}^{\prime}\right)\) be the \(p \times p\) matrix with \((l, m)\) entry \(\operatorname{cov}\left(\tilde{\boldsymbol{\beta}}_{l}(\mathbf{s}), \tilde{\boldsymbol{\beta}_{m}}\left(\mathbf{s}^{\prime}\right)\right)\) and let
\[
C\left(\mathbf{s}, \mathbf{s}^{\prime}\right)_{l m}=\rho\left(\mathbf{s}-\mathbf{s}^{\prime} ; \phi\right) \tau_{l m},
\]
where \(\rho\) is a valid scalar correlation function in two dimensions and \(T_{p \times p}\) is such that \((T)_{l m}=\tau_{l m}\) is positive-definite symmet-
ric. In other words, \(T\) is the covariance matrix associated with an observation vector at any spatial location, and \(\rho\) captures the attenuation in association across space. If we collect the set of \(\rho\left(\mathbf{s}-\mathbf{s}^{\prime} ; \phi\right)\) into an \(n \times n\) matrix \(H(\phi)\) as in the previous section, then, with \(\otimes\) denoting the Kronecker product, the distribution of \(\beta\) is
\[
\tilde{\boldsymbol{\beta}} \sim N\left(\mathbf{1}_{n \times 1} \otimes \boldsymbol{\mu}_{\beta}, H(\phi) \otimes T\right) .
\]

As in the previous section, if \(\tilde{\boldsymbol{\beta}}=\boldsymbol{\beta}+\mathbf{1}_{n \times 1} \otimes \boldsymbol{\mu}_{\beta}\), then we can write (12) as
\[
Y(\mathbf{s})=\mathbf{X}^{T}(\mathbf{s}) \mu_{\beta}+\mathbf{X}^{T}(\mathbf{s}) \beta(\mathbf{s})+\epsilon(\mathbf{s})
\]

In (15), the total error in the regression model is partitioned into \(p+1\) pieces each with an obvious interpretation. Following Section 2, using (12) and (14), we can integrate over \(\beta\) to obtain
\[
\begin{aligned}
& L\left(\boldsymbol{\mu}_{\beta}, \tau^{2}, T, \phi ; \mathbf{y}\right) \\
& \quad=\left|X(H(\phi) \otimes T) X^{T}+\tau^{2} I\right|^{-\frac{1}{2}} \\
& \quad \times \exp \left\{-\frac{1}{2}\left(\mathbf{y}-X\left(\mathbf{1} \otimes \boldsymbol{\mu}_{\beta}\right)\right)^{T}\left(X(H(\phi) \otimes T) X^{T}+\tau^{2} I\right)^{-1}\right. \\
& \left.\quad \times\left(\mathbf{y}-X\left(\mathbf{1} \otimes \boldsymbol{\mu}_{\beta}\right)\right)\right\}
\end{aligned}
\]

The possibly daunting form (16) still involves only \(n \times n\) matrices.

The Bayesian model is completed with a prior \(f\left(\boldsymbol{\mu}_{\beta}\right.\), \(\left.\tau^{2}, T, \phi\right)\) which we assume to take the product form \(f\left(\mu_{\beta}\right) f\left(\tau^{2}\right) f(T) f(\phi)\). Later, these components are normal, inverse gamma, inverse Wishart and gamma, and gamma, where \(\phi=(\gamma, v)\) under the Matèrn correlation function. The model using (16) and such a prior is readily fitted through a sliced Gibbs sampler similar to that described at the end of the previous section.

With regard to prediction analogous to (5), \(f(\tilde{\boldsymbol{\beta}} \mid \mathbf{y})\) can be sampled one-for-one with the posterior samples from \(f\left(\boldsymbol{\mu}_{\beta}, \tau^{2}, T, \phi \mid \mathbf{y}\right)\) using \(f\left(\tilde{\boldsymbol{\beta}} \mid \boldsymbol{\mu}_{\beta}, \tau^{2}, T, \phi, y\right)\), which is \(N(B \mathbf{b}, B)\), where \(B=\left(X^{T} X / \tau^{2}+H^{-1}(\phi) \otimes T^{-1}\right)^{-1}\) and \(\mathbf{b}=X^{T} \mathbf{y} / \tau^{2}+\left(H^{-1}(\phi) \otimes T^{-1}\right)\left(\mathbf{1} \otimes \boldsymbol{\mu}_{\beta}\right)\). \(B\) is \(n p \times n p\), but for sampling \(\tilde{\boldsymbol{\beta}}\), only a Cholesky decomposition of \(B\) is needed and only for the retained posterior samples. Prediction at a new location, say \(\mathbf{s}_{\text {new }}\), requires, analogous to (6), \(f\left(\tilde{\boldsymbol{\beta}}\left(\mathbf{s}_{\text {new }}\right) \mid\right. \left.\tilde{\boldsymbol{\beta}}, \boldsymbol{\mu}_{\beta}, \tau^{2}, T, \phi\right)\). Defining \(\mathbf{h}_{\text {new }}(\phi)\) to be the \(n \times 1\) vector with \(i\) th row entry \(\rho\left(\mathbf{s}_{i}-\mathbf{s}_{\text {new }} ; \phi\right)\), this distribution is normal with mean \(\boldsymbol{\mu}_{\beta}+\left(\mathbf{h}_{\text {new }}^{T}(\phi) \otimes T\right)\left(H^{-1}(\phi) \otimes T^{-1}\right)\left(\tilde{\boldsymbol{\beta}}-\mathbf{1}_{n x 1} \otimes \boldsymbol{\mu}_{\beta}\right)= \boldsymbol{\mu}_{\beta}+\left(\mathbf{h}_{\text {new }}^{T}(\phi) H^{-1}(\phi) \otimes I\right)\left(\tilde{\boldsymbol{\beta}}-\mathbf{1}_{n x 1} \otimes \boldsymbol{\mu}_{\beta}\right)\) and covariance matrix \(T-\left(\mathbf{h}_{\text {new }}^{T}(\phi) \otimes T\right)\left(H^{-1}(\phi) \otimes T^{-1}\right)\left(\mathbf{h}_{\text {new }}(\phi) \otimes T\right)= \left(I-\mathbf{h}_{\text {new }}^{T}(\phi) H^{-1}(\phi) \mathbf{h}_{\text {new }}(\phi)\right) T\). Finally, the predictive distribution for \(Y\left(\mathbf{s}_{\text {new }}\right), f\left(Y\left(\mathbf{s}_{\text {new }}\right) \mid \mathbf{y}\right)\) is sampled analogously to (7).

We conclude this section by noting an extension of (12) when we have repeated measurements at location \(\mathbf{s}\). That is, suppose that we have
\[
Y(\mathbf{s}, l)=\mathbf{X}^{T}(\mathbf{s}, l) \beta(\mathbf{s})+\epsilon(\mathbf{s}, l)
\]
where \(l=1, \ldots, L_{\mathbf{s}}\) with \(L_{\mathbf{s}}\) the number of measurements at \(\mathbf{s}\) and the \(\epsilon(\mathbf{s}, l)\) still white noise. As an illustration, in the real estate context that we mentioned earlier, \(\mathbf{s}\) might denote the location for an apartment block and \(l\) indexes apartments in this
block that have sold, with the \(l\) th apartment having characteristics \(\mathbf{X}(\mathbf{s}, l)\). Suppose further that \(\mathbf{Z}(\mathbf{s})\) denotes an \(r \times 1\) vector of site level characteristics. For an apartment block, these characteristics might include amenities provided or distance to the central business district. Then (17) can be extended to a multilevel model in the sense of Goldstein (1995) or Raudenbush and Bryk (2002). In particular, we can write
\[
\boldsymbol{\beta}(\mathbf{s})=\left(\begin{array}{c}
\mathbf{Z}^{T}(\mathbf{s}) \boldsymbol{\gamma}_{1} \\
\vdots \\
\mathbf{Z}^{T}(\mathbf{s}) \boldsymbol{\gamma}_{p}
\end{array}\right)+\mathbf{W}(\mathbf{s}) .
\]

In (18), \(\boldsymbol{\gamma}_{j}\) is an \(r \times 1\) vector associated with \(\tilde{\beta}_{j}(\mathbf{s})\) and \(\mathbf{W}(\mathbf{s})\) is a mean \(\mathbf{0}\) multivariate Gaussian spatial process as, for example, earlier. In (18), if the \(\mathbf{W}(\mathbf{s})\) were independent, then we would have a usual multilevel model specification. In the case where \(\mathbf{Z}(\mathbf{s})\) is a scalar capturing just an intercept, we return to the initial model of this section.

\section*{4. SPATIALLY VARYING COEFFICIENTS MODELS WITH SPATIO-TEMPORAL DATA}

A natural extension of the modeling of the previous sections is to the case in which data are correlated at spatial locations across time. Such data frequently arise in ecological, environmental, and meteorological settings. If we assume that time is discretized to a finite set of equally spaced points on a scale, then we can conceptualize a time series of spatial processes that are observed only at the spatial locations \(\mathbf{s}_{1}, \ldots, \mathbf{s}_{n}\).

Adopting a general notation that parallels (10), let
\[
Y(\mathbf{s}, t)=\mathbf{X}^{T}(\mathbf{s}, t) \tilde{\beta}(\mathbf{s}, t)+\epsilon(\mathbf{s}, t), \quad t=1,2, \ldots, M .
\]

That is, we introduce spatio-temporally varying intercepts and spatio-temporally varying slopes. Alternatively, if we write \(\beta(\mathbf{s}, t)=\boldsymbol{\beta}(\mathbf{s}, t)+\boldsymbol{\mu}_{\beta}\), then we are partitioning the total error into \(p+1\) spatio-temporal intercept pieces including \(\epsilon(\mathbf{s}, t)\), each with an obvious interpretation. So we continue to assume that the \(\epsilon(\mathbf{s}, t)\) are iid \(N\left(0, \tau^{2}\right)\), but we need to specify a model for \(\tilde{\beta}(\mathbf{s}, t)\). Regardless, (19) defines a nonstationary process with \(E(Y(\mathbf{s}, t))=\mathbf{X}^{T}(\mathbf{s}, t) \tilde{\boldsymbol{\beta}}(s, t), \operatorname{var}(Y(\mathbf{s}, t))= \mathbf{X}^{T}(\mathbf{s}, t) \Sigma_{\tilde{\boldsymbol{\beta}}(\mathbf{s}, t)} \mathbf{X}(\mathbf{s}, t)+\tau^{2}\), and \(\operatorname{cov}\left(Y(\mathbf{s}, t), Y\left(\mathbf{s}^{\prime}, t^{\prime}\right)\right)= \mathbf{X}^{T}(\mathbf{s}, t) \Sigma_{\tilde{\boldsymbol{\beta}}(\mathbf{s}, t), \tilde{\boldsymbol{\beta}}\left(\mathbf{s}^{\prime}, t^{\prime}\right)} \mathbf{X}\left(\mathbf{s}^{\prime}, t^{\prime}\right)\).

We propose four models for \(\beta(\mathbf{s}, t)\). Paralleling the customary longitudinal data modeling assumption when the time series are usually short, we could set
model 1: \(\beta(\mathbf{s}, t)=\beta(\mathbf{s})\),
where \(\beta\) (s) is modeled as in the previous sections. Model 1 can be viewed as a local linear growth curve model.

Next, we have
model 2: \(\beta(\mathbf{s}, t)=\beta(\mathbf{s})+\alpha(t)\),
where \(\beta(\mathbf{s})\) is again as in model 1. In modeling \(\boldsymbol{\alpha}(t)\), we examine two possibilities. The first possibility treats the \(\alpha_{k}(t)\) as time dummy variables, taking this set of \(p M\) variables to be a priori independent and identically distributed. The second possibility models the \(\boldsymbol{\alpha}(t)\) as a random walk or autoregressive process. The components could be assumed to be independent across \(k\), but for greater generality, we take them to be dependent, using an analog of (13), replacing \(\mathbf{s}\) with \(t\) and \(\rho\), now a valid
correlation function in one dimension. Such additivity in space and time has been discussed in the case of usual modeling of the error structure by, for example, Gelfand, Ecker, Knight, and Sirmans (2003).

In model 3 we consider an analog of the nested-effects areal unit specification of Waller, Carlin, Xia, and Gelfand (1997). (see also Gelfand et al. 2002). In particular, we have
model 3: \(\beta(\mathbf{s}, t)=\boldsymbol{\beta}^{(t)}(\mathbf{s})\).
Here we have spatially varying coefficient processes nested within time. The processes are assumed to be independent across \(t\) (essentially time dummy processes) and permit temporal evolution of the coefficient process. Following Section 3, the process \(\boldsymbol{\beta}^{(t)}(\mathbf{s})\) would be mean 0 , second-order stationary Gaussian with cross-covariance specification at time \(t\), \(C^{(t)}\left(\mathbf{s}, \mathbf{s}^{\prime}\right)\), where \(\left(C^{(t)}\left(\mathbf{s}, \mathbf{s}^{\prime}\right)\right)_{l m}=\rho\left(\mathbf{s}-\mathbf{s}^{\prime} ; \phi^{(t)}\right) \tau_{l m}^{(t)}\). We have specified model 3 with a common \(\mu_{\beta}\) across time, which enables some comparability with the other models that we have proposed. However, we can increase flexibility by replacing \(\mu_{\beta}\) with \(\mu_{\beta}^{(t)}\).

Finally, model 4 proposes a separable covariance specification in space and time, extending work of Gelfand, Zhu, and Carlin (2001):
model 4: \(\beta(\mathbf{s}, t)\) such that \(\Sigma_{\left[\boldsymbol{\beta}(\mathbf{s}, t), \boldsymbol{\beta}\left(\mathbf{s}^{\prime}, t^{\prime}\right)\right]}=\rho^{(1)}\left(\mathbf{s}-\mathbf{s}^{\prime} ; \phi\right) \times \rho^{(2)}\left(t-t^{\prime} ; \gamma\right) T\),
where \(\rho^{(1)}\) is a valid two-dimensional correlation function, \(\rho^{(2)}\) is a valid one-dimensional choice, and \(T\) is positive-definite symmetric. Here \(\rho^{(1)}\) obtains spatial association as in the earlier sections, which is attenuated across time by \(\rho^{(2)}\). The resulting covariance matrix for the full vector \(\beta\), blocked by site and time within site, has the convenient form \(H_{2}(\gamma) \otimes H_{1}(\phi) \otimes T\).

In each of the aforementioned models, we can marginalize over \(\beta(\mathbf{s}, t)\) as we did in the earlier sections. Depending on the model, it may be more computationally convenient to block the data by site or by time. We omit the details, noting only that with \(n\) sites and \(T\) time points, the resulting likelihood will involve the determinant and inverse of an \(n T \times n T\) matrix.

Note that all of the foregoing modeling can be applied to the case of cross-sectional data in which the set of observed locations varies with \(t\). This is the case with, for instance, our real estate data. We observe a selling price only at the time of a transaction. With \(n_{t}\) locations in year \(t\), the likelihood for all but model 3 will involve a \(\sum n_{t} \times \sum n_{t}\) matrix.

\section*{5. MODEL COMPARISON}

In either the purely spatial case or the spatio-temporal case, we may seek to compare models. For instance, in the spatial case from Sections 2 and 3, we can consider models in which only a subset of the coefficients vary spatially, where the coefficient processes are independent, or where the multivariate dependence version is assumed. In the spatio-temporal case, we can consider models 1-4 (with possible submodels in the cases of models 2 and 3).

We use the posterior predictive loss approach of Gelfand and Ghosh (1998). Illustrating in the spatio-temporal context, for each ( \(\mathbf{s}, t\) ) in (19), let \(Y_{\text {new }}(\mathbf{s}, t)\) be a new observation obtained at that location and time for the given \(\mathbf{X}(\mathbf{s}, t)\). If under a particular model, \(\mu(\mathbf{s}, t)\) and \(\sigma^{2}(\mathbf{s}, t)\) denote the mean and variance

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/fc6b12fd-a585-4486-9ab4-cef281bdb867-07.jpg?height=879&width=767&top_left_y=202&top_left_x=190}
\captionsetup{labelformat=empty}
\caption{Figure 1. Locations Sampled Within the Parish of Baton Rouge for the Static Spatial Models.}
\end{figure}
of the predictive distribution of \(Y_{\text {new }}(\mathbf{s}, t)\) given \(Y_{\text {obs }}\) (the set of all observed \(Y\) 's), then the criterion becomes
\[
D_{k}=\sum_{(\mathbf{s}, t)}\left(Y_{o b s}(\mathbf{s}, t)-\mu(\mathbf{s}, t)\right)^{2}+k \sum_{(\mathbf{s}, t)} \sigma^{2}(\mathbf{s}, t) .
\]

In (20), the first term is a goodness-of-fit component ( \(G\) ), and the second term is a penalty for model complexity component \((P)\). The constant \(k\) weights these components and is often set to 1 . The model yielding the smallest value of (20) is chosen.

\section*{6. AN EXAMPLE}

We draw our data from a database of real estate transactions in Baton Rouge, LA, during the 8 -year period 19851992. In particular, we focus on modeling the log selling price of single family homes. In the literature it is customary to work with log selling price in order to achieve better approximate normality. A range of house characteristics are available. We use four of the most common choices: age of house, square feet of living area, square feet of other area (e.g., garages, carports, storage) and number of bathrooms. For the static spatial case, a sample of 237 transactions was drawn from 1992. Figure 1 shows the parish of Baton Rouge and the locations contained in an encompassing rectangle within the parish.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1. Values of Posterior Predictive Model Choice Criterion for Two-Dimensional Models (intercept process is always included)}
\begin{tabular}{lccc}
\hline \hline Model & \(G\) & \(P\) & \(D\) \\
\hline Living area & 69.87 & 46.24 & 116.11 \\
Age & 74.52 & 44.58 & 119.10 \\
Other area & 70.24 & 49.87 & 120.11 \\
Bathrooms & 78.02 & 52.93 & 130.95 \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2. Values of Posterior Predictive Model Choice Criterion for Three-Dimensional Models (intercept process is always included)}
\begin{tabular}{lccc}
\hline \hline Model & \(G\) & \(P\) & \(D\) \\
\hline Age, living area & 61.38 & 47.83 & 109.21 \\
Age, other area & 63.80 & 48.45 & 112.25 \\
Age, bathrooms & 67.25 & 48.66 & 114.91 \\
Living area, other area & 72.35 & 50.75 & 123.10 \\
Living area, bathrooms & 78.47 & 49.28 & 127.75 \\
Other area, bathrooms & 74.58 & 47.44 & 122.02 \\
\hline
\end{tabular}
\end{table}

We fitted the following collection of models. In all cases the correlation function is from the Matern class, that is, \(\rho\left(\mathbf{s}-\mathbf{s}^{\prime} ; \phi\right) \propto\left(\gamma\left\|\mathbf{s}-\mathbf{s}^{\prime}\right\|\right)^{v} K_{v}\left(\gamma\left\|\mathbf{s}-\mathbf{s}^{\prime}\right\|\right)\). We used priors that are fairly noninformative and comparable across models as sensible. We began with a spatial-varying intercept and one spatially varying slope coefficient. (The remaining coefficients do not vary). The intercept and coefficient processes follow a two-dimensional version of (13). There are four such models. The model comparison results are given in Table 1. The model with a spatially varying living area coefficient is best here. We then introduce two spatially varying slope coefficient processes along with a spatially varying intercept using a three-dimensional version of (13). There are six models here. Table 2 shows that the model with spatially varying age and living area is best. Finally, we allow five spatially varying processes; an intercept and four coefficients. We tried a model with five independent processes along with a five-dimensional version of (13). From Table 3 the fivedimensional model using (13) is far superior, and the independence model is clearly worst, supporting our earlier intuition.

The prior specification used for the five-dimensional dependent process model are as follows. We take vague \(N\left(\mathbf{0}, 10^{5} I\right)\) for \(\mu_{\beta}\); a five-dimensional inverted wishart, \(\operatorname{IW}(5, \operatorname{diag}(.001))\), for \(T\), and inverted gamma, \(\operatorname{IG}(2,1)\), for \(\tau^{2}\) (mean 1 , infinite variance). For the Màtern correlation function parameters \(\phi\) and \(\nu\), we assume gamma priors \(G(2, .1)\) (with mean 20 and variance 200). For all of the models, three parallel chains were run to assess convergence. Satisfactory mixing was obtained within 3,000 iterations for all the models; another 2,000 samples were generated and retained for posterior inference. The resulting posterior inference summary is provided in Table 4. We note a significant negative overall age coefficient with significant positive overall coefficients for the other three covariates. These are as expected. The contribution to spatial variability from the components of \(\beta\) is captured through the diagonal elements of the \(T\) matrix scaled by the corresponding covariates after the discussion at the end of Section 2. We see that the spatial intercept process contributes most to the error

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3. Values of Posterior Predictive Model Choice Criterion (over all models)}
\begin{tabular}{lccc}
\hline \hline Model & \(G\) & \(P\) & \(D\) \\
\hline Five-dimensional model & 42.21 & 36.01 & 78.22 \\
Three-dimensional model (best) & 61.38 & 47.83 & 109.21 \\
Two-dimensional model (best) & 69.87 & 46.24 & 116.11 \\
Independent process model & 94.36 & 59.34 & 153.70 \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4. Inference Summary for the Five-Dimensional Multivariate Spatially Varying Coefficients Model}
\begin{tabular}{|l|l|l|l|}
\hline Parameter & 2.5\% & 50\% & 97.5\% \\
\hline \(\beta_{0}\) (intercept) & 9.908 & 9.917 & 9.928 \\
\hline \(\beta_{1}\) (age) & -. 008 & -. 005 & -. 002 \\
\hline \(\beta_{2}\) (living area) & . 283 & . 341 & . 401 \\
\hline \(\beta_{3}\) (other area) & . 133 & . 313 & . 497 \\
\hline \(\beta_{4}\) (bathrooms) & . 183 & . 292 & . 401 \\
\hline \(T_{11}\) & . 167 & . 322 & . 514 \\
\hline \(\bar{x}_{1}^{2} T_{22}\) & . 029 & . 046 & . 063 \\
\hline \(\bar{x}_{2}^{2} T_{33}\) & . 013 & . 028 & . 047 \\
\hline \(\bar{x}_{3}^{2} T_{44}\) & . 034 & . 045 & . 066 \\
\hline \(\bar{x}_{4}^{2} T_{55}\) & . 151 & . 183 & . 232 \\
\hline \(T_{12} / \sqrt{T_{11} T_{22}}\) & -. 219 & -. 203 & -. 184 \\
\hline \(T_{13} / \sqrt{T_{11} T_{33}}\) & -. 205 & -. 186 & -. 167 \\
\hline \(T_{14} / \sqrt{T_{11} T_{44}}\) & . 213 & . 234 & . 257 \\
\hline \(T_{15} / \sqrt{T_{11} T_{55}}\) & -. 647 & -. 583 & -. 534 \\
\hline \(T_{23} / \sqrt{T_{22} T_{33}}\) & -. 008 & . 011 & . 030 \\
\hline \(T_{24} / \sqrt{T_{22} T_{44}}\) & . 061 & . 077 & . 098 \\
\hline \(T_{25} / \sqrt{T_{22} T_{55}}\) & -. 013 & . 018 & . 054 \\
\hline \(T_{34} / \sqrt{T_{33} T_{44}}\) & -. 885 & -. 839 & -. 789 \\
\hline \(T_{35} / \sqrt{T_{33} T_{55}}\) & -. 614 & -. 560 & -. 507 \\
\hline \(T_{45} / \sqrt{T_{44} T_{55}}\) & . 173 & . 232 & . 301 \\
\hline \(\phi\) (decay parameter) & . 51 & 1.14 & 2.32 \\
\hline \(v\) (smoothness parameter) & . 91 & 1.47 & 2.87 \\
\hline range (in km) & 2.05 & 4.17 & 9.32 \\
\hline \(\tau^{2}\) & . 033 & . 049 & . 077 \\
\hline
\end{tabular}
\end{table}
variability with, perhaps surprisingly, the "bathrooms" process second. Clearly, spatial variability overwhelms the pure error variability \(\left(\tau^{2}\right)\), showing the importance of the spatial model. The dependence between the processes is evident in the posterior correlation between the components. In fact, under (13), it is straightforward to calculate that \(\operatorname{cov}\left(\tilde{\beta}_{l}(\mathbf{s}), \tilde{\beta}_{m}(\mathbf{s}+\mathbf{h})\right) / \sqrt{\operatorname{cov}\left(\tilde{\beta}_{l}(\mathbf{s}), \tilde{\beta}_{l}(\mathbf{s}+\mathbf{h})\right) \operatorname{cov}\left(\tilde{\beta}_{m}(\mathbf{s}), \tilde{\beta}_{m}(\mathbf{s}+\mathbf{h})\right)}=T_{l m} / \sqrt{T_{l l} T_{m m}}\), regardless of \(\mathbf{h}\). We find the anticipated negative association between the intercept process and the slope processes (apart from that with the "other area" process). Under the Matérn correlation function, by inverting \(\rho(\cdot ; \phi)=.05\) for a given value of the decay parameter \(\gamma\) and the smoothing parameter \(\nu\), we obtain the range, that is, the distance beyond which spatial association becomes negligible. Posterior samples of ( \(\gamma, \nu\) ) produce posterior samples for the range. The resulting posterior median is roughly 4 km over a somewhat sprawling parish, which measures roughly \(22 \mathrm{~km} \times 33 \mathrm{~km}\). The smoothness parameter suggests processes with mean squared differentiable realizations \((v>1)\). Figure 2 shows the posterior mean spatial surfaces for each of the processes. The contour plots are evidently quite different. Table 5 considers a sample of 20 holdout sites on which the model can be validated; in fact, the validation is done with all of the models from Table 3. The entries in bold indicate validation failures. We find that for the first model, 19 of 20 predictive intervals contain the true value while for the independence model, 18 of 20 do so. More importantly, notice how much shorter these intervals are using the five-dimensional model.

Turning to the dynamic models proposed in Section 5, we returned to the Baton Rouge database, drawing a sample of 120 transactions at distinct spatial locations for the years 19891992. We compare models 1-4. In particular, we have two versions of model 2; 2a has the \(\boldsymbol{\alpha}(t)\) as four iid time dummies, and 2 b uses the multivariate temporal process model for \(\boldsymbol{\alpha}(t)\). We also have two versions of model 3 ; 3a has a common \(\boldsymbol{\mu}_{\beta}\) across \(t\), whereas 3 b uses \(\boldsymbol{\mu}_{\beta}^{(t)}\). In all cases, we used the fivedimensional spatially varying coefficient model for \(\beta\) 's. Table 6 gives the results. Model 3, in which space is nested within time, turns out to be the best, with model 4 following closely behind. Finally Table 7 summarizes the posterior inference summary for model 3b. The overall coefficients ( \(\mu_{\beta}^{(t)}\) ) do not change much over time; however, there is some indication that spatial range does change over time.

\section*{7. THE GENERALIZED LINEAR MODEL SETTING}

We briefly consider a generalized linear model version of (12), replacing the Gaussian first stage with
\[
f\left(y\left(\mathbf{s}_{i}\right) \mid \theta\left(\mathbf{s}_{i}\right)\right)=h\left(y\left(\mathbf{s}_{i}\right)\right) \exp \left(\theta\left(\mathbf{s}_{i}\right) y\left(\mathbf{s}_{i}\right)-b\left(\theta\left(\mathbf{s}_{i}\right)\right)\right)
\]
where, using a canonical link, \(\theta\left(\mathbf{s}_{i}\right)=\mathbf{X}^{T}\left(\mathbf{s}_{i}\right) \tilde{\boldsymbol{\beta}}\left(\mathbf{s}_{i}\right)\). The specification generates the models of Diggle et al. (1998). In (19), we could include a dispersion parameter with little additional complication.

The resulting first-stage likelihood becomes
\[
L(\tilde{\boldsymbol{\beta}}: \mathbf{y})=\exp \left\{\sum y\left(\mathbf{s}_{i}\right) \mathbf{X}^{T}\left(\mathbf{s}_{i}\right) \tilde{\boldsymbol{\beta}}\left(\mathbf{s}_{i}\right)-b\left(\mathbf{X}^{T}\left(\mathbf{s}_{i}\right) \tilde{\boldsymbol{\beta}}\left(\mathbf{s}_{i}\right)\right)\right\} .
\]

Taking the prior on \(\tilde{\beta}\) in (14), the Bayesian model is completely specified with a prior on \(\phi, T\), and \(\mu_{\beta}\).

This model can be fitted using a conceptually straightforward Markov chain Monte Carlo algorithm in the form of a Gibbs sampler, which would update the components of \(\mu_{\beta}\) and \(\tilde{\beta}\) using adaptive rejection sampling (Gilks and Wild 1992). With an inverse Wishart prior on \(T\), the resulting full conditional of \(T\) is again inverse Wishart. Updating \(\phi\) is usually very awkward, because it enters in the Kronecker form in (14). Metropolis updates are hard to design but offer perhaps the best possibility. Also problematic is the repeated componentwise updating of \(\tilde{\boldsymbol{\beta}}\). This hierarchically centered parameterization (Gelfand, Sahu, and Carlin 1995, 1996) is preferable to working with \(\mu_{\beta}\) and \(\boldsymbol{\beta}\), but the algorithm still runs very slowly with autocorrelation problems.

An alternative to componentwise updating is to introduce blocked Metropolis updating through a Langevin diffusion (see, e.g., Roberts, Gelman, and Gilks 1997; Christensen, Möller, and Waagepetersen 2000). Updating the entire \(n p \times 1\) vector \(\tilde{\beta}\) seems unrealistic. Blocking \(\tilde{\beta}\) by component process is essentially computationally intractable, and blocking \(\tilde{\beta}\) by site creates strong autocorrelation problems. Clearly, more work is needed to provide efficient model-fitting strategies for these models.
[Received August 2001. Revised July 2002.]

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/fc6b12fd-a585-4486-9ab4-cef281bdb867-09.jpg?height=2348&width=1567&top_left_y=227&top_left_x=247}
\captionsetup{labelformat=empty}
\caption{Figure 2. Mean Posterior Spatial Surfaces for the 5-D SVC Model: (a) Intercept Process, (b) Age Process, (c) Living Area Process, (d) Other Area Process, (e) Bathrooms Process.}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 5. Cross-Validatory Posterior Predictive Intervals for 20 Locations in Baton Rouge}
\begin{tabular}{|l|l|l|l|l|}
\hline True value & Five-dimensional & Three-dimensional (best) & Two-dimensional (best) & Independent \\
\hline 11.18 & (10.85, 11.53) & (10.74, 11.68) & (10.56, 11.52) & (10.44, 11.79) \\
\hline 10.86 & (10.54, 11.22) & (10.47, 11.31) & (10.52, 11.35) & (10.71, 12.08) \\
\hline 11.16 & (11.05, 11.68) & (11.03, 11.96) & (10.97, 11.98) & (10.68, 12.06) \\
\hline 11.46 & (11.39, 11.78) & (11.28, 12.11) & (11.35, 12.26) & (11.51, 12.34) \\
\hline 12.29 & (12.01, 12.72) & (12.09, 12.94) & (12.06, 13.01) & (11.87, 13.24) \\
\hline 10.18 & (9.81, 10.59) & (9.80, 10.60) & (9.73, 10.88) & (9.68, 10.75) \\
\hline 11.26 & (11.03, 11.51) & (11.17, 11.76) & (11.08, 11.72) & (11.21, 12.13) \\
\hline 10.97 & (10.82, 11.42) & (11.00, 11.97) & (10.77, 11.35) & (10.71, 11.52) \\
\hline 11.55 & (11.27, 11.91) & (11.21, 11.86) & (11.28, 12.01) & (11.27, 12.08) \\
\hline 12.22 & (12.11, 12.76) & (12.20, 12.84) & (12.15, 13.05) & (12.29, 13.25) \\
\hline 12.05 & (11.81, 12.55) & (11.85, 12.53) & (11.85, 12.51) & (11.78, 12.49) \\
\hline 11.68 & (11.50, 11.89) & (11.47, 11.92) & (11.49, 11.90) & (11.53, 12.18) \\
\hline 11.95 & (11.60, 12.22) & (11.72, 12.24) & (11.64, 12.18) & (11.52, 12.21) \\
\hline 12.22 & (11.77, 12.35) & (11.62, 12.27) & (11.89, 12.84) & (11.74, 12.82) \\
\hline 11.34 & (11.02, 11.76) & (11.23, 11.86) & (11.22, 11.90) & (11.26, 11.85) \\
\hline 11.05 & (11.10, 11.75) & (10.96, 11.78) & (11.01, 11.81) & (11.12, 12.09) \\
\hline 10.53 & (10.04, 10.86) & (10.05, 10.89) & (10.07, 10.86) & (10.05, 10.88) \\
\hline 11.42 & (11.26, 11.64) & (11.18, 11.82) & (11.16, 11.88) & (11.18, 12.05) \\
\hline 10.68 & (10.58, 11.29) & (10.76, 11.53) & (10.73, 11.78) & (10.57, 11.58) \\
\hline 12.16 & (11.66, 12.19) & (11.85, 12.62) & (11.94, 12.89) & (11.78, 12.82) \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 6. Model Choice Criteria for Various Spatio-Temporal Process Models (see Sec. 6)}
\begin{tabular}{|l|l|l|l|}
\hline Model & G & \(P\) & D \\
\hline \multicolumn{4}{|l|}{Independent process} \\
\hline Model 1 & 88.58 & 56.15 & 144.73 \\
\hline Model 2a & 77.79 & 50.65 & 128.44 \\
\hline Model 2b & 74.68 & 50.38 & 125.06 \\
\hline Model 3a & 59.46 & 48.55 & 108.01 \\
\hline Model 3b & 57.09 & 48.41 & 105.50 \\
\hline Model 4 & 53.55 & 52.98 & 106.53 \\
\hline \multicolumn{4}{|l|}{Dependent process} \\
\hline Model 1 & 54.54 & 29.11 & 83.65 \\
\hline Model 2a & 47.92 & 26.95 & 74.87 \\
\hline Model 2b & 43.38 & 29.10 & 72.48 \\
\hline Model 3a & 43.74 & 20.63 & 64.37 \\
\hline Model 3b & 42.35 & 21.04 & 63.39 \\
\hline Model 4 & 37.84 & 26.47 & 64.31 \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 7. Inference Summary for the 5-Dimension Multivariate Spatio-Temporal Model 3 (space nested within time)}
\begin{tabular}{|l|l|l|l|l|}
\hline Parameter & 1989 & 1990 & 1991 & 1992 \\
\hline \(\beta_{0}\) (intercept) & 9.697 (9.438, 9.956) & 9.164 (8.782, 9.613) & 9.56 (9.237, 9.860) & 9.716 (9.445, 9.980) \\
\hline \(\beta_{1}\) (age) & -. 005 (-.007, -.001) & -. 004 (-.007, -.002) & -. 004 (-.007, -.002) & -. 005 (-.008, -.002) \\
\hline \(\beta_{2}\) (living area) & . 383 (.314, .490) & . 357 (.252, .460) & . 402 (.232, .574) & . 348 (.278, .416) \\
\hline \(\beta_{3}\) (other area) & . 246 (.068, .410) & . 248 (.075, .422) & . 305 (.137, .486) & . 326 (.146, .515) \\
\hline \(\beta_{4}\) (bathrooms) & . 301 (.195, .407) & . 296 (.193, .395) & . 311 (.197, .420) & . 304 (.201, .416) \\
\hline \(\phi\) (decay parameter) & 1.63 (.68, 3.52) & 1.77 (.85, 3.15) & 1.05 (.42, 2.33) & 1.17 (.57, 2.38) \\
\hline \(\nu\) (smoothness parameter) & 1.48 (.92, 3.05) & 1.51 (.93, 3.05) & 1.54 (.86, 3.16) & 1.47 (.88, 2.92) \\
\hline range (in km) & 2.95 (1.40, 7.04) & 2.80 (1.55, 5.65) & 4.55 (2.07, 11.31) & 4.11 (2.03, 8.4) \\
\hline \(\tau^{2}\) & \multicolumn{4}{|c|}{. 038 (.024, .082)} \\
\hline
\end{tabular}
\end{table}

\section*{REFERENCES}

Agarwal, D. K., and Gelfand, A. E. (2001), "Slice Gibbs Sampling for Simulation-Based Fitting of Spatial Models," technical report, University of Connecticut, Dept. of Statistics.
Agarwal, D. K., Gelfand, A. E., Sirmans, C. F., and Thibadeau, T. G. (2003), "Nonstationary Spatial House Price Models," Journal of Statistical Planning and Inference,
Assunçao, J. J, Gamerman, D., and Assunçao, R. M. (1999), "Regional Differences in Factor Productivities of Brazilian Agriculture: A Space-Varying Parameter Approach," technical report, Universidade Federal do Rio de Janeiro, Statistical Laboratory.
Banerjee, S., and Gelfand, A. E. (2002), "Prediction, Interpolation and Regression for Spatially Misaligned Data," Sankhya, Ser. A, 64, 227-245.
- (2003), "On Smoothness Properties of Spatial Processes," Journal of Multivariate Analysis, to appear.
Cressie, N. A. C. (1993), Statistics for Spatial Data (2nd ed.), New York: Wiley.
Christensen, O. F., Moller, J., and Waagepetersen, R. (2000), "Analysis of Spatial Data Using Generalized Linear Mixed Models and Langevin-Type Markov Chain Monte Carlo," Research Report R-00-2009, Aalborg University.
Diggle, P. J., Tawn, J. A., and Moyeed, R. A. (1998), "Model-Based Geostatistics" (with discussion), Applied Statistics, 47, 299-350.
Gamerman, D., Moreira, A. R. B., and Rue, H. (2002), "Space-Varying Regression Models," technical report, Universidade Federal do Rio de Janeiro, Statistical Laboratory.
Gelfand, A. E., Ecker, M. D., Knight, J. R., and Sirmans, C. F. (2003), "The Dynamics of Location In Home Price," Journal of Real Estate Finance and Economics, to appear.
Gelfand, A. E., and Ghosh, S. K. (1998), "Model Choice: A Minimum Posterior Predictive Loss Approach," Biometrika, 85, 1-11.

Gelfand, A. E., Sahu, S. K., and Carlin, B. P. (1995), "Efficient Parametrization for Normal Linear Mixed Effects Models," Biometrika, 82, 479-488.
(1996), "Efficient Parametrization for Generalized Linear Mixed Models," in Bayesian Statistics 5, eds. J. Bernardo et al., Oxford, U.K.: Clarendon Press, pp. 165-180.
Gelfand, A. E., Zhu, L., and Carlin, B. P. (2001), "On the Change of Support Problem for Spatio-Temporal Data," Biostatistics, 2, 31-45.
Gilks, W. R., and Wild, P. (1992), "Adaptive Rejection Sampling for Gibbs Sampling," Journal of the Royal Statistical Society C, 41, 337-348.
Goldstein, H. (1995), Multilevel Statistical Models (2nd ed.), London: Arnold.
Kent, J. T. (1989), "Continuity Properties for Random Fields," Annals of Probability, 17, 1432-1440.
Luo, Z., and Wahba, G. (1998), "Spatio-Temporal Analogues of Temperature Using Smoothing Spline ANOVA," Journal of Climatology, 11, 18-28.
Mardia, K. V., and Goodall, C. (1993), "Spatio-Temporal Analyses of Multivariate Environmental Monitoring Data," in Multivariate Environmental Statistics, eds. G. P. Patil and C. R. Rao, Amsterdam: Elsevier, 347-386.
Neal, R. (2002), "Slice Sampling," The Annals of Statistics, December issue.
Raudenbush, S. W., and Bryk, A. S. (2002), Hierarchical Linear Models (2nd ed.), Newbury Park, CA: Sage.
Roberts, G., Gelman, A., and Gilks, W. (1997), "Weak Convergence and Optimal Scaling of Random Walk Metropolis Algorithms," Annals of Applied Probability, 7, 110-120.
Stein, M. L. (1999a), Interpolation of Spatial Data: Some Theory for Kriging, New York: Springer.
- (1999b), "Predicting Random Fields With Increasingly Dense Observations," Annals of Applied Probability, 9, 242-273.
Waller, L., Carlin, B. P., Xia, H., and Gelfand, A. E. (1997), "Hierarchical Spatio-Temporal Mapping of Disease Rates," Journal of the American Statistical Association, 92, 607-617.