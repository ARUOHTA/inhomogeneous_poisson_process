\title{
Preferential sampling and Bayesian geostatistics: Statistical modeling and examples
}

\author{
Lorenzo Cecconi, \({ }^{1}\) Laura Grisotto, \({ }^{1}\) Dolores Catelan, \({ }^{1,2}\) Corrado Lagazio, \({ }^{3}\) Veronica Berrocal \({ }^{4}\) and Annibale Biggeri \({ }^{\mathbf{1} \boldsymbol{,} \mathbf{2}}\)
}

\begin{abstract}
Preferential sampling refers to any situation in which the spatial process and the sampling locations are not stochastically independent. In this paper, we present two examples of geostatistical analysis in which the usual assumption of stochastic independence between the point process and the measurement process is violated. To account for preferential sampling, we specify a flexible and general Bayesian geostatistical model that includes a shared spatial random component. We apply the proposed model to two different case studies that allow us to highlight three different modeling and inferential aspects of geostatistical modeling under preferential sampling: (1) continuous or finite spatial sampling frame; (2) underlying causal model and relevant covariates; and (3) inferential goals related to mean prediction surface or prediction uncertainty.
\end{abstract}

\section*{Keywords}

Preferential sampling, Bayesian geostatistics, shared component model

\section*{I Introduction}

The purpose in geostatistical analysis is to estimate the continuous spatial distribution of a georeferenced variable. \({ }^{1,2}\) Typically, inference is based on data, e.g. a set of measurements taken at a finite number of locations. Commonly, the choice of the sample locations is considered part of the study design and the analysis is performed conditionally on these sampling constraints. The simplest spatial design encountered in geostatistics consists of sample locations uniformly

\footnotetext{
\({ }^{1}\) Department of Statistics, Computer Science, Applications "G. Parenti", University of Florence, Italy
\({ }^{2}\) Biostatistics Unit, Institute for Cancer Prevention and Research, Tuscany Region, Italy
\({ }^{3}\) Department of Economics, University of Genoa, Italy
\({ }^{4}\) Department of Biostatistics, University of Michigan, USA
Corresponding author:
Lorenzo Cecconi, Department of Statistics, Computer Science, Applications "G. Parenti", University of Florence, Viale Morgagni, 59, 50134 Florence, Italy.
Email: cecconi@disia.unifi.it
}
distributed over the observed region-as in a regular lattice design, see Diggle and Lophaven \({ }^{3}\) for a review of sampling designs in a Bayesian framework. However, there are situations in which the researcher cannot choose the sample locations by design and must rely on available locations assuming that the sampling locations are stochastically independent of the sampling process. Although often used, in many applications, such assumption may be unrealistic. Preferential sampling refers to situations where there is a stochastic relationship between data and locations. \({ }^{4}\) Most specifically, preferential sampling refers to any situation in which the spatial process and the sampling locations are not stochastically independent. We note that the definition of preferential sampling involves a stochastic dependence, as opposed to a functional dependence: in particular, any dependence on a common set of explanatory variables does not constitute preferential sampling in a strict sense. Lee et al. \({ }^{5}\) expanded the definition of preferential sampling by including dependence between spatial process and sample locations through common covariates. In any case, this distinction between functional and stochastic dependence has no consequence as far as we include the informative covariates in the model. As stated in Pati et al. \({ }^{6}\) "... accounting for informative sampling is only necessary when there is an association between the spatial surface of interest and the sampling density that cannot be explained by the [common] spatial covariates."

In this paper, we present two examples of geostatistical analysis in which the usual assumption of stochastic independence between the point process and the measurement process is violated.

Diggle et al., \({ }^{4}\) Pati et al. \({ }^{6}\) and Gelfand et al. \({ }^{7}\) showed that if preferential sampling is ignored, geostatistical inferences and predictions could be biased. Working in an air pollution context, Lee et al. \({ }^{5}\) illustrated that preferential sampling can also have a large impact on the validity and reliability of the health effect estimate.

To account for preferential sampling, Diggle et al. \({ }^{4}\) proposed a class of models that postulate a shared spatial latent process. For parameter estimation, the authors suggested using maximum likelihood estimation with estimates obtained via Monte Carlo methods. A full Bayesian formulation was instead proposed by Pati et al. \({ }^{6}\) Briefly, they proposed a joint model for the point process that accounts for the sampling locations and for the point-referenced spatial process that describes the spatial intensity of the georeferenced variable. The log intensity of the spatial point process of the sampling locations is included as a term in the model for the pointreferenced spatial process. Assuming a finite population of sampling locations, Shaddick and Zidek \({ }^{8}\) extended the model to the spatio-temporal context and took advantage of variations of the HorwitzThompson estimator to provide unbiased parameter estimates when the data are preferentially sampled. A similar idea was pursued by Cecconi et al. \({ }^{9}\) who, working on a finite locations scenario, proposed a two-step modeling strategy: first, a model for the sampling locations probabilities; then, a weighted Bayesian geostatistical model on the georeferenced measurements. In Cecconi et al., \({ }^{9}\) the adjustment for preferential sampling is carried out by introducing weights in the geostatistical model with weights given by the sampling locations probabilities. With a different inferential goal, Grisotto et al. \({ }^{10}\) developed a Bayesian geostatistical model to account for preferential sampling when the main interest is in statistical integration of different sources of information and uncertainty assessment. In particular, in Grisotto et al., \({ }^{10}\) the integration of information refers to the integration of information arising from deterministic model outputsi.e. mathematical and chemico-physical models that represent the spatial process in consideration-adjusting for what Lee et al. \({ }^{5}\) called "preferential sampling through informative covariates," while the stochastic component of the model is shown to have an effect on prediction uncertainty.

Many examples of preferential sampling presented in the literature are in the context of air pollution exposure assessment. In fact, air quality monitors are seldom randomly distributed in
space and are often more commonly located where the level of pollution is higher. \({ }^{8,10,11}\) Other examples of preferential sampling arise in veterinary parasitology where the location of the sampled farms \({ }^{12}\) is not independent of the spatial distribution of parasites. Related preferential sampling problems can also be seen when considering residence of people \({ }^{13}\) or community surveys. \({ }^{14}\) In the first set of examples of preferential sampling, the bias in statistical inference arises as a consequence of the dependence between the spatial process and the probability that a given location be sampled whereas in the latter examples on residence of people or community surveys, the bias in the inference may be a function of the location itself or a function of the characteristics of the subject resident at that location.

In this paper, we present two examples of geostatistical analysis in which the usual assumption of stochastic independence between the point process of the sampling locations and the measurement process is violated. To handle the preferential sampling, we specify and model both processes in a unified Bayesian model formulation that includes a shared spatial random component, resulting in a flexible and general model specification that allows us to address also continuous and finite population locations over space. We present applications of the model to two case studies, one rooted in veterinary parasitology and involving a discrete specification of the point process, and one in exposure assessment with a continuous point process for the sampling locations. In the first case study, the data refer to prevalence of parasitic infection in the Campania region (Italy) and the interest is in predicting probabilities of infection for the entire Campania region. Only a finite number of sampling points is available and of this only a subset can be sampled. In this context, preferential sampling arises because of the farmers' selective reporting. In the second example, \(\mathrm{PM}_{10}\) data from the air quality network of the Environmental Protection Agency of Lombardy region (Italy) are considered. In this application, the sampling points can be any point in the region and preferential sampling arises as a consequence of the policy strategies in locating monitoring stations. Last but not least, we use these two case studies to highlight the effect of preferential sampling on two distinct goals of a geostatistical analysis: namely, prediction of the mean spatial process, our main interest in the first application, and prediction uncertainty, the focus of the second application. A simulation study completes the paper, showing the validity and flexibility of our Bayesian model formulation.

\section*{2 Methods}

Following Diggle et al., \({ }^{2}\) we claim that a geostatistical model is specified by three different processes:
- a field process \(S\)
- a measurement process \(Y\)
- a point process \(P\)

Let \(\{S(x) \mid x \in D\}\) be a field process in the domain \(D \subset \Re^{2}\) and \(y_{i}=Y\left(x_{i}\right)\) denote the measurement made at locations \(x_{i}\), for \(i=1, \ldots, n\), with the sample locations \(x_{i}\) realizations of the point process \(P\).

The measurement process \(Y\) is obviously dependent on the field process \(S\). Under preferential sampling, also the point process \(P\) depends on the field process \(S\).

\subsection*{2.1 General shared model for preferential sampling}

To account for preferential sampling, we propose a general spatial-shared component model as in Held et al. \({ }^{15}\)

Specifically, we model the sampling process as an inhomogeneous point process whose intensity function \(\lambda(x)\) is expressed as a function of the covariates, a spatial-specific component and a spatialshared component. Under a fixed design or a complete random design, the point process of the sampling locations is typically a homogenous Poisson process with constant intensity \(\lambda\). On the contrary, if preferential sampling occurred, the sampling design is dependent on the field process \(S\) and an inhomogeneous point process must be specified.

In turn, the field process \(S\) is modeled as a spatial process with a mean expressed as a function of covariates-the same included in the model for the intensity \(\lambda(x)\)-plus a spatial-specific component and a spatial-shared component.

Finally, the measurement process is modeled to be a noisy version of the field process \(S\). In summary, the model we propose is as follows
\[
\begin{aligned}
X(x) & \sim P P(\lambda(x)) \\
\log (\lambda(x)) & =\alpha^{\prime}+\beta^{\prime} Z(x)+\eta^{\prime}(x)+\delta^{-1} \xi(x) \\
Y(x) & \sim N\left(S(x), \tau^{2}\right) \\
S(x) & \sim \operatorname{GP}\left[\mu(x), \sigma^{2}, \rho(d)\right] \\
\mu(x) & =\alpha^{\prime \prime}+\beta^{\prime \prime} Z(x)+\eta^{\prime \prime}(x)+\delta \xi(x)
\end{aligned}
\]
where \(P P\) is a point process with spatially varying intensity function \(\lambda(x), Y(x)\) is the measurement process, and \(S(x)\) is a Gaussian process (GP) with mean \(\mu(x)\) and covariance matrix induced by the correlation function \(\rho(d)\) that depends on the distance \(d\) between pairs of points and by the spatial variance parameter \(\sigma^{2}\).

The measurement errors, that in geostatistics are typically assumed to be addictive, possibly on a transformed scale, are assumed to be mean-zero random variables with variance \(\tau^{2}\), also called the nugget variance.

The \(Z(x)\) term in equation (1) denotes the available covariates which can be relevant in explaining the mean surface of the field process \(S(x)\) of interest and/or the spatial intensity of the point process.

To account for the dependence between the sampling process \(P\) and the field process \(S\), we introduce in our model an underlying latent spatial process \(\xi(x)\) which induces a correlation between the two processes. In our examples, this shared latent component is specified to be the sum of spatially non-structured random terms \(u(x)\) and spatially structured random terms \(v(x)\), following the suggestion of Besag et al. \({ }^{16}\)

The importance of the shared component in the two processes is controlled by the tuning parameter \(\delta>0\). If in equation (1), the specific component \(\eta^{\prime}(x)=0\), the model can be interpreted as a geostatistical model with an unmeasured covariate for which a surrogate (the location intensity \(\lambda(x))\) is available. \({ }^{6}\)

Here, for matter of generality, we specify a model that includes both the shared spatial process \(\xi(x)\) and the specific spatial random components \(\eta^{\prime}(x), \eta^{\prime \prime}(x)\). However, we observe that the inclusion of the specific spatial components in equation (1) causes problems with identifiability for all the latent spatial fields.

\subsection*{2.2 A simple model specification}

We now show a simplification of the general model introduced in Section 2.1. Let us assume that the spatial process \(S\) is a GP with an appropriate Matérn covariance function.

For the point process, \(P\) of the sampling location, we can distinguish two different possibilities, which we present and discuss more in detail in the two case studies in Section 3:
- \(P\) is an inhomogeneous point process whose intensity function depends on location \(x\) and whose "events," the sampling points, can occur at any point in the region;
- \(P\) is a discrete point process with a finite set of eligible locations where "events" can occur, e.g. with a finite set of eligible sampling points. This means that the number of sampling points is finite and a fixed, known list of locations are available as potential sampling points.

All known covariates which can at least partially explain the common spatial pattern in the sampling points process \(P\) and the field process \(S\) are included as predictors in both modeling equations. Alternatively, if they do not influence the point process but are relevant for the mean surface of interest of the field process \(S\), they should be specified in the equation for the latter only. Finally, if relevant covariates are not available for both spatial processes, then the modeling equations in equation (1) would include only a specific spatial random component and a shared component.

For matter of simplicity, on the domain of interest, we superimpose a fine grid with a degree of fineness that depends on both the spatial distribution of the relevant covariates and the nature of the domain. We employ the same grid for the point process of the sampling locations and the spatial process.

Let \(x_{j}\) be the coordinates of the \(j\) th grid cell centroid and \(X\left(x_{j}\right)\) be the number of observations located in the \(j\) th cell \((j=1, \ldots, J)\). We take the counts \(X\left(x_{j}\right)\) as a realization of the inhomogeneous point process for the sampling locations, and we use them to infer upon the spatially varying intensity function \(\lambda(x)\) over the grid covering the domain, e.g. we estimate \(\lambda\left(x_{j}\right)\), for \(j=1, \ldots, J\).

We employ the same grid to predict the surface of interest \(S(x)\).
Specifically, let \(Y\left(x_{j(i)}\right)\) be the observed value of the spatial process at the \(i\) th sample location \(x_{j(i)}\) for \(i=1, \ldots, n\), falling within the \(j\) th grid cell, for some \(j \in\{1, \ldots, J\}\). We assume \(S\left(x_{j(i)}\right)\) to follow a multivariate Gaussian distribution with a mean vector that depends on covariates \(\mathrm{Z}\left(x_{j(i)}\right)\) and a covariance matrix induced by the exponential correlation function which postulates an exponential decay in the correlation between any two points. We choose the parameters of the exponential covariance function so that the correlation between two sites located at the maximum observed inter-site distance is almost null, while the correlation between two sites located at the minimum observed inter-site distance is very high and close to \(1 .^{17}\)

In summary, given the discrete sampling locations, the model in equation (1) becomes
\[
\begin{aligned}
& X\left(x_{j}\right) \sim P P\left(\lambda\left(x_{j}\right)\right) \\
& \log \left(\lambda\left(x_{j}\right)\right)=\alpha^{\prime}+\beta^{\prime} Z\left(x_{j}\right)+\eta^{\prime}\left(x_{j}\right)+\delta^{-1}\left(u\left(x_{j}\right)+v\left(x_{j}\right)\right) \\
& Y\left(x_{j(i)}\right) \sim N\left(S\left(x_{j(i)}\right), \tau_{Y}\right) \\
& S\left(x_{j(i)}\right) \sim \operatorname{SpatialExp}\left[\mu\left(x_{j(i)}\right) ; \psi\left(x_{j(i)}\right)=(f(\varphi, \gamma) ; \sigma)\right] \\
& \mu\left(x_{j(i)}\right)=\alpha^{\prime \prime}+\beta^{\prime \prime} Z\left(x_{j(i)}\right)+\eta^{\prime \prime}\left(x_{j(i)}\right)+\delta\left(u\left(x_{j(i)}\right)+v\left(x_{j(i)}\right)\right) \\
& \alpha^{\prime}, \alpha^{\prime \prime} \sim \operatorname{flat}() \\
& \beta^{\prime}, \beta^{\prime \prime} \sim \operatorname{Normal}\left(0, \tau_{\beta}\right) \\
& u\left(x_{j(i)}\right) \sim \operatorname{Normal}\left(0, \tau_{u}\right)
\end{aligned}
\]
\[
\begin{aligned}
& v\left(x_{j(i)}\right) \sim \operatorname{CAR}\left(\bar{v}_{j \in i}, \tau_{v}\right) \\
& \log (\delta) \sim \operatorname{TN}\left(0, \tau_{\delta}\right) \\
& (\varphi, \gamma) \sim \text { informative priors }
\end{aligned}
\]
where SpatialExp is the notation we use to indicate a GP with an exponential covariance function, and \(u\left(x_{j(i)}\right)\) and \(v\left(x_{j(i)}\right)\) indicate, respectively, the heterogeneity and clustering random terms as in the Besag et al. \({ }^{16}\) model. The clustering random terms are provided with a Conditional autoregressive (CAR) prior, i.e. with a normal improper conditional autoregressive distribution with a binary adjacency matrix. Finally, as prior distribution on the model parameters we use an improper uniform on \(\alpha^{\prime}, \alpha^{\prime \prime}\) and a distribution symmetric around zero for the tuning term \(\delta\). To allow the shared component to vary by a constant factor, we choose to specify a lognormal distribution on \(\delta\), allowing any value to be equally likely a priori as its reciprocal.

Predictive distributions (and their summaries, e.g. mean, mode and standard deviation) of the measurement process \(Y\) at the set of locations \(x_{j}(j=1, \ldots, J)\), e.g. the centroids of the grid, are derived using standard Bayesian geostatistical formulations; more specifically using the equation
\[
[\tilde{\mathbf{y}} \mid \mathbf{y} ; Z, v, v ; \tilde{Z}, \tilde{v}, \tilde{v}]=\int[\tilde{\mathbf{y}} \mid \mathbf{y}, \Omega ; Z, v, v ; \tilde{Z}, \tilde{v}, \tilde{v}][\Omega \mid \mathbf{y} ; Z, v, v] \mathrm{d} \Omega
\]
where \(\Omega\) is the set of parameters in the mean and covariance function \(\left(\alpha^{\prime} \beta^{\prime \prime} \delta ; \phi \gamma \sigma\right), \tilde{Z}, \tilde{v}, \tilde{v}\) are the covariate and shared spatial component at the prediction points, \({ }^{17}\) and the posterior distribution of \(\Omega\) derives from the preferential sampling model in equation (2).

\section*{3 Case studies}

\subsection*{3.1 Preferential sampling in veterinary parasitological surveillance: Sheep farms in the Campania region}

In the first case study, we have information on parasitic infection in 89 farms, sampled from a total population of 8794 farms (Figure 1). The data come from two different programs on parasitological surveillance conducted in the years 2013 and 2014 in the Campania region (Southern Italy). Here, preferential sampling is due to the farms' selective reporting (see Cecconi et al. \({ }^{9}\) for details). The outcome is the presence or absence of a parasitic infection-the liver fluke Fasciola hepatica-among the farms and the aim is prediction of infection probability over the Campania region, taking into account preferential sampling.

\subsection*{3.1.1 Model}

To tailor the model in equation (2) to this specific application, we modified the model as follows. First, we recorded and georeferenced all members of the farm population in the Campania region, overimposing a \(10 \times 10 \mathrm{~km}\) grid on the region for a total of 184 cells. In each cell, we counted the population denominator, e.g. the total number of farms in the Campania region, and we denote with \(n\left(x_{j}\right)\), the total number of farms belonging to the \(j\) th cell with centroid \(x_{j}(j=1, \ldots, 184)\). Conversely, we indicate with \(X\left(x_{j}\right)\), the number of sampled farm in the \(j\) th cell. Given the type of data available, we specify a binomial likelihood for the point process at the grid cell. In this example, we omit to consider any covariate, and adjustment for preferential sampling due the stochastic dependence

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-07.jpg?height=630&width=1329&top_left_y=230&top_left_x=244}
\captionsetup{labelformat=empty}
\caption{Figure I. Spatial distribution of the 89 analyzed farms in our surveillance dataset (left) and of the total population of farms in the Campania region in years 2013-2014 (right) with the overimposed \(10 \times 10 \mathrm{~km}\) grid.}
\end{figure}
between the point process for the farm location and the spatial process for the infection probability is achieved through the introduction only of the shared latent process.

Thus, the model we use for this first application is
\[
\begin{aligned}
& X\left(x_{j}\right) \sim \operatorname{Binomial}\left(n\left(x_{j}\right), p\left(x_{j}\right)\right) \\
& \operatorname{logit}\left(p\left(x_{j}\right)\right)=\alpha^{\prime}+\eta^{\prime}\left(x_{j}\right)+\delta^{-1}\left(u\left(x_{j}\right)+v\left(x_{j}\right)\right) \\
& Y\left(x_{j(i)}\right) \sim N\left(\tilde{S}\left(x_{j(i)}\right), \tau_{Y}\right) \\
& S\left(x_{j(i)}\right) \sim \operatorname{SpatialExp}\left[\mu\left(x_{j(i)}\right) ; \psi\left(x_{j(i)}\right)=(f(\varphi, \gamma) ; \sigma)\right] \\
& \mu\left(x_{j(i)}\right)=\alpha^{\prime \prime}+\eta^{\prime \prime}\left(x_{j(i)}\right)+\delta\left(u\left(x_{j(i)}\right)+v\left(x_{j(i)}\right)\right) \\
& \\
& \alpha^{\prime}, \alpha^{\prime \prime} \sim \operatorname{flat}() \\
& u\left(x_{j(i)}\right) \sim \operatorname{Normal}\left(0, \tau_{u}\right) \\
& v\left(x_{j(i)}\right) \sim \operatorname{CAR}\left(\bar{v}_{j \in i}, \tau_{v}\right) \\
& \log (\delta) \sim N\left(0, \tau_{\delta}\right) \\
& (\varphi, \gamma) \sim \text { informative priors }
\end{aligned}
\]

WinBugs code for Case Study 1 is reported in the appendix.

\subsection*{3.1.2 Results}

In Figure 1, we report the geographical distribution of the 89 farms in our dataset and of the entire population of farms in the Campania Region in years 2013-2014. The 89 analyzed farms are clearly not a representative sample and they are not homogenously spatially distributed.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-08.jpg?height=1050&width=1149&top_left_y=230&top_left_x=378}
\captionsetup{labelformat=empty}
\caption{Figure 2. Map of posterior probability that a grid cell is in excess with respect to the observed mean prevalence of F. hepatica, in Campania region. Years 2013-2014.}
\end{figure}

The observed prevalence of \(F\). hepatica is low ( \(7.9 \%, 90 \%\) CI \(3.7 ; 14.3 \%\) ) consistently with the climatic characteristics of the region. The map of posterior probability that a grid cell is in excess with respect to the observed mean prevalence is presented in Figure 2. The figure highlights some areas in the North-West and marshland areas within the Southern Province of Salerno with higher probabilities.

The posterior mean of the shared component \(\xi\left(x_{j}\right)\) (Figure 3) is not "flat" but exhibits a strong spatial structure suggesting the presence of selective farms reporting. It is also clear that the association with the map of prevalence of infection.

\subsection*{3.2 Preferential sampling in environmental assessment: uncertainty in PM \(_{10}\) predicted concentration in Lombardy}

The second application that we present is in the context of air pollution exposure assessment and here the main interest is in the statistical integration of different sources of information and prediction uncertainty quantification. The outcome is the annual concentrations of \(\mathrm{PM}_{10}\) in the Lombardy region (Italy). As in the previous example, we overimposed a \(4 \times 4 \mathrm{~km}\) grid on the region for a total of 1679

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-09.jpg?height=1059&width=1163&top_left_y=230&top_left_x=327}
\captionsetup{labelformat=empty}
\caption{Figure 3. Spatial distribution of the shared random term \(\xi\left(x_{j}\right)=u\left(x_{j}\right)+v\left(x_{j}\right)\). F. hepatica, Campania Region. Years 2013-2014.}
\end{figure}
cell. In year 2007, the outcome was observed in 58 air quality monitoring stations (shown in Figure 4), part of the regional air quality monitoring network of the Regional Environmental Protection Agency (ARPA; Milan, Italy). Predictions of \(\mathrm{PM}_{10}\) annual average concentrations were also available from a Eulerian photochemical model, \({ }^{18}\) which was run on the same \(4 \times 4 \mathrm{~km}\) grid covering the Lombardy region.

\subsection*{3.2.I Model}

Differently from the previous application, for the point process of the monitor location, we specify a Negative Binomial point process with the monitor spatial intensity modeled as a function of covariates and a shared spatial component.

Following Grisotto et al., \({ }^{10}\) we integrate the information from the deterministic model output in the geostatistical model for the spatial process \(S\) by including it as a covariate in the model. Additionally, the deterministic model output was also introduced as a covariate in the model for the monitor location intensity \(\lambda(x)\)-in the spirit of what Lee et al. \({ }^{5}\) called preferential sampling through informative covariates. Residual spatial dependence between the spatial process and the spatial location process is still present and is accounted for by the shared spatial component.

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-10.jpg?height=1086&width=1167&top_left_y=228&top_left_x=368}
\captionsetup{labelformat=empty}
\caption{Figure 4. Geographical location of the \(58 \mathrm{PM}_{10}\) monitoring stations reporting annual average \(\mathrm{PM}_{10}\) concentration in year 2007 and the \(4 \times 4 \mathrm{~km}\) grid overimposed on the Lombardy region.}
\end{figure}

More specifically, the model that we use for this application is
\[
\begin{aligned}
& X\left(x_{j}\right) \sim \operatorname{NegBin}\left(p\left(x_{j}\right), r\right) \\
& p\left(x_{j}\right)=\left(r /\left(r+\lambda\left(x_{j}\right)\right)\right) \\
& \log \left(\lambda\left(x_{j}\right)\right)=\alpha^{\prime}+\beta^{\prime} Z\left(x_{j}\right)+\eta^{\prime}\left(x_{j}\right)+\delta^{-1} v\left(x_{j}\right) \\
& Y\left(x_{j(i)}\right) \sim N\left(S\left(x_{j(i)}\right), \tau_{Y}\right) \\
& S\left(x_{j(i)}\right) \sim \operatorname{SpatialExp}\left[\mu\left(x_{j(i)}\right) ; \psi\left(x_{j(i)}\right)=(f(\varphi, \gamma) ; \sigma)\right] \\
& \mu\left(x_{j(i)}\right)=\alpha^{\prime \prime}+\beta^{\prime \prime} Z\left(x_{j(i)}\right)+\eta^{\prime \prime}\left(x_{j(i)}\right)+\delta v\left(x_{j(i)}\right) \\
& \alpha^{\prime}, \alpha^{\prime \prime} \sim \operatorname{flat}() \\
& \beta^{\prime}, \beta^{\prime \prime} \sim \operatorname{Normal}\left(0, \tau_{\beta}\right) \\
& r \sim \operatorname{Gamma}(0.001,0.001) \\
& v\left(x_{j(i)}\right) \sim \operatorname{CAR}\left(\bar{v}_{j \in i}, \tau_{v}\right) \\
& \log (\delta) \sim N\left(0, \tau_{\delta}\right)
\end{aligned}
\]

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-11.jpg?height=1088&width=1163&top_left_y=228&top_left_x=327}
\captionsetup{labelformat=empty}
\caption{Figure 5. Posterior predicted \(\mathrm{PM}_{10}\) concentration for the Lombardy region in year 2007 as yielded by the shared component geostatistical model.}
\end{figure}
where \(x_{j}\) is the coordinates of the \(j\) th grid cell centroid, \(X\left(x_{j}\right)\) is the number of air quality monitoring station located in the \(j\) th cell \((j \in\{1, \ldots, 1679\}), \lambda\left(x_{j}\right)\) is the continuous monitor spatial intensity, \(r\) is the dispersion parameter of the Negative Binomial model, \(Z\left(x_{j}\right)\) is the output from the deterministic model in grid cell \(j\), and \(v\left(x_{j}\right)\) is the shared spatial process, with tuning parameter \(\delta\). Additionally, \(Y\left(x_{j(i)}\right)\) denotes the \(\mathrm{PM}_{10}\) annual average concentration observed at the \(i\) th sample \(x_{j(i)}\), for \(i=1, \ldots, 58\), located in the \(j\) th grid cell, for some \(j \in\{1, \ldots, 1679\}\), which we assume depends on spatial covariates \(\mathrm{Z}\left(x_{j(i)}\right)\), specific \(\eta^{\prime \prime}\left(x_{j(i)}\right)\) and shared spatial random processes \(v\left(x_{j(i)}\right)\). WinBugs code for Case Study 2 is reported in the appendix.

\subsection*{3.2.2 Results}

As Figure 4 shows, the distribution of monitors is clearly not homogeneous over the region: more monitors are located in the North-West part of the domain where the city of Milan and other larger urban center are located. The annual average \(\mathrm{PM}_{10}\) concentration predicted by our geostatistical model using equation (3) is shown in Figure 5 (and the standard deviations in Figure 6). We observe that the spatial distribution of the predicted annual average \(\mathrm{PM}_{10}\) concentration is almost identical to that produced by the Eulerian model, thus suggesting that the integration of the deterministic model output with the information from the air quality monitoring network did not result in a

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-12.jpg?height=1038&width=1114&top_left_y=230&top_left_x=395}
\captionsetup{labelformat=empty}
\caption{Figure 6. Standard deviations of posterior predicted \(\mathrm{PM}_{10}\) concentration for the Lombardy region in year 2007 as yielded by the shared component geostatistical model.}
\end{figure}
predicted concentration surface which is different from the original one yielded by the Eulerian model. However, differently from the deterministic numerical model, our Bayesian hierarchical model formulation allows us to obtain a measure of the prediction uncertainty through an appropriate propagation of uncertainty, which is not possible in the deterministic model.

Figure 7 shows a map of the posterior mean of the shared spatially structured component, which corresponds to the residual spatial variability not explained by the covariates. In this example, our main interest is in showing the effect of preferential sampling on the predictions' standard deviations. A map of the ratio between these standard deviations obtained when accounting vs. not accounting for preferential sampling is shown in Figure 8. As the figure illustrates, we estimate a consistently greater standard deviation (ratio greater than 1 ) in the areas not properly covered by the air quality network when accounting for preferential sampling.

\subsection*{3.3 Simulation study}

To evaluate the proposed modeling approach, we ran a simulation study. We built a realistic frame using the same design of our case study 2 on air pollutant monitors. Specifically, we defined a grid of 137 cells covering the spatial domain, while the spatial process was observed at a fixed set of 58 locations. To assess how the strength of preferential sampling affects inference, we considered three

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-13.jpg?height=1046&width=1115&top_left_y=228&top_left_x=350}
\captionsetup{labelformat=empty}
\caption{Figure 7. Spatial distribution of the posterior mean of the shared random component \(\nu\left(\mathrm{x}_{\mathrm{j}}\right)\). Lombardy region. Year 2007.}
\end{figure}
different preferential sampling scenarios: high, low and no preferential sampling. The degree of preferential sampling was varied by assigning different values to the tuning parameter \(\delta\) that appears in the equation for the spatial intensity of the point process: \(\log \left(\lambda_{i}\right)=\alpha+(1 / \delta) \times b_{i}\). Specifically, we used values of \(\delta\) equal to \(0.4,1.0\) and a very large value \((10,000)\). Indeed, no preferential sampling is obtained when \(\delta \rightarrow \infty\).

To determine the locations of the monitors in a preferential sampling setting, we proceeded as follows. First, we selected grid cells with probabilities proportional to the spatial intensity at the grid cells, and then we sampled the exact locations of the sampled monitors uniformly within each cell. For each sampled monitor, we then generated a \(\log \left(\mathrm{PM}_{10}\right)\) value sampling from a normal density with a standard deviation consistent with that of the real data in case study 2 and with a mean proportional to \(b_{i}\), where the spatial components \(b_{i}\) were assumed to be known (see Figure 9).

For each preferential sampling scenario, we simulated 100 datasets, to which we fit our shared model for preferential sampling using an MCMC algorithm that we ran for 100,000 iterations, with a burn-in of 10,000 iterations. As our simulations are designed to reproduce the realistic situation of an air quality monitoring network where the monitors' locations are strongly spatially structured, under these conditions we expect simple geostatistical models to fail to recover the pollutant concentration surface. To assess whether this was indeed the case, we compared the performance of the shared model for preferential sampling with that of a Bayesian kriging model in terms of (i)

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-14.jpg?height=1042&width=1113&top_left_y=230&top_left_x=397}
\captionsetup{labelformat=empty}
\caption{Figure 8. Ratio of the posterior standard deviations of the predictions of annual average \(\mathrm{PM}_{10}\) concentration obtained when accounting vs. not accounting for preferential sampling. Lombardy region, year 2007.}
\end{figure}
bias of the predictions, calculated as the mean absolute difference between the predicted values and the true ones, averaged across the simulated datasets; (ii) average variance of the predictions, averaged across simulations; and (iii) mean squared error between predicted values and true values, averaged across the simulated datasets. Table 1 shows the results of our simulation study.

As Table 1 indicates, Bayesian kriging behaves poorly in presence of preferential sampling, as expected. In epidemiological applications, the use of covariates may help reduce the bias, but since preferential sampling is defined on the residual association between the spatial process and the covariates, not including the covariates in our simulations does not alter the results. Notice that, as indicated by the results under scenario 3, our model is able to capture even small deviations from pure random sampling. Indeed, the strong spatial structure of the spatial process may induce a small association even with a spatial random sample of locations.

\section*{4 Discussion}

In this paper, we have presented a general Bayesian hierarchical spatial model that accounts for preferential sampling by introducing a latent shared spatial component that is included both in the model for the point process of the sampling locations and in the model for the point-referenced process underlying the measured geo-referenced variable. Our model allows us to address three

\begin{figure}
\includegraphics[max width=\textwidth]{https://cdn.mathpix.com/cropped/f8386915-40ad-486a-8227-0ccad23d3249-15.jpg?height=972&width=1114&top_left_y=230&top_left_x=351}
\captionsetup{labelformat=empty}
\caption{Figure 9. True \(\log \mathrm{PM}_{10}\) concentrations used for the simulation study (mean 3.47, SD 0.19).}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table I. Bias, variance, and mean squared error (MSE) for a Bayesian kriging model and for the proposed shared model that accounts for preferential sampling under the three preferential sampling scenarios considered in the three simulated datasets: high preferential sampling, low preferential sampling, and no preferential sampling, respectively.}
\begin{tabular}{|l|l|l|l|l|}
\hline \multicolumn{2}{|c|}{} & Scenario I: high preferential sampling & Scenario 2: low preferential sampling & Scenario 3: no preferential sampling \\
\hline \multirow[t]{3}{*}{Bayesian kriging} & Bias & 0.1292323 & 0.0727762 & 0.05184599 \\
\hline & Variance & 0.005180182 & 0.004671143 & 0.004124283 \\
\hline & MSE & 0.02188116 & 0.009967519 & 0.006812289 \\
\hline \multirow[t]{3}{*}{Shared model for preferential sampling} & Bias & 0.03004569 & 0.04551036 & 0.04653664 \\
\hline & Variance & 0.00451783 & 0.004164269 & 0.004082967 \\
\hline & MSE & 0.005420574 & 0.006235463 & 0.006248626 \\
\hline
\end{tabular}
\end{table}
different issues that can arise and are of interest when modeling a spatial process under preferential sampling: (1) continuous or finite spatial sampling locations population; (2) underlying causal model and relevant covariates; and (3) two different inferential goals-mean prediction surface or prediction uncertainty. The two original examples that we discuss in the paper highlight each of these issue, respectively.

The typical approach used in geostatistics to account for preferential sampling is to specify an inhomogeneous Poisson point process for the measurement locations. Such specification allows to
account for the dependence between the spatial intensity of the location point process and the measurement spatial process, which otherwise would affect the statistical inference. In many applications, the region of interest may present barriers or areas where measurements cannot be taken-e.g. lakes or rivers, inaccessible areas, etc.-and in some circumstances, the set of potential sampling locations may be finite. For example, in veterinary surveillance, the primary sampling unit is the farm and the number of farms in a region is finite. Shaddick and Zidek \({ }^{8}\) presented a similar finite population context in their example on air pollution monitoring stations in UK. In this instance, conditioning on the measurement location being usable, a Bernoulli model for the sampling probability is more appropriate than an inhomogeneous Poisson process. In our first example, we overimpose a fine grid over the spatial domain (the Campania region) and, having counted the population denominator (e.g. total number of farms) per grid cell, we model the sampling probability for each grid cell using a Binomial model. Notice that we can do that since we have a very large total population of 8794 farms (Figure 1).

In our second example, we use an inhomogeneous point process for the air pollution monitoring station locations. For computational ease, since we want to integrate information from a deterministic model output with the monitoring data, we report the observational data and the numerical model output on the same fine grid. As we have only 58 monitoring stations over 1679 grid cells, resulting in a large excess of zero counts, for the monitoring stations counts per grid cell we use a Negative Binomial model. In a previous paper that employed the same data \({ }^{10}\) we have instead used a Poisson model. The results that we obtain using the two different model specifications are very similar. In the Negative Binomial model, we specify the dispersion parameter in the likelihood while we introduce the shared spatial component in the log linear predictor. Conversely, in the Poisson model of Grisotto et al., \({ }^{10}\) the log linear predictor contains both the heterogeneity and the clustering random terms and, since they are not identifiable, they are considered jointly in the shared component (see for a similar modeling approach, Jürgens et al. \({ }^{19}\) ).

Using the extended definition of preferential sampling by Lee et al., \({ }^{5}\) a model that aims to account for preferential sampling ought to include also the underlying causal model and relevant covariates. Therefore, important covariates explaining both processes must be specified in the joint model-as stated by Diggle et al. \({ }^{4}\) and Pati et al. \({ }^{6}\) Other covariates may be important only for the location process or only for the spatial process, but not for both. It is still left to investigate what is the gain in efficiency or the reduction in potential bias achieved by including other covariates in the model equations. In principle, covariates in the spatial process model equation should improve efficiency, while covariates in the location process model equation may reduce efficiency or even produce bias in case of a collider (see for example a similar discussion in the propensity score literature by Brookhart et al. \({ }^{20}\) ). It is important to observe that if the systematic component of the model is correctly specified, the dependence between the sampling location process and the spatial process depends on the underlying latent residual shared spatial component.

In our example on veterinary parasitology, we did not consider covariates, as only remote sensing covariates were available and we did not expect them to influence farmers' behavior. On the contrary, in our second example, the covariate, represented by the deterministic model predictions of air pollutant concentrations, is expected to be very accurate, and since these predictions are produced by the same agency who is responsible for the monitoring stations locations (the local Environmental Protection Agency), we expect the sampling location process to be strongly influenced by this covariate. The residual shared spatial component is still important and emphasizes the role of other socio-political factors in locating the air pollution monitors.

Most of the literature on preferential sampling has focused mainly on addressing the potential prediction and estimation bias due to the dependence between the sampling and the measurement
processes, see for example Diggle et al. \({ }^{4}\) where the bias in the estimates of geostatistical parameters is discussed. In our second example, we have shown that in exposure assessment of air pollutants, preferential sampling has an effect also on prediction uncertainty. This is a result that, besides an initial attempt by Gelfand et al., \({ }^{7}\) none of the previous literature on preferential sampling has explored fully, and it is the main point of interest and conclusion of our second case study. A large body of literature on air pollution exposure assessment has shown that when trying to obtain an estimate of the spatial distribution of some pollutants, statistical models that leverage the information contained in deterministic or dispersion models are superior to simple kriging models or land use regression approaches (for an example on air quality model, see CMAQ: http://www.epa.gov/asmdnerl/CMAQ). These geostatistical models that combine deterministic model outputs with monitoring data to derive estimates of the spatial distribution of air pollution are also used to obtain uncertainty estimates, see for example Bayesian melding \({ }^{21}\) or downscaling. \({ }^{22}\) However, none of these previous modeling efforts have assessed how the preferential location of the monitors affects the uncertainty in the predictions. Here, we explored this issue more in depth and we showed how an adjustment for preferential sampling can be introduced in this kind of models. To quantify the effect of preferential sampling on the prediction uncertainty, we used the simple ratio of standard deviations obtained when accounting vs. non accounting for preferential sampling. We believe that more work has to be devoted to this issue to better understand how to evaluate different standard deviation surfaces.

To approximate joint posterior densities and joint predictive densities, we used an MCMC approach working on a fine grid superimposed over our spatial domain of interest. In both our case studies, we fitted the models for the locations counts on the same grid used for estimation and prediction of the spatial process. Although advantageous and convenient in our applications, this can be limiting in situation where the data do not consent this reduction. Diggle et al. \({ }^{23}\) extended geostatistical approaches to log-Gaussian Cox processes. This extension allows to perform inference that is independent of the particular partition of the region of interest and is particularly relevant when we cannot consider the grid to be fine on a subject-specific field, or when we must rely on a priori defined areas-for example, administrative units.

\section*{5 Conclusion}

Preferential sampling refers to any situation in which the spatial process and the sampling locations are not stochastically independent. Taking advantage of two examples of geostatistical analysis in which the usual assumption of stochastic independence between the point process and the measurement process is violated, we highlighted three different modeling and inferential aspects of geostatistical modeling under preferential sampling: (1) continuous or finite spatial sampling frame; (2) underlying causal model and relevant covariates; and (3) inferential goals related to mean prediction surface or prediction uncertainty. We specified a flexible and general Bayesian geostatistical model that includes a shared spatial random component. A simulation study showed that our model is able to capture even small deviations from pure random sampling.

\section*{Acknowledgements}

We thank the Air Quality Unit of the Lombardy Regional Agency for the Environment and the Health Directorate of the Region of Lombardy Government for providing the data that made this work possible. We thank Laura Rinaldi and Giuseppe Cringoli of the University of Naples Federico II for parasitological
data, Dario Consonni and Pier Alberto Bertazzi of the University of Milan for environmental epidemiology example, and all of them for help in interpreting results.

\section*{Declaration of conflicting interests}

The author(s) declared no potential conflicts of interest with respect to the research, authorship, and/or publication of this article.

\section*{Funding}

The author(s) disclosed receipt of the following financial support for the research, authorship, and/or publication of this article: This work was supported by the Italian Ministry of University and Scientific Research. The research leading to these results has received funding also from the European Union Seventh Framework Programme FP7-KBBE-2011-5 under grant agreement no 288975 and from the Regional Government of Lombardy (ESSIA project).

\section*{References}
1. Cressie N. Statistics for spatial data. London: Chapman and Hall, 2002.
2. Diggle P, Tawn JA and Moyeed RA. Model-based geostatistics. J Royal Stat Soc C 1998; 47: 299-350.
3. Diggle P and Lophaven S . Bayesian geostatistical design. Scand J Stat 2006; 33: 53-64.
4. Diggle P, Menezes R and Su TL. Geostatistical inference under preferential sampling. J R Stat Soc C 2010; 59: 91-232.
5. Lee A, Szpiro A, Kim SY, et al. Impact of preferential sampling on exposure prediction and health effect inference in the context of air pollution epidemiology. Environmetrics 2015; 26: 255-267.
6. Pati D, Reich BJ and Dunson DB. Bayesian geostatistical modelling with informative sampling locations. Biometrika 2011; 98: 35-48.
7. Gelfand AE, Sahu SK and Holland DM. On the effect of preferential sampling in spatial prediction. Environmetrics 2012; 23: 565-578.
8. Shaddick G and Zidek JV. Unbiasing estimates from preferentially sampled spatial data. Spat Stat 2015; 9: 43.
9. Cecconi L, Biggeri A, Grisotto L, et al. Preferential sampling in veterinary parasitological surveillance. Geospat Health 2015; 11: 412.
10. Grisotto L, Consonni D, Cecconi L, et al. Geostatistical integration and uncertainty in pollutant concentration surface under preferential sampling. Geospat Health 2016; 11: 426.
11. Guttorp P and Sampson P. Discussion of Geostatistical inference under preferential sampling by Diggle PJ, Menezes R and Su T. J Royal Stat Soc C 2010; 59: 191-232.
12. Rinaldi L, Biggeri A, Musella V, et al. Sheep and Fasciola hepatica in Europe: the GLOWORM experience. Geospat Health 2015; 9: 309-317.
13. Vicedo-Cabrera AM, Biggeri A, Grisotto L, et al. A Bayesian kriging model for estimating residential exposure
to air pollution of children living in a high-risk area in Italy. Geospat Health 2013; 8: 87-95.
14. Giorgi E, Sesay SS, Terlouw DJ, et al. Combining data from multiple spatially referenced prevalence surveys using generalized linear geostatistical models. J Royal Stat Soc A 2014; 178: 445-464.
15. Held L, Natario I, Fenton SE, et al. Towards joint disease mapping. Stat Methods Med Res 2005; 14: 61-82.
16. Besag J, York J and Mollié A. Bayesian image restoration, with two applications in spatial statistics. Ann Inst Stat Math 1991; 43: 1-59.
17. Banerjee S, Carlin BP and Gelfand AE. Hierarchical modeling and analysis for spatial data. London: Chapman and Hall, 2006.
18. Silibello C, Calori G, Brusasca G, et al. Modelling of \(\mathrm{PM}_{10}\) concentrations over Milano urban area using two aerosol modules. Environ Model Softw 2008; 23: 333-343.
19. Jürgens V, Ess S, Phuleria HC, et al. Bayesian spatiotemporal modelling of tobacco-related cancer mortality in Switzerland. Geospat Health 2013; 7: 219-236.
20. Brookhart MA, Schneeweiss S, Rothman KJ, et al. Variable selection for propensity score models. Am J Epidemiol 2006; 163: 1149-1156.
21. Fuentes M and Raftery AE. Model evaluation and spatial interpolation by Bayesian combination of observations with outputs from numerical models. Biometrics 2005; 61: 36-45.
22. Berrocal VJ, Gelfand AE and Holland DM. A spatiotemporal downscaler for output from numerical models. \(J\) Agric Biol Environ Stat 2010; 15: 176-197.
23. Diggle PJ, Moraga P, Rowlingson B, et al. Spatial and spatio-temporal Log-Gaussian Cox processes: extending the geostatistical paradigm. Stat Sci 2013; 28: 542-563.

\section*{Appendix 1: Winbugs Code}

\section*{Case study I (Campania region)}
```
model {
    # p: sampling intensity
    # camp[i]: number of sampled farms in the i-th cell
    # pop[i]: number of farms (population) in the i-th cell
    for (i in 1:Ncell) {
        camp[i] ~ dbin(p[i],pop[i])
        logit(p[i]) <- alpha+(1/delta)*(b[i]+c[i]);
        c[i] ~ dnorm(0,tau.c)
        }
    # ICAR prior distribution
    b[1:Ncell]~car.normal(adj[], weights[], num[], tau.b)
    for (k in 1:sumNumNeigh){
        weights[k] <-1
        }
    alpha~dflat()
    tau.b~dgamma(0.001, 0.001)
    tau.c~dgamma(0.001, 0.001)
        # Spatially structured multivariate normal likelihood
        b.clu[1:Nfarms] ~ spatial.exp(mu[], x[], y[], tau.bclu,phi, 1)
    #exponential correlation function
    for (i in 1:Nfarms) {
    obs01[i]~dbern(q[i])
    logit(q[i])<-alphaq+b.clu[i]+delta*(b[ind[i]]+c[ind[i]])
    mu[i] <-0
    }
        # Priors
        alphaq ~ dflat()
        logdelta ~ dnorm(0, 5.9)
        delta<-exp(logdelta)
        tau.bclu ~ dunif(0.1,100)
        # priors for spatial.exp parameters
        phi <-0.02
        # Spatial prediction
        # Single site prediction
    for (j in 1:Ncell){
        b.clu.pred[j] ~ spatial.unipred(mu.pred[j],x.pred[j],y.pred[j],b.clu[])
        }
        for(j in 1:Ncell) {
        mu.pred[j] <-0
        z.pred[j]<-alphaq+b.clu.pred[j]+delta*(b[j]+c[j])
        obs01.pred[j]<-exp(z.pred[j])/(1+exp(z.pred[j]))
        conta[j] <-step(obs01.pred[j]-0.08)
```

```
# mean infection probability: 0.08
    }
}
```


\section*{Case study 2 (Lombardy region)}
```
model {
    # ncentr[i]: number of monitors in the i-th cell
    # lpm10arpa.pred[i]: predictions from deterministic model
    # idcell[i]: cell id which belongs the i-th monitor
    for (i in 1:Ncell) {
        ncentr[i] ~ dnegbin(p[i], r)
        p[i]<-(r/(r+lambda[i]))
        log(lambda[i])<-alfa+(1/delta)*b[i]+betapmint*lpm10arpa.pred[i]
        }
    b[1:Ncell]~car.normal(adj[], weights[], num[], tau.b)
    for (k in 1:sumNumNeigh){
        weights[k] <-1
        }
    alfa~dflat()
    tau.b~dgamma(0.001, 0.001)
    sigma<-sqrt(1/tau)
    r~ dgamma(0.001, 0.001)
    #kriging
        # Spatially structured multivariate normal likelihood
        lpm10[1:Nmonitor] ~ spatial.exp(mu[], x[], y[], tau, phi, 1) #exponential
correlation function
        for(i in 1:Nmonitor) {
        mu[i] <- alpha+beta*lpm10arpa[i]+delta*b[id_cell[i]]
        }
        # Priors
        alpha ~ dflat()
        beta ~ dnorm(0, 0.001)
        logdelta ~ dnorm(0, 5.9)
        delta<-exp(logdelta)
    tau.d<-1/0.017
        tau ~ dgamma(0.001, 0.001)
        sigma2 <- 1/tau
        betapmint~dnorm(0, 0.001)
        # priors for spatial.exp parameters
        phi <-0.2
        # Spatial prediction. Single site prediction
        for(j in 1:Ncell) {
        beta_pred[j]<-alpha+beta*lpm10arpa.pred[j]+ delta*b[j]
        lpm10.pred[j] ~ spatial.unipred(beta_pred[j], x.pred[j], y.pred[j], lpm10[])
        }
    }
```
