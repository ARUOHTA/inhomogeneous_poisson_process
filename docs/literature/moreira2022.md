# ANALYSIS OF PRESENCE-ONLY DATA VIA EXACT BAYES, WITH MODEL AND EFFECTS IDENTIFICATION 

By Guido A. Moreira ${ }^{1, \mathrm{a}}$ and Dani Gamerman ${ }^{2, \mathrm{~b}}$<br>${ }^{1}$ Centro de Biologia Molecular e Ambiental, Universidade do Minho, ${ }^{\mathrm{a}}$ guidoalber@gmail.com<br>${ }^{2}$ Instituto de Matemática, Universidade Federal do Rio de Janeiro, ${ }^{\mathrm{b}}$ dani@im.ufrj.br


#### Abstract

This paper provides an exact modeling approach for the analysis of presence-only ecological data. Our proposal is also based on frequently used inhomogeneous Poisson processes but does not rely on model approximations, unlike other approaches. Exactness is achieved via a data augmentation scheme. One of the augmented processes can be interpreted as the unobserved occurrences of the relevant species, and its posterior distribution can be used to make predictions of the species over the region of study beyond the observer bias. The data augmentation also leads to a natural Gibbs sampler to make Bayesian inference through MCMC. The proposal shows better performance than the currently standard method based on Poisson process with intensity function depending log-linearly on the covariates. Additionally, an identification problem that arises in the traditional model does not seem to affect our proposal in the analyses of real ecological data.


1. Introduction. Species distribution modeling (SDM) is an important area in ecology which aims to explain and predict the occurrences of one or multiple species across a region. Its applications range from conservation and reserve planning, evolution, epidemiology, invasive-species management and other fields (Phillips, Anderson and Schapire (2006)).

Presence-only data is a problem that arises from opportunistic sampling. Instead of a planned and randomized survey, the data consists of a collection of locations where the species has been sighted. The data is usually originated from citizen science programs or from institutional herbaria and museums.

This means that presence-only data configures preferential sampling in the sense that the data collection process is related to what is being modeled. In other words, a species distribution model, based on presence-only data, can be misleading since the data is usually collected close to where the species is more accessible or observable. Consequently, the data might carry the misleading information that the species more often occurs where the data is observable.

Nevertheless, this type of data is the only information available for many species (Elith and Leathwick (2007)). Furthermore, presence-only data has become more abundant over the years (Warton and Shepherd (2010)) and, therefore, can be a valuable source of information if treated adequately. It was shown by Diggle, Menezes and Su (2010) that the bias introduced by preferential sampling in a spatial statistical model is reduced if the model accounts for preferentiality patterns.

This paper benefits from the discussions in Warton and Shepherd (2010), Fithian and Hastie (2013), Renner and Warton (2013), Dorazio (2014) and Fithian et al. (2015). These authors point out that it makes more sense to model presence-only data through a thinned point process. For this purpose the aforementioned works use an inhomogeneous Poisson process (IPP) for the species occurrences using a log-linear intensity function with respect to

[^0]![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-02.jpg?height=590&width=575&top_left_y=166&top_left_x=595)
FIG. 1. 230 Eucalyptus presence locations over the Greater Blue Mountains World Heritage Area near Sydney, Australia.

intensity covariates. The occurrences are then thinned via logistic regression on observability covariates to describe the observed dataset. This approach will be referred to as traditional IPP, in contrast to our proposal.

One problem that arises from the traditional IPP is that variables in intensity and observability sets cannot be highly correlated between them, as their respective slope parameters become unidentified. Our proposal uses a modification of the model based on the work of Goncalves and Gamerman (2018) to address this issue. However, their model uses spatially varying effects and does not have to deal with an occurrence detection process.

In the proposed model there is empirical evidence which implies that issues regarding model parameters identification have been reduced as a byproduct of performing exact evaluation of the likelihood function. Even more relevant than this, evidence supports that model identification is verified; the intensity and observability functions are uncorrelated a posteriori. Thus, the model was able to separate the two aspects of the data.

An additional contrast with the traditional model, as hinted above, is the use of exact inference in the sense that no approximation is applied to the model. This is achieved through a data augmentation technique with latent point processes. Inference on the posterior distribution is carried out by MCMC in a very straightforward manner.

Two datasets are studied in this paper. One is an eucalyptus (Eucalyptus Sparsifolia) dataset in the Greater Blue Mountains World Heritage Area near Sydney, Australia, and a surrounding 100 km buffer zone, whose 230 presence locations are seen in Figure 1. The data is available in Supplementary Material B (Moreira and Gamerman (2022a)).

This dataset has been analyzed in Renner et al. (2015) and can be used to compare the predictive performance of different models. The comparison is based on a separate independent dataset which was sampled by a systematic unbiased survey. The goal of the analysis is to identify the environmental occurrence propensity for an abundant and broadly distributed species across the studied region (Renner et al. (2015)). For that purpose, key environmental features are included as explanatory variables. All records are the result of incidental sightings since 1972. The database is available where the observed point locations were known to a one km degree of accuracy.

The second analysis is with a dataset of an angico tree species (Anadenanthera Colubrina). It is a tree particular of the semiarid region of Brazil. Its presence locations can be seen in Figure 2. The data has been collected and organized by Oliveira-Filho (2017). It is available in the Supplementary Material B (Moreira and Gamerman (2022a)). In Mazzochini et al.

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-03.jpg?height=593&width=1358&top_left_y=167&top_left_x=115)
FIG. 2. Angico tree presence locations. Left: in most of South America (studied region highlighted with darker grey). Right: in the studied region.

(2019), the data has been organized as presences locations in a one km cell size. This was done for 771 woody species in the semiarid region in Brazil. Their work consisted on a joint analysis with SDMs for relating species diversity with ecological stability.

We focus the analysis on only one species. Although this tree has been sighted in a very wide area in the continent, it is studied in the context of a semi-arid region. Consequently only the sightings in this area are considered in the analysis. Only 140 presence locations of the original 549 are used.

No independent dataset exists that can be used to construct model comparison. However, this study requires the use of the same covariates for intensity and observability, which is not possible under the traditional IPP model. Thus, it allows testing the claim that intensity and observability effects for the same covariate might be identified in our approach.

This paper presents a brief review of the existing methods in Section 2. The model is presented in a general form in Section 3 and with a specific parametric form for the intensity function and thinning in Section 4, where some analytical properties are discussed. The inference procedure and prediction are detailed in Section 5. In this section a discussion on model comparison is presented. A probabilistically well-defined method for such prediction that benefits from the discussion in Gelfand and Schliep (2018) is proposed. A few simulation studies are shown in Section 6 to test model features in a controlled environment, and in Section 7 the model is applied to the two aforementioned datasets. The proposed model is compared to the traditional Poisson process model in terms of prediction with the eucalyptus dataset. Issues of parameter and model identification are shown to be reduced with the angico tree dataset. Finally, some concluding remarks are made in Section 8, along with a discussion on current and future challenges of presence-only data modeling.
2. Brief review of the literature. The existing methods for presence-only data have different approaches to the problem and they are evidenced below. Our proposal addresses specific issues with some of them, and this is pointed out as well.

The pseudo-absence model of Pearce and Boyce (2006) chooses random locations in the region and assumes they are absences. Then, the presence and pseudo-absence locations are treated by a logistic regression based on intensity covariates. The main problem with this is that there is no guarantee that the absence locations are actual absences which may lead to undesired bias in the estimates. Also, treating presence-only data as presence-absence locations potentially skews the analysis based on the chosen lattice size (Gelfand and Shirota (2019)).

Our proposal is a point process which makes no such assumptions and correctly treats the data at the continuous point level.

The numerical solution of the MaxEnt method (Phillips, Anderson and Schapire (2006), Phillips, Dudík and Schapire (2004), Phillips et al. (2009)) and the maximum likelihood estimate for the traditional IPP are mathematically equivalent, as proved by Warton and Shepherd (2010), Fithian and Hastie (2013) and Renner and Warton (2013). The latter uses an approximation to the integral of the intensity function over the region which implicitly assumes that it is constant by parts. The unthinned version of the traditional IPP does not include an observability process, which can lead to undesired estimates bias, since it does not deal with preferential sampling. The thinned version, on the other hand, has the same theoretical foundation of the proposed model. However, since observation probabilities are usually small, the model is effectively the same as adding the observability covariates in the intensity set based on the approximation

$$
\begin{equation*}
e^{a} \frac{e^{b}}{1+e^{b}} \approx e^{a+b} \tag{1}
\end{equation*}
$$

when $b$ is negative and $|b|$ is large. Thus, the method usually fits an unthinned Poisson process, adding the observability covariates in the intensity set. Thus, there should not be highly correlated covariates between the intensity and observability sets, as their effects would not be identifiable. This also applies to the intercepts. That is to say that the two linear predictors from the intensity and observability processes are prone to be (approximately) merged into one predictor. Thus, correlated covariates appear in a single linear predictor, leading to the same previously mentioned identification problems.

Method MaxLike (Royle et al. (2012)) promises to estimate species prevalence, defined as the presence probability in a randomly chosen lattice cell of the regularly discretized region. A problem of their approach is that they assume that the cells with observed presences are a uniform random sample of all region cells with presences which ignores preferential sampling. Another problem is that the adequate estimation of the model depends on a logistic relationship between the observations and parameters. Estimates suffer greatly from a deviation of this premise, as discussed by Hastie and Fithian (2013). Our proposal not only incorporates preferential sampling in the modeling, but deviations from the parametric form do not seem to heavily impact the result, as verified in Section 6.

It is also worth noting that Hefley et al. (2013) approaches the issue of preferential sampling using the missingness perspective. Instead of using covariates to model the sampling bias, they use expert knowledge to subjectively elicit species presence probabilities over the region. As such, they provide an independent source of information, which is required to correct Not-Missing-At-Random sampling, as described by Rubin (1976) and Little and Rubin (1991). This is an important detail, considering that presence-only data does not have information about the sampling mechanism. This approach was also applied in Hefley et al. (2015). However, the inclusion of observability covariates in the model provides an alternative approach to treat sampling bias.

Another important work worth mentioning is Gelfand and Schliep (2018) who discuss the probabilistic aspects of modeling presence-only data with various models in detail. They propose an approach inspired in Diggle, Menezes and Su (2010) that considers preferential sampling to model the nonideal data collection in the observations.
3. The inhomogeneous Poisson process model. The inhomogeneous Poisson process (IPP) is a model for point patterns (Cressie (1993), Diggle (2014)), that is, a random set of points over a region, say $\mathcal{D}$, for which both the amount and the locations of the points is random. The model has a positive functional parameter called intensity function $\lambda(s), s \in \mathcal{D}$
whose only restriction is that $\int_{\mathcal{D}} \lambda(s) d s<\infty$. A Poisson process whose intensity function is constant is called a homogeneous Poisson process (HPP) or completely random process.

In addition to the pairwise independence of counts in disjoint subregions, the IPP is characterized by the fact that, for any Borelian sub-region $A \subset \mathcal{D}$, the number of points in $A$ follows a Poisson distribution with rate $\int_{A} \lambda(s) d s$.

The thinning operation (Cressie (1993)) can be used to relate the species occurrences with the presence-only data. The operation is advantageous since it has a direct interpretation in the data. Namely, there is a probability inherent to every location $s \in \mathcal{D}$, say $p(s)$, which is related to its observability, that an eventual occurrence is observed. Additionally, the thinning of a Poisson process with intensity function $\lambda(\cdot)$ is also a Poisson process with intensity function $p(\cdot) \lambda(\cdot)$.

A problem of the IPP model is with respect to its likelihood function. If $X$ is the presenceonly observed locations, then the likelihood function resulting from the model is

$$
\begin{equation*}
L_{x}(\lambda(\cdot), p(\cdot))=e^{-\int_{\mathcal{D}} \lambda(s) p(s) d s} \prod_{s \in x} \lambda(s) p(s) \tag{2}
\end{equation*}
$$

In most practical problems the integral $\int_{\mathcal{D}} \lambda(s) p(s) d s$ cannot be solved analytically. Even if the intensity function is assumed to be known at every location $s \in \mathcal{D}$, its functional form over the region is usually not known. Most solutions apply some sort of numerical approximation for this integral.

Our proposal uses the idea from Adams, Murray and MacKay (2009) and Gonçalves and Gamerman (2018). Those works did not deal with an occurrence detection process. This paper extends them to include a process to represent the unobserved occurrences. In essence, the intensity function $\lambda(\cdot)$ is assumed to be upper bounded and is parameterized so that $\lambda(s)=q(s) \lambda^{*}, s \in \mathcal{D}$, where $q: \mathcal{D} \rightarrow(0,1)$. In this parameterization, $\lambda^{*}$ represents an upper bound for the intensity function. Thus, the presence-only process $X$ has intensity function $q(\cdot) p(\cdot) \lambda^{*}$.

The proposed model adds two latent processes $X^{\prime}$ and $U$, which are conditionally independent between them and with respect to the observed points $X$, given the parameters. The process $X^{\prime}$ represents all unobserved occurrences and has intensity function $q(\cdot)(1-p(\cdot)) \lambda^{*}$. The points in $U$ can be viewed as an additional set of points needed to evaluate the analytically exact likelihood function because it is modeled with intensity $(1-q(\cdot)) \lambda^{*}$. A reviewer noted that "in terms of the experiment, $U$ can be thought of as a hypothetical set of points that might have been sampled and the species was not present if there had been uniform sampling across the entire region."

More explicitly, the augmented model is

$$
\begin{align*}
X & \sim \operatorname{IPP}\left(q(\cdot) p(\cdot) \lambda^{*}\right),  \tag{3}\\
X^{\prime} & \sim \operatorname{IPP}\left(q(\cdot)(1-p(\cdot)) \lambda^{*}\right),  \tag{4}\\
U & \sim \operatorname{IPP}\left((1-q(\cdot)) \lambda^{*}\right) . \tag{5}
\end{align*}
$$

By writing the likelihood function resulting from this joint model, the integral becomes trivial, since $q(s) p(s)+q(s)(1-p(s))+(1-q(s))=1, \forall s \in \mathcal{D}$,

$$
\begin{align*}
& L_{x}\left(q(\cdot), p(\cdot), \lambda^{*}, x^{\prime}, u\right) \\
& \quad=L_{x}\left(q(\cdot), p(\cdot), \lambda^{*}\right) f\left(x^{\prime} \mid q(\cdot), p(\cdot), \lambda^{*}\right) f\left(u \mid q(\cdot), \lambda^{*}\right) \\
& \quad=\frac{e^{-\int_{\mathcal{D}} \lambda^{*} d s}}{n_{x^{\prime}}!n_{u}!} \prod_{s \in x} q(s) p(s) \lambda^{*} \prod_{s \in x^{\prime}} q(s)(1-p(s)) \lambda^{*} \prod_{s \in u}(1-q(s)) \lambda^{*}  \tag{6}\\
& \quad=\frac{e^{-\lambda^{*}|\mathcal{D}|}}{n_{x^{\prime}}!n_{u}!} \lambda^{* n_{x}+n_{x^{\prime}}+n_{u}} \prod_{s \in x} q(s) p(s) \prod_{s \in x^{\prime}} q(s)(1-p(s)) \prod_{s \in u}(1-q(s)),
\end{align*}
$$

where $n_{x^{\prime}}$ and $n_{u}$ are, respectively, the number of points of $X^{\prime}$ and $U$ in $\mathcal{D}$, and $|\mathcal{D}|$ is the Lebesgue measure of $\mathcal{D}$, namely, its area. Note that the troublesome integral disappears in the augmented likelihood.

An issue about the introduction of $\lambda^{*}$ in the model must be addressed, however. It is, strictly speaking, an upper bound for the occurrences intensity function in the region. However, if a value $c$ is adequate for this parameter, then any other value greater than $c$ is also adequate with an appropriately lower value for $q(\cdot)$. This can be seen as an identification problem for this parameter.

It is argued in Gonçalves and Gamerman (2018) that some prior information should be given to increase the model's identifiability. Taking on a parametric form for $q(\cdot)$ will also help with this issue, yet if there is any information about the species' maximum rate of occurrence in a given subregion, then it should be included in the prior distribution. Such information can be procured from a specialist, for example.
4. Modeling the intensity function and retaining probabilities. A desired feature of a species distribution model is to relate the occurrence of the species with environmental or geographical explanatory variables. In the traditional model this is achieved by modeling the log of the intensity function as a linear regression of these variables, say $Z(\cdot) \beta$, where $\beta$ is a parameter vector to be estimated and the intensity covariates set $Z(s)$ is assumed to be known at every location $s \in \mathcal{D}$.

Since the observability is also inherent to each location, the retaining probabilities of the thinning operation is also usually modeled as a logistic-linear function of a regression on another set of explanatory variables, say $W(\cdot) \delta$. The parameter vector $\delta$ is also to be estimated, and $W(s)$ is, analogously to $Z(s)$, a set of observability covariates assumed to be known at every location $s \in \mathcal{D}$.

As a direct consequence, the traditional IPP model for presence-only data deals with an intensity function for the data, such as

$$
\begin{equation*}
\lambda(s) p(s)=e^{Z(s) \beta} \frac{e^{W(s) \delta}}{1+e^{W(s) \delta}}, \quad s \in \mathcal{D} \tag{7}
\end{equation*}
$$

However, Warton and Shepherd (2010), Fithian and Hastie (2013), Renner and Warton (2013), Dorazio (2014) and Fithian et al. (2015) argue that, since the observability is commonly low, equation (7) is approximately $\exp \{Z(s) \beta+W(s) \delta\}$. These authors claim that variables in $Z(\cdot)$ should not be highly correlated to variables in $W(\cdot)$. They also argue that there cannot be an intercept for both $\beta$ and $\delta$ vectors. This is particularly problematic if a covariate that explains species occurrence also happens to explain observability.

For our proposal, however, the set of covariates $Z(\cdot)$ relates to the intensity function through $q(\cdot)$ which must vary between 0 and 1 . Any cumulative distribution function can be used, including the logistic. Therefore, the aforementioned identification issue only happens when both $q(\cdot)$ and $p(\cdot)$ are simultaneously low.

Note that only the detection probabilities $p(\cdot)$ need to be low for there to be identification problems with the traditional log-linear link for the intensity function. However, the species occurs more often where the intensity function is in its higher values and is more commonly registered when the observability is in the higher values of its spectrum. This means that, for the proposed model, it is unlikely that there will be identification issues, unless the species is rare and the detection probabilities are low everywhere, unlike the traditional model where the problem arises whenever there is only improbable detection.

Thus, equations (3), (4), (5) and (6) are completed with the parameterization

$$
\begin{equation*}
q(s)=\frac{e^{Z(s) \beta}}{1+e^{Z(s) \beta}} \quad \text { and } \tag{8}
\end{equation*}
$$

$$
\begin{equation*}
p(s)=\frac{e^{W(s) \delta}}{1+e^{W(s) \delta}} \tag{9}
\end{equation*}
$$

As an added feature, the parametric forms of equations (8) and (9) relate our proposal to the traditional IPP model which can be verified in Proposition 1.

Proposition 1. Let $Z(\cdot) \beta=\beta_{0}+\tilde{Z}(\cdot) \tilde{\beta}$, where $\tilde{Z}(\cdot) \tilde{\beta}$ represents a linear regression with no intercept. Also, let $\beta_{0}$ and $\lambda^{*}$ move jointly by the relationship

$$
\begin{equation*}
\frac{e^{\beta_{0}}}{1+e^{\beta_{0}}} \lambda^{*}=c \tag{10}
\end{equation*}
$$

for some positive constant $c<\lambda^{*}$. Then,

$$
\begin{equation*}
\lim _{\substack{\lambda^{*} \rightarrow \infty \\ \beta_{0} \rightarrow-\infty}} \frac{e^{\beta_{0}+\tilde{Z}(\cdot) \tilde{\beta}}}{1+e^{\beta_{0}+\tilde{Z}(\cdot) \tilde{\beta}}} \lambda^{*}=e^{\beta_{0}^{\prime}+\tilde{Z}(\cdot) \tilde{\beta}} \tag{11}
\end{equation*}
$$

where $e^{\beta_{0}^{\prime}}=c$.

Proof of Proposition 1 is found in the Supplementary Material A, Section 1 (Moreira and Gamerman (2022b)). Proposition 1 shows that our proposal approaches the traditional model when $\lambda^{*}$ grows unbounded.
5. Inference and prediction. Inference on the proposed model is based on the posterior distribution, obtained via the Bayes rule,

$$
\begin{equation*}
\pi(\Theta \mid x) \propto L_{x}(\Theta) \pi(\Theta) \tag{12}
\end{equation*}
$$

where $\Theta=\left(\beta, \delta, \lambda^{*}, x^{\prime}, u\right)$. The resulting posterior distribution is not known in closed form, so approximation via Markov chain Monte Carlo (MCMC) (Gamerman and Lopes (2006)) is performed. The prior distribution could be chosen to introduce as little information as possible. The prior $\pi\left(\beta, \delta, \lambda^{*}, x^{\prime}, u\right)$ is factorized as $f\left(x^{\prime} \mid \beta, \delta, \lambda^{*}\right) f\left(u \mid \beta, \lambda^{*}\right) \pi(\beta) \pi(\delta) \pi\left(\lambda^{*}\right)$. This choice leads to a complete Gibbs sampler described below.

For $X^{\prime}$ and $U$, the priors are already in the model's likelihood function in equation (6). Their respective full conditional distributions are also IPP. They can be sampled jointly following the algorithm from Lewis and Shedler (1979).

For the full conditional sampling of both $\beta$ and $\delta$, a sampling scheme is built based on logistic regression. By removing constant terms relative to each of these parameter vectors from the likelihood function, the form of a logistic regression is found. The locations $x \cup x^{\prime}$ configure the successes for $\beta$ while $u$ are the failures. Regarding $\delta$, this happens for the sets $x$ and $x^{\prime}$, respectively. Since they are independent a priori, their respective full conditionals are the posterior of such a regression model. In order for a Gibbs sampler to be used, the data augmentation scheme from Polson, Scott and Windle (2013) is employed. This method achieves a conditional likelihood for a weighted linear regression model, given the augmented data. Any shrinkage inducing or model selection prior for linear regression can be used as $\pi(\beta)$ and $\pi(\delta)$ for analogous effect.

Finally, if a Gamma prior for $\lambda^{*}$ is used, then its full conditional is also Gamma, and so this choice is made. In particular, if $\lambda^{*} \sim \operatorname{Gamma}(a, b)$, then its full conditional is $\operatorname{Gamma}(a+ \left.n_{x}+n_{x^{\prime}}+n_{u}, b+|\mathcal{D}|\right)$.

The full algorithm is provided with a pseudo-code in Algorithm 1. Furthermore, this algorithm has been implemented in C++ under the R package bayesPO (Moreira (2021)).

```
Algorithm 1 MCMC procedure
    Initialize $\lambda^{*}, \beta, \delta$ as iteration $t=0$.
    for $t$ from 1 to the chosen number of iterations do
        Sample $X^{\prime}$ and $U$ :
        Sample $\mathcal{Y} \sim \operatorname{Poisson}\left(\lambda^{*}|\mathcal{D}|\right)$.
        Spread $\mathcal{Y}$ points uniformly through $\mathcal{D} .^{\star \ddagger}$
        for each point $s$ spread in the last step do
            Calculate $q(s)$ and sample $\mathcal{U} \sim \operatorname{Unif}(0,1)$.
            if $\mathcal{U}>q(s)$ then
                Assign $s$ to $U$.
            else
                Calculate $p(s)$.
                if $\mathcal{U}>q(s) p(s)$ then
                    Assign $s$ to $X^{\prime}$.
                else
                    Discard $s$.
                end if
            end if
        end for
        Sample $\lambda^{*} \sim \operatorname{Gamma}\left(a+n_{x}+n_{x^{\prime}}+n_{u}, b+|\mathcal{D}|\right)$
        Sample $\beta$ :
        Define a $\left(n_{x}+n_{x^{\prime}}\right)$-dimensional vector filled with 1 , associated with corresponding in-
        tensity covariates from $X$ and $X^{\prime}$.
        Define a $n_{u}$-dimensional vector filled with 0 , associated with corresponding intensity
        covariates from $U$.
        Sample $\omega_{b} \mid \cdot \sim P G(Z(s) \beta) \cdot{ }^{\dagger}$
        Sample $\beta \mid \cdot \sim \mathcal{N}\left(m_{b}, V_{b}\right) .^{\dagger}$
        Sample $\delta$ :
        Define a $n_{x}$-dimensional vector filled with 1 , associated with corresponding observabil-
        ity covariates from $X$.
        Define a $n_{X^{\prime}}$-dimensional vector filled with 0 , associated with corresponding observabil-
        ity covariates from $X^{\prime}$.
        Sample $\omega_{d} \mid \cdot \sim P G(W(s) \delta) .^{\dagger}$
        Sample $\delta \mid \cdot \sim \mathcal{N}\left(m_{d}, V_{d}\right) .^{\dagger}$
        Store $X^{\prime}, U, \lambda^{*}, \beta$ and $\delta$ as iteration $t$.
    end for
    * If $\mathcal{D}$ is very irregular and it is hard to spread points uniformly in it, this can be
    equivalently achieved by choosing an enveloping rectangle $\mathcal{A}$ around $\mathcal{D}$, sampling $\mathcal{Y}^{\prime} \sim$
    Poisson $\left(\lambda^{*}|\mathcal{A}|\right)$, spreading $\mathcal{Y}^{\prime}$ points uniformly in $\mathcal{A}$ then discarding the points that fall
    outside of $\mathcal{D}$.
    ${ }^{\ddagger}$ If $\mathcal{D}$ is stored in the computer with identically sized pixels, this step can be achieved by
    uniformly sampling the pixels with replacement.
    ${ }^{\dagger} P G$ stands for the Pólya-Gamma distribution. Its definition, sampling and relation to
    logistic regression are described in Polson, Scott and Windle (2013). The adequate values
    for $m_{b}, m_{d}, V_{b}$ and $V_{d}$ are, therefore, explicitly acquired from their discussion.
```

Only a few modeling options are available at the time of writing. Other options will be made available in the future. Supplementary Material B (Moreira and Gamerman (2022a)) includes the package code required for estimation of the models.

A point needs to be made about computation time. Due to the sampling of the latent processes $X^{\prime}$ and $U$, which have random sizes, each iteration of the MCMC requires a computer to evaluate a random number of mathematical operations. In fact, higher values of $\lambda^{*}$ in the chain increases the computation time of that iteration. This implies that the total runtime for estimating the model is random. For this reason, it is recommended to use an initial value for the chain close to the posterior mode which is fast to compute, as it improves the time to achieve Markov chain convergence and the expected computing time.
5.1. Prediction. Prediction in the present-only context is not selfevident. The most natural way to make predictions with Poisson processes in the continuous scale is by evaluating the intensity function over the studied region. However, that might not be useful for ecologists in the sense that they usually want the predicted unobserved occurrences in a specific discretization of the region. Another use for prediction is comparing different models in their predictive capabilities, but this can be also tricky. Both cases require careful consideration when making the actual prediction.

The prediction of unobserved species occurrences in any specific subset of the region of interest $\mathcal{D}$ is straightforward for our proposed model. The component that represents those very occurrences is $X^{\prime}$, defined in equation (4). This process is sampled in the MCMC procedure, and these realizations can be seen as a sample from the predictive distribution of the infinite-dimensional process of unobserved occurrences. The posterior presence probability in an arbitrary subset of $\mathcal{D}$, such as a lattice cell, is approximated via MCMC by the proportion of samples where $X^{\prime}$ points occur at least once in it. This is done for both real data exercises in Section 7. Analogous reasoning can be applied to obtain the posterior probability of 1,2 or more occurrences in such a subregion.

On the other hand, comparing predictive models is not straightforward in this setting. Excluding a portion of the data to reserve for a test set is problematic since the occurrences exclusion pattern is part of what is being modeled. The best way to measure accuracy is by means of an independent dataset obtained by a randomized survey. Ideally, surveyors go to each previously selected location and verify species occurrence or lack thereof. Such dataset is known as presence-absence data.

It is still problematic to predict presence-absence data using a presence-only model. As discussed in Gelfand and Shirota (2019), a point pattern is not comparable with a dataset composed of zeros and ones. The points made in Gelfand and Schliep (2018) are referred for a deeper and more thorough discussion on point processes prediction.

Still, an argument is made here since binary response surveys check a small area for the species. By being an area, a direct relationship with point processes counts can be made. In the particular case of Poisson processes, the presence-absence dataset is composed of $m$ small areas $\Delta_{i}, i=1, \ldots, m$. For each one of these areas, a single value for each covariate is given, and it is assumed that the covariates are approximately constant in the area. This means that the intensity function is approximately constant in each $\Delta_{i}$ and are denoted $\lambda_{i}, i=1, \ldots, m$. So, the presence probability can be evaluated as

$$
\begin{equation*}
P\left(\text { presence at } \Delta_{i} \mid \Theta\right) \simeq \frac{\lambda_{i}\left|\Delta_{i}\right|}{1+\lambda_{i}\left|\Delta_{i}\right|} \tag{13}
\end{equation*}
$$

for any given parameter $\Theta$. The approximation becomes equality when the covariates are exactly constant in every $\Delta_{i}$. The full derivation of this formula is given in the Supplementary Material A, Section 2 (Moreira and Gamerman (2022b)).

Model accuracy can be assessed with this derivation. A common metric for comparing models in ecology is the AUC, that is, area under ROC Curve (Fletcher Jr. et al. (2019), Journé et al. (2020)), where ROC stands for receiving operating curve (Fielding and Bell (1997)).

Table 1
Chosen values for the parameters to simulate data with the complementary log-log link
| Scenario | $\beta_{0}$ | $\beta_{1}$ | $\delta_{0}$ | $\delta_{1}$ |
| :--- | ---: | ---: | ---: | ---: |
| 1 | -1 | 2 | 3 | 4 |
| 2 | 1 | 2 | 3 | 4 |
| 3 | 1 | 2 | -3 | 4 |


Ideally, the AUC ranges between 0.5 and 1 , but it can be lower than 0.5 . A lower value is undesirable, as 0.5 is the AUC of a random classifier. The highest value of 1 represents the perfect classifier.

Furthermore, note that the AUC is a deterministic function of an estimated model. Therefore a sample from its posterior distribution can be easily obtained. In particular, the posterior mean AUC could be used for a pointwise comparison between the models. A reviewer has pointed out that AUC is not a local and proper scoring rule given the difference between its theoretical and empirical realized values (Byrne (2016)). This means that its use is not to be encouraged. The discussion of an ideal choice on prediction metrics is beyond the scope of this manuscript.
6. Simulated data analysis. As previously seen, the proposed approach of data augmentation to a Poisson process, combined with a reparameterization of the intensity function, may help reduce identification problems with presence-only methods. This is probably achieved by the use of a parametric form for $q(\cdot)$ in equation (6).

As with any parametric choice, it has the potential weakness of displaying poor performance when the data generating process has deviating patterns relative to it. For the proposal, any function taking values in ( 0,1 ) is possible for $q(\cdot)$ and $p(\cdot)$ in equation (6). The logistic link has been explored in a simulated study to test its robustness under a different link function. This is particularly important in the presence-only context, given that Hastie and Fithian (2013) discuss the shortcomings of the MaxLike method regarding this very problem.

In particular, another function in the same domain is tested for $q(\cdot)$. The inverse of the complementary log-log link, given by equation (14),

$$
\begin{equation*}
\operatorname{clog} \log ^{-1}(x)=1-e^{-e^{x}} \tag{14}
\end{equation*}
$$

was chosen. It has a few differences with respect to the logistic one, including not being symmetric around 0 .

Additionally, data from the traditional IPP model were generated, and the proposal was used for the fit. Finally, a simulation study was performed on sensitivity with respect to the choice of prior parameters for $\lambda^{*}$. All simulations used $\mathcal{D}$ as the unit square.
6.1. Data generated with the cloglog link. Data was generated from the proposed model using the cloglog inverse link, defined in equation (14) for $q(\cdot)$, but a model with the logistic link was fitted. The generating value for $\lambda^{*}$ was chosen as 1000 . A single covariate was used for the intensity set and also for the observability set. In both cases the covariate values were simulated independently and identically distributed (i.i.d.) from a standard normal distribution. Three sets of values were chosen for the $\beta$ and $\delta$ vectors and are displayed in Table 1.

Additionally the study was replicated 100 times for each scenario. Table 2 shows a summary for the number generated in each scenario.

Table 2
Summary of number of presence-only points locations drawn from the simulations using the complementary log-log link for $q(\cdot)$
| Statistic | Scenario |  |  |
| :--- | :--- | :--- | :--- |
|  | 1 | 2 | 3 |
| Min. | 261.00 | 281.00 | 285.00 |
| 1st Qu. | 301.75 | 300.75 | 308.50 |
| Mean | 314.07 | 313.68 | 317.14 |
| 3rd Qu. | 326.50 | 325.25 | 327.00 |
| Max. | 361.00 | 367.00 | 362.00 |


Since the objective of the analysis is to correctly estimate the intensity function of the species occurrence, the true intensity as a function of the covariate value is compared with the estimated ones. This comparison can be plotted since the intensity function only varies for different values of the single intensity covariate in the set. Results for the first scenario are displayed in Figure 3. The other two results can be found in the Supplementary Material A, Section 3 (Moreira and Gamerman (2022b)).

The estimated intensity functions are close to the true function. In terms of relative bias, the largest differences interestingly happen at the lowest values of the function. This is due to the fact that the denominator is very small. The relative biases stabilize under 0.13 in absolute value when the function is larger. Much of bias at large values of the function happens due to the variability in the estimation of $\lambda^{*}$ which is further explored in Section 6.3. In other words, misspecification of the parametric form of the intensity function does not heavily impact model estimation.
6.2. Data generated from the traditional IPP model. The result from Proposition 1 shows that the proposal is related to the traditional IPP model, that is, the model where $\log \lambda(s)= Z(s) \beta$.

For the purpose of checking if the proposal can be adequate if the traditional IPP is the true model, data was generated from the traditional model and fitted with the proposal using the logistic link. Note that, this time, there is no true value for $\lambda^{*}$. As before, a single covariate was used to both the intensity and observability sets, and they were generated i.i.d. from the standard normal distribution. Additionally, it is better to compare the true intensity function

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-11.jpg?height=392&width=1184&top_left_y=1798&top_left_x=203)
FIG. 3. Estimation and relative bias of the intensity function for the first scenario simulated with the complementary log-log link: (a) Point estimation. The solid black line represents the true function and the gray lines represent the different mean posterior intensity functions, each estimated from a replicated dataset. (b) The relative bias for each of the estimated intensity functions from the different datasets. Black lines represent a $75 \%$ pointwise envelope for the bias at each value of the covariate.

Table 3
Parameters chosen for the generating log-linear model
| $\beta_{0}$ | $\beta_{1}$ | $\delta_{0}$ | $\delta_{1}$ |
| :---: | :---: | :---: | :---: |
| 5 | 2 | 2 | -2 |


Table 4
Number of presence-only points in each dataset drawn the log-linear model
| Min. | 1st Qu. | Mean | 3rd Qu. | Max. |
| :--- | :---: | :---: | :---: | :---: |
| 787.00 | 826.25 | 850.20 | 870.00 | 917.00 |


with its estimator, viewed as a function of the intensity covariate. The chosen parameter values are displayed in Table 3.

The summary of presence points quantity for the different datasets is seen in Table 4.
The prior for $\lambda^{*}$ was set so that it has prior variance equal to $10^{20}$ in order to allow it to be very large. The result of the estimation can be seen in Figure 4. The same comments made about Figure 3 apply here.

As expected, the model was able to properly recover the traditional log-linear intensity. Not only was it able to filter the observability process but also to describe the true intensity function with low relative bias.
6.3. Prior for $\lambda^{*}$ sensitivity. Sensitivity analysis was performed in order to justify the use of certain hyperparameters for the Gamma prior for $\lambda^{*}$, as described in Section 5.

In terms of prior information for $\lambda^{*}$, it is important to remember that it may have some identification issues with $q(\cdot)$ or, in case of the logistic link proposed in Section 4, with either the intercept $\beta_{0}$ or $\delta_{0}$. This means that, indirectly, an informative prior for these parameters leads to prior information for $\lambda^{*}$ as well.

For the simulation a value for 10,000 was chosen for $\lambda^{*}$ and the logistic link was used for both $q(\cdot)$ and $p(\cdot)$. The prior hyperparameters were the same for $a$ and $b$, parameterized so that the prior expected value is 1 the variance is $\frac{1}{b}$.

As with the other simulated studies, the simulations are done with i.i.d. standard normal covariates. Only one covariate was used in the intensity and observability sets each.

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-12.jpg?height=390&width=1181&top_left_y=1800&top_left_x=296)
FIG. 4. Estimation and relative bias of the intensity function for the first scenario simulated with the traditional log-linear IPP: (a) Point estimation. The solid black line represents the true function, and the gray lines represent the different mean posterior intensity functions, each estimated from a replicated dataset. (b) The relative bias for each of the estimated intensity functions from the different datasets. Black lines represent a $75 \%$ pointwise envelope for the bias at each value of the covariate.

The results are displayed in the Supplementary Material A, Section 3 (Moreira and Gamerman (2022b)). It can be seen that the posterior mean values are estimated around the true $\lambda^{*}$, albeit with some high variability. This can be confirmed by the relative bias, which is mostly under 0.1 in absolute value. The only exception is when the prior variance is very small ( $1 \%$ of the true value) in which case the parameter is under-estimated. Given that this does not happen even when the variance is $10 \%$ of the true value, the recommendation is to set the prior variance to be at least equal to the number of observed points in the sample when there is no information about the species abundance.
6.4. Simulation conclusion. The simulated results reveal robust and stable behavior by the model. In a controlled environment it has the capability of adequately estimating data when the true intensity is the traditional log-linear IPP and also when there is misspecification of the intensity link. There is also some liberty in defining prior information for $\lambda^{*}$, as various levels of prior uncertainty lead to posteriors with equivalently low bias.
7. Real data analysis. Two datasets are analyzed in this section, and the proposed model's features are tested. Since the species distribution model is primarily a predictive model, its prediction capabilities are tested. Then, an instance of model identification is exemplified. In the analyses, all covariates have been standardized before model fitting.
7.1. Making and comparing prediction on Eucalyptus Sparsifolia data. The model is applied to the eucalyptus data that can be seen in Figure 1. The covariates used were the same from Renner et al. (2015), namely, number of fires since 1943, minimum and maximum annual temperature, annual rainfall and a categorical soil variable. Observability covariates were distance to main roads and distance to urban areas.

The proposed model has been estimated using independent Normal prior distributions for the parameter vectors $\beta$ and $\delta$. The variance for each individual element of these vectors have been set to 10 . This value has been chosen since slope parameters should seldom be higher than 3 in absolute value in a logistic scale when multiplying a standardized covariate. This is in accordance with the discussion in Gelman, Simpson and Betancourt (2017) regarding prior information. Parameters for the $\lambda^{*}$ Gamma prior distribution were 0.0001 and 0.0001 to represent little prior information.

A summary of the marginal posterior of the parameters is found in Table 5. The results are consistent with already known aspects of these data (Renner et al. (2015)). Namely, the species is associated with low nutrient soils and has a strong association with the minimum annual temperature quadratic effect. Differently from other models, however, the results show a relatively important quadratic rain effect.

It is worth noting that the observability covariates show impact on the data. The distance to roads seems to have more effect in the quadratic term, while the distance to urban centers strongly affect the observations linearly.

In addition to the individual parameter values, Figure 5 shows the prediction of the latent process $X^{\prime}$ which represents all unobserved occurrences of the model.

For the purpose of model comparison, an independent presence-absence dataset must be used for prediction. One such dataset is available in the studied region and can be seen in Figure 6.

After the proposed model is estimated, a sample of the posterior is available from the MCMC procedure. This sample was used to calculate a posterior sample for the AUC, using probability of presence given in equation (13). The traditional IPP can be fitted using R function ppm from package spatstat (Baddeley, Rubak and Turner (2015)), but there is a problem with this.

Table 5
Posterior summary for the estimation of the proposed model for the Eucalyptus Sparsifolia data. Covariates were standardized prior to estimation
| Covariate | Parameter | Mean | Std. Deviation | $q_{2.5}$ | Median | $q_{97.5}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Intercept | $\beta_{0}$ | -6.47 | 1.04 | -8.42 | -6.49 | -4.67 |
| Fire ${ }^{2}$ | $\beta_{1}$ | -0.26 | 0.32 | -0.88 | -0.24 | 0.31 |
| Fire | $\beta_{2}$ | 5.47 | 2.07 | 1.78 | 5.25 | 9.82 |
| $\mathrm{MinT}^{2}$ | $\beta_{3}$ | -3.17 | 1.03 | -5.31 | -3.11 | -1.34 |
| MinT | $\beta_{4}$ | 1.74 | 2.33 | -2.60 | 1.57 | 6.56 |
| MaxT ${ }^{2}$ | $\beta_{5}$ | -1.92 | 1.85 | -5.42 | -1.95 | 1.63 |
| MaxT | $\beta_{6}$ | 0.42 | 1.99 | -3.47 | 0.41 | 4.32 |
| Rain ${ }^{2}$ | $\beta_{7}$ | -2.07 | 0.47 | -3.06 | -2.09 | -1.10 |
| Rain | $\beta_{8}$ | -1.42 | 2.02 | -5.48 | -1.43 | 2.41 |
| Fire x MinT | $\beta_{9}$ | -4.81 | 0.81 | -6.58 | -4.77 | -3.48 |
| MinT x MaxT | $\beta_{10}$ | 2.72 | 1.42 | -0.10 | 2.68 | 5.55 |
| MaxT x Rain | $\beta_{11}$ | 5.62 | 1.79 | 2.30 | 5.58 | 9.23 |
| Fire x MaxT | $\beta_{12}$ | -0.40 | 1.60 | -3.85 | -0.31 | 2.59 |
| MinT x Rain | $\beta_{13}$ | -1.65 | 0.82 | -3.22 | -1.65 | -0.10 |
| Fire x Rain | $\beta_{14}$ | -0.45 | 0.26 | -1.00 | -0.44 | 0.04 |
| Soil 1 | $\beta_{15}$ | 1.49 | 0.97 | -0.04 | 1.52 | 3.34 |
| Soil 2 | $\beta_{16}$ | 3.09 | 0.96 | 1.59 | 3.08 | 4.95 |
| Soil 3 | $\beta_{17}$ | 1.73 | 0.96 | 0.13 | 1.62 | 3.55 |
| Soil 4 | $\beta_{18}$ | 1.17 | 1.06 | -0.75 | 1.22 | 3.26 |
| Soil 5 | $\beta_{19}$ | 1.70 | 1.05 | -0.26 | 1.72 | 3.75 |
| Intercept | $\delta_{0}$ | -4.23 | 0.39 | -4.95 | -4.23 | -3.49 |
| D.MainRoads ${ }^{2}$ | $\delta_{1}$ | -0.32 | 0.15 | -0.62 | -0.31 | -0.05 |
| D.MainRoads | $\delta_{2}$ | 0.23 | 0.17 | -0.08 | 0.23 | 0.57 |
| D.Urban ${ }^{2}$ | $\delta_{3}$ | 0.91 | 0.81 | -0.72 | 0.96 | 2.40 |
| D.Urban | $\delta_{4}$ | -3.16 | 0.44 | -4.06 | -3.15 | -2.34 |
| D.MainRoads x D.Urban | $\delta_{5}$ | 0.57 | 0.21 | 0.17 | 0.57 | 1.00 |
| - | $\lambda^{*}$ | 6313.62 | 784.95 | 5077.00 | 6340.31 | 7678.95 |


![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-14.jpg?height=585&width=661&top_left_y=1677&top_left_x=554)
FIG. 5. Prediction of the latent $X^{\prime}$ process for the Eucalyptus Sparsifolia data over the study region. The colors represent the frequency at which the posterior draw for the process fell in each cell. The colors can be interpreted as the posterior probability that an unobserved occurrence happens in each cell.

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-15.jpg?height=639&width=622&top_left_y=162&top_left_x=484)
FIG. 6. Eucalyptus Sparsifolia presence-absence locations. Red locations are sampled absences, and black locations are sampled presences.

The scale for the area measure must be the same for the presence-only dataset used to fit the model and the presence-absence dataset used to evaluate the AUC. Since the ppm function does not provide the area measure for the region, it cannot be compared with the area from the presence-absence dataset. Therefore, the maximum likelihood estimator for the traditional model was approximated with R function optim which numerically maximizes a function. The integral from equation (2) was approximated by a very fine discretization $s_{i}, i=1, \ldots, N, N=8,620,092$, of the region $\mathcal{D}$, as is shown in equation (15),

$$
\begin{equation*}
\int_{\mathcal{D}} e^{\beta_{0}+Z(s) \beta+W(s) \delta} d s \approx \frac{|\mathcal{D}|}{N} \sum_{s_{i}=1}^{N} e^{\beta_{0}+Z\left(s_{i}\right) \beta+W\left(s_{i}\right) \delta} \tag{15}
\end{equation*}
$$

With the resulting estimated values, the AUC was calculated using a compatible area measure in equation (13). Finally, the posterior mean AUC has also been plotted in the comparison seen in Figure 7.

Not only is the pointwise estimate of our proposal $(\sim 0.618)$ higher than the traditional model ( $\sim 0.587$ ), but the posterior for the AUC has substantial mass above the traditional model's AUC, namely, 0.908 . So, there is over $90 \%$ probability that the AUC under our proposal is larger than the AUC under the traditional approach. This represents substantial evidence that the proposal has improved this particular metric. Thus, our proposal seems to have provided better prediction than the traditional model in the context of this application.

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-15.jpg?height=428&width=702&top_left_y=1875&top_left_x=443)
FIG. 7. AUC comparison for different models. The thick black vertical line indicates the posterior mean of the AUC under our proposal for the prediction. The value 0.5 is highlighted for comparison.

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-16.jpg?height=436&width=731&top_left_y=160&top_left_x=523)
FIG. 8. Marginal bivariate posterior distribution for the pair of parameters for which identifiability is verified.

7.2. Checking model identification on Anadenanthera Colubrina data. The statistical analysis of the dataset from Figure 2 illustrates model identification discussed in Section 4. The appropriate analysis of this dataset requires the observability covariate to be included in the intensity set. This inclusion is not possible for the traditional model. Nevertheless, our approach should be able to handle this case.

The intensity variables were Altitude, Slope, Annual Mean Temperature, Mean Temperature of Warmest Quarter, Mean Temperature of Coldest Quarter, Annual Precipitation, Precipitation of Wettest Month, Precipitation of Driest Month, Precipitation Seasonality (Coefficient of Variation), Precipitation of Wettest Quarter, Precipitation of Driest Quarter, Precipitation of Warmest Quarter, Precipitation of Coldest Quarter, Mean Diurnal Range (Mean of monthly*(max temp-min temp)), Isothermality (Mean Diurnal Range/Temperature Annual Range) (* 100), Temperature Seasonality (standard deviation *100), Max Temperature of Warmest Month, Min Temperature of Coldest Month, Temperature Annual Range (BIO5BIO6), Mean Temperature of Wettest Quarter, Mean Temperature of Driest Quarter. The observability covariate was the distance to the closest research institute.

There are 21 intensity covariates and one observability covariate. Including the quadratic and interaction covariates in the intensity set would result in a 252 -dimensional parameter vector for the intensity set alone which is too large for a 140 -sized dataset. Therefore, the quadratic effect was included only for the observability set.

The linear version of the observability covariate is included in the intensity set, so there are 22 slope parameters in the $\beta$ vector and two in the $\delta$ vector. Thus, the parameters $\beta_{22}$ and $\delta_{2}$ are not identifiable in the traditional IPP model.

Figure 8 shows that these parameters seem to be identified in our proposal after the MCMC procedure has been performed, as they seem to be nearly independent a posteriori. More relevant than parameter identification is model identification. In particular, functions $q(s)$ and $p(s)$ from equation (6) should be identified for all $s$. For this purpose a posterior sample of these functions can be obtained in all locations $s \in \mathcal{D}$, both in the region and in the presenceonly dataset. Correlations close to zero suggest model identification.

Figure 9 displays the posterior correlations between these functions in all locations. It shows that the posterior correlations between the two logistic functions are generally close to zero. Although some spatial dependence can be observed in these correlations, no clear relation with the presence data can be found. In conclusion, the model seems to be completely identified, even when adding observability covariates in the intensity set.

The summary of the posterior distribution for the parameters can be seen in Supplementary Material A, Section 4 (Moreira and Gamerman (2022b)). The species can be seen to have relation to many covariates. In particular, as expected, it is negatively related to the seasonal precipitation and the precipitation in the dry quarters, as it prefers dry weather.

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-17.jpg?height=366&width=972&top_left_y=160&top_left_x=306)
FIG. 9. Model identification evidenced by the posterior correlations between $p(\cdot)$ and $q(\cdot)$ in every location in the region: Left: displayed on the map; Right: displayed in a smoothed histogram where PO represents the 140 presence-only locations.

The prediction of the unobserved occurrences process $X^{\prime}$ is also available for this dataset and can be viewed in Figure 10.

It is interesting to see that the prediction of unobserved occurrences is low around some regions in Figure 10. These are the locations of the research institutions and the distance to them is the single observability covariate used. This is expected, as the prediction is done for the unobserved locations. This means that the regions close to institutions are highly observable, where unobserved occurrences are not frequent.
8. Discussion. Our proposal shows good prediction results, and there is empirical evidence that suggest that the model can be identified, as shown in the results. Although other datasets exist for model comparison, as explained in Elith et al. (2020), the two sets presented here provide a good overview of the model's capabilities.

It is important to note that the apparent model identification and the possibility to estimate intercepts in the two linear predictors is, among other things, a consequence of the parametric logistic form in equations (8) and (9). A question that directly arises is how the model would fare in a context that deviates from this form. The simulation study in Section 6 addresses this point. The model using a logistic link was able to adequately recover the intensity function from data generated using different links, including the traditional log-linear IPP.

The problem associated with the identifiability of $\lambda^{*}$ is discussed in Gonçalves and Gamerman (2018) and is also noted toward the end of Section 3. MCMC convergence can be slowed due to this issue, and thus a starting value for the chain based on a numerical approximation

![](https://cdn.mathpix.com/cropped/70393354-d898-41c0-9734-7d0b9cf96d8a-17.jpg?height=544&width=606&top_left_y=1720&top_left_x=490)
FIG. 10. Prediction of the latent $X^{\prime}$ process for the Anadenanthera Colubrina data. The colors represent the frequency at which the posterior draw for the process fell in each cell. The colors can be interpreted as the posterior probability that there is at least one unobserved occurrence in the cell.

of the posterior mode is recommended in addition to prior information. A traceplot of the model parameters can be found in the Supplementary Material A, Section 4 (Moreira and Gamerman (2022b)).

Nevertheless, the forms for $q(\cdot)$ and $p(\cdot)$ can be changed to be any version of binary classifier. This can be treated analytically by defining the processes $X, X^{\prime}$ and $U$ to be thinning of a homogeneous process with constant rate $\lambda^{*}$. Computationally, this only implies a change in the Gibbs sampling step of $\beta$ and $\delta$. Thus, if $q(\cdot)$ is changed to a probit link, the sampling scheme will have to be changed into one which samples from a probit regression. An argument can be made toward a nonparametric binary classifier, but the theoretical foundations of such a model become less clear and identification problems between $q(\cdot)$ and $\lambda^{*}$ may arise.

Currently, presence-only methodology development has somewhat stabilized. Instead, multiple-species modeling (Fithian et al. (2015), Gelfand and Schliep (2018), Renner, Louvrier and Gimenez (2019), Renner et al. (2015), Shirota, Gelfand and Banerjee (2019)) and joining presence-only with presence-absence data (Gelfand and Shirota (2019)) have gained attention. Furthermore, as Renner et al. (2015) and Renner and Warton (2013) discuss, it is not uncommon to face a situation where there needs to be further spatial structure in the model beyond that which is present in the covariates. That can be achieved by the introduction of a latent term in the regression predictor governed by a Gaussian process, for example. That is an instance of the advantage of a probabilistic model over an information theory, one like MaxEnt.

Despite the imperfect sampling, presence-only data is available for a large number of species where a systematic survey is not. The development of new methodology for presenceonly data does not mean that systematic surveys are not encouraged. It only means that there is available information that should not be discarded. The analysis of multiple species is an important step forward. The borrowing of information between different species might mean that the sensibility to lackluster information of presence-only methods can be diminished.

The greatest challenge of a presence-only dataset is still the fact that it is lacking in information. The choice of informative covariates for both intensity and observability sets is instrumental to the quality of estimation and prediction. There are good and reliable repositories for environmental and geographical variables in a fine lattice, which is commonly useful for intensity covariates. However, observability covariates are usually not so easily acquired. Obtaining distance to urban centers, for example, may require some expertise with GIS software. Still, it is advised in any presence-only study to form a close partnership with an expert who should be consulted regarding the choice of covariates.

Acknowledgments. The authors thank Dr. Carolina Levis for introducing us to the area of presence-only in ecology, Dr. Fernando Figueiredo for useful references in the area, Professor David Warton for valuable insights and guidance in acquiring data, Dr. Guilherme Mazzochini for his aid in acquiring data and for sharing his ecological expertise, Professor Flávio Gonçalves for his technical comments and Professor Ian Renner for making the Eucalyptus data available. They are grateful for the added value provided by the reviewers. Finally, the authors thank the graduate program in statistics at the Federal University of Rio de Janeiro (UFRJ). This paper is based on the doctoral thesis of the first author developed under the supervision of the second author. The authors thank the hospitality provided by the Department of Statistics, UFMG.

Funding. The first author was funded by research grants from CAPES, Brazil. The second author was funded by grants from CNPq and FAPERJ from Brazil. The support from all these research supporting agencies is gratefully acknowledged by the authors.

## SUPPLEMENTARY MATERIAL

A-Supplementary derivations and results (DOI: 10.1214/21-AOAS1569SUPPA; .pdf). This file includes the proof of Proposition 1, the derivation of presence probability in small areas, additional simulated results and more results from the real data analysis.

B-Code for the package and data (DOI: 10.1214/21-AOAS1569SUPPB; .zip). This part consists of a compressed file which contains the code for the package and for simulating data. The results from Section 6 can be replicated with this code. The file also includes the data analyzed in Section 7.

## REFERENCES

Adams, R. P., Murray, I. and MacKay, D. J. C. (2009). Tractable Nonparametric Bayesian Inference in Poisson Processes with Gaussian Process Intensities. In Proceedings of the 26 th Annual International Conference on Machine Learning. ICML'09 9-16. Association for Computing Machinery, New York, NY, USA.
Baddeley, A., Rubak, E. and Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with $R$. CRC Press/CRC Press, London.
Byrne, S. (2016). A note on the use of empirical AUC for evaluating probabilistic forecasts. Electron. J. Stat. 10 380-393. MR3466187 https://doi.org/10.1214/16-EJS1109
Cressie, N. A. C. (1993). Spatial Point Patterns. Wiley, New York.
Diggle, P. J. (2014). Statistical Analysis of Spatial and Spatio-temporal Point Patterns, 3rd ed. Monographs on Statistics and Applied Probability 128. CRC Press, Boca Raton, FL. MR3113855
Diggle, P. J., Menezes, R. and Su, T. (2010). Geostatistical inference under preferential sampling. J. R. Stat. Soc. Ser. C. Appl. Stat. 59 191-232. MR2744471 https://doi.org/10.1111/j.1467-9876.2009.00701.x
Dorazio, R. M. (2014). Accounting for imperfect detection and survey bias in statistical analysis of presenceonly data. Glob. Ecol. Biogeogr. 23 1472-1484.
Elith, J. and Leathwick, J. (2007). Predicting species distributions from museum and herbarium records using multiresponse models fitted with multivariate adaptive regression splines. Diversity and Distributions $\mathbf{1 3}$ 265-275.
Elith, J., Graham, C., Valavi, R., Abegg, M., Bruce, C., Ford, A., Guisan, A., Hijmans, R., Huettmann, F. et al. (2020). Presence-only and Presence-absence Data for Comparing Species Distribution Modeling Methods. Biodiversity Informatics 1569-80.
Fielding, A. H. and Bell, J. F. (1997). A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24 38-49.
Fithian, W. and Hastie, T. (2013). Finite-sample equivalence in statistical models for presence-only data. Ann. Appl. Stat. 7 1917-1939. MR3161707 https://doi.org/10.1214/13-AOAS667
Fithian, W., Elith, J., Hastie, T. and Keith, D. A. (2015). Bias correction in species distribution models: Pooling survey and collection data for multiple species. Methods Ecol. Evol. 6 424-438. https://doi.org/10. 1111/2041-210X. 12242
Fletcher Jr., R. J., Hefley, T. J., Robertson, E. P., Zuckerberg, B., McCleery, R. A. and Dorazio, R. M. (2019). A practical guide for combining data to model species distributions. Ecology 100 e02710.
Gamerman, D. and Lopes, H. F. (2006). Markov Chain Monte Carlo: Stochastic simulation for Bayesian inference, 2nd ed. Texts in Statistical Science Series. CRC Press/CRC, Boca Raton, FL. MR2260716
Gelfand, A. E. and Schliep, E. M. (2018). Bayesian Inference and Computing for Spatial Point Patterns. NSF-CBMS Regional Conference Series in Probability and Statistics 10. IMS, Beachwood, OH. MR3890052
Gelfand, A. E. and Shirota, S. (2019). Preferential sampling for presence/absence data and for fusion of presence/absence data with presence-only data. Ecol. Monogr. 89 e 01372.
Gelman, A., Simpson, D. and Betancourt, M. (2017). The prior can generally only be understood in the context of the likelihood.
Gonçalves, F. B. and Gamerman, D. (2018). Exact Bayesian inference in spatiotemporal Cox processes driven by multivariate Gaussian processes. J. R. Stat. Soc. Ser. B. Stat. Methodol. 80 157-175. MR3744716 https://doi.org/10.1111/rssb. 12237
Hastie, T. and Fithian, W. (2013). Inference from presence-only data; the ongoing controversy. Ecography 36 864-867. https://doi.org/10.1111/j.1600-0587.2013.00321.x
Hefley, T. J., Tyre, A. J., Baasch, D. M. and Blankenship, E. E. (2013). Nondetection sampling bias in marked presence-only data. Ecol. Evol. 3 5225-5236.

Hefley, T. J., Baasch, D. M., Tyre, A. J. and Blankenship, E. E. (2015). Use of opportunistic sightings and expert knowledge to predict and compare Whooping Crane stopover habitat. Conserv. Biol. 29 1337-1346.
Journé, V., Barnagadd, J.-Y., Bernard, C., Crochet, P.-A. and Morin, X. (2020). Correlative climatic niche models predict real and virtual species distributions equally well. Ecology $\mathbf{1 0 1}$ e02912. https://doi.org/10. 1002/ecy. 2912
Lewis, P. A. W. and Shedler, G. S. (1979). Simulation of nonhomogeneous Poisson processes by thinning. Nav. Res. Logist. Q. 26 403-413. MR0546120 https://doi.org/10.1002/nav. 3800260304
Little, R. J. A. and Rubin, D. (2014). Statistical Analysis with Missing Data, 2nd ed. Wiley, New York.
Mazzochini, G. G., Fonseca, C. R., Costa, G. C., Santos, R. M., Oliveira-Filho, A. T. and Ganade, G. (2019). Plant phylogenetic diversity stabilizes large-scale ecosystem productivity. Glob. Ecol. Biogeogr. 28 1430-1439.
Moreira, G. A. (2021). bayesPO: Bayesian Inference for Presence-Only Data. R package version 0.3.1.
Moreira, G. A. and Gamerman, D. (2022a). Supplement to "Analysis of presence-only data via exact Bayes, with model and effects identification." https://doi.org/10.1214/21-AOAS1569SUPPA
Moreira, G. A. and Gamerman, D. (2022b). Supplement to "Analysis of presence-only data via exact Bayes, with model and effects identification." https://doi.org/10.1214/21-AOAS1569SUPPB
Oliveira-Filho, A. T. (2017). NeoTropTree, Arborea flora of the Neotropical Region: A Database involving biogeography, diversity and consevation. Universidade Federal de Minas Gerais. Available at http://www.neotroptree.info.
Pearce, J. L. and Boyce, M. S. (2006). Modelling distribution and abundance with presence-only data. J. Appl. Ecol. 43405-412.
Phillips, S. J., Anderson, R. P. and Schapire, R. E. (2006). Maximum entropy modeling of species geographic distributions. Ecol. Model. 190 231-259.
Phillips, S. J., DudíK, M. and Schapire, R. E. (2004). A Maximum Entropy Approach to Species Distribution Modeling. In Proceedings of the Twenty-first International Conference on Machine Learning. ICML'04 83. ACM, New York, NY, USA.

Phillips, S. J., Dudík, M., Elith, J., Graham, C. H., Lehmann, A., Leathwick, J. and Ferrier, S. (2009). Sample selection bias and presence-only distribution models: Implications for background and pseudoabsence data. Ecol. Appl. 19 181-197.
Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference for logistic models using PólyaGamma latent variables. J. Amer. Statist. Assoc. 108 1339-1349. MR3174712 https://doi.org/10.1080/ 01621459.2013 .829001

Renner, I. W., Louvrier, J. and Gimenez, O. (2019). Combining multiple data sources in species distribution models while accounting for spatial dependence and overfitting with combined penalised likelihood maximisation. BioRxiv.
Renner, I. W. and Warton, D. I. (2013). Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics 69 274-281. MR3058074 https://doi.org/10.1111/ j.1541-0420.2012.01824.x

Renner, I. W., Elith, J., Baddeley, A., Fithian, W., Hastie, T., Phillips, S. J., Popovic, G. and Warton, D. I. (2015). Point process models for presence-only analysis. Methods Ecol. Evol. 6 366-379.
Royle, J. A., Chandler, R. B., Yackulic, C. and Nichols, J. D. (2012). Likelihood analysis of species occurrence probability from presence-only data for modelling species distributions. Methods Ecol. Evol. 3 545-554.
Rubin, D. B. (1976). Inference and missing data. Biometrika 63 581-592. MR0455196 https://doi.org/10.1093/ biomet/63.3.581
Shirota, S., Gelfand, A. E. and Banerjee, S. (2019). Spatial joint species distribution modeling using Dirichlet processes. Statist. Sinica 29 1127-1154. MR3932512
Warton, D. I. and Shepherd, L. C. (2010). Poisson point process models solve the "pseudo-absence problem" for presence-only data in ecology. Ann. Appl. Stat. 4 1383-1402. MR2758333 https://doi.org/10.1214/ 10-AOAS331


[^0]:    Received November 2020; revised October 2021.
    Key words and phrases. Point process, presence-only, spatial statistics, Bayesian analysis, species distribution model.

