# Classifying the OXPHOS status of skeletal muscle fibres
An R package for the Bayesian classification of myofibres according to OXPHOS protein expression profiles.

__PLEASE READ BEFORE ATTEMPTING TO USE PACAKGE__

## Update R 
Before continuing it is suggested that R is updated, this is imported to be able
to use Rtools and to configure a C++ toolchain. The package was built
on R version 4.3.1 and so at least this version is suggested. The [article](https://www.listendata.com/2015/08/how-to-update-r-software.html)
contains instructions on how to do this. Note that if you are using Rstudio 
(Posit) you may need to set the version of R being used, [this](https://support.posit.co/hc/en-us/articles/200486138-Changing-R-versions-for-the-RStudio-Desktop-IDE) post on the Posit website should help with this. 

For windows users updating R can easily done within R itself, using the `installr` package. 
```{R echo=TRUE}
install.packages("installr")
library("installr")
updateR()
```

## C++ toolchain
C++ toolchain can be created through R, using some helpful packages along the way. 

For Windows users, this is relatively easy and can done through `Rtools`, which
can be installed [here](https://cran.r-project.org/bin/windows/Rtools/). 
Alternatively it can be installed via the `installr::installRtools` function.
If the package `installr` is not installed this should be done first. 
```{R echo=TRUE}
# install.packages("installr")
library("installr")
install.Rtools()
```

For Mac a C++ toolchain can be created using the `macrtools` package, a guideline for
its installation and use can be found 
[here](https://mac.thecoatlessprofessor.com/macrtools/). The toolchain can then
be configured following instructions [here](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac).

For Linux a toolchain can be constructed and configured following the
instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

The RStan GitHub repository should be an up to date place for information and
guides for creating a C++ toolchain and is where the information presented
here is from. The [Rstan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) page will
hopefully guide you if any issues arise. The page is particularly
relevant as this package, `analysis2Dmito`, uses `rstan` to execute inference.

## Dependencies and installation
Once a suitable C++ toolchain has been created and R is updated (if required)
the package dependencies can be installed. Here we also install `devtools`, this
is a package which allows you install packages directly from GitHub. The following
code installs `devtools` as well as dependencies and loads them into your environment.

```{r echo=TRUE}
install.packages(c("data.table", "dplyr", "readr", "tidyr", "plyr", "devtools"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("devtools")
```
Once the dependencies are installed successfully the package itself can be installed. Downloadingthe package may take a couple minutes and produce a series of warning 
messages This is to be expected and as long as the package has a non-zero exit 
status i.e. the package is downloaded successfully, this should not be a problem. 
```{r echo=TRUE}
install_github("jordanbchilds/analysis2Dmito", quiet = TRUE, upgrade="never")
library("analysis2Dmito")
```

## Example Script

There is an example script in the main folder of the repo, `example_analysis.R`.
The output from each model fit is saved and split across five files; posterior 
draws, prior draws, posterior predictive, prior predictive and fibre 
classifications, with the suffixes; "POST", "PRIOR", "POSTPRED", "PRIORPRED" and 
"CLASSIF" respectively. The script also plots three key aspects of the output; 
MCMC output, prior and posterior comparisons and fibre classifications. The 
pipeline used in the example script and plots produced are discussed in more 
detail in the following sections, the scripts itself has some minimal 
comments throughout. 

# Pipeline

## Data

Getting data into the correct form is crucial to be able to use the inference 
and plotting functions which are part of this package. The particular format of 
data used in this package has been chosen for consistency with historical 
datasets generated in our group. The package comes with an example dataset which 
is already in the correct form, so to start we will inspect this. The function 
`get_exampleData` loads the dataset. 

```{r echo=TRUE include=TRUE}
exampleData = get_exampleData()
head(exampleData, 3)
```

This should show the first three rows of the example dataset, here we have 
called it `exampleData`. A description of each column is given in the table 
below.

| Variable Name | Description |
| ------------- | ----------- |
| `value` | A numerical value representing the expression level (average pixel intensity) for a particular fibre within a sample for a particular channel. |
| `sampleID` | A unique string identifying the sample from which an observation was made e.g. `C01`, `C02`, `P01`, `P02`, ... |
|`fibreID` | An integer identifying a single fibre in its sample. This identifier is not unique throughout the whole dataset, it is only unique within the sample it comes from. It should be consistent across channels for the sample. |
| `sbj_type` | A string identifying which samples are control subjects and which are patients. Control samples must be labelled as "control" and patient samples labelled as "patient". |
| `channel` | A string labeling the channel or protein measured. |


__Any dataset used with the functions in this package must have at least these five columns__, other columns are allowed but are not necessary. The data is in a format 
commonly referred to as long form, a helpful function to be able to get data in 
this form is the `tidyr::pivot_longer` function (from the `tidyr` package). For 
information on the `pivot_longer` function and examples see this [blog post](https://tidyr.tidyverse.org/reference/pivot_longer.html).

## Transforming the data

Before moving on to inference it is advisable to explore the data. For good 
results the healthy control data should show a linear relationship and constant 
variance, two assumptions made during linear modelling. The healthy and 
deficient fibres should also appear to be distinct populations although the 
shape of the deficient population is less important. This may be achieved by 
data transformation. The usual transformation seen in literature is the log 
transform, this works well for many datasets although it may not be the best 
for all datasets. For example, data that shows a strong V-shape, where one 
branch of the V is healthy and the other deficient, a log-log transformation 
(logging the expression values twice) may separate the two branches better than 
a log transformation. The log-log transformations requires all raw expression 
values to be greater than 1.0, if this is not the case simply adding 1.0 to all 
expression values would not affect classification. The code snippet below plots 
each patient data with the control data, for the raw, log transformed and log-
log transformed data, for the example dataset.

```{r echo=TRUE include=TRUE}
# the 2Dmito plot x-axis - known for your dataset
mitochan = "VDAC"
# return all channels which are not mitochan
channels = unique( grep(mitochan, exampleData$channel, value=TRUE, invert=TRUE) )

# extract control sample IDs
ctrlIDs = unique( exampleData[exampleData$sbj_type=="control", "sampleID"] )
# extract patient sample IDs
patIDs = unique( exampleData[exampleData$sbj_type=="patient", "sampleID"] )
patIDs = sort(patIDs)

# plot patient data and control data
for( chan in channels ){
  xDat_ctrl = exampleData[exampleData$sbj_type=="control" & exampleData$channel==mitochan, "value"]
  yDat_ctrl = exampleData[exampleData$sbj_type=="control" & exampleData$channel==chan, "value"]

  for( pat in patIDs ){
    png(paste0("ExploratoryAnalysis_",chan,"_",pat,".png"), type = "cairo-png", width=1920, height=1080,pointsize=36)
    op = par(mfrow=c(1,3))
    xDat_pat = exampleData[exampleData$sampleID==pat & exampleData$channel==mitochan, "value"]
    yDat_pat = exampleData[exampleData$sampleID==pat & exampleData$channel==chan, "value"]
    plot( xDat_ctrl, yDat_ctrl, pch=20, col="black",
          xlab=mitochan, ylab=chan, main=pat,
          xlim=range(c(xDat_ctrl, xDat_pat)), ylim=range(c(yDat_ctrl, yDat_pat)) )
    points( xDat_pat, yDat_pat, pch=20, col="green")
    
    plot( log(xDat_ctrl), log(yDat_ctrl),  pch=20, col="black",
          xlab=mitochan, ylab=chan, main=paste("log", pat),
          xlim=log(range(c(xDat_ctrl, xDat_pat))), ylim=log(range(c(yDat_ctrl, yDat_pat))) )
    points( log(xDat_pat), log(yDat_pat), pch=20, col="green")
    
    plot( log(log(xDat_ctrl)), log(log(yDat_ctrl)),  pch=20, col="black",
          xlab=mitochan, ylab=chan, main=paste("log log", pat),
          xlim=log(log(range(c(xDat_ctrl, xDat_pat)))), ylim=log(log(range(c(yDat_ctrl, yDat_pat)))) )
    points( log(log(xDat_pat)), log(log(yDat_pat)), pch=20, col="green")
    par(op)
    dev.off()
  }
}
```
A comparison of the untransformed data and the log and log-log transformations are shown in the figure below, for a specific channel and patient in the example dataset. Fibres from healthy control subjects are shown in black and patient fibres (both healthy and deficient) are shown in green. The log and log-log transforms increase the difference between healthy and deficient populations fibres with the lowest expression levels but reduces the difference in fibres with higher expression levels. Depending on the data, one may be more appropriate than the other. Here we will choose the log transform and the data is transformed at the beginning of the next code snippet. 

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/transformations_ex.png?raw=true)

## Choosing prior beliefs

Prior beliefs can be constructed based on the control data available. The
inference function within the package, `stan_inference`, does this automatically
based upon the data that it is given. Detail about how the prior beliefs are
chosen are described [here](link to paper). In summary, linear models are fitted
to each control subject data independently and the mean value of their
slopes and intercepts are used as the expected values of our prior beliefs. 

## Fit the model

If the data is transformed so that healthy fibres show a linear relationship and is in the correct form then we can now fit the model. Before doing this we must pass the data through the `getData_mats` function which organises the data into matrices to be passed to `stan_inference`. The following code snippet should run the inference for the first patient in `patIDs` and the channel, `chan`, defined above. 

```{r echo=TRUE}

pat = "P01"

dataMats = getData_mats(data=exampleData, 
                        ctrlID=ctrlIDs,
                        channels=c(mitochan, chan),
                        pts=pat, 
                        getIndex=TRUE)

output = stan_inference(dataMats)
```

## Understanding inference output

The inference function outputs several things in a list. Two of which are the `POST` and `PRIOR` matrices, containing prior and posterior draws for each parameter in the model in the form of a matrix, where each column is a different parameter. The parameters outputted from the inference can be seen by the column names of the two matrices (the two should match exactly). 

```{r echo=TRUE}
colnames(output$POST)
colnames(output$PRIOR)
```

We see here that there are several `m` and `c` parameters. This is because each control sample and patient sample has its own slope and intercept. The slopes for the control samples are given first, as such `m[1]`, `m[2]`, `m[3]` and `m[4]` are the slopes for the control samples if there are four controls, these will be in the order of `ctrlIDs` passed to the `getData_mats` function. Passing a single patient sample to the model will result in an additional slope, `m[5]`. The same is true for the intercept parameters `c[1]`,..., `c[5]`. The parameters `m_pred` and `c_pred` seen in both the prior and posterior matrices are the predictive distribution of the slope and intercept respectively. These can also be seen as the our belief about the slope and intercepts for all possible samples, given the data we have. They are the layer in the hierarchy from which value of individual slopes and intercepts for a specific sampler are drawn. 

Comparing prior and posterior distributions is useful as this shows how our parameter beliefs have been updated. Below is a snippet for how to inspect the prior and posterior densities for the proportion of deficiency (the proportion of deficient fibres within a sample). Note that although the smoothed density estimates used to visualise the distribution may extend outside the range [0,1], inspection of the output itself shows that no value is outside this range.

```{r echo=TRUE}
dPost = density(output$POST[,"probdiff"])
dPrior = density( output$PRIOR[,"probdiff"] )
xlims = range(c(dPrior$x, dPost$x))
ylims = range(c(dPrior$y, dPost$y))

plot(dPrior, col=alphaPink(0.7), lwd=2,
    xlim=xlims, ylim=ylims,
     main="Proportion of Deficiency", xlab="")
lines( dPost, col=alphaGreen(0.7), lwd=2)
legend("topright", legend=c("post", "prior"), lty=1, col=c(alphaGreen(1.0), alphaPink(1.0)))
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/pi_postprior_ex.png?raw=true)

The `POSTPRED` and `PRIORPRED` matrices in the output list contain prior and posterior predictive interval - marginalised over parameter uncertainty - for the patient sample passed to the model. The first column in this is called `mitochan` and are the values on the x-axis for the prediction. The remaining six columns are the 2.5\%, 50\% and 97.5\% quantiles of the prediction for the healthy patient fibres and the deficient patient fibres. This is found by calculating the predictive distribution of the protein expression at each value stored in the `mitochan` column and then calculating its quantiles. 

```{r echo=TRUE}
plot(dataMats$pts, pch=20, col=alphaBlack(0.2),
     xlab=paste0("log(", mitochan ,")"), ylab="log(MTCO1)", 
     main="Posterior predictive linear regression")
lines(output$POSTPRED[,"mitochan"], output$POSTPRED[,"medNorm"], 
      lty=1, col=alphaPink(0.9))
lines(output$POSTPRED[,"mitochan"], output$POSTPRED[,"lwrNorm"], 
      lty=2, col=alphaPink(0.9))
lines(output$POSTPRED[,"mitochan"], output$POSTPRED[,"uprNorm"], 
      lty=2, col=alphaPink(0.9))
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/postpred_ex.png?raw=true)

The last item in the list is called `CLASSIF` and is matrix of every posterior classification for every patient fibre. Each column of the matrix includes our posterior belief about the classification of a fibre and the many repeated measures include our uncertainty about that classification. We can find the probability that an individual fibre is deficient by calculating the mean of the column. This can be done using the `colMeans` function, as seen below. The classifications themselves can be plotted using these probabilities and the `classcols` function which converts the probabilities to a colour on a scale between blue and red, where a deficiency probability of 0.0 is strongly blue and a probability of 1.0 is strongly red.  

```{r echo=TRUE}
def_prob = colMeans(output$CLASSIF)
plot(dataMats$pts, pch=20, col=classcols(def_prob),
     xlab=paste0("log(", mitochan ,")"), ylab="log(MTCO1)", 
     main="Posterior fibre classifications")
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/classif_ex.png?raw=true)

## Plotting model output

There are a few plotting functions in the package which are designed to make things easier. One of which is a plot which was not seen in the previous section, it plots the MCMC output of the inference. Checking the MCMC output is important to make sure the inference has worked and that there are no major issues with it. There are a number of metrics which can be used to check for problems with the inference, see [this blog](https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics#:~:text=The%20simplest%20way%20to%20diagnose,results%20on%20all%20the%20chunks.) for more information. 

The diagnostic plotter, `MCMCplot`, in this package plots three things per parameter; trace plot, autocorrelation plot (ACF) and a density estimate of the posterior (and prior if given). The trace plot plots each posterior draw in the order they were generated and the autocorrelation shows the correlation between draws for a range of lags. The density estimate is like that which has been already seen, a kernel density estimate of distribution given a number of draws from it. Using these plots it is important to check some key elements of the output.

1. Has the chain reached a stationary distribution?
2. Are the draws (almost) independent?

To check that the chain has reached a stationary distribution, check the trace plots. An ideal chain would look like a thick straight black line with whiskers coming off it. If the chain is curved at beginning and then flattens out, increase the `warmup` parameter in the `stan_inference` function. A chain that has not reached a stationary distribution will look like a single jagged line (similar to a stock price), again increasing `warmup` may solve the problem. If the chain is jumping between two or more stationary distributions this will be clear as the density plot will be bimodal. A harder issue to solve, however if there are strong reasons to believe that the parameter should be the lower of higher peak then changing the prior to reflect should help. 

Autocorrelated samples can be quickly diagnosed in the ACF plot. A sample of uncorrelated draws will have no slope in the autocorrelations, showing a spike at a lag of zero and then small random correlations there after. Generally a small amount of correlation in successive samples is fine, when samples are extremely autocorrelated the resulting model can suffer. The amount of correlation can be reduced by increasing the `MCMCthin` parameter in the `inference` function. Below is an example of the `MCMCplot` use and results for the same parameter before and after increasing thinning. A larger amount of thinning resulted in the correlations betwen successive draws being massively decreased. Note: increasing thinning will increase inference time. 

```{r echo=TRUE}
MCMCplot(post=output$POST, 
         prior=output$PRIOR,
         nRow=1)
```
An example of an MCMCPlot with autocorrelated samples (increase `MCMCthin` in `inference` function to improve this):
![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/mcmc_m1_ex.png?raw=true)

An example of an MCMCPlot with no autocorrelation problems:
![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/mcmc_m1_thin_ex.png?raw=true)


To be able to visualise the prior and posterior densities for all variables we can use the `postPlot` function. The function plots all the priors, posteriors, and predictives (for hierarchical distributions) as well as the model classifications and posterior predictive for the healthy patient data. Unlike `MCMCplot` this does not automatically create a plotting environment, so we have to define this before calling the function. For this inference there are multiple slopes and intercepts but one is of particular interest, those for the patient sample, as such this one is in a solid colour while the others are more transulcent. 

```{r echo=TRUE}
op = par(mfrow=c(3,3)) # plotting grid: 3x3

postPlot(post=output$POST,
         prior=output$PRIOR,
         postpred=output$POSTPRED,
         dataMats=dataMats,
         classifs=def_prob,
         var.names=c("mu_m", "tau_m", "m", "mu_c", "tau_c", "c", "probdiff", "tau_norm"), 
         mitoPlot_xlab="log(VDAC)",
         mitoPlot_ylab="log(MTCO1)"
         )

par(op) # end plotting grid 
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/postPlot_ex.png?raw=true)

In the above output we see that little is learnt about the parameters `mu_m`, `tau_m`, `mu_c` and `tau_c`. This is due to the strong prior beliefs placed on them, they do, however, still allow for a range of slopes and intercepts to be fitted to different datasets as can be seen in their respective posterior distributionss. Our prior beliefs about the model precision, `tau_norm`, are also strong and so this has only been updated marginally. Our prior beliefs about the proportion of deficient samples were completely uniformed, being flat on the interval [0,1], and have been updated a lot. The prior beliefs placed on parameters may affect the resulting model fitted to the data.

To just see the fibre classification, the `classif_plot` function can be used. The function creates a plot similar to that which was seen before but allows you to plot the posterior predictive without tediously writing the `lines` function. The function also automatically plots the control data in a translucent black colour. 

```{r echo=TRUE}
classif_plot(dataMats=dataMats,
             classifs=def_prob,
             postpred=output$POSTPRED,
             xlab="log(VDAC)",
             ylab="log(MTCO1)")
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/classif_plot_ex.png?raw=true)

# Checking model output
The inference is run by default with a 10,000 iteration burn-in period, it is not
guaranteed that this is enough to effectively explore the parameter space (find 
the posterior distribution). If the output of the analysis appears to be wrong 
this could be the problem. There are a few things which can be done to solve this
issue. 

## Assessing output  

### Multiple chains
Running multiple chains is a good idea in general. The `stan_inference` function
allows this by passing arguments to the sampling function in `rstan`. 
using the `chains` argument and they can be split across cores using 
the `cores` argument, following the style of `rstan` sampling functions. The 
following code snippet will execute five chains, each with a 100,000 iteration 
burn-in period and produce 2,000 draws from the posterior and split the chains
across two cores. 

```{R echo=TRUE}
output = stan_inference(dataMats, chains=5, ncores=2)
```

Once multiple chains have run, we can check have reached the same posterior
distribution using the `rstan::Rhat` function. It is suggested that if the value
of rhat is greater than 1.05 then the chains have NOT converged to the same 
distribution. The function requires a matrix of the draws for a specific parameter. 
More information about the `Rhat` (and other metrics for convergence) can be
found in the function documentation. The code snippet below loops through the 
parameters and calculates its rhat value, storing them in a vector called `rhat_vec`.

```{R echo=TRUE}
post = output$POST
chains = 5
mcmcOut = 2000
rhat_vec = vector("numeric", length=ncol(post))
names(rhat_vec) = colnames(post)
for( param in colnames(post) ){
  param_mat = matrix(post[,param], nrow=mcmcOut, ncol=chains, byrow=FALSE)
  rhat_vec[param] = Rhat(param_mat)
}
print(rhat_vec)
```
If there are no signs that the chains have reached different posterior distributions,
then the chains can be combined and all draws can be considered to be from the 
same distribution. 

### Effective sample size

The effective sample size (ESS) is an estimate of the number of __uncorrelated__
draws from the posterior distribution that the inference has produced. Ideally,
this number would be close to the number of posterior draws as possible (or even
larger). STAN usually produces a large effective sample size but if the chain has 
not converged this will result in a low ESS. 

There are no exact rules to follow as to what is an acceptable ESS. However, here
we can use the large amount of datasets to guide us. By choosing a dataset whose
inference shows no signs of non-convergence, we have a decent idea of what a good 
ESS looks like for the model. From experience inference that has not converged 
have drastically different ESS to those which have, resulting in a noticeable 
divide if the ESS is calculated for all output. 

Usually the ESS is calculated for a single parameter but we can also calculate 
multivariate ESS. This gives an estimated effective sample size for the whole 
inference. Multivariate ESS can be inspected and though of much the same way 
single variate ESS can and a poor Multivariate ESS will stand out amongst a set
of good ones. 

The following snippet uses the `coda` and `mcmcse` packages to calculate the 
single and multivariate ESS for an output from the `stan_inference` function.

```{R echo=TRUE}
install.packages("coda")
install.packages("mcmcse")
library("coda")
library("mcmcse")

(ess = coda::effectiveSize( output$POST ))
(multiESS = mcmcse::multiESS( output$POST ))
```
For the example dataset provided with the package an acceptable multivariate ESS
was at least 95\% of the number of posterior draws and an acceptable univariate
ESS was approximately greater than 70\% of posterior draws. 

## Improving output

### Increasing burn-in

If testing has shown that there are lack of convergence in your model output, be
it from multiple chains or just one, increasing the burn-in period should be the
first step to try. Increasing burn-in gives the inference scheme a longer time 
to explore the parameter space and find the posterior distribution. By default,
when using the `stan_inference` function, the burn-in period is 20,000. It can 
be increased by setting the `warmup` argument of the function, if this is done the
`iter` argument must also be updated as this is the total number of interations 
and therefore must always be greater than `warmup`. 

```{R echo=TRUE}
output = stan_inference(dataMats, warmup=100000, iter=105000)
```

### Changing prior beliefs

If the inference is still not converging to a posterior distribution it is 
possible to change the prior beliefs to make them more precise. 
It is also possible to use truncated prior densities to remove the possibility 
of some values even being suggested. This is not guaranteed to work and is not 
advisable unless there is good reason to do so.




