# Classifying the OXPHOS status of skeletal muscle fibres
An R package for the Bayesian classification of myofibres according to OXPHOS protein expression profiles.

## Required Software
The package uses the statistical computing software JAGS which must be installed before use. This cannot be done through R or RStudio/Posit but there are many online resources to help install JAGS. It can be installed directly though [sourceforge.net](https://sourceforge.net/projects/mcmc-jags/files/). A guide to installing R, RStudio, JAGS and R packages can be found [here](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119287995.app1). 

## Installation
Once JAGS has been installed on your machine, you can install this package and its dependencies within R or Rstudio. 
```{r echo=TRUE}
# install dependencies
install.packages(c("data.table", "dplyr", "readr", "tidyr", "rjags"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("rjags")

install.packages("devtools")
library("devtools")
devtools::install_github("jordanbchilds/analysis2Dmito")

library("analysis2Dmito")
```

## An Example Classification Pipe

### Getting example data
Getting data into the correct form is crucial to be able to use the inference and plotting functions which are part of this package. The particular format of data used in this package have been chosen for consistency with historical datasets generated in our group.  The package comes with an example dataset which is already in the correct form, so to start we will look at this. The function `get_exampleData` loads the example dataset. 
```{r echo=TRUE include=TRUE}
exampleData = get_exampleData()
head(exampleData)
```
This should show the first six rows of the example dataset, here we have called it `exampleData`. The description of the columns are given in the table below.
| Variable Name | Description |
| ------------- | ----------- |
| `sampleID` | A (unique) string identifying of the sample e.g. `C01`, `C02`, `P01`, `P02`, ... |
|`fibreID` | An integer identifying of a single fibre from that sample. This identifier is not unique throughout the whole dataset, it is only unique within the sample it comes from. |
| `sbj_type` | A string identifying which samples are control subjects and which are patients. Control samples must be labelled as "control" and patient samples labelled as "patient". |
| `channel` | A string labelling the channel or protein measured.  |
| `value` | A numerical value representing the expression level (average pixel intensity) for this particular fibre on this particular channel. |

__Any dataset used with the functions in this package must have at least these five columns__, other columns are allowed but are not necessary. The data is in a format commonly referred to as long form, a helpful function to be able to get data in this form is the `tidyr::pivot_longer` function (from the `tidyr` package). For information on the `pivot_longer` function and examples see this [blog post](https://tidyr.tidyverse.org/reference/pivot_longer.html).

### Explore the data

#### Transforming the data

Before moving on to inference it is advisable to explore the data. For good results the healthy control data should show a linear relationship and constant variance, two assumptions made during linear modelling. The healthy and deficient fibres should also appear to be distinct populations although the shape of the deficient population is less important. This may be achieved by data transformation. The usual transformation seen in literature is the log transformation, this works well for many datasets although it may not be the best for all datasets. For example, data that shows a strong V-shape, where one branch of the V is healthy and the other deficient, a log-log transformation (logging the expression values twice) may separate the two branches better than a log transformation. The log-log transformations requires all raw expression values to be greater than 1.0, if this is not the case simply adding 1.0 to all expression values would not affect classification. The code snippet below plots each patient data with the control data, for the raw, log transformed and log-log transformed data, for the example dataset.

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
A comparison of the log and log-log transformations are shown in the figure below, for a specific channel and patient in the example dataset. Healthy control fibres are shown in black and patient fibres (both healthy and deficient) are shown in green. The log-log transform increases the difference between healthy and deficient populations in lowest expression fibres but reduces the difference in fibres with higher expressions. Depending on the data one may be more appropriate than the other. Here we will choose the log transform and the data is transformed at the beginning of the next code snippet. 

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/transformations_ex.png?raw=true)

#### Choosing prior beliefs

Once a transformation has been chosen the expression data can be changed and prior beliefs about parameters can be considered. By fitting linear models to the control samples individually we know what values of model parameters might be likely and can gain an idea the expected value of model parameters. For ease the linear models can be fit in a frequentist setting. The code snippet below fits a linear model to each control sample and each protein, and saves relevant output in three matrices; `slopes`, `intercepts`, and `precisions`.

```{r echo=TRUE}
exampleData$value = log( exampleData$value )

slopes = matrix(NA, nrow=length(channels), ncol=length(ctrlIDs))
rownames(slopes) = channels
colnames(slopes) = ctrlIDs

intercepts = slopes # defines an empty matrix with row and col names as wanted
precisions = slopes

for(chan in channels){
  for(crl in ctrlIDs){
    x = exampleData[exampleData$sampleID==crl & exampleData$channel==mitochan, "value"]
    y = exampleData[exampleData$sampleID==crl & exampleData$channel==chan, "value"]
    df = data.frame(mitochan=x, chan=y)
    
    lnmod = lm(chan~mitochan, data=df)
    
    slopes[chan, crl] = lnmod$coefficients["mitochan"]
    intercepts[chan, crl] = lnmod$coefficients["(Intercept)"]
    precisions[chan, crl] = 1 / summary(lnmod)$sigma^2
  }
}
```

The matrices; `slopes`, `intercepts` and `precisions`, contain the frequentist estimates of the respective parameters for each control sample in each channel. One way to specify prior beliefs would be to set the expected value of a parameter to the mean of the appropriate values from the frequentist fits and give the parameter a low prior variance. Choosing a small prior variance implies that we have a high level of confidence in our prior beliefs, which should be true in this case as we expect that healthy patient fibres resemble heatlhy control subject fibres strongly. Below we calculate the mean slope, intercept and precision for each channel to use as the prior expectations for their respective parameters. We also calculate the variance of the frequentist estimates, to use as the prior expectations for the parameters variance. By splitting the data by channel we can specify a different prior belief for each channel in the dataset, although it is not necessary to do so and the same prior may be set for all channels.

```{r echo=TRUE}
slope_mean = apply(slopes, 1, mean)
inter_mean = apply(intercepts, 1, mean)
prec_mean = apply(precisions, 1, mean)

slope_var = apply(slopes, 1, var)
inter_var = apply(intercepts, 1, var)
prec_var = apply(precisions, 1, var)
```

Before constructing the rest of the prior beliefs, let us inspect the prior beliefs for the slope, 'm'. Is is defined by two parameters; `mu_m` and `tau_m`, both unknown with a normal and gamma prior distribution placed on them respectively. The expected value of `mu_m` can be set to the mean of the frequentist estimates of the slope and its precision chosen by us. Similarly, the expected value of `tau_m` can be chosen by the precision of the frequentist estimates of the slope by setting this to the mode of the prior distribution and its variance chosen by us. For a specified channel, `chan`, the snippet below plots the prior densities for `mu_m`, `tau_m` and `m`. The prior density for `m` is found by repeatedly sampling a value of `mu_m` and `tau_m` from their prior distributions and then sampling a value of `m` from its distribution, using the sampled `mu_m` and `tau_m` as its expectation and precision. By sampling `m` this way we can account for the uncertainty in its governing parameters. The resulting densities are seen in the below figure.

```{r echo=TRUE}
chan = "MTCO1"

# define the mean and precision of mu_m
mean_mu_m = slope_mean[chan] # mean of frequentist slopes
prec_mu_m = 1 / 0.01^2 # chosen by us

# define the shape and rate of tau_m
tau_m_mode = 1/slope_var[chan] # precision of frequentist slopes
tau_m_var = 1 # chosen by us
rate_tau_m = 0.5 * (tau_m_mode + sqrt(tau_m_mode ^ 2 + 4 * tau_m_var)) / tau_m_var
shape_tau_m = 1 + tau_m_mode * rate_tau_m

# plot mu_m prior
curve(dnorm(x, mean_mu_m, 1/sqrt(prec_mu_m)), from=0.5, to=2.5, 
      col=alphaPink(1.0), lwd=2, 
      main="Prior mu_m", xlab="mu_m", ylab="Density")

# plot tau_m prior
curve(dgamma(x, shape_tau_m, rate_tau_m), from=0, to=50, 
      col=alphaPink(1.0), lwd=2, 
      main="Prior tau_m", xlab="tau_m", ylab="Density")

# sample from priors of mu_m and tau_m 
mu_ms = rnorm(1e4, mean_mu_m, 1/sqrt(prec_mu_m))
tau_ms = rgamma(1e4, shape_tau_m, rate_tau_m)

# sample from m prior and plot its density estimate
ms = rnorm(1e4, mu_ms, 1/sqrt(tau_ms))
plot(density( ms ), col=alphaPink(1.0), lwd=2, 
     main="Prior m", xlab="m", ylab="Density")

```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/m_priors_ex.png?raw=true)

The rest of the prior parameters can be inspected in much the same way, although this is not shown here. The remaining parameters must be defined and all parameter values stored in a list to be passed to the inference function. The names of the parameters in the list must be the same as those used in the function otherwise the function will not be able to identify them. If a prior parameter is not chosen and passed to the inference function, a default parameter is used. In practice no parameters need to be passed to the inference function however the default parameters are unlikely to suit all datasets. 

```{r echo=TRUE}
mean_mu_c = inter_mean[chan]
prec_mu_c = 1 / 0.02^2

tau_mode_c = 1/inter_var[chan]
tau_var_c = 0.1
rate_tau_c = 0.5 * (tau_mode_c + sqrt(tau_mode_c ^ 2 + 4 * tau_var_c)) / tau_var_c
shape_tau_c = 1 + tau_mode_c * rate_tau_c

tau_mode = prec_mean[chan]
tau_var = 1
rate_tau = 0.5 * (tau_mode + sqrt(tau_mode ^ 2 + 4 * tau_var)) / tau_var
shape_tau = 1 + tau_mode * rate_tau


paramVals = list(shape_tau=shape_tau, rate_tau=rate_tau, 
                 shape_tau_c=shape_tau_c, rate_tau_c=rate_tau_c, 
                 shape_tau_m=shape_tau_m, rate_tau_m=rate_tau_m,
                 mean_mu_m=mean_mu_m, prec_mu_m=prec_mu_m, 
                 mean_mu_c=mean_mu_c, prec_mu_c=prec_mu_c)
```

### Fit the model

If the data is transformed so that healthy fibres show a linear relationship and is in the correct form then we can now fit the model. Before doing this we must pass the data through the `getData_mats` function which organises the data into matrices to be passed to `rjags`. The following code snippet should run the inference for the first patient in `patIDs` and the channel, `chan`, defined above. We use the parameter values defined in the previous section, without this the function will use a default set of parameters. 

```{r echo=TRUE}

pat = patIDs[1]

dataMats = getData_mats(data=exampleData, 
                        ctrlID=ctrlIDs,
                        channels=c(mitochan, chan),
                        pts=pat, 
                        getIndex=TRUE)

output = inference(dataMats, parameterVals=paramVals)
```

### Understanding inference output

The inference function outputs several things in a list. The first of which are the `POST` and `PRIOR` matrices, containing a list of prior and posterior draws for each parameter in the model in the form of a matrix, where each column is a different parameter. The parameters outputted from the inference can be seen by the column names of the two matrices (which should match exactly). 

```{r echo=TRUE}
colnames(output$POST)
colnames(output$PRIOR)
```

We see here that there are several `m` and `c` parameters. This is because each control sample and patient sample has its own slope and intercept. The slopes for the control samples are given first, as such `m[1]`, `m[2]`, `m[3]` and `m[4]` are the slopes for the control samples, if there are four controls, these will be in the order of `ctrlIDs` passed to the `getData_mats` function. Passing a single patient sample to the model will result in one additional slope, `m[5]`. The same is true for the intercept parameters `c[1]`,..., `c[5]`. Similarly to when we were inspecting the prior distribution for the slope, `m`, we needed to marginalise out the uncertainty from it's governing parameters, `mu_m` and `tau_m`. The inference function has already done this and outputted the resulting draws in the `m_pred` column, in the same way that was done when inspecting the priors before - repeated sampling. Similarly the `c_pred` is the distribution of the intercept. 

Comparing prior and posterior distributions is useful as this shows how our parameter beliefs have been updates. Below is a snippet for how to inspect the prior and posterior densities for the proportion of deficiency. Note that although the densities may extend outside the range [0,1] this is because a density estimate is being plotted, inspection of the output itself shows that no value is outside this range. The resulting plot can be seen in the figure below. 

```{r echo=TRUE}
dPost = density(output$POST[,"probdiff"])
dPrior = density( output$PRIOR[,"probdiff"] )
xlims = range(c(dPrior$x, dPost$x))
ylims = range(c(dPrior$y, dPost$y))

plot(dPrior, col=alphaPink(0.7), lwd=2,
    xlim=xlims, ylim=ylims,
     main="Proportion of Deficiency", xlab="")
lines( dPost, col=alphaGreen(0.7), lwd=2)
legend("topright", legend=c("prior", "post"), lty=1, col=c(alphaGreen(1.0), alphaPink(1.0)))
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/pi_postprior_ex.png?raw=true)

The `POSTPRED` and `PRIORPRED` matrices in the output list contain prior and posterior predictive interval - marginalised over parameter uncertainty - for the patient sample passed to the model. The first column in this is called `mitochan` and are the values on the x-axis for the prediction. The remaining six columns are the 2.5\%, 50\% and 97.5\% quantiles of the prediction for the healthy patient fibres and the deficient patient fibres. This is found by calculating the predictive distribution of the protein expression at each value stored in the `mitochan` column and then calculating its quantiles. 

```{r echo=TRUE}
plot(dataMats$pts, pch=20, col=alphaBlack(0.2),
     xlab=paste0("log(", mitochan ,")"), ylab="log(MTCO1)", 
     main="Posterior Predictive")
lines(output$POSTPRED[,"mitochan"], output$POSTPRED[,"medNorm"], 
      lty=1, col=alphaPink(0.9))
lines(output$POSTPRED[,"mitochan"], output$POSTPRED[,"lwrNorm"], 
      lty=2, col=alphaPink(0.9))
lines(output$POSTPRED[,"mitochan"], output$POSTPRED[,"uprNorm"], 
      lty=2, col=alphaPink(0.9))
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/postpred_ex.png?raw=true)

The last item in the list is called `CLASSIF` and is matrix of every posterior classification for every patient fibre. Each column of the matrix are the posterior classifications of a fibre. We can find the probability that an individual fibre is deficient by calculating the mean of the column. This can be done using the `apply` function, as seen below. The classifications themselves can be plotted using these probabilities and the `classcols` function which converts the the probabilities to a colour on a scale between blue and red, where a deficiency probability of 0.0 is strongly blue and a probability of 1.0 is strongly red.  

```{r echo=TRUE}
def_prob = apply(output$CLASSIF, 2, mean)
plot(dataMats$pts, pch=20, col=classcols(def_prob),
     xlab=paste0("log(", mitochan ,")"), ylab="log(MTCO1)", 
     main="Posterior Predictive")
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/classif_ex.png?raw=true)

### Plotting model output

There are a few plotting functions in the package which are designed to reduce the effort needed by the user. The first of which is a plot which was not seen in the previous section, it plots the MCMC output of the inference. Checking the MCMC output is important to make sure the inference has worked and that there are no major issues with it. There are a number of metrics which can be used, see [this blog](https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics#:~:text=The%20simplest%20way%20to%20diagnose,results%20on%20all%20the%20chunks.) for more information. 

The diagnostic plotter, `MCMCplot`, in this package plots three things per parameter; trace plot, autocorrelation plot (ACF) and a density estimate of the posterior (and prior if given). The trace plot plots each posterior draw in order and the autocorrelation shows the correlation between draws for a range of lags. The density estimate is like that which has been already seen, a kernel density estimate of distribution given a number of draws from it. Using these plots it is important to check some key elements of the output.

1. Has the chain reached a stationary distribution?
2. Are the draws (almost) independent?

To check that the chain has reached a stationary distribution, check the trace plots. An ideal chain would look like a thick straight black line with whiskers coming off it. If the chain is curved at beginning and then flattens out, increase the `MCMCburnin` parameter in the `inference` function. A chain that has not reached a stationary distribution will look like a single jagged line (similar to a stock), again increasing `MCMCburnin` may solve the problem. If the chain is jumping between two or more stationary distributions this will be clear as the density plot will be bimodal. A harder issue to solve, however if there are strong reasons to believe that the parameter should be the lower of higher peak then changing the prior to reflect should help. 

Autocorrelated samples can be quickly diagnosed in the ACF plot. A sample of uncorrelated draws will have no slope in the autocorrelations, showing a spike at a lag of zero and then small random correlations there after. Generally a small amount of correlation in successive samples is fine, when samples are extremely autocorrelated the resulting model can suffer. The amount of correlation can be reduced by increasing the `MCMCthin` parameter in the `inference` function. Below is an example of the `MCMCplot` use and results for the same parameter before and after increasing thinning. A larger amount of thinning resulted in the correlations betwen successive draws being massively decreased. Note: increasing thinning will increase inference time. 

```{r echo=TRUE}
MCMCplot(post=output$POST, 
         prior=output$PRIOR,
         nRow=1)
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/mcmc_m1_ex.png?raw=true)

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
         mitoPlot_ylab="log(NDUFB8)"
         )

par(op) # end plotting grid 
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/postPlot_ex.png?raw=true)

To just see the fibre classification, the `classif_plot` function can be used. The function creates a plot similar to that which was seen before but allows you to plot the posterior predictive without tediously writing the `lines` function. The function also automatically plots the control data in a translucent black colour. 

```{r echo=TRUE}
classif_plot(dataMats=dataMats,
             classifs=def_prob,
             postpred=output$POSTPRED)
```

![alt text](https://github.com/jordanbchilds/analysis2Dmito/blob/main/readme_png/classif_plot_ex.png?raw=true)









