# analysis2Dmito
An R package for the classification of fibres by their protein expression levels

## Required Software
The package uses the statistical computing softer JAGS which must be installed before use. This cannot be done through R or RStudio/Posit but there are many online resources to help install JAGS. It can be installed directly though [sourceforge.net](https://sourceforge.net/projects/mcmc-jags/files/). A guide to installing R, RStudio, JAGS and R packages can be found [here](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119287995.app1). 

## Installation
To install this package and its dependencies run the following R code. 
```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jordanbchilds/analysis2Dmito")

library("analysis2Dmito")

# install dependencies
install.packages(c("data.table", "dplyr", "readr", "tidyr"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
```

## An Example ClassificationPipe

### Getting Example data
Getting data into the correct form is crucial to be able to use the inference and plotting functions which are part of this package. The package comes with an example dataset which is already in the correct form, so to start we will look at this. The function `get_exampleData` loads the example dataset. 
```{r echo=TRUE include=TRUE}
exampleData = get_exampleData()
head(exampleData)
```
This should show the first six rows of the example dataset, here we have called it `exampleData`. The description of the columns are given in the table below.
| Variable Name | Description |
| ------------- | ----------- |
| `sampleID` | The (unique) identifier of the sample e.g. `C01`, `C02`, `P01`, `P02`, ... |
|`fibreID` | The identifier of that fibre for that sample. The identifier is not unique throughout the whole dataset, it is only unique for the sample it comes from. |
| `sbj_type` | Identifies which sample are control subjects and which are patient subjects. Where control samples are labelled as "control" and patient samples are labelled as "patient". |
| `channel` | The channel or protein on which the measurement is made, a string.  |
| `value` | The expression level for this particular fibre on this particular channel. |

__Any dataset used with the functions in this package must have these five columns__, other columns are allowed but are not necessary. The data is commonly referred to as a long form, a helpful function to be able to get data in this form is the `tidyr::pivot_longer` function (from the `tidyr` package). For information on the `pivot_longer` function and examples see this [blog post](https://tidyr.tidyverse.org/reference/pivot_longer.html).

### Explore the data
Before moving on to inference, although not necessary, it is advisable to explore the data. For good results the healthy control data should show a strong linear relationship, if this is not the case the data should be transformed e.g. by log or square root tranformations.

```{r echo=TRUE include=TRUE}

# the 2Dmito plot x-axis - known for your dataset
mitochan = "raw_porin"
# return all channels which are not mitochan
channels = unique( grep(mitochan, exampleData$channel, value=TRUE, invert=TRUE) )

# extract control sample IDs
ctrlIDs = unique( exampleData[exampleData$sbj_type=="control", "sampleID"] )
# extract patient sample IDs
patIDs = unique( exampleData[exampleData$sbj_type=="patient", "sampleID"] )
patIDs = sort(patIDs)

# plot control data
for(crl in ctrlIDs){
    xDat_crl = exampleData[exampleData$sampleID==crl & exampleData$channel==mitochan, "value"]
  for( chan in channels ){
    yDat_crl = exampleData[exampleData$sampleID==crl & exampleData$channel==chan, "value"]
    plot(xDat_crl, yDat_crl, pch=20, xlab=mitochan, ylab=chan, main=crl )
  }
}

# plot patient data and control data
for( chan in channels ){
  xDat_ctrl = exampleData[exampleData$sbj_type=="control" & exampleData$channel==mitochan, "value"]
  yDat_ctrl = exampleData[exampleData$sbj_type=="control" & exampleData$channel==chan, "value"]

  for( pat in patIDs ){
    xDat_pat = exampleData[exampleData$sampleID==pat & exampleData$channel==mitochan, "value"]
    yDat_pat = exampleData[exampleData$sampleID==pat & exampleData$channel==chan, "value"]
    plot( xDat_ctrl, pch=20, col="black",
          yDat_ctrl, xlab=mitochan, ylab=chan, main=pat,
          xlim=range(exampleData$value), ylim=range(exampleData$value) )
    points( xDat_pat, yDat_pat, pch=20, col="green")
  }
}
```
Before choosing values to summarise prior beliefs about parameter values, we transform the data. Here, we choose a log (to base e) transformation and replot the data to examine the relationship between the transformed variables. 
```{r echo=TRUE}

exampleData$value = log(exampleData$value)
exampleData$channel = gsub("raw_", "log_", exampleData$channel)

mitochan = gsub("raw_", "log_", mitochan)
channels = gsub("raw_", "log_", channels)

for( chan in channels ){
  xDat_ctrl = exampleData[exampleData$sbj_type=="control" & exampleData$channel==mitochan, "value"]
  yDat_ctrl = exampleData[exampleData$sbj_type=="control" & exampleData$channel==chan, "value"]
  
  for( pat in patIDs ){
    xDat_pat = exampleData[exampleData$sampleID==pat & exampleData$channel==mitochan, "value"]
    yDat_pat = exampleData[exampleData$sampleID==pat & exampleData$channel==chan, "value"]
    plot( xDat_ctrl, pch=20, col="black",
          yDat_ctrl, xlab=mitochan, ylab=chan, main=pat,
          xlim=range(exampleData$value), ylim=range(exampleData$value) )
    points( xDat_pat, yDat_pat, pch=20, col="green")
  }
}

```
Prior information can choose by examination of the control data. By fitting linear models to the control samples individually we know what are likely values of the parameters of Bayesian model. For ease this can be done in a frequentist setting. The code snippet below fits a linear model each control sample for each protein and saves relavent output. 
```{r echo=TRUE}
slopes = matrix(NA, nrow=length(channels), ncol=length(ctrlIDs))
rownames(slopes) = channels
colnames(slopes) = ctrlIDs

intercepts = slopes
errors = slopes

for(chan in channels){
  for(crl in ctrlIDs){
    x = exampleData[exampleData$sampleID==crl & exampleData$channel==mitochan, "value"]
    y = exampleData[exampleData$sampleID==crl & exampleData$channel==chan, "value"]
    
    df = data.frame(mitochan=x, chan=y)
    
    lnmod = lm(chan~mitochan, data=df)
    
    slopes[chan, crl] = lnmod$coefficients[1]
    intercepts[chan, crl] = lnmod$coefficients[2]
    errors[chan, crl] = summary(lnmod)$sigma
    
    xSyn = seq(min(exampleData$value)-2, max(exampleData$value)+2, length.out=1000)
    df_pred = data.frame(mitochan=xSyn)
    
    pred = predict.lm(lnmod, newdata=df_pred, interval="prediction")
    
    plot(df, 
         pch=20,
         col=alphaBlack(0.1),
         xlab=mitochan, 
         ylab=chan,
         xlim=range(exampleData$value),
         ylim=range(exampleData$value))
    lines(xSyn, pred[,"lwr"], lty=2, col=alphaPink(0.9), lwd=2)
    lines(xSyn, pred[,"upr"], lty=2, col=alphaPink(0.9), lwd=2)
    lines(xSyn, pred[,"fit"], lty=1, col=alphaPink(0.9), lwd=2)
  }
}
```
The matrices; `slopes`, `intercepts` and `errors`, contain the frequentist estimates for each channel and control sample of their name sakes. The for a specific channel and we can specify priors using the mean and variance of these values. The snippet below calculates the means and standard deviations of the slopes, intercepts and errors, as fitted to the control samples, for each channel. 
```{r echo=TRUE}
slope_mean = apply(slopes, 1, mean)
inter_mean = apply(intercepts, 1, mean)
prec_mean = apply(precisions, 1, mean)

slope_sd = apply(slopes, 1, sd)
inter_sd = apply(intercepts, 1, sd)
prec_sd = apply(precisions, 1, sd)
```
We may choose to set the exprected values of our parameters _a priori_ to the means calculated above and their variances to be small, as we are confident in out beliefs. For the channel, `chan`, here specified to be the first channels in the `channels` vector we can calculate our prior parameters as follows. Savingn them in list called `paramVals` allows them to be passed to the inference function. The names of the parameter values in the list have to be the same as their names used in the model discription otherwise the function would not know what parameter you are tryinng to define. 

```{r echo=TRUE}

chan = channels[1]

mean_mu_m = slope_mean[chan]
prec_mu_m = 1 / 0.01^2
mean_mu_c = inter_mean[chan]
prec_mu_c = 1 / 0.02^2
tau_mode_c = 1/inter_sd[chan]^2
tau_var_c = 0.1
rate_tau_c = 0.5 * (tau_mode_c + sqrt(tau_mode_c ^ 2 + 4 * tau_var_c)) / tau_var_c
shape_tau_c = 1 + tau_mode_c * rate_tau_c
tau_mode_m = 1/slope_sd[chan]^2
tau_var_m = 0.1
rate_tau_m = 0.5 * (tau_mode_m + sqrt(tau_mode_m ^ 2 + 4 * tau_var_m)) / tau_var_m
shape_tau_m = 1 + tau_mode_m * rate_tau_m
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
If the data is transformed so that healthy fibres show a linear relationship and in the correct form then we can now fit the model. Before doing this we must pass the data through the `getData_mats` function which organises the data into matrices to be passed to `rjags`. The following code snippet should run the inference for the first patient in `patIDs` defined above. We use the parameter values defined in the previous section, without this the function will use a default set of parameters. 

```{r echo=TRUE}
exampleData$value = log(exampleData$value)
chan = channels[1]
pat = patIDs[1]

dataMats = getData_mats(data=exampleData, 
                        ctrlID=ctrlIDs,
                        channels=c(mitochan, chan),
                        pts=pat, 
                        getIndex=TRUE)

output = inference(dataMats, parameterVals=paramVals)
```

### Understanding inference output
The inference function outputs several things in a list. The first of which are the `post` and `prior` matrices. This contains a list of posterior draws for each parameter in the model in the form of a matrix where each column is a different parameter. 

The `postpred` and `priorpred` objects are two matrices of the prior and posterior predictive interval - marginalised over parameter uncertainty. The first column in this is called `mitochan` and is the values of the x-axis for the prediction. The remaining six columns are the 2.5\%, 50\% and 97.5\% quantiles of the prediction for the healthy, lie control patient fibres and the deficient fibres. This is found by calculating the predicitive distribution of the y variable (protein expression) at each value stored in the `mitochan` column. 

The last item in the list is called `classif` and is matrix of every posterior classification for every patient fibre. Each column of the matrix are the classifications for a different fibre. We can therefor get the average classification and the probability that an individual fibre is deficient by calculating the mean of the column. This can be done using the `apply` function.

```{r}
def_prob = apply(output$CLASSIF, 2, mean)
```

### Plotting model output
There a few functions which help to visualise the output of the inference. The first of which, `MCMCplot()` is used to check some diagnostics of the MCMC output. The function produces trace plots, autocorrelation (ACF) plots and kernel density estimates an MCMC output that is passed to it. A key thing to look for are that successive posteriors are not highly correlated. This is seen in the ACF plots. 
```{r}
MCMCplot(post=output$post, 
         prior=output$prior,
         nRow=3)
```

To be able to visualise the prior and posterior densities for all variables we can use the `postPlot()` function. The function plots all the priors, posteriors, and predictives (from hierarchical distributions) as well as the model classifications and posterior predictive for the healthy patient data. Unlike `MCMCplot()` this does not automatically create a plotting environment, so we have to define this before calling the function.
```{r}
# plotting grid: 3x3
op = par(mfrow=c(3,3))
postPlot(post=output$POST,
         prior=output$PRIOR,
         postpred=output$POSTPRED,
         dataMats=dataMats,
         classifs=def_prob,
         var.names=c("mu_m0", "tau_m0", "m", "mu_c0", "tau_c0", "c", "prebdef", "tau_norm"))

# end plotting grid 
par(op)
```

If you just want to see the classification, then we can use the `classif_plot()` function. 

```{r}
classif_plot(dataMats=dataMats,
             classifs=def_prob,
             postpred=output$POSTPRED)
```












