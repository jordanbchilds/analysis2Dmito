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

__Any dataset used with the functions in this package must have these five columns__, other columns are allowed but are not necessary. The data is commonly referred to as a long form, a helpful function to be able to get data in this form is the `tidyr::pivot_longer` function (from the `tidyr` package). For information on this see [this blog](https://tidyr.tidyverse.org/reference/pivot_longer.html).

### Explore the data
Before moving on to inference, although not necessary, it is advisable to explore the data. For good results the healthy control data should show a strong linear relationship, if this is not the case the data should be transformed e.g. log or sqrt.

```{r echo=TRUE include=TRUE}

# the 2Dmito plot x-axis - known for your dataset
mitochan = "raw_porin"
# return all channels which are not mitochan
channels = unique( grep(mitochan, exampleData$channel, value=TRUE, invert=TRUE) )

# extract control sample IDs
ctrlID = unique( exampleData[exampleData$sbj_type=="control", "sampleID"] )
# extract patient sample IDs
patIDs = unique( exampleData[exampelData$sbj_type=="patient", "sampleID"] )

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
          xlims=range(xDat_ctrl, xDat_pat), ylims=range(yDat_ctrl, yDat_pat) )
    points( xDat_pat, yDat_pat, pch=20, col="green")
  }
}
```

### Fit the model
If the data is transformed so that healthy fibres show a linear relationship and in the correct form then we can now fit the model. Before doing this we must pass the data through the `getData_mats` function which organises the data into matrices to be passed to `rjags`. The following code snippet should run the inference for the first patient in `patIDs` defined above. It will run with the default parameter set, but these can be altered using the `parameterValue` argument to the function. 

```{r echo=TRUE}
exampleData$value = log(exampleData$value)
pat = patIDs[1]

dataMats = getData_mats(data=exampleData, 
                        ctrlID=ctrlID, 
                        pts=pat, 
                        getIndex=TRUE)

output = inference(dataMats)
```

### Understanding inference output
The inference function outputs several things in a list. The first of which are the `post` and `prior` matrices. This contains a list of posterior draws for each parameter in the model in the form of a matrix where each column is a different parameter. 

The `postpred` and `priorpred` objects are two matrices of the prior and posterior predictive interval - marginalised over parameter uncertainty. The first column in this is called `mitochan` and is the values of the x-axis for the prediction. The remaining six columns are the 2.5\%, 50\% and 97.5\% quantiles of the prediction for the healthy, lie control patient fibres and the deficient fibres. This is found by calculating the predicitive distribution of the y variable (protein expression) at each value stored in the `mitochan` column. 

The last item in the list is called `classif` and is matrix of every posterior classification for every patient fibre. Each column of the matrix are the classifications for a different fibre. We can therefor get the average classification and the probability that an individual fibre is deficient by calculating the mean of the column. This can be done using the `apply` function.

```{r}
def_prob = apply(output$classif, 2, mean)
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

postPlot(post=output$post,
         prior=output$prior,
         postpred=output$postpred,
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
             postpred=output$postpred)
```












