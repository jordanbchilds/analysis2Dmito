# analysis2Dmito
An R package for the classification of fibres by their protein expression levels

## Installation
```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jordanbchilds/analysis2Dmito")

library("analysis2Dmtio")
```

## An Example ClassificationPipe

### Getting data
Getting data into the correct form is crucial to be able to use the inference and plotting function which are part of this package. The package comes with andexample dataset which is already in the correct form, so to start we will look at this. 
```{r echo=TRUE include=TRUE}
exampleData = get_exampleData()
head(exampleData)
```

