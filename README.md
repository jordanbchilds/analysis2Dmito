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
Getting data into the correct form is crucial to be able to use the inference and plotting function which are part of this package. The package comes with an example dataset which is already in the correct form, so to start we will look at this. 
```{r echo=TRUE include=TRUE}
exampleData = get_exampleData()
head(exampleData)
```
This should show the first six rows of the example dataset, here we have called it `exampleData`. The description of the columns are given in the table below
| Variable Name | Description |
| ------------- | ----------- |
| `sampleID` | The identifier of the sample e.g. `C01`, `C02`, `P01`, `P02`, ... |
|`fibreID` | The identifier of that fibre for that sample. The identifier is not unique throughout the whole dataset, it is only unique for the sample it comes from |
| `sbj_type` | Identifies which sample are control subjects and which are patient subjects |
| `channel` | The channel or protein on which the measurement is made |
| `value` | The raw expression level for this particular fibre on this particular channel |
| --- | --- |
