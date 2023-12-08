# --- READ THIS ---
# It is advised that before using `analysis2Dmito` package for analysis you read
# the README file available on the GitHub page.

# This example uses the dataset from the paper; A stagewise response to
# mitochondrial dysfunction in mitochondrial DNA maintenance disorders by
# Vincent et al. (https://www.biorxiv.org/content/10.1101/2023.10.24.563545v1).
# A subset of the dataset itself is available to download through this package
# with the function `analysis2Dmito::get_exampleData`.
#
# The example data has already been formatted so that it can be used by functions
# within this package, see GitHub (https://github.com/jordanbchilds/analysis2Dmito)
# for more details.
#
# The measure of mitochondrial mass is manually set, as this varies
# between datasets and cannot be generalised, the remaining channels present in
# the dataset are assumed to be OXPHOS protein measurements on which we wish to
# fit the model. The mitochondrial mass protein channel is the first thing to be
# defined in the script.
#
# If something in this script produces an error check lines which are commented
# out with a single '#', it could that a required package is not installed.
# Lines commented out with a '# ---' are comments to guide you.

# --- measure of mitochondrial mass (x-axis in 2Dmito)
mitochan = "VDAC"

# --- install and load dependencies
# install.packages(c("data.table", "dplyr", "readr", "tidyr", "rjags"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("rstan")
library("parallel")

# install.packages("devtools")
library("devtools")
# --- install analysis2Dmito from GitHub (https://github.com/jordanbchilds/analysis2Dmito)
devtools::install_github("jordanbchilds/analysis2Dmito")
library("analysis2Dmito")

# --- load example data
data_raw = analysis2Dmito::get_exampleData()

# --- transform expression values
data = data_raw
data$value = log( data_raw$value ) # log-transform data

# --- find control and patient subject IDs
sbj = unique(data$sampleID)
ctrlIDs = grep("C", sbj, value=TRUE)
pts = sort( sbj[!(sbj %in% ctrlIDs)] )

# --- all channels present in the data, including mitochondrial mass measure
channels_all = unique(data$channel)
# --- channels of interest (y-axis in 2Dmito) i.e. any channel which is not the
# --- measure of mitochondrial mass
channels = channels_all[channels_all != mitochan]

# --- create a list of datasets to be passed to the inference function
data_list = list()
for (chan in channels) {
  for (pat in pts) {
    root = paste(chan, pat, sep = "_")
    data_list[[root]] = analysis2Dmito::getData_mats(
      data,
      pts = pat,
      channels = c(mitochan, chan),
      ctrlID = ctrlIDs)
  }
}

# ---------------------------- #
# --- RUN INFERENCE & SAVE --- #
# ---------------------------- #

# --- The inference is split across half the number of available coreson your
# --- machine. If your machine is starts to crash this can be changed to a lower
# --- number.
# ncores = 1
ncores = floor( detectCores() / 2 )
# --- The following definitions are the default of the stan_inference function,
# --- but they are needed to be able to check chain convergence post inference
cl  = makeCluster(ncores)
{
  output = parLapply(cl,
                     data_list,
                     analysis2Dmito::stan_inference)
}
stopCluster(cl)
# --- The output of the inference in saved in a list called `output`. Each
# --- the output from a single 2Dmito plot is an element of this list.
# --- The names of the list elements are such that if it is the output for
# --- channel, 'Ch1', and patient, 'P01', then its names is 'Ch1_P01'.

# ------------------- #
# --- SAVE OUTPUT --- #
# ------------------- #

# --- Create a directory (folder) in the current working directory called
# --- "Output" in which the outputs are saved
dir.create("Output")

# --- Save the output using the `analysis2Dmito::list_saver` function.
# --- The output is saved such that for channel, Ch1, and patient, P01, the
# --- posterior files are saved with the file path `./Output/Ch1_P01__POST.txt`
for (root in names(output)) {
  analysis2Dmito::list_saver(output[[root]], file.path("Output", root), rootSep = "_")
}

# ------------------- #
# --- PLOT OUTPUT --- #
# ------------------- #

dir.create(file.path("PDF"), showWarnings = FALSE)

# ------ MCMC plot
pdf(file.path("PDF", "MCMCplot.pdf"), width=13, height=8)
{
  for( root in names(output) ){
    post = output[[root]]$POST
    prior = output[[root]]$PRIOR

    analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100,
                             main_title=root)
  }
}
dev.off()

# ------ postPlot
pdf(file.path("PDF", "postPlot.pdf"), width=13, height=8)
{
  for( chan in channels ){
    for( pat in pts ){
      root = paste0(chan, "_", pat)
      post = output[[root]]$POST
      prior = output[[root]]$PRIOR
      postpred = output[[root]]$POSTPRED
      class = colMeans(output[[root]]$CLASSIF)

      dataMats = analysis2Dmito::getData_mats(data=data,
                                              channels=c(mitochan, chan),
                                              ctrlID=ctrlIDs,
                                              pts=pat)

      op = par(mfrow=c(3,3), mar=c(4.3,4.1,1.5,2), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
      analysis2Dmito::postPlot(post=post, prior=prior,
               postpred=postpred,
               classifs=class,
               dataMats = dataMats,
               var.names=c("mu_m", "tau_m", "tau_norm", "mu_c", "tau_c", "probdiff", "m", "c"),
               mitochan=paste0("log(", mitochan, ")"),
               chan=paste0("log(", chan, ")"),
               pat_id=pat)
      par(op)
    }
  }
}
dev.off()

# ------ Classif plot
pdf(file.path("PDF", "classif.pdf"), width=13, height=8)
{
  op = par( mfrow = c(1, 1), mar = c(6, 6, 6, 3), cex.main = 2, cex.axis = 1.5, cex.lab = 2)
  for (chan in channels) {
    for (pat in pts) {
      root = paste0(chan, "_", pat)

      class = apply(output[[root]]$CLASSIF, 2, mean)
      post = output[[root]]$POST

      pi_est = round(mean(post[, "probdiff"]), 3)

      postpred = output[[root]]$POSTPRED

      dataMats = analysis2Dmito::getData_mats(
        data = data,
        channels = c(mitochan, chan),
        ctrlID = ctrlIDs,
        pts = pat
      )

      analysis2Dmito::classif_plot(
        dataMats = dataMats,
        postpred = NULL,
        classifs = class,
        xlab = mitochan,
        ylab = chan
      )
      title(main = bquote(atop(.(pat), "nPat:" ~ .(nrow(dataMats$pts)) * ",   E(" * pi * "|X)=" * .(pi_est))))
    }
  }
  par(op)
}
dev.off()

# --- Calculating ESS
install.packges("coda")
install.packages("mcmcse")
library("coda")
library("mcmcse")

ess_list = list()
multiESS = list()
for( chan in channels ){
  for( pat in pts ){
    root = paste0(chan, "_", pat)
    post = output[[root]]$POST
    ess_list[[root]] = coda::effectiveSize( post )
    multiESS[[root]] = mcmcse::multiESS( post )

    print(root)
    print(multiESS[[root]])
    print(ess_list[[root]])
  }
}









