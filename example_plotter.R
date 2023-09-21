library("devtools")
#install_github("jordanbchilds/analysis2Dmito", force=TRUE)
library(analysis2Dmito)

# load data in
data_raw = get_exampleData()

# transform data
data = data_raw
data$value = log( data_raw$value )

# find control and patient subject IDs
sbj = unique(data$sampleID)
ctrlIDs = grep("C", sbj, value=TRUE)
pts = sort( sbj[!(sbj %in% ctrlIDs)] )


# measure of mitochondrial mass in the dataset (x-axis in 2Dmito)
mitochan = "VDAC"
# proteins of interest (y-axis in 2Dmito)
channels_all = unique(data$channel)
channels = channels_all[channels_all != mitochan]

# create folder to save output - if folder already exists do nothing
dir.create(file.path("PDF"), showWarnings = FALSE)

# ------ MCMC plot
pdf(file.path("PDF", "MCMCplot.pdf"), width=13, height=8)
{
  for(chan in channels){
    for(pat in pts) {
      outroot = paste(chan, pat, sep="_")
      fn_post = file.path("Output", paste0(outroot, "__POST.txt"))
      if(file.exists(fn_post)){
        post = as.data.frame(fread(fn_post))
        fn_prior = file.path("Output", paste0(outroot, "__PRIOR.txt"))
        prior = as.data.frame(fread(fn_prior))
        analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100)
        title(main=paste(chan, pat), outer=TRUE, line=-1)
      }
    }
  }
}
dev.off()

# ------ postPlot
pdf(file.path("PDF", "postPlot.pdf"), width=13, height=8)
{
  for(chan in channels){
    for( pat in pts ){

      outroot = paste(chan, pat, sep="_")
      fn_post = file.path("Output", paste0(outroot, "__POST.txt"))
      if(file.exists(fn_post)){
        post = as.data.frame(fread(fn_post))
        fn_prior = file.path("Output", paste0(outroot, "__PRIOR.txt"))
        prior = as.data.frame(fread(fn_prior))
        fn_postpred = file.path("Output", paste0(outroot, "__POSTPRED.txt"))
        postpred = as.data.frame(fread(fn_postpred))

        fn_class = file.path("Output", paste0(outroot, "__CLASSIF.txt"))
        class = apply(as.matrix(fread(fn_class)), 2, mean)


        dataMats = analysis2Dmito::getData_mats(data=data,
                                                channels=c(mitochan, chan),
                                                ctrlID=ctrlIDs,
                                                pts=pat)

        op = par(mfrow=c(3,3), mar=c(4,4,3,3), cex.main=2, cex.lab=1.5, cex.axis=1.5)
        analysis2Dmito::postPlot(post=post, prior=prior,
                                 postpred=postpred,
                                 classifs=class,
                                 dataMats = dataMats,
                                 var.names=c("mu_m", "tau_m", "tau_norm", "mu_c", "tau_c", "probdiff", "m", "c"),
                                 mitoPlot_xlab=paste0("log(", mitochan, ")"),
                                 mitoPlot_ylab=paste0("log(", chan, ")"))
        title(main=paste(chan, pat), outer=TRUE, line=-2)
        par(op)
      }
    }
  }
}
dev.off()

# ------ Classif plot
pdf(file.path("PDF", "classif.pdf"), width=8, height=5)
{
  op = par(mfrow=c(1,1), mar=c(6,6,3,3), cex.main=2, cex.axis=1.5, cex.lab=2)
  for(chan in channels){
    for( pat in pts ){
      outroot = paste(chan, pat, sep="_")
      fn_class = file.path("Output", paste0(outroot, "__CLASSIF.txt"))
      if(file.exists(fn_class)){
        class = apply(as.matrix(fread(fn_class)), 2, mean)

        fn_postpred = file.path("Output", paste0(outroot, "__POSTPRED.txt"))
        postpred = as.data.frame(fread(fn_postpred))

        dataMats = analysis2Dmito::getData_mats(data=data,
                                                channels=c(mitochan, chan),
                                                ctrlID=ctrlIDs,
                                                pts=pat)

        analysis2Dmito::classif_plot(dataMats=dataMats,
                                     postpred=NULL,
                                     classifs=class,
                                     xlab=paste0("log(", mitochan, ")"),
                                     ylab=paste0("log(", chan, ")"))
      }
    }
  }
  par(op)
}
dev.off()
