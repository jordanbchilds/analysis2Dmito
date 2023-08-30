# install dependencies
# install.packages(c("data.table", "dplyr", "readr", "tidyr", "rjags"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("rjags")
library("parallel")

# install.packages("devtools")
library("devtools")
devtools::install_github("jordanbchilds/analysis2Dmito", force=TRUE)
library("analysis2Dmito")

data_raw = get_exampleData()

data = data_raw
data$value = log( data_raw$value)

sbj = unique(data$sampleID)
ctrlIDs = grep("C", sbj, value=TRUE)
pts = sort( sbj[!(sbj %in% ctrlIDs)] )


# measure of mass (x-axis in 2Dmito)
mitochan = "VDAC"
# proteins of interest (y-axis in 2Dmito)
channels_all = unique(data$channel)
channels = channels_all[channels_all != mitochan]

slopes = matrix(NA, nrow=length(channels), ncol=length(ctrlIDs))
rownames(slopes) = channels
colnames(slopes) = ctrlIDs

intercepts = slopes # defines an empty matrix with row and col names as wanted
precisions = slopes

for(chan in channels){
  for(crl in ctrlIDs){
    x = data[data$sampleID==crl & data$channel==mitochan, "value"]
    y = data[data$sampleID==crl & data$channel==chan, "value"]
    df = data.frame(mitochan=x, chan=y)
    
    lnmod = lm(chan~mitochan, data=df)
    
    slopes[chan, crl] = lnmod$coefficients["mitochan"]
    intercepts[chan, crl] = lnmod$coefficients["(Intercept)"]
    precisions[chan, crl] = 1 / summary(lnmod)$sigma^2
  }
}

# slope_mean = apply(slopes, 1, mean)
# inter_mean = apply(intercepts, 1, mean)
# prec_mean = apply(precisions, 1, mean)
# 
# slope_var = apply(slopes, 1, var)
# inter_var = apply(intercepts, 1, var)
# prec_var = apply(precisions, 1, var)


data_list = list()
for( chan in channels ){
  for(pat in pts){
    data_list[[paste(chan, pat, sep="_")]] = getData_mats(data, 
                                                          pts=pat, 
                                                          channels=c(mitochan, chan), 
                                                          ctrlID=ctrlIDs)
  }
}
  
mean_mu_m = mean(slopes)
prec_mu_m = 1 / 0.05 ^ 2

mean_mu_c = mean(intercepts)
prec_mu_c = 1 / 0.05 ^ 2

tau_m_mode = 1 / var(c(slopes))
tau_m_var = 1
rate_tau_m = 0.5 * (tau_m_mode + sqrt(tau_m_mode ^ 2 + 4 * tau_m_var)) / tau_m_var
shape_tau_m = 1 + tau_m_mode * rate_tau_m

tau_c_mode = 1 / var(c(intercepts))
tau_c_var = 1
rate_tau_c = 0.5 * (tau_c_mode + sqrt(tau_c_mode ^ 2 + 4 * tau_c_var)) / tau_c_var
shape_tau_c = 1 + tau_c_mode * rate_tau_c

tau_mode = mean(precisions)
tau_var = 1
rate_tau = 0.5 * (tau_mode + sqrt(tau_mode ^ 2 + 4 * tau_var)) / tau_var
shape_tau = 1 + tau_mode * rate_tau

tau_def = 0.001

paramVals = list(
  mean_mu_m = mean_mu_m,
  prec_mu_m = prec_mu_m,
  mean_mu_c = mean_mu_c,
  prec_mu_c = prec_mu_c,
  shape_tau_m = shape_tau_m,
  rate_tau_m = rate_tau_m,
  shape_tau_c = shape_tau_c,
  rate_tau_c = rate_tau_c,
  tau_def = tau_def
)


ncores = detectCores() - 1
cl  = makeCluster(ncores)
{
  output = parLapply(
    cl,
    data_list,
    inference,
    MCMCburnin = 10,
    MCMCout = 1000,
    MCMCthin = 1000,
    parameterVals = paramVals
  )
}
stopCluster(cl)

for( root in names(output)){
  list_saver(output[[root]], file.path("Output", root))
}


