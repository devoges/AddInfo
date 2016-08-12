#' Ext.est
#'
#' This function allows you to perform the external data adpative estimation described in "Using auxilary information in statistical function estimate", ESAIM (2006)
#' @param data Dataset to be used
#' @param func.theta The function you wished to estimate, must return an S by 1 vector
#' @param means A list of mean vectors
#' @param vars A list of covariance matrices (of the estimators, not the  population)
#' @param funcs A list of functions, to return vectors of size of vectors in means
#' @param B Number of bootstrap samples, defaults to 500
#' @param eig.keep Keep eig.keep*100 percent of eigen values when inverting a matrix, defaults to 1
#' @export
#' @references
#' \itemize{
#' \item Tarima S, Pavlov D. Using auxilary information in probability estimation. \emph{ESAIM: Probability and Statistics} 2006; 10: 11-23.
#' \item Tarima S, Slavova S, Fritsch T, Hall, L. Probability estimation when some observations are grouped. \emph{Statistics in Medicine} 2007; 26: 1745-61.
#' }
#' @examples
#' ### Use Survival information from Klein and Moeschberger (2003)
#' deaths = c(1,3,3,4,10,13,13,16,16,24,26,27,
#'      28,30,30,32,41,51,65,67,70,72,73,77,91,93,96,
#'      100,104,157,167,1,3,4,5,5,8,12,13,18,23,26,
#'      27,30,42,56,62,69,104,104,112,129,181)
#' cens = c(61,74,79,80,81,87,87,88,
#'      89,93,97,101,104,108,109,120,131,150,231,240,400,
#'      8,67,76,104,176,231)
#' status = 1 - c(rep(0, length(deaths)), rep(1, length(cens)))
#' times = c(deaths, cens)
#' data = data.frame(times, status)
#' ###For basic KM estimate
#' km = survfit(Surv(times, status) ~ 1, data = data)
#' ###But according to
#' ###http://www.cancerresearchuk.org/cancer-help/type/mouth-cancer/treatment/
#' ###statistics-and-outlook-for-mouth-cancers
#' ###About 50% will be alive after 5 years
#' ###We can add this to our own estimator using this function
#'
#' ###First we need to create what we want to estimate, let's say 32 month survival
#' ###Note the only argument should be a dataframe
#' func.theta = function(d) {
#'      km = survfit(Surv(times, status) ~ 1, data = d)
#'      kmest = stepfun(km$time, c(1, km$surv))
#'      kmest(32)
#' }
#' func.theta(data) # = 0.634 is our estimate of survival based on just our data
#' ###With approx 95% CI: 0.5363 to 0.750
#' ###Next we need to create a means list that is our external estimate (at 5 years)
#' means = list(); means[[1]] = 0.5
#' ###Next we need a variance list, for our estimate at 5 years
#' ###In this example, we'll treat 0.5 as exact, perfect, information
#' vars = list(); vars[[1]] = 0
#' ###Lastly, we need to write a list of functions, that will estimate the
#' ###outside data estimate on our own data
#' funcs = list()
#' funcs[[1]] = function(d) {
#'      km = survfit(Surv(times, status) ~ 1, data = d)
#'      kmest = stepfun(km$time, c(1, km$surv))
#'      kmest(60)
#' }
#' ###Now we can use the function to add this outside information
#' res = Ext.est(data, func.theta, means, vars, funcs, B = 500); t = res[[1]]; v = res[[2]]
#' lower = t - 1.96*sqrt(v); upper = t + 1.96*sqrt(v)
#' ###Now our estimate, using the external data, is 0.56, with interval
#' ###0.51 to 0.61, a much smaller interval than with just our data
#'
#' ###Using exact information is really realistic, let's suppose the 0.5 estimate
#' ###at 5 years was due to 20 patients. All we need to change is the vars
#' vars[[1]] = (0.5*0.5)/20
#' res = Ext.est(data, func.theta, means, vars, funcs, B = 500); t = res[[1]]; v = res[[2]]
#' lower = t - 1.96*sqrt(v); upper = t + 1.96*sqrt(v)
#' ###Now the estimate, using the uncertain external data, is still aobut 0.62,
#' which shows that the external data is having less of an influence, with interval
#' ###0.52 to 0.72
#'
#'
#' ### Let's suppose now you want to jointly estimate the survival at 32 and 60 months,
#' ###and you had external information on 24 and 60 months with survival 0.25 and 0.5 based on 20 patients
#' func.theta = function(d){
#'      km = survfit(Surv(times, status) ~ 1, data = d)
#'      kmest = stepfun(km$time, c(1, km$surv))
#'      c(kmest(32), kmest(60))
#' }
#' means[[1]] = c(0.25, 0.5)
#' vars[[1]] = matrix(c(0.25*0.75/20, 0.5, 0.5, 0.5*0.5/20), 2)
#' funcs[[1]] = function(d) {
#'      km = survfit(Surv(times, status) ~ 1, data = d)
#'      kmest = stepfun(km$time, c(1, km$surv))
#'      c(kmest(24), kmest(60))
#' }
#' Ext.est(data, func.theta, means, vars, funcs, B = 500)
Ext.est = function(data, func.theta, means, vars, funcs, B = 500, eig.keep = 1){
  ###Check lengths of lists
  if(!(length(means) == length(vars) & length(vars) == length(funcs)))
    stop("Means, Vars, and Funcs are not all the same length!")

  ###Create function to get length of data, regardless of 1d or not
  len = function(d){
    if(is.vector(d)){
      return(length(d))
    } else {
      return(nrow(d))
    }
  }

  ###Create I and J
  I = length(means); J = length(unlist(means))

  ###Find theta.hat and create S
  theta.hat = func.theta(data); S = length(theta.hat)

  ###Create function to find beta.hat & find it on data
  find.beta.hat = function(d){
    beta.hat = c()
    for(i in 1:length(funcs)) beta.hat = c(beta.hat, funcs[[i]](d))
    return(beta.hat)
  }
  beta.hat = find.beta.hat(data)

  ###Create beta.tilde vector
  beta.tilde = unlist(means)

  ###Create cov(beta.tilde, beta.tilde), k22pp
  k22pp = as.matrix(bdiag(vars))

  ###Start bootstrap to estimate k11, k12, k22p, and k22
  b.beta.hat = matrix(NA, B, J); b.theta.hat = matrix(NA, B, S)
  for(i in 1:B){
    d = as.data.frame(data)[sample(1:len(data), len(data), TRUE),]
    b.beta.hat[i, ] = find.beta.hat(d)
    b.theta.hat[i, ] = func.theta(d)
  }
  k11 = cov(b.theta.hat); k12 = cov(b.theta.hat, b.beta.hat)
  k22p = cov(b.beta.hat); k22 = k22p + k22pp

  ###Find inverse of K22
  if(length(k22) > 1) {
  decomp = svd(k22); sum.eig = sum(decomp$d)
  if(sum(decomp$d == 0) > 0) {
    warning("There are zero eigen values!")
  }
  which.keep = cumsum(decomp$d)/sum(decomp$d)
  inveigs = c()
  for(i in 1:length(which.keep)) {
    if(which.keep[i] <= eig.keep || i == 1) {
      inveigs[i] = 1/decomp$d[i]
    } else {
      inveigs[i] = 0
    }
  }
  invk22 = decomp$v%*%diag(inveigs)%*%t(decomp$u)
  } else {
    invk22 = 1/k22
  }

  ###Get adaptive estimator of theta and covariance k0
  theta.hat.star = theta.hat - k12%*%invk22%*%(beta.hat - beta.tilde)
  k0 = k11 - k12%*%invk22%*%t(k12)

  out = list()
  out$theta.hat.star = theta.hat.star
  out$k0hat =  k0
  return(out)
}
