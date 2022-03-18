# -------------------------------
# Collection of various functions
# -------------------------------

# Function to clear the console
clc <-function() {
  cat("\f");
}

# Function that read the dataset
read_data <-function(directory) {
  setwd(directory);
  colClass <-c("character", "integer", "character", "integer", "character", "character", "character", "integer", 
                 rep("numeric", 12), "integer", rep("numeric", 20));
  data <-read.xlsx("data.xlsx", sheetName = "DATI", rowIndex = 3:1212, colIndex = 2:42, as.data.frame = TRUE,
                            header = FALSE, colClasses = colClass);
  colnames(data) <-c("lev", "sem", "camp", "cfu", "ssd", "cds", "sex", "age", "pm_i", "pf_i", "pi_i", "pf_i",
                     "pm_cs", "pf_cs", "pi_cs", "pf_cs", "pm_ac", "pf_ac", "pi_ac", "pf_ac", "ns",
                     "mq1", "mq2", "mq3", "mq4", "mq5", "mq6", "mq7", "mq8", "mq9", "mq10",
                     "mq11", "mq12", "mq13", "mq14", "mq15", "mq16", "mq17", "mq18", "mq19", "mq20");
  data$cds <-gsub("\n", "", data$cds, fixed = TRUE);
  save(data, file = "data.rda");
}

# Function that fix the NA to the respective columns' averages
na_fixto_av <-function(data_na) {
  N <-dim(data_na)[1];
  p <-dim(data_na)[2];
  data <-data_na;
  av_vec <-colmean(data_na[complete.cases(data_na), ]);
  row_na <-1:N;
  row_na <-row_na[!complete.cases(data_na)];
  for (i in row_na) {
    col_na <-1:p;
    col_na <-col_na[is.na(data_na[i, ])];
    for (j in col_na) {
      data[i, j] <-av_vec[j];
    }
  }
  return(data);
}

# Function to compute the grand mean over the schede
gr_mean <-function(data_mat, n) {
  sum_mat <-t(data_mat)%*%n;
  return(sum_mat/sum(n));
}

# Function to compute means over columns
colmean <-function(data_mat) {
  return(sapply(data_mat, mean, na.rm = TRUE));
}

# Function to compute sds over columns
colsd <-function(data_mat) {
  return(sapply(data_mat, sd, na.rm = TRUE));
}

# Function to build the graph prior adjacency matrix p20 is a 20x1 column vector, p_arg is a 4x1 column vector
graphPrior_adj <-function(p20, p_arg) {
  graph_p <-matrix(rep(0.5, times = 20*20), nrow = 20, ncol = 20);
  graph_p[, 20] <-t(p20);
  graph_p[20, ] <-p20;
  for (i in 1:20) {
    for (j in 1:20) {
      if ((i >= 1) && (i <= 6) && (j >= 1) && (j <= 6)) { # Insegnamento
        graph_p[i, j] <-p_arg[1];
      }
      if ((i >= 7) && (i <= 13) && (j >= 7) && (j <= 13)) { # Docenza
        graph_p[i, j] <-p_arg[2];
      }
      if ((i >= 14) && (i <= 16) && (j >= 14) && (j <= 16)) { # Attivita didattiche integrative
        graph_p[i, j] <-p_arg[3];
      }
      if ((i >= 17) && (i <= 19) && (j >= 17) && (j <= 19)) { # Infrastrutture per questo insegnamento
        graph_p[i, j] <-p_arg[4];
      }
    }
  }
  diag(graph_p) <-rep(0, times = 20);
  return(graph_p);
}

# Function to extract the number of edges at each iteration from bdgraph outputs
nof_edges_extr <-function(graph_list, graph_id) {
  iter <-length(graph_id);
  nof_edges <-rep(0, times = iter);
  for (i in 1:iter) {
    nof_edges[i] <-length(gregexpr(pattern = "1", graph_list[graph_id[i]])[[1]]);
  }
  return(nof_edges);
}

# Function to convert stan output to coda
stan2coda <- function(fit) 
{
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

# Function to compute LPLM from a stanfit object
LPML <-function(fit, par_name = "lik") {
  lik <-as.matrix(extract(fit, pars = par_name)$lik);
  N <-dim(lik)[1]; Ns <-dim(lik)[2];
  CPO <-rep(0, times = N);
  for (i in 1:Ns) {
    CPO[i] <-1/(mean(1./lik[, i]));
  }
  LPLM <-(1/N)*sum(log(CPO));
  return(LPLM);
}

# Function to compute WAIC from a stanfit object
WAIC <-function(fit, par_name = "lik") {
  lik <-as.matrix(extract(fit, pars = par_name)$lik);
  N <-dim(lik)[2];
  lppd <-sum(log(apply(lik, MARGIN = 2, FUN = mean)));
  ro <-sum(apply(log(lik), MARGIN = 2, FUN = var));
  WAIC <--(2/N)*(lppd - ro);
  return(WAIC)
}

# Function to compute 2logBF of hypothesis H0: theta <= q vs H1: theta > q with same prior i.e BF = posterior odds
BF01sp <-function(theta, q) {
  post0 <-length(which(theta <= q))/length(theta);
  post1 <-1-post0;
  return(round(2*log(post0 / post1), 2));
}

# Function to compute the mean of a skew-multivariate-normal distribution
skMean <-function(mu, w, a, R) {
  d <-(1/sqrt(1+quad.form(R, a)))*R%*%a;
  return(mu + sqrt(2/pi)*(diag(w)*d));
}

# Monte Carlo estimate volume out of 1-4 hypercube support from multivariate-skew-normal density samples
mtc_cdf <-function(x) {
  N <-dim(x)[1];
  p <-dim(x)[2];
  c1 <-c4 <-0;
  for (i in 1:N) {
    for (j in 1:p) {
      if (x[i, j] < 1) {
        c1 <-c1+1;
        break;
      }
      if (x[i, j] > 4) {
        c4 <-c4+1;
        break;
      }
    }
  }
  return(c(c1/N, c4/N));
}
