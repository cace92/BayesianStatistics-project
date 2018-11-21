#------------------------
# Intro analysis of data
#------------------------

# Initialization--------------------------------------------------
library(xlsx);
library(plotrix);
library(gridExtra);

source("functions.r");
# If RData object has not been created yet, please run next function
#read_data("C:/Users/Cace/Desktop/Project/Scripts");
load("data.rda");
N <-dim(data)[1];
p <-dim(data)[2];


# Dataset analysis---------------------------------------------
attach(data);
t_small <-ttheme_default(base_size = 5);
# Frequency table: level
lev_t <-tableGrob(data.frame(table(lev)));
# Frequency table: semester
sem_t <-tableGrob(data.frame(table(sem)));
# Frequency table: campus
camp_t <-tableGrob(data.frame(table(camp)));
# Frequency table: cfu
cfu_t <-tableGrob(data.frame(table(cfu)));
# Frequency table: sex
sex_t <-tableGrob(data.frame(table(sex)));
x11();
hcombine <-gtable_combine(lev_t, sem_t, camp_t, cfu_t, sex_t, along = 1);
grid.arrange(hcombine);
# Histogram: age
x11();
hist(age, main = "Age frequencies");
age_stat <-data.frame(mean = signif(mean(age), digits = 4), sd = signif(sd(age), digits = 3),
                      min = min(age), max = max(age));
grid.table(age_stat);
# Frequency table: ssd
x11();
ssd_t1 <-tableGrob(data.frame(table(ssd)[1:21]), theme = t_small);
ssd_t2 <-tableGrob(data.frame(table(ssd)[22:42]), theme = t_small);
ssd_t3 <-tableGrob(data.frame(table(ssd)[43:63]), theme = t_small);
ssd_t4 <-tableGrob(data.frame(table(ssd)[64:84]), theme = t_small);
ssd_t5 <-tableGrob(data.frame(table(ssd)[85:98]), theme = t_small);
hcombine <-gtable_combine(ssd_t1, ssd_t2, ssd_t3, ssd_t4, ssd_t5, along = 1);
grid.arrange(hcombine);
# Frequency table: cds
x11();
cds_t1 <-tableGrob(data.frame(table(cds)[1:28]), theme = t_small);
cds_t2 <-tableGrob(data.frame(table(cds)[29:57]), theme = t_small);
hcombine <-gtable_combine(cds_t1, cds_t2, along = 1);
grid.arrange(hcombine);

# Coping with NA--------------------------------------------------
detach(data);
row_na <-!complete.cases(data);
data_na <-data[row_na,];
attach(data_na);
# Frequency table of na per cds
x11();
grid.table(data.frame(table(cds)));
# Plot of na values
mat_na <-data.frame(matrix(as.numeric(is.na(data_na)), nrow = dim(data_na)[1], ncol = dim(data_na)[2]));
colnames(mat_na) <-colnames(data_na);
x11()
color2D.matplot(mat_na, axes = FALSE, xlab = "Variables", ylab = "Row containing NA val");
boxed.labels(0:(dim(mat_na)[2]-1)+0.5, rep(dim(mat_na)[1]+1, dim(mat_na)[2]), labels = colnames(mat_na),
             border = FALSE, srt = 90, cex = 0.9, adj = 0.1);
detach(data_na);
