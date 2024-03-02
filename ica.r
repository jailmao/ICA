#install.packages("tidyverse") 
#uncomment if you dont have it installed
library(tidyverse)

n225 <- read.csv("N225.csv", stringsAsFactors=FALSE)
fchi <- read.csv("FCHI.csv", stringsAsFactors=FALSE)

n225_r <- function(t) {
    return (log(as.numeric(n225[t+1,2]))-log(as.numeric(n225[t,2])))
}
fchi_r <- function(t) {
    return (log(as.numeric(fchi[t+1,2]))-log(as.numeric(fchi[t,2])))
}
numOfWeeks = nrow(n225)
returnsDataSet <- data.frame(matrix(ncol=3, nrow=0, dimnames=list(NULL, c("Week", "n225_returns", "fchi_returns"))))

for (i in 1:numOfWeeks-1){
    returnsDataSet[i,1]=i
    returnsDataSet[i,2]=n225_r(i)
    returnsDataSet[i,3]=fchi_r(i)
}


ggplot(data = returnsDataSet, aes(x = Week)) +
    geom_point(aes(y = n225_returns, color = "N225 Returns"), size = 1) +
    geom_point(aes(y = fchi_returns, color = "FCHI Returns"), size = 1) +
    labs(x = "Week", y = "log returns", color = "Returns") +
    scale_color_manual(values = c("blue", "red")) +
    theme(legend.position = "top")
print(head(returnsDataSet))
