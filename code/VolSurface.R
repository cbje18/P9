rm(list = ls())
graphics.off()

library(quantmod)
library(dplyr)
library(plotly)


#First get options data from Yahoo
ticker <- "MSFT"
expiration_dates <- c("2023-09-29", "2023-10-06", "2023-10-13", "2023-10-20", "2023-10-27",
                      "2023-11-03", "2023-11-17", "2023-12-15")

options_data <- getOptionChain(ticker, Exp = expiration_dates)

options_data_matrix <- NULL
for(i in 1:length(expiration_dates)){
  options_data_matrix <- rbind(options_data_matrix, options_data[[i]]$calls)
}

#Set current stock price
S_t <- 312.2

options_data_matrix$log_moneyness <- log(S_t/options_data_matrix$Strike)

options_data_matrix$TimeToExp <- as.numeric(as.Date(options_data_matrix$Expiration) - as.Date(Sys.Date()))/365.25

#Create plot
fig <- plot_ly(x = options_data_matrix$log_moneyness, 
               y = options_data_matrix$TimeToExp,
               z = options_data_matrix$IV,
               type = 'mesh3d')
fig

