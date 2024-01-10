rm(list = ls())

library(lubridate)

setwd("C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/GITHUB_EXAMPLE")

# loading functions

functions_path <- list.files(path = "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/GITHUB_EXAMPLE/Functions/",
                             pattern = ".R", full.names = TRUE)

sapply(functions_path,source)

# data generation setting

library(readr)

simulation_plan <- read_csv("Data/simulation_values.csv")

true_lat <- simulation_plan$eq_lat

true_lon <- simulation_plan$eq_lon

true_depth <- simulation_plan$eq_depth

true_time <- dmy_hms(paste0(simulation_plan$eq_origin,sep = " ","00:00:00"))

base_time <- dmy_hms(simulation_plan$eqn_time) - 120 

true_time <- as.numeric(difftime(true_time, base_time, units = "secs"))

true_alpha <- simulation_plan$alpha

true_pi <- simulation_plan$pi

true_vP <- simulation_plan$v_P

true_invlambda <- simulation_plan$inverse_lambda

# data importation

data <- read_csv("Data/sim_data.csv")

# z0 prior 

sigmalon <- 1
sigmalat <- 1

# d0 prior 

s <- 3.939867
r <- 0.2437115

# alpha prior

a <- 0.5
b <- 0.5

# vP prior 

c <- 69.94574
d <- 10.94473

# parallel tempering variables

d_log_temp <- rep(1,9)

L <- length(d_log_temp) + 1

beta <- inv_temperatures(d_log_temp,L)

# MCMC setting

G <- 5000
burnin <- 5000
thin <- 5

eta_1 <- array(dim = c(4,4,L))

for (l in 1:L) {
  eta_1[,,l] <- diag(c(0.1,0.1,10,1))
}

eta_2 <- rep(0.1,L)

eta_3 <- rep(0.1,L)

eta_4 <- rep(0.1,L)

eta_5 <- rep(0.1,L)

eta_6 <- rep(0.1,L)

s_1 <- rep(0.1,L)

s_2 <- rep(0.1,L)

s_3 <- rep(0.1,L)

s_4 <- rep(0.1,L)

s_5 <- rep(0.1,L)

s_6 <- rep(0.1,L)


# EQN estimates

event_lat <- simulation_plan$eqn_lat
event_lon <- simulation_plan$eqn_lon
event_time <- simulation_plan$eqn_time

n <- dim(data)[1]

data$delta <- rep(1,n)
data$delta[which(data$time=="NaT")] <- 0

data$time[which(data$time=="NaT")] <- event_time
data$time <- as.numeric(difftime(dmy_hms(data$time), base_time, units = "secs"))
data <- data[which(data$time>=0),]

# initial values for the MCMC

init_lat <- event_lat
init_lon <- event_lon
init_t0 <- as.numeric(difftime(dmy_hms(event_time), base_time, units = "secs"))

# z0 prior mean 

latstar <- mean(data$lat[which(data$delta==1)])
lonstar <- mean(data$lon[which(data$delta==1)])

# initial values

set.seed(2024)

th0 <- c(runif(1,(init_lat - 1),(init_lat + 1)), 
         runif(1,(init_lon - 1),(init_lon + 1)),
         runif(1,5,100),
         runif(1,0,init_t0),
         runif(1),
         runif(1),
         runif(1,3,7),
         runif(1,0,3.5),
         runif(1,800,80000))

library(tictoc)

tic()

samples <- mh_surv_quake_2(G = G, burnin = burnin, thin = thin, th0 = th0,
                           eta_1 = eta_1, eta_2 = eta_2, eta_3 = eta_3, eta_4 = eta_4, eta_5 = eta_5, eta_6 = eta_6,
                           s_1 = s_1, s_2 = s_2, s_3 = s_3, s_4 = s_4, s_5 = s_5, s_6 = s_6,
                           y = data$time, lats = data$lat, lons = data$lon, delta = data$delta,
                           latstar = latstar, lonstar = lonstar, 
                           sigmalat = sigmalat, sigmalon = sigmalon, 
                           r = r, s = s, a = a, b = b, c = c, d = d,
                           beta = beta, d_log_temp = d_log_temp)
toc()

colnames(samples) <- c("lat0","lon0","d0","t0","alpha","pi","vP","tau","invlambda")

# posterior estimates

pdens <- list()
hdiR <- list()
hdi95 <- list()

library(HDInterval)

for (p in 1:9) {
  pdens[[p]] <- density(samples[,p,1])
  hdiR[[p]] <- hdi(pdens[[p]],allowSplit = TRUE,credMass = 0.0001)
  hdi95[[p]] <- hdi(pdens[[p]],allowSplit=TRUE,credMass = 0.95)
}

post_estimates <- c(median(hdiR[[1]][1,]),
                    median(hdiR[[2]][1,]),
                    median(hdiR[[3]][1,]),
                    median(hdiR[[4]][1,]),
                    median(hdiR[[5]][1,]),
                    median(hdiR[[6]][1,]),
                    median(hdiR[[7]][1,]),
                    median(hdiR[[8]][1,]),
                    median(hdiR[[9]][1,]))

post_err <- coord_dist(post_estimates[1],post_estimates[2],true_lat,true_lon)
detect_err <- coord_dist(event_lat,event_lon,true_lat,true_lon)

# chain diagnostics

library(coda)

samples.mcmc <- mcmc.list(mcmc(samples[,,1]))

library(ggmcmc)
library(latex2exp)

samples.ggs <- ggs(samples.mcmc, keep_original_order = TRUE)

true_values <- c(true_lat,true_lon,true_depth,true_time,
                 true_alpha,true_pi,true_vP,0.67,true_invlambda)

data_hline <- data.frame(Parameter = as.factor(c("lat0","lon0","d0","t0","alpha","pi","vP","tau","invlambda")),
                         true_values = true_values)

samples.ggs$Parameter <- factor(samples.ggs$Parameter,
                                levels = c("lat0","lon0","d0","t0","alpha","pi","vP","tau","invlambda"),
                                labels = c("lat[0]","lon[0]","d[0]","t[0]","alpha","pi","v[P]","tau","lambda^-1"))

data_hline$Parameter <- factor(data_hline$Parameter,
                               levels = c("lat0","lon0","d0","t0","alpha","pi","vP","tau","invlambda"),
                               labels = c("lat[0]","lon[0]","d[0]","t[0]","alpha","pi","v[P]","tau","lambda^-1"))

x11()
ggs_traceplot(samples.ggs) +
  facet_wrap(~ Parameter, ncol = 3, scales = "free", labeller = label_parsed) +
  geom_hline(data = data_hline, aes(yintercept = true_values), color = 2) +
  theme_bw() + 
  ylab("Value")

x11()
ggs_density(samples.ggs) +
  facet_wrap(~ Parameter, ncol = 3, scales = "free", labeller = label_parsed) +
  geom_vline(data = data_hline, aes(xintercept = true_values), color = 2) +
  theme_bw() +  
  ylab("Density")

x11()
ggs_autocorrelation(samples.ggs) +
  facet_wrap(~ Parameter, ncol = 3, scales = "free", labeller = label_parsed) +
  theme_bw() 

# results visualization

library(leaflet)

icon.true <- makeIcon(iconUrl = "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/epicenter.png",
                      iconWidth = 30, iconHeight = 30)
icon.eqn <- makeIcon(iconUrl = "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/eqn_estimate.png",
                     iconWidth = 30, iconHeight = 30)
icon.post <- makeIcon(iconUrl = "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/post_estimate_4.png",
                      iconWidth = 30, iconHeight = 30)
icon.no_triggered <- makeIcon(iconUrl = "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/not_triggered.png",
                              iconWidth = 12, iconHeight = 12)

library(leaflegend)

threshold <- vector(length = sum(data$delta))

aux_time <- data$time[which(data$delta==1)]

for (i in 1:sum(data$delta)) {
  if (aux_time[i] < 25){
    threshold[i] <- "y < 25"
  } else if (aux_time[i] >= 25 & aux_time[i] < 50){
    threshold[i] <- "25 \u2264 y < 50"
  } else if (aux_time[i] >= 50 & aux_time[i] < 75){
    threshold[i] <- "50 \u2264 y < 75"
  } else if (aux_time[i] >= 75 & aux_time[i] < 100){
    threshold[i] <- "75 \u2264 y < 100"
  } else if (aux_time[i] >= 100){
    threshold[i] <- "100 \u2264 y"
  }
}

threshold <- factor(threshold)
threshold <- relevel(threshold,"75 \u2264 y < 100")
threshold <- relevel(threshold,"50 \u2264 y < 75")
threshold <- relevel(threshold,"25 \u2264 y < 50")
threshold <- relevel(threshold,"y < 25")
threshold

my_colors <- RColorBrewer::brewer.pal(9, "Greens")[5:9]
pal <- colorFactor(my_colors, domain = levels(threshold),ordered = TRUE)

leaflet(data = data,options = list(zoomControl = FALSE,attributionControl=FALSE)) %>%
  addTiles() %>%
  setView(mean(data$lon), mean(data$lat), zoom = 7) %>%
  addMapPane("no_triggers", zIndex = 410) %>%
  addMapPane("triggers", zIndex = 420) %>%
  addMapPane("EQN", zIndex = 430) %>%
  addMapPane("true", zIndex = 440) %>%
  addMapPane("SEQM", zIndex = 450) %>%
  addSymbols(data = data[which(data$delta==0),], lat = ~lat, lng = ~lon, color = "#E14F77FF",
             shape = "plus", width = 12, strokeWidth = 0.1,
             options = pathOptions(pane = "no_triggers")) %>%
  addCircleMarkers(data = data[which(data$delta==1),],
                   radius = 6,
                   color = ~pal(threshold),
                   stroke = FALSE, fillOpacity = 1,
                   options = pathOptions(pane = "triggers")) %>%
  addMarkers(lat = true_lat, lng = true_lon, icon = icon.true,
             options = pathOptions(pane = "true")) %>%
  addMarkers(lat = event_lat, lng = event_lon, icon = icon.eqn,
             options = pathOptions(pane = "EQN")) %>%
  addMarkers(lat = post_estimates[1], lng = post_estimates[2], icon = icon.post,
             options = pathOptions(pane = "SEQM")) %>%
  addLegendFactor(pal = pal, values = ~factor(threshold), opacity = 1,
                  shape = 'circle', title = "Trigger time",
                  position ="topright", width = 15, height = 15) %>%
  addLegendImage(images = c("C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/epicenter.png",
                            "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/eqn_estimate.png",
                            "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/post_estimate_4.png",
                            "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/not_triggered.png"),
                 labels = c("True epicentre","EQN estimate","SEQM estimate","Active smartphone"),
                 labelStyle = "font-size: 14px; vertical-align: middle;",
                 position = "topleft", width = c(20,20,20,15), height = c(20,20,20,15)) %>%
  addScaleBar(position = "bottomright")

# with the posterior plotted

z0.post <- as.data.frame(samples[,1:2,1])
names(z0.post)[1] <- "lat"
names(z0.post)[2] <- "lon"

library(KernSmooth)
kde <- bkde2D(z0.post[ , c(2,1)],
              bandwidth = c(bw.nrd0(z0.post[,2]), bw.nrd0(z0.post[,1])),
              gridsize = c(2500,2500))

CL <- contourLines(kde$x1 , kde$x2 , kde$fhat, nlevels = 50)

## EXTRACT CONTOUR LINE LEVELS
LEVS <- as.factor(sapply(CL, `[[`, "level"))
NLEV <- length(levels(LEVS))

library(sp)
## CONVERT CONTOUR LINES TO POLYGONS
pgons <- lapply(1:length(CL),
                function(i) Polygons(list(Polygon(cbind(CL[[i]]$x, CL[[i]]$y))), ID=i))
spgons = SpatialPolygons(pgons)

pal2 <- colorNumeric(heat.colors(NLEV, NULL)[LEVS], domain = kde$fhat)

leaflet(data = data,options = list(zoomControl = FALSE,attributionControl=FALSE)) %>%
  addTiles() %>%
  setView(post_estimates[2],post_estimates[1], zoom = 9) %>%
  addMapPane("posterior points", zIndex = 400) %>%
  addMapPane("posterior", zIndex = 410) %>%
  addMapPane("no_triggers", zIndex = 420) %>%
  addMapPane("triggers", zIndex = 430) %>%
  addMapPane("EQN", zIndex = 440) %>%
  addMapPane("true", zIndex = 450) %>%
  addMapPane("SEQM", zIndex = 460) %>%
  addCircleMarkers(data = z0.post,
                   radius = 6,
                   stroke = FALSE, fillOpacity = 0.1,
                   options = pathOptions(pane = "posterior points")) %>%
  addPolygons(data = spgons, color = heat.colors(NLEV, NULL)[LEVS],
              options = pathOptions(pane = "posterior")) %>%
  addSymbols(data = data[which(data$delta==0),], lat = ~lat, lng = ~lon, color = "#E14F77FF",
             shape = "plus", width = 12, strokeWidth = 0.1,
             options = pathOptions(pane = "no_triggers")) %>%
  addCircleMarkers(data = data[which(data$delta==1),],
                   radius = 6,
                   color = ~pal(threshold),
                   stroke = FALSE, fillOpacity = 1,
                   options = pathOptions(pane = "triggers")) %>%
  addMarkers(lat = true_lat, lng = true_lon, icon = icon.true,
             options = pathOptions(pane = "true")) %>%
  addMarkers(lat = event_lat, lng = event_lon, icon = icon.eqn,
             options = pathOptions(pane = "EQN")) %>%
  addMarkers(lat = post_estimates[1], lng = post_estimates[2], icon = icon.post,
             options = pathOptions(pane = "SEQM")) %>%
  addLegend(pal = pal2,
            value = kde$fhat,
            title = "Density") %>%
  addLegendImage(images = c("C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/epicenter.png",
                            "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/eqn_estimate.png",
                            "C:/Dati/Dottorato/Aiello_Argiento_Finazzi_Paci/Implementazione/MWG/Fast/California/Images/post_estimate_4.png"),
                 labels = c("True epicentre","EQN estimate","SEQM estimate"),
                 labelStyle = "font-size: 14px; vertical-align: middle;",
                 position = "topleft", width = 20, height = 20) %>%
  addScaleBar(position = "bottomright")

# Posterior predictive checks

time <- 1:120

chain_lat0 <- samples[,1,1]
chain_lon0 <- samples[,2,1]
chain_d0 <- samples[,3,1]
chain_t0 <- samples[,4,1]
chain_alpha <- samples[,5,1]
chain_pi <- samples[,6,1]
chain_vP <- samples[,7,1]
chain_tau <- samples[,8,1]
chain_invlambda <- samples[,9,1]

chain_surv <- list()

for (g in 1:G) {
  
  if ((g / G * 100) %in% seq(10, 100, by = 10)) {
    cat("### Progress:", g / G * 100, " % \n")
  }
  
  mu <- mus(data$lat, data$lon, chain_lat0[g], chain_lon0[g], chain_t0[g], chain_d0[g], chain_vP[g])
  
  mu_P <- mu[, 1]
  mu_S <- mu[, 2]
  
  aux <- sapply(1:dim(data)[1], function(i) {
    exp(-time / chain_invlambda[g]) *
      (chain_pi[g] + (1 - chain_pi[g]) *
         (chain_alpha[g] * pnorm(time, mu_P[i], chain_tau[g], lower.tail = FALSE) +
            (1 - chain_alpha[g]) * pnorm(time, mu_S[i], chain_tau[g], lower.tail = FALSE)))
  })
  
  chain_surv[[g]] <- colMeans(t(aux))
  
}

chain_surv_sub_mean <- do.call(rbind,chain_surv)

# Function to calculate summary statistics for a column
aux_quantiles <- function(column) {
  c(quantile(column, 0.500),
    quantile(column, 0.025),
    quantile(column, 0.975))
}

# Apply the function to each column of the matrix
surv_quantiles <- t(apply(chain_surv_sub_mean, 2, aux_quantiles))
surv_quantiles <- as.data.frame(cbind(time,surv_quantiles))

colnames(surv_quantiles) <- c("time","median","lb","ub")

library(ggsurvfit)

survfit2(Surv(data$time, data$delta) ~ 1) %>%
  ggsurvfit() +
  labs(
    x = "Seconds",
    y = "Standardised survival function"
  ) +
  geom_line(data = surv_quantiles, aes(x=time,y=median), color = 2) +
  geom_vline(aes(xintercept = true_time), color = 3) +
  geom_ribbon(data = surv_quantiles, aes(x=time,ymin=lb,ymax=ub),
              color = 2, fill = 2, alpha = 0.33) +
  scale_x_continuous(breaks = seq(0,120, by = 20))

  
