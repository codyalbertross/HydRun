


# Baseflow separation -----------------------------------------------------


#' Baseflow separation
#'
#' \code{separate.baseflow} is used to separate a stream hydrograph into baseflow and stormflow components.
#' @param hydrograph A stream hydrograph data frame (format is 2 columns: datetime and discharge).
#' @param filter The filter coefficient used in baseflow separation - range 0.9 - 0.95.
#' @param passes The number of filter passes.
#' @return A list containing stormflow and baseflow of \code{hydrograph}.
#' @examples 
#' hydrograph_separated <- separate.baseflow(hydrograph, 0.9, 3)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @references Fuka, D. R., Walter, M. T., Archibald, J. A., Steenhuis, T. S., Easton, Z. M., Fuka, M. D., & KeepSource, T. R. U. E. (2014). Package EcoHydRology.
#' @export

separate.baseflow <- function(hydrograph, filter, passes) {
  
  n <- dim(hydrograph)[1]  # dimensions of initial hydrograph  
  datetime <- hydrograph[,1]  # time vector of initial hydrograh 
  hydrograph <- hydrograph[,2]  # initial hydrograph
  baseflow <- rep(NA, n)  # empty baseflow vector
  baseflow_p <- hydrograph # baseflow from previous pass, baseflow_p is the initial hydrograph to start
  
  for (j in 1:passes) {
    
    #set forward and backward pass
    if ((j%%2) == 1) {
      sidx <- 1
      eidx <- n
      incr <- 1
    }else {
      sidx <- n
      eidx <- 1
      incr <- -1
    }
    
    #set the initial value for baseflow
    baseflow[sidx] <- baseflow_p[sidx]
    
    for (i in seq(from = (sidx+incr),to = eidx, by = incr)) {
      tmp <- filter*baseflow[i-incr] + (1-filter)*(baseflow_p[i]+baseflow_p[i-incr])/2
      baseflow[i] <- min(tmp, baseflow_p[i], na.rm = TRUE)
    }
    
    baseflow_p <- baseflow
  }
  
  # Calculate stormflow from the input hydrograph and the calculated baseflow
  stormflow <- hydrograph - baseflow
  
  # Return a list containing two data frames, one of stormflow, the other of baseflow
  return(list(stormflow = data.frame(datetime,stormflow), baseflow = data.frame(datetime,baseflow)))
  
}



# Extract runoff events ---------------------------------------------------

#' Extract runoff events
#'
#' Extracts runoff events from a stormflow timeseries.
#' @param stormflow A stream stormflow data frame (format is 2 columns: datetime and stormflow).
#' @param min_diff The minimum difference between the start (or end) of runoff and the runoff peak.
#' @param return_ratio Determines where runoff terminates: runoff is terminated when it declines below a dynamic threshold (peak discharge)*\code{return_ratio}.
#' @param b_slope The beginning slope (slope threshold used to cut flat head).
#' @param e_slope The ending slope (slope threshold to end event).
#' @param sc The smoothing coefficient, which determines the number of passes applied to \code{stormflow}.
#' @param min_dur The minimum duration of a runoff event, expressed as the number of timesteps.
#' @param dyslp The dynamic slope threshold used to cut the flat head and end the runoff event.
#' @return A list containing the runoff events, \code{RunoffEvents}, and the number of runoff events, \code{nRunoffEvent}. 
#' @examples 
#' runoff_events <- extract.runoff(stormflow, 1.1, 0.001, 0.001, 0.35, 1, 4, 0.001)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

extract.runoff <- function(stormflow, min_diff, return_ratio, b_slope, e_slope, sc, min_dur, dyslp){
  
  
  return.constant <- min_diff/3
  
  #Step 1: Hydrograph Smoothing
  hydrograph <- stormflow
  hydrograph[,2] <- smooth.curve(hydrograph[,2], sc)
  
  #Step 2: Identify and extract turning points (TP)
  #Identify local max and min (i.e., peaks and valleys)
  TP <- find.tp(hydrograph[,2])  #col1: index, col2: label (valley = 0, peak = 1)
  TP <- cbind(TP, hydrograph[TP[,1],2]) #col3: discharge
  
  #The first element in hydrograph is considered a valley point if it is very low (< average Q / 10)
  if ((TP[1,2] == 1) & (hydrograph[1,2] < (mean(hydrograph[,2],na.rm=TRUE)/10))) {
    TP <- rbind(c(1, 0, hydrograph[1,2]), TP)
  }
  
  #The last element in hydrograph is considered a valley point if it is very low (< average Q / 10)
  if ((TP[dim(TP)[1],2] == 1) & (hydrograph[dim(hydrograph)[1],2] < (mean(hydrograph[,2],na.rm=TRUE)/10))) {
    TP <- rbind(TP, c(dim(hydrograph)[1], 0, hydrograph[dim(hydrograph)[1],2]))
  }
  
  #Remove incomplete event(s) at the beginning and end
  nTP <- dim(TP)[2]
  while (TP[1,2] == 1) {
    TP <- matrix(TP[-1,], ncol = nTP)
  }
  
  while (TP[dim(TP)[1],2] == 1) {
    TP <- matrix(TP[-dim(TP)[1],], ncol = nTP)
  }
  
  #Step 3: Identify the start and end points of runoff events
  #Get differences between peak and valley
  TP <- cbind(TP, c(diff(TP[,3]),0))
  
  #identify start and end points of event flow
  i <- 1
  c <- 1
  st <- vector()
  ed <- vector()
  isComplete <- 1
  nInflect <- dim(TP)[1]
  
  while (i < (nInflect-1)) {
    j <- 1 
    d <- TP[i, 4] + TP[i+j, 4]
    
    while (d > max((return_ratio * max(abs(TP[i:(i+j), 4]), na.rm = TRUE)), return.constant, na.rm = TRUE)) {
      
      if ((i + j) < (nInflect - 1)) {
        j <- j + 1 
        d <- d + TP[i+j, 4]
      }else {
        isComplete <- 0
        break
      }
    }
    
    
    st[c] <- i
    ed[c] <- i + j 
    i <- i + j + 1 
    c <- c + 1
  }
  
  
  
  if (isComplete == 0) {
    st <- st[1:(length(st) - 1)]; 
    ed <- ed[1:(length(ed) - 1)]; 
  }
  
  #Step 4: Extract runoff events and put each of them into a list
  nEvent <- 0
  RunoffEvents <- list()
  
  for (i in 1:length(st)) {
    
    datetime <- stormflow[TP[st[i], 1]:TP[ed[i]+1, 1],1]
    event_smooth <- hydrograph[TP[st[i], 1]:TP[ed[i]+1, 1],2]
    event_stormflow <- stormflow[TP[st[i], 1]:TP[ed[i]+1, 1],2]
    
    event_flow <- data.frame(datetime, event_stormflow, event_smooth)
    temp1 <- diff(event_flow[,3])
    
    #Select runoff events whose peak exceeds threshold (min_diff) 
    if ((max(event_flow[,2],na.rm=TRUE)-event_flow[1,2]) > min_diff &
        (max(event_flow[,2],na.rm=TRUE)-event_flow[dim(event_flow)[1],2]) > min_diff) {
      
      #dynamic slope threshold
      dyslp <- dyslp * range(event_flow[,2], na.rm=TRUE) 
      
      #shorten runoff event by removing flat head and end
      #check slope on smoothed flow data
      while ((length(temp1) > 0) & (temp1[1] < min(c(b_slope,dyslp), na.rm = TRUE))) {
        event_flow <- event_flow[2:dim(event_flow)[1],]
        temp1 <- temp1[2:length(temp1)]
      }
      
      while((length(temp1) > 0) & (temp1[length(temp1)] > -min(c(e_slope,dyslp),na.rm=TRUE))) {
        event_flow <- event_flow[1:(dim(event_flow)[1]-1),]
        temp1 <- temp1[1:(length(temp1)-1)]
      }
      
      #check slope on original flow data
      if (length(temp1) > 0) {
        temp2 <- diff(event_flow[,2])
        
        while ((length(temp2) > 0) & (temp2[1] <= min(c(b_slope,dyslp),na.rm=TRUE))) {
          
          event_flow <- event_flow[2:dim(event_flow)[1],]
          temp2 <- temp2[2:length(temp2)]
        }
        
        while ((length(temp2) > 0) & (temp2[length(temp2)] >= -min(c(e_slope, dyslp),na.rm=TRUE))) {
          
          event_flow <- event_flow[1:(dim(event_flow)[1]-1),]
          temp2 <- temp2[1:(length(temp2)-1)]
        }
        
        #select runoff events whose duration exceeds threshold min_dur
        if (length(temp2) > min_dur) {
          nEvent <- nEvent + 1
          RunoffEvents[[i]] <- event_flow[,1:2]
        }
      }
      
    } 
  }
  
  #remove runoff Events
  if(length(RunoffEvents)>0) {
    RunoffEvents <- RunoffEvents[-which(sapply(RunoffEvents,is.null))]
    names(RunoffEvents) <- paste("Event",1:length(RunoffEvents),sep='')
    nRunoffEvent <- length(RunoffEvents)
  } else {
    RunoffEvents <- NA
    nRunoffEvent <- NA
  }
  
  return(list(RunoffEvents = RunoffEvents, nRunoffEvent = nRunoffEvent))
  
}


#' Smooth runoff response data
#'
#' Smooths runoff response data for runoff event delineation using a moving average filter (window is 3).
#' @param X The response data vector.
#' @param Pass The number of filter passes used on \code{X}.
#' @return A smoothed response vector.
#' @examples 
#' smoothed_response <- smooth.curve(response, 3)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

smooth.curve <- function(X, Pass) {
  
  Y <- X
  
  if (Pass > 0) {
    
    for (i in 1:Pass) {
      for (j in 2:(length(X)-1)) {
        Y[j] <- (Y[j-1] + 2*Y[j] + Y[j+1])/4
      }
    }
  }
  
  return(Y)
  
}

#' Identify turning points in a line
#'
#' Identifies turning/inflection points in a line.
#' @param X The line being assessed for turning points.
#' @return A data frame (format is 2 columns: index of turning points in \code{X} and labels for peaks (1) and valleys (0)).
#' @examples 
#' turning_points <- find.tp(response)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

find.tp <- function(X) {
  
  FOD <- diff(X) #first order derivative
  sDiff <- diff(sign(FOD)) #sign difference
  idx <- which(sDiff != 0) #find the index of the turning points
  
  #Turning points
  c1 <- idx + 1
  c2 <- sDiff[idx]
  
  #Mark local minimums with 0
  c2[c2 > 0] <- 0
  
  #Mark local maximums with 1
  c2[c2 < 0] <- 1
  
  return(cbind(c1,c2))
  
}


# Extract precipitation events --------------------------------------------

#' Extract precipitation events
#'
#' Extracts precipitation events from a precipitation timeseries.
#' @param precip A precipitation data frame (format is 2 columns: datetime and precipitation).
#' @param MAG The maximum allowable gap for an intermittent rainfall, expressed in hours.
#' @param Thr The minimum amount of precipitation.
#' @return A list containing the precipitation events, \code{PrecipEvents} the number of precipitation events, \code{nPrecipEvent}, and the portion of NAs in a rainfall event, \code{pNA}. 
#' @examples 
#' precipitation_events <- extract.precip(precip, 24, 0)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

extract.precip <- function(precip, MAG, Thr){
  
  #Replace NAs with zero
  oprecip <- precip
  precip[is.na(precip[,2]),2] <- 0
  
  #Convert the unit of MAG from hours to the number of time steps
  MAG <- round((MAG/2) / (get.mode(as.numeric(diff(precip[,1]), units = "days")*24))) 
  
  if (MAG < 1) { MAG <- 0 }
  
  mat <- matrix(data = 0, nrow = dim(precip)[1], ncol = 2*MAG+1)
  mat[,(MAG+1)] <- precip[,2]
  
  for (i in 1:MAG) {
    mat[1:(dim(mat)[1]-i), (MAG+1-i)] <- precip[(1+i):dim(precip)[1], 2]
    mat[(1+i):dim(mat)[1], (MAG+1+i)] <- precip[1:(dim(precip)[1]-i), 2]
  }
  
  tmp <- rowSums(mat)
  idx1 <- as.numeric(tmp != 0)
  
  tmp <- diff(idx1)
  tmp <- c(0,tmp,0)
  
  #Complete the possibly incomplete event at the beginning and end
  tmp2 <- tmp[tmp != 0]
  
  if (tmp2[1] == -1) {tmp[1] <- 1}
  if(tmp2[length(tmp2)] == 1) { tmp[length(tmp)] <- -1}
  
  st <- which(tmp == 1)
  ed <- which(tmp == -1)
  lens <- ed-st
  
  precip2 <- precip[as.logical(idx1),]
  PrecipEvents <- list()
  
  for (i in 1:length(lens)) {
    PrecipEvents[[i]] <- precip2[1:lens[i],]
    precip2 <- precip2[-(1:lens[i]),]
  }
  
  PrecipEvents <- lapply(PrecipEvents, function(x) x[(MAG+1):(dim(x)[1]-MAG),])
  
  oprecip2 <- oprecip[as.logical(idx1),]
  oPrecipEvents <- list()
  
  for (i in 1:length(lens)) {
    oPrecipEvents[[i]] <- oprecip2[1:lens[i],]
    oprecip2 <- oprecip2[-(1:lens[i]),]
  }
  
  pNA <- sapply(oPrecipEvents, function(x) sum(is.na(x[,2]))/dim(x)[1])
  
  idx2 <- sapply(PrecipEvents, function(x) sum(x[,2])>Thr)
  
  PrecipEvents <- PrecipEvents[idx2]
  names(PrecipEvents) <- paste("Event",1:length(PrecipEvents),sep='')
  
  pNA <- pNA[as.logical(idx2)]
  
  nEvent <- length(PrecipEvents)
  
  return(list(PrecipEvents = PrecipEvents, nPrecipEvent = nEvent, Precip_pNA = pNA))
}


#' Calculate mode of vector
#'
#' Calculate mode of vector.
#' @param v A vector to calculate the mode of.
#' @return The mode
#' @examples 
#' vector_mode <- get.mode(V)
#' @export

get.mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# Match precipitation and runoff events -----------------------------------

#' Match precipitation and runoff events
#'
#' \code{match.precip.runoff} is used to match precipitation-runoff events.
#' @param PrecipEvents A list containing precipitation events, where each event is a data frame (format is 2 columns: datetime and precipitation).
#' @param RunoffEvents A list containing runoff events, where each event is a data frame (format is 2 columns: datetime and stormflow).
#' @param n The search window. The left edge of the window is n hours before the start of the runoff event.
#' @return A list containing matched precipitation and runoff events.
#' @examples 
#' events <- match.precip.runoff(PrecipEvents, RunoffEvents, n = 96)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @references Fuka, D. R., Walter, M. T., Archibald, J. A., Steenhuis, T. S., Easton, Z. M., Fuka, M. D., & KeepSource, T. R. U. E. (2014). Package EcoHydRology.
#' @export

match.precip.runoff <- function(PrecipEvents, RunoffEvents, n) {
  
  
  stRainEvent <- bind_rows(lapply(PrecipEvents, function(x) x[1,]))[,1]
  edRainEvent <- bind_rows(lapply(PrecipEvents, function(x) x[dim(x)[1],]))[,1]
  
  nRunoffEvent <- length(RunoffEvents)
  
  nMRE <- rep(NA, nRunoffEvent)                           #number of matched precipitation events
  idxMRE <- matrix(NA, nrow = nRunoffEvent, ncol = 20)    #index of matched precipitation events
  
  RR_events <- list()
  
  for (i in 1:nRunoffEvent) {
    
    curr_runoff <- RunoffEvents[[i]]
    curr_itc <- compute.hydro.ITC(curr_runoff)
    
    swb <- curr_itc$start - hours(n)                #the beginning of the search window
    
    tmp <- smooth.curve(curr_runoff[,2], 2)
    infl <- find.tp(tmp)
    idxPeak <- infl[infl[,2] == 1, 1]
    swe <- curr_runoff[idxPeak[length(idxPeak)], 1]      #the end of search window (the last peak of runoff event)
    
    if (i > 1) {
      pre_runoff <- RunoffEvents[[i-1]]
      pre_runoffITC <- compute.hydro.ITC(pre_runoff)
      swb <- max(c(curr_itc$start - hours(n), pre_runoffITC$centroid), na.rm = TRUE)
    }
    
    #Find the precipitation event that overlapped with the runoff event
    col_1 <- as.numeric(sign(stRainEvent - swb) + sign(stRainEvent - swe))
    col_2 <- as.numeric(sign(edRainEvent - swb) + sign(edRainEvent - swe))
    
    idx <- col_1 == 0 | col_2 == 0 | (col_1 * col_2) < 0
    tmp <- which(idx)
    nMRE[i] <- length(tmp)
    
    if (nMRE[i] > 0) { idxMRE[i, 1:nMRE[i]] <- tmp}
    
    curr_precipitation <- PrecipEvents[idx]
    curr_precipitation <- bind_rows(curr_precipitation)
    
    if (length(curr_precipitation) > 1) {
      RR_events[[i]] <- list(precipitation = curr_precipitation, runoff = curr_runoff)
    }else {
      RR_events[[i]] <- list(precipitation = NA, runoff = curr_runoff)
    }
    
  }
  
  return(RR_events)
}


# Event timing characteristics --------------------------------------------

#' Event hydrograph timing characteristics
#'
#' \code{compute.hydro.ITC} is used to calculate timing characteristics of a runoff event storm hydrograph.
#' @param RunoffEvent A data frame containing a runoff event (format is 2 columns: datetime and stormflow).
#' @return A list containing the time of runoff event \code{start}, \code{end}, \code{peak}, and \code{centroid}.
#' @return The unit of time characteristics is in hours.
#' @examples 
#' events <- compute.hydro.ITC(RunoffEvent)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

compute.hydro.ITC = function(RunoffEvent) {
  
  #Find the highest peak
  iPeak <- which(RunoffEvent[,2] == max(RunoffEvent[,2]))
  
  #Identify the time instants for runoff event
  ITC.start <- RunoffEvent[1,1]
  ITC.end <- RunoffEvent[dim(RunoffEvent)[1],1]
  ITC.peak <- RunoffEvent[iPeak,1]
  ITC.centroid <- sum(as.numeric(RunoffEvent[,1])*RunoffEvent[,2], na.rm = TRUE)/sum(RunoffEvent[,2])
  ITC.centroid <- as.POSIXct(ITC.centroid, origin = '1970-01-01', tz = "UTC")
  
  return(list(start = ITC.start, end = ITC.end, peak = ITC.peak, centroid = ITC.centroid))
}

#' Event hyetograph timing characteristics
#'
#' \code{compute.hyeto.ITC} is used to calculate timing characteristics of a precipitation event hyetograph.
#' @param PrecipEvent A data frame containing a precipitation event  event (format is 2 columns: datetime and precipitation).
#' @return A list containing the time of runoff event \code{start}, \code{end}, and \code{centroid}.
#' @return The unit of time characteristics is in hours.
#' @examples 
#' events <- compute.hyeto.ITC(PrecipEvent)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

compute.hyeto.ITC = function(PrecipEvent) {
  
  if (all(is.na(PrecipEvent))) {
    ITC.start <- NA
    ITC.end <- NA
    ITC.centroid <- NA
    
  }else if (sum(PrecipEvent[,2], na.rm=TRUE) == 0) {
    ITC.start <- NA
    ITC.end <- NA
    ITC.centroid <- NA
  }else{
    ITC.start <- PrecipEvent[1,1]
    ITC.end <- PrecipEvent[dim(PrecipEvent)[1],1]
    ITC.centroid <- sum(as.numeric(PrecipEvent[,1])*PrecipEvent[,2], na.rm = TRUE)/sum(PrecipEvent[,2])
    ITC.centroid <- as.POSIXct(ITC.centroid, origin = '1970-01-01', tz = "UTC")
  }
  
  return(list(start = ITC.start, end = ITC.end, centroid = ITC.centroid))
  
}


#' Rainfall-runoff event timing characteristics
#'
#' \code{compute.TC} is used to calculate timing characteristics of a rainfall-runoff event.
#' @param PrecipEvent A data frame containing a precipitation event  event (format is 2 columns: datetime and precipitation).
#' @param RunoffEvent A data frame containing a runoff event  event (format is 2 columns: datetime and stormflow).
#' @return A list containing the rainfall duration \code{Tw}, response lag \code{TLR}, 
#' time of rise \code{Tr}, lag-to-peak \code{TLP}, centroid lag-to-peak \code{TLPC}, centroid lag \code{TLC},
#' time base \code{Tb}.
#' @return The unit of time characteristics is in hours.
#' @examples 
#' events <- compute.TC(PrecipEvent, RunoffEvent)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

compute.TC = function(PrecipEvent, RunoffEvent) {
  
  rainITC <- compute.hyeto.ITC(PrecipEvent) 
  runoffITC <- compute.hydro.ITC(RunoffEvent) 
  
  TimeCharac = list()
  
  TimeCharac[[1]] <- as.numeric(difftime(rainITC$end, rainITC$start, units = "hours"))
  TimeCharac[[2]] <- as.numeric(difftime(runoffITC$start, rainITC$start, units = "hours"))
  TimeCharac[[3]] <- as.numeric(difftime(runoffITC$peak, runoffITC$start, units = "hours"))
  TimeCharac[[4]] <- as.numeric(difftime(runoffITC$peak, rainITC$start, units = "hours"))
  TimeCharac[[5]] <- as.numeric(difftime(runoffITC$peak, rainITC$centroid, units = "hours"))
  TimeCharac[[6]] <- as.numeric(difftime(runoffITC$centroid, rainITC$centroid, units = "hours"))
  TimeCharac[[7]] <- as.numeric(difftime(runoffITC$end, runoffITC$start, units = "hours"))
  TimeCharac[[8]] <- as.numeric(difftime(runoffITC$end, rainITC$end, units = "hours"))
  
  names(TimeCharac) <- c("Tw","TLR","Tr","TLP","TLPC","TLC","Tb","Tc")
  
  
  return(TimeCharac)
}

# Other event characteristics ---------------------------------------------

#' Event runoff ratio
#'
#' \code{compute.RR} is used to calculate the runoff ratio of a rainfall-runoff event
#' @param PrecipEvent A data frame containing a precipitation event  event (format is 2 columns: datetime and precipitation).
#' @param RunoffEvent A data frame containing a runoff event  event (format is 2 columns: datetime and stormflow).
#' @param drainageArea The drainage area of the catchment.
#' @return A list containing the event runoff ratio \code{runoffRatio}, precipitation volume \code{precip_vol}, 
#' and runoff volume \code{runoff_vol}.
#' @return The unit of the returned volumes is in m^3
#' @examples 
#' events <- compute.hydro.ITC(RunoffEvent)
#' @note This is an R implementation of a function from the MATLAB toolbox HydRun.
#' @note It is assumed that the input precipitation is measured in mm, runoff in m^3/s, and drainage area in km^2
#' @references Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.
#' @export

compute.RR = function(PrecipEvent, RunoffEvent, drainageArea) {
  
  if (!all(is.na(RunoffEvent))) {
    
    interval <- as.numeric(difftime(RunoffEvent[2,1], RunoffEvent[1,1], units = "sec"))
    runoffVol <- sum(RunoffEvent[,2] * interval, na.rm=TRUE)
    
  }else {
    runoffVol <- NA
  }
  
  #calculate total volume of rainfall in m^3
  
  if (!all(is.na(PrecipEvent))) {
    precipVol <- sum(PrecipEvent[,2], na.rm = TRUE) * drainageArea * 1000
    
  }else {
    precipVol <- NA
  }
  
  #calculate the runoff ratio 
  runoffRatio = runoffVol/precipVol
  
  return(list(RR = runoffRatio, precip_vol = precipVol, runoff_vol = runoffVol))
}