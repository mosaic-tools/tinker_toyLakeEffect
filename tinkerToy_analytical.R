# tinkerToy
# script containing support function for running simple, equilibrium model of scaling in hydrologic network
# augmented with functions to think about stream order from NHD output and add lakes to a stream network
# mosaics working group

# this package relies on information about each reach in a network where a reach is a stream segment, lake, or reservoir between confluences
# information needed per reach includes:
#   reachNumber: numerical code for reach
#   reachID: id placeholder for lakeCat and streamCat COMIDs or what have you
#   length: length of reach (m)
#   width: width of reach (m)
#   depth: depth of reach (m)
#   localArea: local contributing area (m2)
#   localPET: water yield (precip - evapotranspiration) for local contributing area (m d-1)
#   localCin: carbon concentration in water exiting local contributing area (gC m-3)
#   d: a organic carbon decay rate (d-1)
#   connect: the other reaches that directly flow into this reach; if none use 0, if multiple separate with commas
# this information is passed to the function as a dataframe with a row per reach and the above columns with those names


sortNetwork_HD <- function(df) {
  library(dplyr)
  df = df %>% arrange(connect) #could do this outside of this function. But speeds up loop considerably!
  
  for (i in 1:nrow(df)) {
    # print(i)
    repeat{
      if (df$connect[i] == 0) {break} #first order streams at top
      #check for upstream reaches lower in the table. If none break loop
      if (!any(unlist(strsplit(df$connect[i],",")) %in% df$reachID[i:nrow(df)])) {break} 
      df = bind_rows(slice(df,-i), slice(df,i)) #put row at bottom, and repeat
    }
  }
  print("The network sorted!")
  return(df)
}

sortNetwork<-function(network, counterMax=200,verbose=FALSE){
  # arguments
  # network - setup dataframe for network as described above
  # counterMax - number of iterations that can be used to sort the network for cumulative flow calculation
  
  # add check to make sure all contributing reachIDs are in the network
  allContributing=unlist(strsplit(network$connect,","))
  missingReaches=allContributing[!allContributing%in%c(0,network$reachID)]
  if(length(missingReaches)>0){
    print(paste("There are",length(missingReaches),"contributing reaches that are not included in your network!"))
    # could add a step that automatically deleted these missingReaches from the connect field
  }
  
  ####sort network to allow accumulation of flow across the network
  network=network[order(network$connect),]
  sorted=network
  sorted$reachID=0
  sorted$connect="0"
  # fail safe limit on number of attempts to sort the network
  counter=1
  # move headwater systems to the top
  Nsorted=sum(network$connect=="0")
  sorted[1:Nsorted,]=network[network$connect==0,]
  if(verbose){
    print("iteration - remaining to sort")
  }
  while((Nsorted<nrow(network) & counter<counterMax)){
    # move reaches that are filled from previously sorted reaches into sorted
    remaining=network[!(network$reachID%in%sorted$reachID),]
    if(verbose){
      print(paste(counter,nrow(remaining),sep=" - "))
    }
    for(i in 1:nrow(remaining)){
      if(all(unlist(strsplit(remaining$connect[i],","))%in%sorted$reachID)){
        sorted[(Nsorted+1),]=remaining[i,]
        Nsorted=Nsorted+1
      }
    }
    # increment fail safe counter
    counter=counter+1
  }
  
  if(counter==counterMax){
    print("The network did not sort in the maximum allowed iterations. If your network is very large, increase counterMax. If it isn't, check to make sure you have your connections right.")
  }else{
    print("The network sorted!")
    return(sorted)
  }
}

 
solveNetwork<-function(network){ 
  #arguments
  # network - a SORTED network data frame
  
  # would be good to add check to be sure network is sorted and if not throw an error
  
  
    #Define number of reaches
    nReaches <- nrow(network)
    #Define length of each reach (m)
    length <- network$length
    #Define width  of each reach (m)
    width <- network$width
    #Define depth of each reach (m)
    depth <- network$depth
    #Define volume of each reach (m3)
    V <- length*width*depth
    #Define area of each reach (m2)
    A <- length*width
    #Define local watershed inflows (m3 d-1)
    qL <- network$localArea*network$localPET
    #Define local watershed carbon concentrations (g m-3)
    cL <- network$localCin
    #Define decay rates (d-1)
    d <- network$d
    #Define network inflows  (m3 d-1)
    qN <- rep(0,nReaches)
    #Define carbon loads from upstream network (g d-1)
    qcN <- rep(0,nReaches)
    
    C=numeric(nrow(network))
    for(i in 1:nrow(network)){
      #calculate flows from upstream network
      donors=which(network$reachID%in%unlist(strsplit(network$connect[i],",")))
      if(length(donors)>0){
        qN[i]=sum(qL[donors])+sum(qN[donors])
        qcN[i]=sum((qL[donors]+qN[donors])*C[donors]/V[donors])
      }

      #state variable C is mass of carbon in the reach, g C
      #q is flow, m3 d-1
      #c is concentration of carbon, g m-3
      #L subscript denotes local inputs from the immediate watershed of a reach
      #N subscript denotes network inputs from all of the upstream watershed of a reach, not including immediate watershed
      #d is decay rate, d-1
      C[i]=(qL[i]*cL[i]+qcN[i])/(d[i]+(qL[i]+qN[i])/V[i])
    }
    
    #Calculate concentrations in each reach
    CConc <- C/V
    
    #What is the total carbon export from the bottom of the network? (g d-1)
    COut <- C[nReaches]/V[nReaches]*(qN[nReaches]+qL[nReaches])
    
    #What was the total carbon input?
    CIn <- sum(cL*qL)
    
    #Fraction of C processed
    CLost <- (1-COut/CIn)
    
    #For each reach, what was carbon export?
    COutReach <- C/V*(qN+qL)
    
    #For each reach, what was the carbon input?
    CInReach <- qL*cL + qcN
  
    #For each reach, what was the fraction of C processed/lost
    CLostReach <- (1-COutReach/CInReach)
  
    #Residence times
    resTimeReach <- V/(qL+qN)
  
    resTime <- sum(V)/(qL[nReaches]+qN[nReaches])  # NOT SURE THIS WORKS
  
    #Respiration rate in each reach (g m-2 d-1)
    respReach <- d*C/A
  
    #C spiraling length (m)
    spiralC <- (qN+qL)*(C/V)/(respReach*width)
  
    #C uptake velocity (m d-1)
    vF <- respReach/(C/V)
  
    #Summary
    reachSummary <- data.frame(length,width,depth,V,qL,qN,resTimeReach,d,cL,CInReach,COutReach,CLostReach,C=C,CConc=CConc,respReach,spiralC,vF)
    networkSummary <- list(V=sum(V),resTime=resTime,CIn=CIn,COut=COut,CLost=CLost)
  
  
    return(list(network=network,C=C,CConc=CConc,reachSummary=reachSummary,networkSummary=networkSummary))
}


solveNetwork_wcsed<-function(network){ 
  #arguments
  # network - a SORTED network data frame
  
  # would be good to add check to be sure network is sorted and if not throw an error
  
  
  #Define number of reaches
  nReaches <- nrow(network)
  #Define length of each reach (m)
  length <- network$length
  #Define width  of each reach (m)
  width <- network$width
  #Define depth of each reach (m)
  depth <- network$depth
  #Define volume of each reach (m3)
  V <- length*width*depth
  #Define area of each reach (m2)
  A <- length*width
  #Define local watershed inflows (m3 d-1)
  qL <- network$localArea*network$localPET
  #Define local watershed carbon concentrations (g m-3)
  cL <- network$localCin
  #Define network inflows  (m3 d-1)
  qN <- rep(0,nReaches)
  #Define carbon loads from upstream network (g d-1)
  qcN <- rep(0,nReaches)
  #Define reach hydraulic load
  Hl=(qL+qN)/A
  #Define aggregate water column and sediment decay rates (d-1)
  d <- network$d+network$vF/depth
  
  C=numeric(nrow(network))
  for(i in 1:nrow(network)){
    #calculate flows from upstream network
    donors=which(network$reachID%in%unlist(strsplit(network$connect[i],",")))
    if(length(donors)>0){
      qN[i]=sum(qL[donors])+sum(qN[donors])
      qcN[i]=sum((qL[donors]+qN[donors])*C[donors]/V[donors])
    }
    
    #state variable C is mass of carbon in the reach, g C
    #q is flow, m3 d-1
    #c is concentration of carbon, g m-3
    #L subscript denotes local inputs from the immediate watershed of a reach
    #N subscript denotes network inputs from all of the upstream watershed of a reach, not including immediate watershed
    #d is decay rate, d-1
    C[i]=(qL[i]*cL[i]+qcN[i])/(d[i]+(qL[i]+qN[i])/V[i])
  }
  
  #Calculate concentrations in each reach
  CConc <- C/V
  
  #What is the total carbon export from the bottom of the network? (g d-1)
  COut <- C[nReaches]/V[nReaches]*(qN[nReaches]+qL[nReaches])
  
  #What was the total carbon input?
  CIn <- sum(cL*qL)
  
  #Fraction of C processed
  CLost <- (1-COut/CIn)
  
  #For each reach, what was carbon export?
  COutReach <- C/V*(qN+qL)
  
  #For each reach, what was the carbon input?
  CInReach <- qL*cL + qcN
  
  #For each reach, what was the fraction of C processed/lost
  CLostReach <- (1-COutReach/CInReach)
  
  #Residence times
  resTimeReach <- V/(qL+qN)
  
  resTime <- sum(V)/(qL[nReaches]+qN[nReaches])  # NOT SURE THIS WORKS
  
  #Respiration rate in each reach (g m-2 d-1)
  respReach <- d*C/A
  
  #C spiraling length (m)
  spiralC <- (qN+qL)*(C/V)/(respReach*width)
  
  #C uptake velocity (m d-1)
  vF <- respReach/(C/V)
  
  #Summary
  reachSummary <- data.frame(length,width,depth,V,qL,qN,resTimeReach,d,cL,CInReach,COutReach,CLostReach,C=C,CConc=CConc,respReach,spiralC,vF)
  networkSummary <- list(V=sum(V),resTime=resTime,CIn=CIn,COut=COut,CLost=CLost)
  
  
  return(list(network=network,C=C,CConc=CConc,reachSummary=reachSummary,networkSummary=networkSummary))
}

# function to generate count of reaches of each order in a network
networkOrderCount<-function(network){
  return(table(c(network$order,1:max(network$order)))-1)
}

# function to add random lakes to a network
addRandomLake<-function(network,order,area,depth,lake_d=0.005,verbose=FALSE){
  # ARGUMENTS
  # network:  network to add lakes to
  # order: vector of order streams you'd like to convert to lakes
  # area: desired area (m2) of lakes you'll be adding
  # depth: desired depth (m) of lakes you'll be adding
  # lake_d: decay rate of solute in a lake (d-1)
  # verbose: print out changes made
  
  # make sure width and depth are the same length as order, if not repeat what is supplied to make same length
  if(length(area)<length(order)){
    area=rep(area,length.out=length(order))
  }
  if(length(depth)<length(order)){
    depth=rep(depth,length.out=length(order))
  }
  # new network dataframe to be returned with lakes added
  newNet=network
  
  # summary of network (number of reaches of each order)
  networkOrderCount=table(c(network$order,1:max(c(network$order,order))))-1
  
  # number of reaches of various orders to be changed to lakes
  toChangeCount=table(c(order,1:max(c(network$order,order))))-1
  
  # make random changes
  for(i in 1:length(toChangeCount)){
    # if we're changing any order i reaches to lakes
    if(toChangeCount[i]>0){
      # if there are enough order i reaches to make the desired number of lakes
      if(toChangeCount[i]<networkOrderCount[i]){
        # randomly select desired number of reaches of order i to turn into lakes by setting width, depth, and d to values passed as arguments
        toChange=sample((1:nrow(network))[network$order==names(toChangeCount)[i]],toChangeCount[i],replace=FALSE)
        newNet$width[toChange]=area[order==names(toChangeCount)[i]]/network$length[toChange]
        newNet$depth[toChange]=depth[order==names(toChangeCount)[i]]
        newNet$d[toChange]=lake_d
        newNet$lake[toChange]=1
        if(verbose){
          print(paste(toChangeCount[i]," of ",networkOrderCount[i]," order ",names(networkOrderCount)[i]," changed to lakes.",sep=""))
        }
      # if there are NOT enough order i reaches to make the desired number of lakes, change all of them and print warning if verbose is TRUE 
      }else{
        toChange=(1:nrow(network))[network$order==names(toChangeCount)[i]]
        newNet$width[toChange]=area[order==names(toChangeCount)[i]]/network$length[toChange]
        newNet$depth[toChange]=depth[order==names(toChangeCount)[i]]
        newNet$d[toChange]=lake_d
        newNet$lake[toChange]=1
        
        print(paste("ONLY ",length(toChange)," of ",networkOrderCount[i]," order ",names(networkOrderCount)[i]," changed to lakes.",sep=""))
      }
    # if we are NOT changing any order i reaches to lakes  
    }else{
      if(verbose){
        print(paste("0 of ",networkOrderCount[i]," order ",names(networkOrderCount)[i]," changed to lakes.",sep=""))
      }
    }
  }
  return(newNet)
}
