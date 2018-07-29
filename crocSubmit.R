#author Jiaqi Pan
#date£º2017-9-30


#only need to predict current state, backwards algorithm is not needed
#And given current state, there is no need to memorize the markov chain before


#find the state after one transition
transition=function(state,edges){
  #find transition matrix
  tranM = matrix(0,nrow = 40,ncol = 40)
  for(i in 1:nrow(tranM)){
    nextNode = getOptions(i,edges)
    tranM[i,nextNode] = 1/length(nextNode)
  }
  #matrix multiplication
  state = state %*% tranM
  return(state)
  
}

#find the state after multiplied by emission
emission=function(state,probs,readings){
  #use readings and probs find out the probablity of in one hole.
  emisProb<- c(rep(0,40))
  silinityP <- 0
  phosphateP <- 0
  nitrogenP <- 0
  
  #the emission of HMM is the readings, calculate the possibility of generate this readings in different holes
  for(i in 1:40){
    salinityP = dnorm(readings[1], mean=probs$salinity[i,1], sd=probs$salinity[i,2])
    phosphateP = dnorm(readings[2], mean=probs$phosphate[i,1], sd=probs$phosphate[i,2])
    nitrogenP = dnorm(readings[3], mean=probs$nitrogen[i,1], sd=probs$nitrogen[i,2])
    emisProb[i] = salinityP * phosphateP * nitrogenP
  }
  
  #multiply by state and then figure out current croc position probablity
  
  return(state*emisProb)
  
}

initialState= function(readings,state,probs){
  #initial state itself has its emission
  state<-c(rep(1,40))
  state = emission(state,probs,readings)
  state = normalizeState(state)
  return(state)
}

modifyState= function(state,positions){

  #************swedes emit information*************
  if(!is.na(positions[1])){
    if(positions[1] < 0){
      state = rep(0,40)
      state[positions[1]*(-1)] = 1;
    }else{
      state[positions[1]] = 0;
    }
  }  
  if(!is.na(positions[2])){
    if(positions[2] < 0){
      state = rep(0,40)
      state[positions[2]*(-1)] = 1;
    }
    else{
      state[positions[2]] = 0;
    }
  }
 
  return(state)
}

normalizeState=function(state){
  sum = sum(state[1:length(state)])
  for(i in 1: length(state)){
    state[i] = state[i]/sum
  }
  return(state)
}

findMostPossible=function(state){
  #find the most probablity place of croc and confidence for the inference
  likelyCroc = which.max(state)
  confidence = state[likelyCroc]
  return(c(likelyCroc,confidence))
}


findPath= function(positions,croc,edges){
  #find the shortest path to croc 
  #Because all the edges have same weight,so choose breadth first search
  visited<-c()
  start<-list(id = positions[3],parent = NULL)
  queue<-list()
  queue[[1]] = start;
  chosen = queue[[1]]
  
  while(chosen$id != croc[1]){
    visited[length(visited)+1] = chosen$id
    neighbor = c(edges[which(edges[,1]==chosen$id),2],edges[which(edges[,2]==chosen$id),1])
    for(i in neighbor){
      #if i not in visited
      if(i %in% visited == FALSE){
        #put in queue
        queue[[length(queue)+1]] = list(id = i,parent = chosen)
      }
    }
    queue[[1]]=NULL
    chosen = queue[[1]]
  }
  #now chosen$id = croc, backtrack to next step;
  p<-list()
  p = chosen
  path<-c()
  
  #if croc is where I am, p$parent = NULL
  if(!is.null(p$parent)){
    while(p$parent$id != start$id){
      path = cbind(p$id,path)
      p = p$parent
    }
  }
  #when path length is 1, that is, croc is next to me, path is empty
  #path[1] is the next next edge to traverse
  if(length(path) >= 1){
    return(c(p$id,path[1]))
  }
  else{
    return(c(p$id,0))
  }
}

HMMprocess=function(moveInfo,readings,positions,edges,probs,state){
#state is the current state
#efore search the croc, have readings from croc to infer where croc is now 
  state = transition(state,edges)
  state = emission(state,probs,readings)
  state = modifyState(state,positions)
  state = normalizeState(state)
  toGo = findMostPossible(state)
  
  nextStep = findPath(positions,toGo[1],edges)
  return(list(state,nextStep))

}


#' @export
MarkovWC=function(moveInfo,readings,positions,edges,probs) {
#use moveInfo$mem[[1]] to store last state
#moveInfo$mem[[2]] to store last position searched
  
  if(length(moveInfo$mem)==0){
    state <- c()
    state = initialState(readings,state,probs)
    moveInfo$mem[[1]] = state
    croc = findMostPossible(state)
    nextStep = findPath(positions,croc,edges)
    moveInfo$moves=nextStep
    
  }else{
    stateStep = HMMprocess(moveInfo,readings,positions,edges,probs,moveInfo$mem[[1]])
    moveInfo$mem[[1]] = stateStep[[1]]
    moveInfo$moves= stateStep[[2]]
    
  }
  return(moveInfo)
}
