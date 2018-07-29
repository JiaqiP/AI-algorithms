 
#Group 16
#shortest_total_route

generateHeuristics = function (offset,packages,goal){
  #use manhantten distance
  h = matrix(nrow = 10, ncol = 10)
  for (i in 1:10){
    for(j in 1:10){
      h[i,j]=(abs(i-packages[goal,1+offset]) + abs(j-packages[goal,2+offset]))
    }
  }
  return (h)
}

#returns the total cost and the next node
astar = function(offset,goal,roads,car,packages){
  
  visited <- list()
  frontier <- list()
  neighbours <- list()
  
  #calculate heuristic
  h = generateHeuristics(offset,packages,goal)
  
  #add the starting point first
  start = list(x = car$x, y = car$y, parent = NULL, f = h[car$x,car$y], g = 0, heuristic = h[car$x,car$y])
  frontier[[1]]<-start
  
  while (length(frontier)!=0 ){
    #take a node with shortest f from frontier
    frontier = frontier[order(sapply(frontier,'[[',4))]    #sort according to f value
    
    #delete it from frontier add it to the visited
    chosen <- frontier[[1]]
    visited[[length(visited)+1]] <- chosen
    frontier[[1]] = NULL
    
    
    #if goal is in visited then terminate
    if((chosen$x == packages[goal,offset+1]) && (chosen$y == packages[goal,offset+2]))
      break;
    
    #expand the node, set its neighbour(not in visited) into frontier
    
    
    #find neighbours, calculate all the attributes of the neighbour
    if(chosen$x-1 > 0){
      g = roads$hroads[chosen$y,chosen$x-1] + chosen$g
      heuristic = h[chosen$x-1,chosen$y]
      f = g + heuristic
      neighbours[[length(neighbours)+1]] = list(x = chosen$x-1, y = chosen$y, parent = chosen, f = f, g = g, heuristic = heuristic )
    }
    
    if(chosen$x+1 <= 10){
      g = roads$hroads[chosen$y,chosen$x] + chosen$g
      heuristic = h[chosen$x+1,chosen$y]
      f = g + heuristic
      neighbours[[length(neighbours)+1]] = list(x = chosen$x+1, y = chosen$y, parent = chosen, f = f, g = g, heuristic = heuristic )
    }
    
    if(chosen$y-1 > 0){
      g = roads$vroads[chosen$y-1,chosen$x] + chosen$g
      heuristic = h[chosen$x,chosen$y-1]
      f = g + heuristic
      neighbours[[length(neighbours)+1]] = list(x = chosen$x, y = chosen$y-1, parent = chosen, f = f, g = g, heuristic = heuristic )
    }
    
    if(chosen$y+1 <= 10){
      g = roads$vroads[chosen$y,chosen$x] + chosen$g
      heuristic = h[chosen$x,chosen$y+1]
      f = g + heuristic
      neighbours[[length(neighbours)+1]] = list(x = chosen$x, y = chosen$y+1, parent = chosen, f = f, g = g, heuristic = heuristic )
    }
    
    
    #for each of the neighbouring nodes,check if it is already in visited or frontier
    
    for(node in neighbours){
      exist = 0
      #check if it is visited
      for(j in visited){
        if((node$x == j$x) && (node$y == j$y)){
          exist = 1;
          break;
        }
      }
      for(i in frontier){
        #check if it is in frontier
        if((node$x == i$x) && (node$y == i$y)){
          exist = 1
          #if has a lower cost, replace
          if(node$f < i$f)
            node = i
          break;
        }
      }
      #else not in, add to frontier
      if(exist == 0)
        frontier[[length(frontier)+1]] = node;
    }
    
  }
  
  #find the next step by backtracking to the start 
  #now chosen is the goal
  p = chosen
  if(!identical(p,start)){
    while(!identical(p$parent, start)){
      p = p$parent
    }
  }
  return (c(chosen$f,p))
  
}

#generate all the possibilities of delivery order
###REFERANCE: the permutations function is from StackOverflow######
permutations <- function( x, prefix = c() )
{
  if(length(x) == 0 ) return(prefix)
  do.call(rbind, sapply(1:length(x), FUN = function(idx) permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
}

findCost <- function(order,packages,car){
  cost = abs(packages[order[1],1]-car$x)+abs(packages[order[1],2]-car$y)
  #package delivery point i to i+1 package pick up point
  for(i in 1:(length(order)-1)){
    cost = cost + abs(packages[order[i],3]-packages[order[i+1],1]) + abs(packages[order[i],4]-packages[order[i+1],2])
  }
  return(cost)
}


#find best order for delivery
generateOrder <- function(packages,car){
  toPick = which(packages[,5]==0)
  if(length(toPick) <= 1)
  {
    return(c(toPick))
  }
  
  allOrder = permutations(toPick)
  
  costVector = c()
  for(i in 1:nrow(allOrder)){
    costVector = rbind(costVector, findCost(allOrder[i,],packages,car))
  }
  allOrder = cbind(allOrder,costVector)
  decision = allOrder[which.min(allOrder[,ncol(allOrder)]),]
  return(decision)
}

AStarDM = function(roads, car, packages){
  nextMove=0
  toGo=0
  offset=0
  p <- list()
  
  if (car$load==0) {
    order = generateOrder(packages,car)
    toGo = order[1]
  } 
  else {
    if(toGo != car$load)
      toGo=car$load  
    offset=2
  }
  
  p = astar(offset,toGo,roads,car,packages)
  
  if (car$x<p$x) {nextMove=6}
  else if (car$x>p$x) {nextMove=4}
  else if (car$y<p$y) {nextMove=8}
  else if (car$y>p$y) {nextMove=2}
  else {nextMove=5}
  
  
  car$nextMove=nextMove
  car$mem=list()
  return (car)
  
}
