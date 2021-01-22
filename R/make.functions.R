make.exp.decay<-function(t=1,s=1){
  if( !(is.numeric(t) && is.numeric(s)) ){
    stop("t and s both must be numbers")
  }
  if(t<=0){
    stop("t must be positive")
  }
  myfunction<-function(distance){
    return(exp( -((distance^s)/t) ))
  }
  return(myfunction)
}

make.linear.decay<-function(threshold){
  if( !is.numeric(threshold) ){
    stop("threshold must be a number")
  }
  if(threshold<=0){
    stop("threshold must be positive")
  }
  myfunction<-function(distance){
    if(distance>threshold){
      return(0)
    } else{
      return(1-(distance/threshold))
    }
  }
}

make.stepwise<-function(rates,distances){
  if( !(is.numeric(rates) && is.numeric(distances))){
    stop("rates and distances both must be vectors of numbers")
  }
  if(length(rates) != length(distances) ){
    stop("the length of the rate vector should be equal to the length of the distance vector")
  }
  if(sum(rates<0)!=0){
    stop("Rates should all be nonnegative")
  }
  if(sum(distances<=0)!=0){
    stop("Distances should all be positive")
  }

  myfunction<-function(distance){
    x<-which(distances>=distance)
    if(length(x)==0){
      stop(paste("Rates are defined on the range 0 to ",max(distances),". You gave a distance of ",distance,sep = "" ))
    }
    return(rates[ min(which(distances>=distance))])
  }
  return(myfunction)
}

make.polynomial.decay<-function(threshold,degree=1){
  if( !is.numeric(threshold) ){
    stop("threshold must be a number")
  }
  if(threshold<=0){
    stop("threshold must be positive")
  }
  if(degree<=0){
    stop("degree must be positive")
  }
  myfunction<-function(distance){
    if(distance>threshold){
      return(0)
    } else{
      return(1- ((distance/threshold)^degree) )
    }
  }
}


