rm( list=ls())

# Define Variables about Simulation
iterations <- 4*10^6 # Times to read ensemble
itperread <- 100 # Number of steps per ensemble read
b.move <- 0.2 # Probability of displacement
b.add <- 0.5 # Probability of addition
delx <- c(5, 5, 5)  # Average change in position [=] A

graphene <- function( n ){
  i = 0:(n-1)
  x1 = 123 + i*246
  x2 = i*246
  
  j = 0:n
  y1 = j*(142+71)
  y2 = 71+j*(142+71)
  
  x = c()
  for (i in 1:(n+1)){
    if( i%%2 == 1){
      x = c( x, x1,x2)
    }else{
      x = c( x, x2,x1)
    }
  }
  x = x/100
  
  y = c()
  for (i in 1:(n+1)){
    y= c( y, rep(y1[i],n), rep(y2[i],n))
  }
  y = y/100
  z <- rep( 0, length(x))
  X.c <- matrix( c(x,y,z), ncol = 3)
  
  return(X.c)
}  


# Define Variables about System

X.c <- graphene( 17 )

Space <- c(max(X.c[,1]), max(X.c[,2]), 70) # Area being worked with A
V <- Space[1]*Space[2]*Space[3] # Volume of adsorbtion spaced

# Matrix of positions in x,y,z for starting hydrogen particles 
X.h <- Space*matrix( c(0.5,0.5,0.6,
                       0.1,0.5,0.6,
                       0.2,0.7,0.4), nrow=3,ncol=3, byrow=TRUE) 
# Lennard-Jones Parameters for h2-h2
s.hh <- 2.958 # Sigma A 
e.hh <- 0#36.7 #E/kb K

# Lennard-Jones Parameters for h2-C
s.hc <- 3.216 # m
e.hc <- 41.924 #K/kb K

# Lennard-Jones Parameters for Ar-Ar
s.aa <- 3.4
e.aa <- 120.0

# Lennard-Jones Parameters for C-C (graphene)
s.cc <- 3.4
e.cc <- 28.0

# Lennard-Jones Parameters for Ar-C
s.ac <- (s.aa + s.cc)/2
e.ac <- sqrt( e.aa*e.cc)

gamma.A <- -0.54
gamma.R <- 0.38
  
s.min <- 3.60984

s.xx <- s.aa
e.xx <- e.aa
s.xc <- s.ac
e.xc <- e.ac


library(emdbook)
#P <- lseq( 1e-5, 1e3, 100) #Pa
P <- 10000
T <- 77
V <- Space[1]*Space[2]*Space[3] # A
kb <- 1.38064852*10^(-23) #J/k
lilkb <- 1.38064852



dist <- function( object, others, Space){
  # Distance in Angstroms
  dx <- ifelse( abs(object[1] - others[,1]) >= Space[1]/2,
                abs(object[1] - others[,1])- Space[1]/2,abs(object[1] - others[,1])  )
  dy <- ifelse( abs(object[2] - others[,2]) >= Space[2]/2,
                abs(object[2] - others[,2])- Space[2]/2,abs(object[2] - others[,2])  )
  dz <- ifelse( abs(object[3] - others[,3]) >= Space[3]/2,
                abs(object[3] - others[,3])- Space[3]/2,abs(object[3] - others[,3])  )
  dist <- sqrt( dx^2 + dy^2 + dz^2 )
  angle <- acos( (c(0,0,1)%*%t(matrix(c(dx,dy,dz),ncol=3)))/dist )
  
  return(matrix(c(dist,angle), ncol=2, byrow = FALSE))
}
LJ.c <- function( dist , epsilon , sigma){
  # Potential Energy/kb in K
  U <- 4*epsilon*( (sigma/dist[1])^12*(1 + gamma.R*(1 - 6/5*cos(dist[2])^2)) -
                   (sigma/dist[1])^6*(1 + gamma.A*(1 - 3/2*cos(dist[2])^2)) )
  return(U)
}
LJ <- function( dist , epsilon , sigma, Space){
  # Potential Energy/kb in K
  U <- 4*epsilon*( (sigma/dist)^12 - (sigma/dist)^6 )
  return(U)
}
move <- function( X.h , X.c , delx , Space , T){
  
  p <- 0
  if ( nrow(X.h) >= 3 ){
    id.o <- as.integer(runif(1)*nrow(X.h)) + 1 # Choose a particle to move
    
    U.oh <- LJ( dist(X.h[id.o,], X.h[-id.o,], Space), e.xx, s.xx) # K
    U.oc <- LJ.c( dist(X.h[id.o,], X.c, Space), e.xc, s.xc ) # K
    U.o <- sum(U.oh) + sum(U.oc)  # K
    
    X.n <- (X.h[id.o,] + runif(3, 0, 1)*delx)%%Space #Potential new Location
    
    U.nh <- LJ( dist(X.n, X.h[-id.o,],Space), e.xx, s.xx) # K
    U.nc <- LJ.c( dist(X.n, X.c,Space), e.xc, s.xc ) # K
    U.n <- sum(U.nh) + sum(U.nc) # K
    
    DU <- U.n - U.o # Difference in potential caused by move
    p <- min( c( 1, exp( DU/T) ) )
  }
  
  if( runif(1) <= p){
    X.h[id.o,] <- X.n
  }
  
  return( X.h )
}
add <- function( X.h , X.c , Space, T , V, P , lilkb){
  
  X.n <- runif(3)*Space
  
  U.h <- LJ( dist( X.n, X.h,Space ), e.xx, s.xx)
  U.c <- LJ.c( dist( X.n, X.c,Space ), e.xc, s.xc )
  DU <- sum(U.h) + sum(U.c)
  
  p <- min( c(1 , 10^(-7)*V*P/lilkb/T/(nrow(X.h)+1)*exp( -DU/T ) ) )
  if( runif(1) <= p){
    X.h <- rbind( X.h, as.vector(X.n))
  }
  
  return( X.h )
}
remove <- function( X.h , X.c , Space , T , V, P , lilkb){
  p <- 0
  if ( nrow(X.h) >= 3 ){
    
    id.r <- as.integer(runif(1)*nrow(X.h)) + 1
    
    U.h <- LJ( dist( X.h[id.r,], X.h[-id.r,],Space ), e.xx, s.xx)
    U.c <- LJ.c( dist( X.h[id.r,], X.c,Space ), e.xc, s.xc )
    DU <- -sum(U.h) - sum(U.c)
    
    p <- min ( c(1 , 10^7*nrow(X.h)*lilkb*T/V/P*exp( -DU/T ) ))
  }
  
  if( runif(1) <= p){
    X.h <- X.h[-id.r,]
  }
  
  return( X.h )
}

N<- array( data = NA, dim =  c(iterations,length(P),length(T) ) )
theta <- array( data = NA, dim =  c(iterations,length(P),length(T) ) )
Time <- array( data = NA, dim =  c(iterations,length(P),length(T) ) )

#N <- rep( NA, iterations)
#theta <- rep( NA, iterations)
#N.p <- rep( NA, length(P))
#theta.p <- rep( NA, length(P))
#Time <- rep(NA, iterations)
#Time.p <- rep(NA, length(P))
k = 1
#for( k in 1:length(T)){

for (p in 1:length(P)){
  t0.p <- Sys.time()
  X.h <- Space*matrix( c(0.5,0.5,0.6,
                         0.1,0.5,0.6,
                         0.2,0.7,0.4), nrow=3,ncol=3, byrow=TRUE) 
  
  for ( i in 1:iterations ){
    
    t0 <- Sys.time()
    if( runif(1) <= b.move){
      X.h <- move( X.h , X.c , delx , Space , T[k])
    }else if( runif(1) <= b.add){
      X.h <- add( X.h , X.c , Space , T[k] , V, P[p] , lilkb)
    }else{
      X.h <- remove( X.h , X.c , Space , T[k] , V, P[p] , lilkb)
    }
    N[i,p,k] <- nrow(X.h) 
    theta[i,p,k] <- sum((X.h[,3] >= 5 + s.min) |(X.h[,3] <= 5 - s.min))/
      length((X.h[,3] >= 5 + s.min) |(X.h[,3] <= 5 - s.min))
    t <- Sys.time()
    Time[i,p,k] <- t - t0
    
    if (i%%(iterations/10) == 0){
      #print( paste( 100*i/iterations, '% of the way there' ) )
      #print( paste( Time[i], 'per iteration'))
    }
  }
  
  #N.p[p] <- nrow(X.h)
  #theta.p[p] <- sum((X.h[,3] >= 5 + s.min) |(X.h[,3] <= 5 - s.min))/
              #length((X.h[,3] >= 5 + s.min) |(X.h[,3] <= 5 - s.min))
  #t.p <- Sys.time()
  #Time.p[p] <- t.p - t0.p
  
  print( paste( P[p] , 'Pressure Completed'))

}
#}
library(reshape2)
#N.m <- melt(N)
#theta.m <- melt(theta)
#Time.m <- melt(Time)

#library(scatterplot3d)
#library(rgl)
#plot3d( X.h[,1],X.h[,2],X.h[,3])
#plot3d( X.c[,1],X.c[,2],X.c[,3])
library(ggplot2)


write.csv( N, file = "C:/Users/lehmkudc/Desktop/N.csv")
write.csv( theta, file = "C:/Users/lehmkudc/Desktop/theta.csv")
write.csv( Time, file = "C:/Users/lehmkudc/Desktop/Time.csv")
write.csv( P, file = "C:/Users/lehmkudc/Desktop/P.csv")


x1 <- 1:iterations

#ggplot(data = N.m) + geom_point( aes(x= Var1, y=value, col = Var2))
#ggplot(data = theta.m) + geom_point( aes(x= Var1, y=value, col = Var2))


#x1 <- 1:length(P)
#plot( x1,N, type = 'l')
#plot( x1, theta, type = 'l', ylim = c(0.2,1))
#plot( N, Time)
#plot( N, theta)
# Do stuff after the simulation is over

#plot( N.p , theta.p)

#plot( P, theta.p)

#plot( P, N.p)
