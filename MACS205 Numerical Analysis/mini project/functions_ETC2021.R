
## Interpolation polynomiale
dividif=function(x,y){
  ## Computes the divided differences (coefficients on the Newton basis) for
  ##  Lagrange interpolation.
  ##
  ## @title dividif: Newton's Divided differences
  ## @param x a vector containing the interpolation nodes 
  ## @param y a vector of same size as x: values of the interpolated function at
  ##          the nodes
  ## @return a vector of same size as x: the divided differences
  ##         \eqn{f_[x_0, ... x_k]} of order 'length(x) -1'.
  
  n = length(x) -1 ## n: degree of Lagrange polynomial.
  Tmat = matrix(ncol = n+1, nrow = n+1)
  Tmat[,1]  = y ## initialization of the vector of divided differences: 
  if(n ==0) {return(diag(Tmat))} 
  for (j in 2:(n+1) ) {
    Tmat[j : (n+1), j ] = (Tmat[j : (n+1), (j-1)] - Tmat[(j-1) : n, (j-1)]) / (x[j:(n+1)] - x[1 : (n+2-j)]) 
  }
  return(diag(Tmat))
}

interpolDividif=function(x,y,z){
  ## Efficient Lagrange interpolation using Horner's method with  
  ## Newton basis for evaluation
  ## @param x : vector containing the interpolation nodes 
  ## @param y : vector of same size as x: values of the interpolated
  ##            function at the nodes
  ## @param z : vector of points where the  interpolating polynomial
  ##            needs to be evaluated. 
  ## @return  : vector of same size as z: the value of the
  ##            interpolating polynomial at points z.
  
  a = dividif(x,y)
  
  return(hornerNewton(a,x,z))
}

affineTransfo = function(a,b,u){
  return(0.5*(a+b)+0.5*(b-a)*u)
}

interpolLagrange =function(n, a, b, neval, noeuds, FUN, Plot){
  ## Generic Lagrange interpolation, with equidistant or Chebyshev nodes. 
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param a : left end-point of the interval
  ## @param b : right end-point of the interval
  ## @param neval :number of evaluation points (a regular grid will be
  ## used on [a,b]
  ## @param nodes :string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## using Chebyshev nodes.
  ## @param FUN: the function to be interpolated 
  ## @param Plot : logical. Setting 'Plot' to TRUE produces a plot
  ## showing the graph of
  ## the true functions and its interpolation.  
  ## @return : vector of size 'neval': the values of the Lagrange
  ## polynomial on an equi-distant grid.
  
  if (noeuds == "equi"){
    x =c()
    for(i in 0:n){
      x[length(x)+1] <- a+i*(b-a)/n
    }
  }
  else if (noeuds == "cheby"){
    x =  c()
    for(k in 0:n){
      x[length(x)+1] <- cos((k+0.5)*pi/(n+1))
    }
    x=affineTransfo(a,b,x)
  }
  else{stop("the nodes must be either 'equi' or 'cheby'") }
  z=seq(a, b, length.out=neval)
  f=interpolDividif(x , y = FUN(x) , z = z)
  
  
  if( Plot ){
    if (noeuds == "equi"){ methodName = " equidistant "}
    else {   methodName = " Chebyshev "}
    
    plot(z, sapply(z,FUN), type="l", ylim=range(c(y,f)) )
    title(main = paste("Lagrange interpolation with ",
                       toString(n), methodName,
                       " nodes", sep=""))
    lines(z,f, col = 'blue') 
    
    legend('topright', legend=c('function','interpolation'),
           col = c('black','blue'), lwd=1)
    
  }
  return(f)              
}


piecewiseInterpol=function(n,nInt,a,b,neval, nodes, FUN, Plot){
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param nInt :  the number of sub-intervals
  ## @param a, b : endpoints of the interval
  ## @param neval : the number of points on the interpolating grid (on
  ## each subinterval)
  ## @param nodes : string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## chebyshev nodes.
  ## @param FUN the function to be interpolated
  ## @param Plot : logical. Should the result be plotted ?
  ## @return : a matrix with 2 rows and neval * nInt -neval + 1:
  ## values of the interpolated funtion on a regular grid (first row)
  ## and the corresponding abscissas (second row).
  
  intEndPoints = seq(a,b,length.out = nInt+1)
  f = c()
  z = c()
  for (m in 1:nInt){
    A = intEndPoints[m]; B = intEndPoints[m+1] 
    
    fm = interpolLagrange(n, A, B, neval, nodes, FUN, Plot) 
    zm = seq(A, B, length.out=neval)
    
    if( m >= 2 && nodes == "equi"){
      ## remove first element of zm, fm to avoid
      ## duplicate values of the  interpolating vector
      
      zm = zm[zm!=zm[1]]
      fm = fm[fm!=fm[1]]
    }
    z = c(z,zm)
    f = c(f,fm)
  }
  
  
  if (Plot == 1){
    if (nodes == "equi") {methodName = " equidistant "}
    else  {methodName = " Chebyshev "}
    methodName = "equidistant"
    
    plot(z, sapply(z,FUN),type="l")
    title(main = paste("Piecewise  Lagrange  interpolation with ",
                       toString(n+1), methodName, " nodes  on ",
                       toString(nInt), " Intervals", sep=""))
    lines(z,f, col='red', lwd=2)
    legend('topright', legend = c('function','interpolation'),
           lwd=c(1,2), col=c('black','red'))
  }
  return(rbind(f,z))
}

findBestPiecewiseEqui=function(a,b,neval,FUN){
  list_n = c()
  list_M = c()
  erreur=c()
  for(n in 1:149){
    M = 149%/%n
    f = piecewiseInterpol(n,M,a,b,(neval-1+M)%/%M,"equi",FUN,FALSE)
    erreur[length(erreur)+1]=max(abs(FUN(f[2,])-f[1,]))
    list_n = c(list_n,n)
    list_M = c(list_M,M)
  }
  Mn = rbind(list_M,list_n)
  ind = which.min(erreur)
  n = Mn[2,]
  plot(n, erreur,type="l")
  title(main = paste("Error with different n"))
  return(c(Mn[,ind], min(erreur)))
}


findBestPiecewiseCheby=function(a,b,neval,FUN){
  list_n = c()
  list_M = c()
  erreur=c()
  for(n in 1:149){
    M = 149%/%n
    f = piecewiseInterpol(n,M,a,b,(neval-1+M)%/%M,"cheby",FUN,FALSE)
    erreur[length(erreur)+1]=max(abs(FUN(f[2,])-f[1,]))
    list_n = c(list_n,n)
    list_M = c(list_M,M)
  }
  Mn = rbind(list_M,list_n)
  ind = which.min(erreur)
  n = Mn[2,]
  plot(n, erreur,type="l")
  title(main = paste("Error with different n"))
  return(c(Mn[,ind], min(erreur)))
}





## Méthodes de quadrature


## Méthode des trapèzes
trapezeInt =function(FUN,a,b,M){
  ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : number of intervals (each of size (b-a)/M)
  ##' @return: the value of the composite trapezoidal quadrature. 
  x = seq(a,b, length.out= M+1)
  y = sapply(x, FUN)
  h = (b-a)/M 
  w = rep(1,M+1)
  w[1] <- 0.5
  w[M+1] <- 0.5
  q = h*(sum(y*w))
  
  return(q)
}

## Division du pas des trapèzes par 2
refineTrapeze=function(FUN,a,b,M,q){
  ##' refinement of the subdivision step: incremental method
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : initial number of intervals (each of size (b-a)/M)
  ##'  having been used to compute q
  ##' @param  q : the value of the trapezoidal  quadrature method
  ##'  of stepsize (b-a)/M
  ##' @return : the value of the quadrature for a stepsize h' = h/2
  h = (b-a)/M
  x = seq(a+h*0.5, b-h*0.5, length.out= M) ## complete here : 
  ##  x : a vector of size M :
  ##     the additional abscissas where 'fun' must be evaluated.
  y = sapply(x, FUN)
  Q = 0.5*q+0.5*h*sum(y)##  complete here : a function of y, h and q. 
  return(Q)
}

## Méthode de Simpson
simpsonInt = function(FUN,a,b,M){
  ##' Simpson integration via trapeze rule
  ##' uses the fact that 
  ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
  qtrapeze = trapezeInt(FUN,a,b,M)
  qrefined = refineTrapeze(FUN,a,b,M,qtrapeze)
  q =  (4/3)*(qrefined-(1/4)*qtrapeze) 
  return(q)
}

## Evaluation de l'erreur a posteriori
evalErrSimpson=function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4* qth2 )
  q = simps_h2  
  E = (simps_h-q)/15 
  return(c(E,simps_h))
}


erreurRelative=function(a,b,M){
  temp = evalErrSimpson(densite,a,b,M/2)
  err1 = temp[1]
  err2 = temp[2]-pgamma(b,shape=K,scale=THETA)
  return(abs(err1-err2)/err2)
}



## Méthode de Romberg
romberg =function(FUN,n,a,b,M){## methode de Romberg avec n etapes
  ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
  ## pas initial h = (b-a)/M
  h= (b-a)/M 
  A = matrix(ncol = n+1, nrow = n+1)
  A[1,1] = trapezeInt(FUN,a,b,M);
  Mc = M
  ## initialisation des differences divisees
  for( i in 2:(n+1)){
    A[i,1] = refineTrapeze( FUN,a,b, Mc, q= A[i-1,1])
    Mc = 2*Mc 
  }
  delta = 1/4;
  lx = log(delta)*(0:n)
  x = exp(lx);
  for (j in 2:(n+1)){
    A[j : (n+1), j ] = (A[j : (n+1), (j-1)] - (x[j:(n+1)])*A[(j-1) : n, (j-1)]) / (1 - x[j:(n+1)])
  }
  return(diag(A))
}