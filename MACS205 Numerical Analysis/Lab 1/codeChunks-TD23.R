


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
    print(n)
    Tmat = matrix(ncol = n+1, nrow = n+1)
    Tmat[,1]  = y ## initialization of the vector of divided differences: 
    if(n ==0) {return(diag(Tmat))} 
    for (j in 2:(n+1) ) {
      Tmat[j : (n+1), j ] = (Tmat[j : (n+1), (j-1)] - Tmat[(j-1) : n, (j-1)]) / (x[j:(n+1)] - x[1 : (n+2-j)]) 
    }
    return(diag(Tmat))
}
nodes = c(0, 1, 2, 3)
om0 = function(x){1}
om1 = function(x){x- nodes[1]}
om2 = function(x){(x- nodes[1])*(x- nodes[2])}
om3 = function(x){(x- nodes[1])*(x- nodes[2])*(x- nodes[3])}
A <- c(0.1,0.2, -0.1 , -0.3)
divtest <- function(x) {
  return( A[1] * om0(x) + A[2]*om1(x) + A[3]*om2(x)+ A[4]*om3(x) )
}
divA <- dividif(x= nodes, y = sapply(nodes, divtest) )
divA
A


hornerNewton = function(a,x,z){
    ## Horner's method: Evaluates  a polynom P at points z, given
    ## nodes x and the coefficients a of P in Newton's basis
    ##
    ## @param a : vector: the  coefficients of the polynomial in
    ##           Newton's basis
    ## @param x : the interpolation nodes. 
    ## @param z : vector of points where the polynom needs to be
    ##            evaluated. 
    ## @return  : a vector of same size as z: the value of the
    ##            polynomial at points z.
    ## 
    n <- length(x) - 1 ## degree of the Lagrange poynomial 
    if( (n < 0) || (length(a) != (n+1)) )
    {
        stop('at least one interpolating point is needed,
              a and x should have same length')
    }
    f <- a[n+1]*rep(1,length(z))
    if(n >= 1){
        for( k in 1:n){
           f <- f*(z-x[n+1-k])+ a[n+1-k] ## complete the code 
        }
    }
    return(f)
}

hornerNewton(a=A, x = nodes, z = c(0,0.5,2,3))
sapply( c(0,0.5,2,3), divtest)
plot(sapply( c(0,0.5,2,3), divtest))


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

n <- 11
x <-  seq(-1,1, length.out =n) #ensemble de points équidistants 

myfun = function(x){cos(5*x/(2+x))*1/(1+(1*x)^2)}

interpolDividif(x , y = myfun(x) , z = seq(-1,1, length.out=100))
z=seq(-1,1, length.out=100)
plot(z, myfun(z), col="red")
lines(z,interpolDividif(x , y = myfun(x) , z = z), type='l')


funRunge=function(x){ 1/(1+(5*x)^2)}
xx <- seq(-1, 1, length.out =n)
plot(z, funRunge(z), type='l', col="green")
lines(z, interpolDividif(x, funRunge(x), z= z), type='l')

affineTransfo = function(a,b,u){
  return(0.5*(a+b)+0.5*(b-a)*u)
}


  interpolLagrange =function(n, a, b, neval, nodes = 'equi', FUN, Plot){
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
      
      if (nodes == "equi"){
          x =c()
          for(i in 0:(n-1)){
            x[length(x)+1] <- a+i*(b-a)/(n-1)
          }
              }
      else if (nodes == "cheby"){
        x =  c()
        for(k in 0:(n-1)){
          x[length(x)+1] <- cos((k+0.5)*pi/n)
        }
        x=affineTransfo(a,b,x)
                      }
      else{stop("the nodes must be either 'equi' or 'cheby'") }
    print(x)
      z=seq(a, b, length.out=neval)
      f=interpolDividif(x , y = FUN(x) , z = z)

      
      if( Plot ){
          if (nodes == "equi"){ methodName = " equidistant "}
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
  
interpolLagrange(10,-1,1,100,"equi",funRunge, TRUE)

piecewiseInterpol=function(n,nInt,a,b,neval, nodes = "equi", FUN, Plot){
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
                
                if( m >= 2){
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

piecewiseInterpol(4,10,-1,1,80,"equi",funRunge,TRUE)

  
