# Goal: Sample small-shape gamma random variables via accept-reject
# Paper: arXiv:1302.1884
rgamss <- function(n, shape, scale=1, do.log=TRUE) 
{  
	a <- shape
	if(a > 0.2)
	{
		oo <- ifelse(do.log,log(rgamma(n, a, 1 / scale)),rgamma(n, a, 1 / scale))
	}
	else 
	{
		e1 <- 2.71828182845905
		L <- 1 / a - 1
		w <- a / e1 / (1 - a)
		ww <- 1 / (1 + w)
		eta <- function(z)
		{ 
			ifelse(z >= 0, exp(-z), w * L * exp(L * z)) 
		}
		h <- function(z)
		{
			exp(-z - exp(-z / a))
		} 
		rh <- function(a)
		{            
			repeat 
			{        
				U <- runif(1)
				z <- ifelse(U <= ww, -log(U / ww), log(runif(1)) / L)
				if(h(z) / eta(z) > runif(1)) 
					return(z)        
			}
		}
		Z <- numeric(n)
		for(i in 1:n) 
			Z[i] <- rh(a)
		o <- log(scale) - Z / a
		if(!do.log) 
		{      
			oo <- exp(o)
			if(any(oo == 0)) 
			{        
				oo <- o
				warning("Output given on log-scale since shape is small")        
			}      
		} 
		else 
			oo <- o    
	}
	return(oo)
}