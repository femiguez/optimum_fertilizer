#' Function to estimate the economic optimum N rate for any model
#' 
#' Current (2021-11-24) price of corn: 5.7 $/bu
#' one bushel is 25.4 kg. Therefore, one kilogram is worth 0.22 $/kg
#' 
#' Current price (2021-11-24) of 1 US ton of nitrogen fertilizer 1,111 dollars.
#' Therefore, 1,111 / 907 kg. (1 US ton equals 907kg) or 1.23/0.22 
#' 
#' A ratio for today would be about 5.6.
#' 
#' @name eonr
#' @param model a linear (quadratic) or nonlinear model
#' @param ratio the ratio of fertilizer to grain price
#' @param x.name name of the covariate if the 
#' default method fails (optional)
#' @param max.rate The maximum rate to consider as a possible EONR.
#' Useful when the model compute an Agronomic Optimum Nitrogen Rate.
#' @param data data.frame in case the function is used inside another
#' function
#' @param newdata optional newdata object with the levels of specific
#' factors. Useful for computing the prediction for specific treatment/level
#' combinations
#' @param level prediction level (required for models of class lme or nlme)
#' @param warning whether to print warnings when EONR is at some boundary
#' @param method method used for finding the EONR. Methods \sQuote{optimize},
#' \sQuote{L-BFGS-B}, \sQuote{DEoptim} are derivative-based. Method \sQuote{brute-force}
#' searches the range of N rates in the data and finds the one that maximizes profit
#' based on yield and nitrogen prices. The resolution is given by argument \sQuote{step}.
#' @param itermax maximum number of iterations (only used for DEoptim method)
#' @param tol tolerance (only used for \sQuote{optimize} method)
#' @param y.price yield price. If missing, the ratio will be used
#' @param n.price nitrogen price. If missing, the ratio will be used
#' @param value whether to return the EONR value or the full table.
#' @export


eonr <- function(model, ratio = 5.6, x.name, max.rate, 
                 data, step = 1, newdata = NULL, level = 0,
                 warning = TRUE, 
                 method = c("brute-force", "optimize", "Brent", "L-BFGS-B", "DEoptim"),
                 itermax = 200, tol = .Machine$double.eps^0.25,
                 y.price, n.price,
                 value = c("eonr", "table")){
 
  method <- match.arg(method)
  
  value <- match.arg(value)
  
  if(method != "brute-force"){
    if(!requireNamespace("numDeriv", quietly = TRUE)){
      warning("The numDeriv package is required for this function")
      return(NULL)
    }
    
    if(!requireNamespace("DEoptim", quietly = TRUE)){
      warning("The DEoptim packege is required when using this method")
      return(NULL)
    }    
  }

  if(missing(x.name)){
    ## Guessing the name of the covariate
    resp.var <- all.vars(formula(model))[1]
    av <- all.vars(formula(model))
    x.name <- setdiff(av, resp.var)[1]
  }
  
  if(missing(data)){
    data <- eval(model$call$data, envir = environment(formula(model)))
  }
  
  if(method == "brute-force"){
    
    if(missing(y.price) && missing(n.price)){
      y.price <- 1
      n.price <- ratio
    }
      
    nrates <- seq(from = min(data[[x.name]], na.rm = TRUE), 
                  to = max(data[[x.name]], na.rm = TRUE),
                  by = step)
    xdat <- data.frame(nrate = nrates)
    names(xdat) <- x.name
    yields <- predict(model, newdata = xdat)
    profit.table <- data.frame(nrate = nrates, yield = yields)
    profit.table$yield.revenue <- y.price * profit.table$yield
    profit.table$nrate.cost <- n.price * profit.table$nrate
    profit.table$profit <- profit.table$yield.revenue - profit.table$nrate.cost
    
    if(value == "eonr"){
      ans <- profit.table[which.max(profit.table$profit), "nrate", drop = FALSE]
      if(nrow(ans) > 1){
        ans <- mean(ans[,"nrate"])
      }else{
        ans <- ans[,"nrate"]
      }
    }else{
      ans <- profit.table
    }
  }
  
  if(method != "brute-force"){
    ## Define the gradient function
    f <- function(cf, ratio, squared = TRUE){
      f_prd <- function(x){
        xdat <- data.frame(x = x)
        names(xdat) <- x.name    
        if(!is.null(newdata)){
          xdat <- cbind(xdat, newdata)
        }
        
        ## Check if the response function actually declines
        min.yield <- predict(model, newdata = xdat[which.min(xdat[[x.name]]), , drop = FALSE], level = level)
        max.yield <- predict(model, newdata = xdat[which.max(xdat[[x.name]]), , drop = FALSE], level = level)
        if(max.yield < min.yield){
          if(warning) warning("Yield declines with input. Lowest input level chosen as EONR")
          return(min(xdat[[x.name]]))
        }
        
        if(inherits(model, c("lme", "nlme", "gnls", "gls"))){
          res <- predict_nlme(model, newdata = xdat, level = level)
        }
        if(inherits(model, c("nls", "lm"))){
          res <- predict_nls(model, newdata = xdat)
        }
        res
      } 
      gprd <- numDeriv::grad(f_prd, cf)
      
      if(squared){
        ans <- (gprd - ratio)^2
      }else{
        ans <- gprd - ratio
      }
      return(ans)
    }
    
    interval <- range(data[[x.name]])
    if(!missing(max.rate)) interval[2] <- max.rate
    
    input <- seq(interval[1], interval[2], by = step)
    
    trgt <- f(input, ratio = ratio, squared = FALSE)
    
    dat <- data.frame(input = input, trgt = trgt)
    
    print(dat)
    
    if(all(dat$trgt < 0)){
      if(warning) warning("All targets (gradient - ratio) are negative so the EONR is zero")
      ans <- 0
    }
    
    if(all(dat$trgt > 0)){
      ## If it is a linear model, the slope is constant
      if(isTRUE(all.equal(dat$trgt - mean(dat$trgt), rep(0, length(dat$trgt))))){
        ## If the slope is greater than the ratio, the max input rate should be used
        if(mean(dat$trgt) > ratio){
          ans <- max(dat$input)
        }else{
          ## if the slope is less than the ratio then the lowest input should be used
          ans <- min(dat$input)
        }
      }else{
        ## Is this even ever used?
        dat2 <- dat[all.equal(max(dat$trgt), dat$trgt), ]
        ans <- min(dat2$input)      
      }
      if(warning) warning("All targets (gradient - ratio) are positive (> 0) so the EONR is determined
                        based on the slope (if linear) or the lowest rate which maximizes yield (if not linear)")
    }
    
    if(max(dat$trgt) > 0 && min(dat$trgt) < 0){
      ## Find the minimum target
      wch.min.trgt <- which.min(abs(dat$trgt)) ## First guess is when the absolute value of the target is closer to zero
      interval <- c(dat[wch.min.trgt, 1] - 5, dat[wch.min.trgt, 1] + 5)
      
      if(method == "optimize"){
        ans <- stats::optimize(f, interval = interval, ratio = ratio, tol = tol)$minimum    
      }
      
      if(method == "Brent"){
        ans <- stats::optim(mean(interval), f, 
                            lower = interval[1],
                            upper = interval[2],
                            method = "Brent", ratio = ratio)$par    
      }
      
      if(method == "L-BFGS-B"){
        ans <- stats::optim(mean(interval), f, 
                            lower = interval[1],
                            upper = interval[2],
                            method = "L-BFGS-B", ratio = ratio)$par    
      }
      
      if(method == "DEoptim"){
        ans <- DEoptim::DEoptim(f, 
                                lower = interval[1],
                                upper = interval[2],
                                ratio = ratio, 
                                control = list(trace = FALSE, itermax = itermax))$optim$bestmem[[1]]  
      }
      
    }    
  }

  return(ans)
}

## plotting function to visualize the gradient - ratio
plot_eonr <- function(model, ratio = 5.6, x.name, max.rate, data, step = 1, newdata = NULL, level = 0){
  
  if(!requireNamespace("numDeriv", quietly = TRUE)){
    warning("The numDeriv is required for this function")
    return(NULL)
  }

  if(!requireNamespace("ggplot2", quietly = TRUE)){
    warning("The ggplot2 packege is required for this function")
    return(NULL)
  }
    
  if(missing(x.name)){
    ## Guessing the name of the covariate
    resp.var <- all.vars(formula(model))[1]
    av <- all.vars(formula(model))
    x.name <- setdiff(av, resp.var)[1]
  }  
  
  ## Define the gradient function
  f <- function(cf, ratio){
    f_prd <- function(x){
      xdat <- data.frame(x = x)
      names(xdat) <- x.name    
      if(!is.null(newdata)){
        xdat <- cbind(xdat, newdata)
      }
      if(inherits(model, c("lme", "nlme", "gnls", "gls"))){
        res <- predict_nlme(model, newdata = xdat, level = level)
      }
      if(inherits(model, c("nls", "lm"))){
        res <- predict_nls(model, newdata = xdat, level = level)
      }
      res
    } 
    gprd <- numDeriv::grad(f_prd, cf)
    ans <- gprd - ratio
  }
  
  if(missing(data)){
    data <- eval(model$call$data, envir = environment(formula(model)))
  }
  
  interval <- range(data[[x.name]])
  if(!missing(max.rate)) interval[2] <- max.rate
  
  input <- seq(interval[1], interval[2], by = step)
  
  trgt <- f(input, ratio = ratio)
  
  dat <- data.frame(input = input, ans = trgt)

  gp1 <- ggplot(data = dat, aes(x = input, y = ans)) + 
               geom_point() + 
               geom_line()
  print(gp1)
}


