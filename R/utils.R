conversion_copula_tau = function(copula_par, copula_family){
  if(copula_family == "frank"){
    return(
      copula::tau(copula = copula::frankCopula(copula_par))
    )
  }
  else if(copula_family == "gaussian"){
    #convert to correct scale
    correlation_scale = (exp(copula_par) - 1)/(exp(copula_par) + 1)
    return(
      copula::tau(copula = copula::normalCopula(correlation_scale))
    )
  }
  else if(copula_family == "clayton"){
    return(
      copula::tau(copula = copula::claytonCopula(copula_par))
    )
  }
  else if(copula_family == "gumbel"){
    return(
      copula::tau(copula = copula::gumbelCopula(copula_par))
    )
  }
}
