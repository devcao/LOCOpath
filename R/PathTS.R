
#' This function wraps up ApproxPath.TS and ExactPath.Ts
#' @description This function wraps up ApproxPath.TS and ExactPath.Ts
#' @param exact: boolean, using exact calculation or not
#' @param ...: args passed to ExactPath.TS and ApproxPath.TS 
#' @return a vector of test statistic for each test specified
#' @export
Path.TS <- function(exact = TRUE,...){

	if(exact){
		return(ExactPath.TS(...))

	}else{
		return(ApproxPath.TS(...))

	}

}


#' This function wraps up ApproxPath.TS.Para and ExactPath.Ts.Para
#' @description This function could be called in parallel 
#' @param exact: boolean, using exact calculation or not
#' @param ...: args passed to ExactPath.TS and ApproxPath.TS 
#' @return a vector of test statistic for each test specified
#' @export
Path.TS.Para <- function(exact = TRUE,...){

	if(exact){
		return(ExactPath.TS.Para(...))

	}else{
		return(ApproxPath.TS.Para(...))

	}

}

#' This function helps to do ExactPath.TS in parallel
#' @description This function could be called in parallel 
#' @param mat: A matrix which contains X and Y
#' @param ...: args passed to ExactPath.TS 
#' @return a vector of test statistic for each test specified
#' @export
ExactPath.TS.Para <- function(mat, ...){
# Calculate PATH statistic exactly, could run this in parallel

	n = nrow(mat)
	p = ncol(mat) - 1
	X = mat[,1:p] 
	Y = mat[,p+1] 

	return(ExactPath.TS(X,Y,...))

}


#' This function helps to do ApproxPath.TS in parallel
#' @description This function could be called in parallel 
#' @param mat: A matrix which contains X and Y
#' @param ...: args passed to ApproxPath.TS 
#' @return a vector of test statistic for each test specified
#' @export
ApproxPath.TS.Para <- function(mat, ...){
# Calculate PATH statistic using interpolating, could run this in parallel

	n = nrow(mat)
	p = ncol(mat) - 1
	X = mat[,1:p] 
	Y = mat[,p+1] 

	return(ApproxPath.TS(X,Y,...))

}



