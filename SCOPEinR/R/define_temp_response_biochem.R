#' get the variables needed in temp_response_biochem
#'
#' @param getTDP is TRUE get a list with TDP
#'
#' @return get a list with all parameters for TDP
#' @export
#'
#' @examples
#' TDP <- define_temp_response_biochem(getTDP = TRUE)
define_temp_response_biochem<-function(getTDP=T){

# get a list with all parameters
TDP<-list()
TDP$delHaV = 65330;
TDP$delSV	= 485;
TDP$delHdV = 149250;

TDP$delHaJ = 43540;
TDP$delSJ	= 495;
TDP$delHdJ	= 152040;

TDP$delHaP = 53100;
TDP$delSP	= 490;
TDP$delHdP	= 150650;

TDP$delHaR	= 46390;
TDP$delSR	= 490;
TDP$delHdR	= 150650;

TDP$delHaKc	= 79430;
TDP$delHaKo	= 36380;
TDP$delHaT	= 37830;

TDP$Q10	= 2;
TDP$s1	= 0.3;
TDP$s2	= 313.15;
TDP$s3	= 0.2;
TDP$s4	= 288.15;
TDP$s5	= 1.3;
TDP$s6	= 328.15;
return(TDP)
}
