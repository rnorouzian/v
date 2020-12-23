
Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   Welcome to 'EDP 380C.12 Survey of Multivariate Methods'.
   Programs developed by Reza Norouzian, Copyright (C) 2019-present"

message(Break, notice, Break)







need <- c('car','psych','reshape','tidyverse','lme4','nlme','MASS','CCA','matrixcalc',
          'parallel','rela','gplots','ICSNP','mvtnorm','mvnormtest','normtest',
          'normwhn.test','nortest','biotools','effects','ez','yacca')

not.have <- need[!(need %in% installed.packages()[,"Package"])]
if(length(not.have)) install.packages(not.have)


suppressWarnings(
suppressMessages({ 
  
  for(i in need){
    library(i, character.only = TRUE)
  }
}))
