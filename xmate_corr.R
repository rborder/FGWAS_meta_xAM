
################################## INPUT ###################################
## ff: female mates
## mm: male mates

################################## OUTPUT ##################################
## output_r: cross-mate correlations
## output_se: standard errors
## output_n : sample sizes
############################################################################

ordinal_vars <- c("Ever.smoker", "Hypertension",  "Hypertension.2", 
                  "ADHD", "Allergic.rhinitis", "Asthma", "Morning.Person", 
                  "Myopia", 
                  "Migraine", "Depression")

continuous_vars <- setdiff(names(ff[-(1:6)]),ordinal_vars)
continuous_vars_male <- setdiff(continuous_vars,c("Age.at.first.birth..women.",
 "Age.at.menarche"))


## winsorize continuous vars within sex
ff[continuous_vars] <- lapply(ff[continuous_vars], DescTools::Winsorize, 
                              probs = c(.005,.995),na.rm=T)
mm[continuous_vars_male] <- lapply(mm[continuous_vars_male], 
                                   DescTools::Winsorize, probs = c(.005,.995),
                                   na.rm=T)

## adjust continuous vars for age and age**2 within sex
regress_out <- function(y,x) {
    if (all(is.na(y))) {
        return(y)
    } else {
       mm <- lm(y~x+I(x**2))
        return(y - predict(mm, newdata=data.frame(x=x)))
    }
}
ff[continuous_vars] <- lapply(ff[continuous_vars], regress_out, x=ff$age)
mm[continuous_vars_male] <- lapply(mm[continuous_vars_male], 
                                   regress_out, x=ff$age)




output_r <- output_se <- output_n <- 
    matrix(NA, nrow=length(c(ordinal_vars, continuous_vars)), 
           ncol=length(c(ordinal_vars, continuous_vars_male)),
           dimnames=list(c(ordinal_vars, continuous_vars),
                         c(ordinal_vars, continuous_vars_male)))


partial_pseudo_r_one_way <- function(data0, ind) {
    augmented <- (amod<-lrm(y ~ x + cov, data=data0[ind,]))$stats['R2']
    compact <- lrm(y ~ cov, data=data0[ind,])$stats['R2']
    return(c(r=sqrt(augmented-compact)*sign(amod$coef[2])))
} 

partial_pseudo_r_two_way <- function(data0, ind) {
    augmented1 <- (amod1 <- lrm(y1 ~ y2 + cov1 + cov2 ,
                                data=data0[ind,]))$stats['R2']
    augmented2 <- (amod2 <- lrm(y2 ~ y1 + cov1 + cov2 ,
                                data=data0[ind,]))$stats['R2']
    compact1 <- lrm(y1 ~ cov1 + cov2 , data=data0[ind,])$stats['R2']
    compact2 <- lrm(y2 ~ cov1 + cov2 , data=data0[ind,])$stats['R2']
    return(c(r=(sqrt(augmented1-compact1)*sign(amod2$coef[2]) + 
                sqrt(augmented2-compact2)*sign(amod2$coef[2]))/2))
} 


library(rms)
library(boot)

for (f_trait in c(ordinal_vars, continuous_vars)) {
    for (m_trait in c(ordinal_vars, continuous_vars_male)) {
        message(f_trait,' ',m_trait)
        try(
            if (f_trait %in% continuous_vars & m_trait %in% continuous_vars) {
            tmp <- polycor::hetcor(ff[[f_trait]],mm[[m_trait]])
            output_r[f_trait,m_trait] <- tmp$correlations[1,2]
            output_se[f_trait,m_trait] <- tmp$std.errors[1,2]
            output_n[f_trait,m_trait] <- tmp$n
        } else if (f_trait %in% ordinal_vars & m_trait %in% continuous_vars) {
            tmp <- na.omit(data.frame(y=ff[[f_trait]], x=mm[[m_trait]], 
                                      cov=ff$age))
            bb <- boot(tmp, R=100, partial_pseudo_r_one_way,parallel='multicore')
            output_r[f_trait,m_trait] <-bb$t0
            output_se[f_trait,m_trait] <- sd(bb$t)
            output_n[f_trait,m_trait] <- nrow(tmp)
        } else if (m_trait %in% ordinal_vars & f_trait %in% continuous_vars) {
            tmp <- na.omit(data.frame(y=mm[[m_trait]], x=ff[[f_trait]], 
                                      cov=mm$age))
            bb <- boot(tmp, R=100, partial_pseudo_r_one_way,parallel='multicore')
            output_r[f_trait,m_trait] <-bb$t0
            output_se[f_trait,m_trait] <- sd(bb$t)
            output_n[f_trait,m_trait] <- nrow(tmp)
        } else if (m_trait %in% ordinal_vars & f_trait %in% ordinal_vars) {
            tmp <- na.omit(data.frame(y1=mm[[m_trait]], y2=ff[[f_trait]], 
                cov1=mm$age, cov2=ff$age))
            bb <- boot(tmp, R=100, partial_pseudo_r_two_way,parallel='multicore')
            output_r[f_trait,m_trait] <-bb$t0
            output_se[f_trait,m_trait] <- sd(bb$t)
            output_n[f_trait,m_trait] <- nrow(tmp)           
        }
        )
        message('  r = ',output_r[f_trait,m_trait])
        message('  se = ',output_se[f_trait,m_trait])
        message('  n = ',output_n[f_trait,m_trait])
    }
}

