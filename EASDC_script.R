# This is the EASDC script file
# 
# retrieving data
easdc_data <- read.csv(file="easdc_data.csv", header=TRUE, sep=",")
# View(easdc_data)
#
# preliminary checks
{
# checking class bias
#
table(easdc_data$regime_change)    # regime changes
table(easdc_data$rc_elections)     # regime changes through elections
table(easdc_data$rc_el_lost)       # regime changes through lost elections
table(easdc_data$rc_el_no_run)     # regime changes through elections in which the incumbent does not run
table(easdc_data$rc_coups)         # regime changes through coups
table(easdc_data$rc_uprisings)     # regime changes through uprisings
table(easdc_data$default_lag)      # sovereign debt crises (full, lagged)
table(easdc_data$default)          # sovereign debt crises (first 5 years only, lagged)
table(easdc_data$el_any_no_lag)    # direct executive or legislative elections (electoral years)
table(easdc_data$nea_no_lag)       # non-electoral autocratic country-years
table(easdc_data$ea_no_lag)        # all electoral autocratic country-years
table(easdc_data$ncea_no_lag)      # non-competitive electoral autocratic country-years
table(easdc_data$cea_no_lag)       # competitive electoral autocratic country-years
}
#
# descriptive statistics
{
library(psych)
des_easdc_data <- easdc_data[- c(1, 2, 3, 4, 20, 21, 26, 27, 29, 30, 33, 34, 36, 37, 40, 41, 43, 45, 47, 48, 49, 51, 52, 53, 54, 55, 57, 59, 61)]
names(des_easdc_data) <- c("duration", "spell", "regime_change", "prevrc", "democratization", "aa_transition", "rc_elections", "rc_el_lost", "rc_el_no_run", "rc_uprisings", "rc_coups", "party", "military", "personal", "monarchy", "default", "el_exec", "el_leg", "el_any", "ea", "ea_exp", "nea", "cea", "cea_exp", "ncea", "gdppc", "ln_gdppc", "gdppcgr", "oilgas", "vdem_index", "polity2", "polity2_avg", "europe", "latam", "mideast", "africa", "eastasia", "southasia", "cgv_civilian", "gwf_civilian", "civil_war", "1940s", "1950s", "1960s", "1970s", "1980s", "1990s", "2000s", "cold_war", "low_income", "lower_middle_income", "upper_middle_income", "high_income")
des = describe(des_easdc_data, fast=TRUE, omit=TRUE)
des <- des[-c(2), -c(1, 7, 8)]
print(des, digits=3)
#
# creating table 5
library(stargazer)
stargazer(des, title="Table 5. Descriptive statistics", summary=FALSE, type="html", out="easdc_table5.doc")
}
#
#
# PART 1: Elections as electoral institutions and autocratic regime stability during sovereign debt crises
#
{
library("stats")
library("pglm")
library("survival")
}
#
# baseline models: debt crises, electoral regime types, and autocratic regime stability
{
function1 <- function() {
  #
  # non-electoral autocracies
  pooled_nea <- glm (regime_change ~ default*nea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_nea <- pglm (regime_change ~ default*nea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_nea <- clogit (regime_change ~ default*nea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # non-competitive electoral autocracies
  pooled_ncea <- glm (regime_change ~ default*ncea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_ncea <- pglm (regime_change ~ default*ncea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_ncea <- clogit (regime_change ~ default*ncea + ln_gdppc + gdppcgr  + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # competitive electoral autocracies
  pooled_cea <- glm (regime_change ~ default*cea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_cea <- pglm (regime_change ~ default*cea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_cea <- clogit (regime_change ~ default*cea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  myres <- list(pooled_nea, re_nea, fe_nea, pooled_ncea, re_ncea, fe_ncea, pooled_cea, re_cea, fe_cea)
  return(myres)
}
easdc_data2 <- easdc_data
myres1 <- function1()
# creating table 1
{
library(texreg)
 extract.pglm <- function (model, include.nobs = TRUE, include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$estimate)
  coefficients <- s$estimate[, 1]
  standard.errors <- s$estimate[, 2]
  significance <- s$estimate[, 4]
  loglik.value <- s$loglik
  n <- nrow(model$model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, loglik.value)
    gof.names <- c(gof.names, "Log-Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                     se = standard.errors, pvalues = significance, gof.names = gof.names, 
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)
}
setMethod("extract", signature = className("maxLik", "maxLik"), definition = extract.pglm)
#
}
htmlreg(list(myres1[[1]], myres1[[2]], myres1[[3]], myres1[[4]], myres1[[5]], myres1[[6]], myres1[[7]], myres1[[8]], myres1[[9]]), file="easdc_table1.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 1. Debt crises, electoral regime types, and regime change in autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
}
#
# checking for the effects of different electoral regime types on regime change
{
function2 <- function() {
  #
  # non-electoral autocracies
  pooled_nea <- glm (regime_change ~ nea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_nea <- pglm (regime_change ~ nea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_nea <- clogit (regime_change ~ nea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # non-competitive electoral autocracies
  pooled_ncea <- glm (regime_change ~ ncea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_ncea <- pglm (regime_change ~ ncea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_ncea <- clogit (regime_change ~ ncea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # competitive electoral autocracies
  pooled_cea <- glm (regime_change ~ cea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_cea <- pglm (regime_change ~ cea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_cea <- clogit (regime_change ~ cea + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  myres <- list(pooled_nea, re_nea, fe_nea, pooled_ncea, re_ncea, fe_ncea, pooled_cea, re_cea, fe_cea)
  return(myres)
}
myres2 <- function2()
# creating table 6
htmlreg(list(myres2[[1]], myres2[[2]], myres2[[3]], myres2[[4]], myres2[[5]], myres2[[6]], myres2[[7]], myres2[[8]], myres2[[9]]), file="easdc_table6.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 6. Electoral regime types and regime change in autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
}
#
# checking for the effects of debt crises on autocratic regime stability in different electoral regime types
{
function3 <- function() {
  #
  # all autocracies
  pooled_all <- glm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_all <- pglm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_all <- clogit (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # non-electoral autocracies
  pooled_nea <- glm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2 [ which(easdc_data2$nea!='0'), ])
  re_nea <- pglm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2 [ which(easdc_data2$nea!='0'), ])
  fe_nea <- clogit (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2 [ which(easdc_data2$nea!='0'), ], method=c("efron"))
  #
  # non-competitive electoral autocracies
  pooled_ncea <- glm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2 [ which(easdc_data2$ncea!='0'), ])
  re_ncea <- pglm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2 [ which(easdc_data2$ncea!='0'), ])
  fe_ncea <- clogit (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2 [ which(easdc_data2$ncea!='0'), ], method=c("efron"))
  #
  # competitive electoral autocracies
  pooled_cea <- glm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2 [ which(easdc_data2$cea!='0'), ])
  re_cea <- pglm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2 [ which(easdc_data2$cea!='0'), ])
  fe_cea <- clogit (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2 [ which(easdc_data2$cea!='0'), ], method=c("efron"))
  #
  myres <- list(pooled_all, re_all, fe_all, pooled_nea, re_nea, fe_nea, pooled_ncea, re_ncea, fe_ncea, pooled_cea, re_cea, fe_cea)
  return(myres)
}
myres3 <- function3()
# creating table 7
htmlreg(list(myres3[[1]], myres3[[2]], myres3[[3]], myres3[[4]], myres3[[5]], myres3[[6]], myres3[[7]], myres3[[8]], myres3[[9]], myres3[[10]], myres3[[11]], myres3[[12]]), file="easdc_table7.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table 7. Debt crises and regime change in different electoral types of autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# checking for the effects of experience with elections
{
  function4 <- function() {
    #
    # all electoral autocracies
    pooled_ea <- glm (regime_change ~ default*ea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2[ which(easdc_data2$ea!='0'), ])
    re_ea <- pglm (regime_change ~ default*ea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2[ which(easdc_data2$ea!='0'), ])
    fe_ea <- clogit (regime_change ~ default*ea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2[ which(easdc_data2$ea!='0'), ], method=c("efron"))
    #
    # non-competitive electoral autocracies
    pooled_ncea <- glm (regime_change ~ default*ea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2[ which(easdc_data2$ncea!='0'), ])
    re_ncea <- pglm (regime_change ~ default*ea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2[ which(easdc_data2$ncea!='0'), ])
    fe_ncea <- clogit (regime_change ~ default*ea_exp + ln_gdppc + gdppcgr  + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2[ which(easdc_data2$ncea!='0'), ], method=c("efron"))
    #
    # competitive electoral autocracies
    pooled_cea <- glm (regime_change ~ default*cea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2[ which(easdc_data2$cea!='0'), ])
    re_cea <- pglm (regime_change ~ default*cea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2[ which(easdc_data2$cea!='0'), ])
    fe_cea <- clogit (regime_change ~ default*cea_exp + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2[ which(easdc_data2$cea!='0'), ], method=c("efron"))
    #
    names(pooled_cea$coefficients) = names(pooled_ncea$coefficients)
    names(re_cea$estimate) = names(re_ncea$estimate)
    names(fe_cea$coefficients) = names(fe_ncea$coefficients)
    myres <- list(pooled_ea, re_ea, fe_ea, pooled_ncea, re_ncea, fe_ncea, pooled_cea, re_cea, fe_cea)
    return(myres)
  }
  easdc_data2 <- easdc_data
  myres4 <- function4()
  # creating table 2  
  htmlreg(list(myres4[[1]], myres4[[2]], myres4[[3]], myres4[[4]], myres4[[5]], myres4[[6]], myres4[[7]], myres4[[8]], myres4[[9]]), file="easdc_table2.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 2. Debt crises and experience with elections in electoral autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13))
}
#
#
# ROBUSTNESS TESTS FOR PART 1
#
# Cold War period only
{
easdc_data2 <- easdc_data [ which(easdc_data$cold_war=='1'), ]
myres1_cw <- function1()
# myres2_cw <- function2()
# myres3_cw <- function3()
# creating tables A1a, A1b, A1c
htmlreg(list(myres1_cw[[1]], myres1_cw[[2]], myres1_cw[[3]], myres1_cw[[4]], myres1_cw[[5]], myres1_cw[[6]], myres1_cw[[7]], myres1_cw[[8]], myres1_cw[[9]]), file="easdc_tableA1.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A1. Debt crises, electoral regime types, and regime change in autocracies: Cold War period</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_cw[[1]], myres2_cw[[2]], myres2_cw[[3]], myres2_cw[[4]], myres2_cw[[5]], myres2_cw[[6]], myres2_cw[[7]], myres2_cw[[8]], myres2_cw[[9]]), file="easdc_tableA1b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A1b. Electoral regime types and regime change in autocracies: Cold War period</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_cw[[1]], myres3_cw[[2]], myres3_cw[[3]], myres3_cw[[4]], myres3_cw[[5]], myres3_cw[[6]], myres3_cw[[7]], myres3_cw[[8]], myres3_cw[[9]], myres3_cw[[10]], myres3_cw[[11]], myres3_cw[[12]]), file="easdc_tableA1c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A1c. Debt crises and regime change in different electoral types of autocracies: Cold War period</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# post-Cold War period only
{
easdc_data2 <- easdc_data [ which(easdc_data$cold_war=='0'), ]
myres1_pcw <- function1()
# myres2_pcw <- function2()
# myres3_pcw <- function3()
# creating tables A2a, A2b, A2c
htmlreg(list(myres1_pcw[[1]], myres1_pcw[[2]], myres1_pcw[[3]], myres1_pcw[[4]], myres1_pcw[[5]], myres1_pcw[[6]], myres1_pcw[[7]], myres1_pcw[[8]], myres1_pcw[[9]]), file="easdc_tableA2.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A2. Debt crises, electoral regime types, and regime change in autocracies: post-Cold War period</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_pcw[[1]], myres2_pcw[[2]], myres2_pcw[[3]], myres2_pcw[[4]], myres2_pcw[[5]], myres2_pcw[[6]], myres2_pcw[[7]], myres2_pcw[[8]], myres2_pcw[[9]]), file="easdc_tableA2b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A2b. Electoral regime types and regime change in autocracies: post-Cold War period</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_pcw[[1]], myres3_pcw[[2]], myres3_pcw[[3]], myres3_pcw[[4]], myres3_pcw[[5]], myres3_pcw[[6]], myres3_pcw[[7]], myres3_pcw[[8]], myres3_pcw[[9]], myres3_pcw[[10]], myres3_pcw[[11]], myres3_pcw[[12]]), file="easdc_tableA2c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A2c. Debt crises and regime change in different electoral types of autocracies: post-Cold War period</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# excluding the 1980s
{
easdc_data2 <- easdc_data [ which(easdc_data$X1980s=='0'), ]
myres1_80 <- function1()
# myres2_80 <- function2()
# myres3_80 <- function3()
# creating tables A3a, A3b, A3c
htmlreg(list(myres1_80[[1]], myres1_80[[2]], myres1_80[[3]], myres1_80[[4]], myres1_80[[5]], myres1_80[[6]], myres1_80[[7]], myres1_80[[8]], myres1_80[[9]]), file="easdc_tableA3.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A3. Debt crises, electoral regime types, and regime change in autocracies: 1980s excluded</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_80[[1]], myres2_80[[2]], myres2_80[[3]], myres2_80[[4]], myres2_80[[5]], myres2_80[[6]], myres2_80[[7]], myres2_80[[8]], myres2_80[[9]]), file="easdc_tableA3b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A3b. Electoral regime types and regime change in autocracies: 1980s excluded</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_80[[1]], myres3_80[[2]], myres3_80[[3]], myres3_80[[4]], myres3_80[[5]], myres3_80[[6]], myres3_80[[7]], myres3_80[[8]], myres3_80[[9]], myres3_80[[10]], myres3_80[[11]], myres3_80[[12]]), file="easdc_tableA3c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A3c. Debt crises and regime change in different electoral types of autocracies: 1980s excluded</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# excluding Latin America
{
easdc_data2 <- easdc_data [ which(easdc_data$latam=='0'), ]
myres1_la <- function1()
# myres2_la <- function2()
# myres3_la <- function3()
# creating tables A4a, A4b, A4c
htmlreg(list(myres1_la[[1]], myres1_la[[2]], myres1_la[[3]], myres1_la[[4]], myres1_la[[5]], myres1_la[[6]], myres1_la[[7]], myres1_la[[8]], myres1_la[[9]]), file="easdc_tableA4.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A4. Debt crises, electoral regime types, and regime change in autocracies: Latin America excluded</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_la[[1]], myres2_la[[2]], myres2_la[[3]], myres2_la[[4]], myres2_la[[5]], myres2_la[[6]], myres2_la[[7]], myres2_la[[8]], myres2_la[[9]]), file="easdc_tableA4b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A4b. Electoral regime types and regime change in autocracies: Latin America excluded</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_la[[1]], myres3_la[[2]], myres3_la[[3]], myres3_la[[4]], myres3_la[[5]], myres3_la[[6]], myres3_la[[7]], myres3_la[[8]], myres3_la[[9]], myres3_la[[10]], myres3_la[[11]], myres3_la[[12]]), file="easdc_tableA4c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A4c. Debt crises and regime change in different electoral types of autocracies: Latin America excluded</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# civilian autocracies only
{
easdc_data2 <- easdc_data [ which(easdc_data$gwf_civilian=='1'), ]
myres1_civ <- function1()
# myres2_civ <- function2()
# myres3_civ <- function3()
# creating tables A5a, A5b, A5c
htmlreg(list(myres1_civ[[1]], myres1_civ[[2]], myres1_civ[[3]], myres1_civ[[4]], myres1_civ[[5]], myres1_civ[[6]], myres1_civ[[7]], myres1_civ[[8]], myres1_civ[[9]]), file="easdc_tableA5.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A5. Debt crises, electoral regime types, and regime change in autocracies: civilian autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_civ[[1]], myres2_civ[[2]], myres2_civ[[3]], myres2_civ[[4]], myres2_civ[[5]], myres2_civ[[6]], myres2_civ[[7]], myres2_civ[[8]], myres2_civ[[9]]), file="easdc_tableA5b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A5b. Electoral regime types and regime change in autocracies: civilian autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_civ[[1]], myres3_civ[[2]], myres3_civ[[3]], myres3_civ[[4]], myres3_civ[[5]], myres3_civ[[6]], myres3_civ[[7]], myres3_civ[[8]], myres3_civ[[9]], myres3_civ[[10]], myres3_civ[[11]], myres3_civ[[12]]), file="easdc_tableA5c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A5c. Debt crises and regime change in different electoral types of autocracies: civilian autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# military autocracies only
{
easdc_data2 <- easdc_data [ which(easdc_data$cgv_civilian=='0'), ]
myres1_mil <- function1()
#  myres2_mil <- function2()
# myres3_mil <- function3()
# creating tables A6a, A6b, A6c
htmlreg(list(myres1_mil[[1]], myres1_mil[[2]], myres1_mil[[3]], myres1_mil[[4]], myres1_mil[[5]], myres1_mil[[6]], myres1_mil[[7]], myres1_mil[[8]], myres1_mil[[9]]), file="easdc_tableA6.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A6. Debt crises, electoral regime types, and regime change in autocracies: military autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_mil[[1]], myres2_mil[[2]], myres2_mil[[3]], myres2_mil[[4]], myres2_mil[[5]], myres2_mil[[6]], myres2_mil[[7]], myres2_mil[[8]], myres2_mil[[9]]), file="easdc_tableA6b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A6b. Electoral regime types and regime change in autocracies: military autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_mil[[1]], myres3_mil[[2]], myres3_mil[[3]], myres3_mil[[4]], myres3_mil[[5]], myres3_mil[[6]], myres3_mil[[7]], myres3_mil[[8]], myres3_mil[[9]], myres3_mil[[10]], myres3_mil[[11]], myres3_mil[[12]]), file="easdc_tableA6c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A6c. Debt crises and regime change in different electoral types of autocracies: military autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# low and lower-middle income countries only
{
easdc_data2 <- easdc_data [ which(easdc_data$high_income!='1' & easdc_data$upper_middle_income!='1'), ]
myres1_low <- function1()
# myres2_low <- function2()
# myres3_low <- function3()
# creating tables A7a, A7b, A7c
htmlreg(list(myres1_low[[1]], myres1_low[[2]], myres1_low[[3]], myres1_low[[4]], myres1_low[[5]], myres1_low[[6]], myres1_low[[7]], myres1_low[[8]], myres1_low[[9]]), file="easdc_tableA7.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A7. Debt crises, electoral regime types, and regime change in autocracies: low and lower-middle income autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_low[[1]], myres2_low[[2]], myres2_low[[3]], myres2_low[[4]], myres2_low[[5]], myres2_low[[6]], myres2_low[[7]], myres2_low[[8]], myres2_low[[9]]), file="easdc_tableA7b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A7b. Electoral regime types and regime change in autocracies: low and lower-middle income autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_low[[1]], myres3_low[[2]], myres3_low[[3]], myres3_low[[4]], myres3_low[[5]], myres3_low[[6]], myres3_low[[7]], myres3_low[[8]], myres3_low[[9]], myres3_low[[10]], myres3_low[[11]], myres3_low[[12]]), file="easdc_tableA7c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A7c. Debt crises and regime change in different electoral types of autocracies: low and lower-middle income autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# high and upper-middle income countries only
{
easdc_data2 <- easdc_data [ which(easdc_data$low_income!='1' & easdc_data$lower_middle_income!='1'), ]
myres1_up <- function1()
# myres2_up <- function2()
# myres3_up <- function3()
# creating tables A8a, A8b, A8c
htmlreg(list(myres1_up[[1]], myres1_up[[2]], myres1_up[[3]], myres1_up[[4]], myres1_up[[5]], myres1_up[[6]], myres1_up[[7]], myres1_up[[8]], myres1_up[[9]]), file="easdc_tableA8.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A8. Debt crises, electoral regime types, and regime change in autocracies: high and upper-middle income autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13, 15, 16, 17, 18))
# htmlreg(list(myres2_up[[1]], myres2_up[[2]], myres2_up[[3]], myres2_up[[4]], myres2_up[[5]], myres2_up[[6]], myres2_up[[7]], myres2_up[[8]], myres2_up[[9]]), file="easdc_tableA8b.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table A8b. Electoral regime types and regime change in autocracies: high and upper-middle income autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 13, 14))
# htmlreg(list(myres3_up[[1]], myres3_up[[2]], myres3_up[[3]], myres3_up[[4]], myres3_up[[5]], myres3_up[[6]], myres3_up[[7]], myres3_up[[8]], myres3_up[[9]], myres3_up[[10]], myres3_up[[11]], myres3_up[[12]]), file="easdc_tableA8c.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)"), caption.above=TRUE, caption="<b>Table A8c. Debt crises and regime change in different electoral types of autocracies: high and upper-middle income autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# creating figure 4 (robustness tests: summary)
#
{
library("jtools")
library("ggplot2")
figure4a = plot_summs (myres1_cw[[1]], myres1_pcw[[1]], myres1_80[[1]], myres1_la[[1]], myres1_civ[[1]], myres1_mil[[1]], myres1_low[[1]], myres1_up[[1]], scale=TRUE, legend.title="Models", model.names=c("Cold War period", "post-Cold War period", "1980s excluded", "Latin America excluded", "civilian autocracies", "military autocracies", "low and lower-middle income autocracies", "high and upper-middle income autocracies"), colors=c("#d78ed2", "#49b8fc", "#fe7b01", "#18d799", "#ff2c83", "#0433ff", "#e5d203", "#999999"), ci_level= .90, coefs=c("default"="default", "nea"="nea", "default:nea"="default:nea")) +
  labs(title="non-electoral autocracies") +
  theme(plot.title=element_text(size=10, hjust=0.5))
# figure4a
figure4b = plot_summs (myres1_cw[[4]], myres1_pcw[[4]], myres1_80[[4]], myres1_la[[4]], myres1_civ[[4]], myres1_mil[[4]], myres1_low[[4]], myres1_up[[4]], scale=TRUE, legend.title="Models", model.names=c("Cold War period", "post-Cold War period", "1980s excluded", "Latin America excluded", "civilian autocracies", "military autocracies", "low and lower-middle income autocracies", "high and upper-middle income autocracies"), colors=c("#d78ed2", "#49b8fc", "#fe7b01", "#18d799", "#ff2c83", "#0433ff", "#e5d203", "#999999"), ci_level= .90, coefs=c("default"="default", "ncea"="ncea", "default:ncea"="default:ncea")) +
  labs(title="non-competitive electoral autocracies") +
  theme(plot.title=element_text(size=10, hjust=0.5))
# figure4b
figure4c = plot_summs (myres1_cw[[7]], myres1_pcw[[7]], myres1_80[[7]], myres1_la[[7]], myres1_civ[[7]], myres1_mil[[7]], myres1_low[[7]], myres1_up[[7]], scale=TRUE, legend.title="Models", model.names=c("Cold War period", "post-Cold War period", "1980s excluded", "Latin America excluded", "civilian autocracies", "military autocracies", "low and lower-middle income autocracies", "high and upper-middle income autocracies"), colors=c("#d78ed2", "#49b8fc", "#fe7b01", "#18d799", "#ff2c83", "#0433ff", "#e5d203", "#999999"), ci_level= .90, coefs=c("default"="default", "cea"="cea", "default:cea"="default:cea")) +
  labs(title="competitive electoral autocracies") +
  theme(plot.title=element_text(size=10, hjust=0.5))
# figure4c
#
library("ggpubr")
figure4 = ggarrange(figure4a, figure4b, figure4c, ncol=3, nrow=1, common.legend=TRUE, legend="bottom") +
  labs(title="Figure 4: Robustness tests for non-electoral, non-competitive electoral, and competitive electoral autocracies") +
  theme(plot.title=element_text(size=12, face="bold"))
figure4
}
# 
#
# PART 2: Elections as electoral events and autocratic regime survival during sovereign debt crises
#
# checking for the type of regime change during sovereign debt crises in competitive electoral autocracies
#
{
  function5 <- function() {
  #
  # regime change via elections
  pooled_el <- glm (rc_elections ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_el <- pglm (rc_elections ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_el <- clogit (rc_elections ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # regime change via military coups
  pooled_c <- glm (rc_coups ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_c <- pglm (rc_coups ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_c <- clogit (rc_coups ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # regime change via popular uprisings
  pooled_up <- glm (rc_uprisings ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_up <- pglm (rc_uprisings ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_up <- clogit (rc_uprisings ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + military + party + personal + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  myres <- list(pooled_el, re_el, fe_el, pooled_c, re_c, fe_c, pooled_up, re_up, fe_up)    
  return(myres)
}
  easdc_data2 <- easdc_data [ which(easdc_data$cea!='0'), ]
  myres5_cea <- function5()
  # creating table 3 
  htmlreg(list(myres5_cea[[1]], myres5_cea[[2]], myres5_cea[[3]], myres5_cea[[4]], myres5_cea[[5]], myres5_cea[[6]], myres5_cea[[7]], myres5_cea[[8]], myres5_cea[[9]]), file="easdc_table3.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 3. Debt crises and types of regime change in competitive electoral autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# checking for the effects of elections on election-triggered regime change during sovereign debt crises in competitive electoral autocracies
{
function6 <- function() {
  #
  # debt crises and election-triggered regime change
  pooled_rc1 <- glm (regime_change ~ default*any_el + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
  re_rc1 <- pglm (regime_change ~ default*any_el + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
  fe_rc1 <- clogit (regime_change ~ default*any_el + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
  #
  # the efect of elections as events on regime change during sovereign debt crises
  pooled_rc2 <- glm (regime_change ~ any_el + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2[ which(easdc_data2$default=='1'), ])
  re_rc2 <- pglm (regime_change ~ any_el + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2[ which(easdc_data2$default=='1'), ])
  fe_rc2 <- clogit (regime_change ~ any_el + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2[ which(easdc_data2$default=='1'), ], method=c("efron"))
  #
  # the effect of sovereign debt crises on regime change in electoral years
  pooled_rc3 <- glm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), data=easdc_data2[ which(easdc_data2$el_any_no_lag=='1'), ])
  re_rc3 <- pglm (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2[ which(easdc_data2$el_any_no_lag=='1'), ])
  fe_rc3 <- clogit (regime_change ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + personal + military + duration + prevrc + strata(ccode), data=easdc_data2[ which(easdc_data2$el_any_no_lag=='1'), ], method=c("efron"))
  #
  myres <- list(pooled_rc1, re_rc1, fe_rc1, pooled_rc2, re_rc2, fe_rc2, pooled_rc3, re_rc3, fe_rc3)
  return(myres)
}
easdc_data2 <- easdc_data [which(easdc_data$cea!='0'), ]
myres6_cea <- function6()
# creating table 4
htmlreg(list(myres6_cea[[1]], myres6_cea[[2]], myres6_cea[[3]], myres6_cea[[4]], myres6_cea[[5]], myres6_cea[[6]], myres6_cea[[7]], myres6_cea[[8]], myres6_cea[[9]]), file="easdc_table4.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 4. Debt crises and election-triggered regime change in competitive electoral autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13))
}
#
#
# ADDITIONAL TESTS FOR PART 2
#
# checking for the type of regime change during sovereign debt crises in competitive electoral autocracies, controls for military regime type dropped
{
  function7 <- function() {
    #
    # regime change through elections
    pooled_el <- glm (rc_elections ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
    re_el <- pglm (rc_elections ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
    fe_el <- clogit (rc_elections ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
    #
    # regime change through coups
    pooled_c <- glm (rc_coups ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
    re_c <- pglm (rc_coups ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
    fe_c <- clogit (rc_coups ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
    #
    # regime change through uprisings
    pooled_up <- glm (rc_uprisings ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc, family=binomial(link="probit"), data=easdc_data2)
    re_up <- pglm (rc_uprisings ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc, family=binomial(link="probit"), model="random", index=c("ccode", "year"), data=easdc_data2)
    fe_up <- clogit (rc_uprisings ~ default + ln_gdppc + gdppcgr + oilgas + polity2_avg + party + duration + prevrc + strata(ccode), data=easdc_data2, method=c("efron"))
    #
    myres <- list(pooled_el, re_el, fe_el, pooled_c, re_c, fe_c, pooled_up, re_up, fe_up)    
    return(myres)
  }
  easdc_data2 <- easdc_data [ which(easdc_data$cea!='0'), ]
  myres7_cea <- function7()
  # creating table 8 
  htmlreg(list(myres7_cea[[1]], myres7_cea[[2]], myres7_cea[[3]], myres7_cea[[4]], myres7_cea[[5]], myres7_cea[[6]], myres7_cea[[7]], myres7_cea[[8]], myres7_cea[[9]]), file="easdc_table8.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 8. Debt crises and types of regime change in competitive electoral autocracies, controls for military regime type dropped</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 2))
}
#
# checking for the type of regime change during sovereign debt crises in all electoral autocracies
{
  easdc_data2 <- easdc_data [ which(easdc_data$ea!='0'), ]
  myres5 <- function5()
  # creating table 9  
  htmlreg(list(myres5[[1]], myres5[[2]], myres5[[3]], myres5[[4]], myres5[[5]], myres5[[6]], myres5[[7]], myres5[[8]], myres5[[9]]), file="easdc_table9.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 9. Debt crises and types of regime change in electoral autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2))
}
#
# checking for the effects of elections on election-triggered regime change during sovereign debt crises in all electoral autocracies
{
  easdc_data2 <- easdc_data [ which(easdc_data$ea!='0'), ]
  myres6 <- function6()
  # creating table 10 
  htmlreg(list(myres6[[1]], myres6[[2]], myres6[[3]], myres6[[4]], myres6[[5]], myres6[[6]], myres6[[7]], myres6[[8]], myres6[[9]]), file="easdc_table10.doc", doctype=TRUE, digits=3, custom.model.names=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"), caption.above=TRUE, caption="<b>Table 10. Debt crises and election-triggered regime change in electoral autocracies</b>", stars=c(0.01, 0.05, 0.1), custom.note="Notes: Clustered standard errors in brackets. <br> Significance levels: &#42;p&lt;0.1, &#42;&#42;p&lt;0.05, &#42;&#42;&#42;p&lt;0.01", inline.css=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE, single.row=TRUE, reorder.coef=c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 2, 3, 13))
}
#
#