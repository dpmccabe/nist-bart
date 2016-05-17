gc()
options(java.parameters = "-Xmx10g")

sapply(unlist(strsplit("gdata dplyr plyr randomForest caret ROCR ResourceSelection bartMachine gtools caret doParallel", " ")), function(pkg) {
  if (!is.element(pkg, installed.packages()[,1])) install.packages(pkg, dep = T)
  library(pkg, character.only = T, quietly = T)
})

cl = makePSOCKcluster(8)
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

set.seed(20160513)

rm(list = ls())
load("~/Dropbox/2016-1 (Spring)/Bayesian (2530)/project2/nist.rdata")
load("~/Dropbox/2016-1 (Spring)/Bayesian (2530)/project2/nistcc.rdata")
load("~/Dropbox/2016-1 (Spring)/Bayesian (2530)/project2/bt_caret_out.rdata")
load("~/Dropbox/2016-1 (Spring)/Bayesian (2530)/project2/bt_fin.rdata")

####################################################################### utilities
######################################################################
explode = function(s) {
  return(unlist(strsplit(s, " ")))
}

####################################################################### data cleanup
######################################################################
nist[nist[, "CKUP_AGE"] %in% c(77, 99), "CKUP_AGE"] = NA

fac_cols = explode("PDAT2 ASTHMA CKUP_11_12 CKUP_LAST RISK_EVER RISK_HH RISK_NOW CPOX_HAD HPVI_ANY_REC HPVI_INTENTR IMM_ANY MEN_ANY_REC NOSCHOOLR TET_ANY_REC VISITS AGEGRP_M_I CEN_REG CHILDNM EDUC1 INCPOV1 I_HISP_K LANGUAGE MARITAL2 MOBIL_I RACE_K RENT_OWN SEX FACILITY VFC_ORDER P_U13HPV P_UTDHPV P_UTDHPV2 P_U13HPV3 P_UTDHPV11 P_UTDHPV12 P_UTDHPV13 P_UTDHPV3C TIS_INS_1 TIS_INS_2 TIS_INS_3 TIS_INS_3A TIS_INS_4_5 TIS_INS_6 TIS_INS_11")
fac_na_vals = c("DON'T KNOW", "IN ERROR", "REFUSED", "UNKNOWN", "UNKNOWN/DON'T KNOW")

for (col in fac_cols) {
  nist[nist[, col] %in% fac_na_vals, col] = NA
  nist[, col] = factor(nist[, col], ordered = F)
}

nist = droplevels(nist)
str(nist[, fac_cols])
summary(nist[, fac_cols])

####################################################################### insurance indicator variables
######################################################################
ins_cols = explode("TIS_INS_1 TIS_INS_2 TIS_INS_3 TIS_INS_3A TIS_INS_4_5 TIS_INS_6 TIS_INS_11")
n_ins_present = rowSums(nist[, ins_cols] == "YES" | nist[, ins_cols] == "NO", na.rm  = T)
table(n_ins_present)

nist$has_ins_info = n_ins_present > 0

summary(nist[, ins_cols])
nist[, ins_cols][is.na(nist[, ins_cols])] <- as.factor("NO")
summary(nist[, ins_cols])

####################################################################### vax indicator variables
######################################################################
vax_cols = explode("P_U13HPV P_UTDHPV P_UTDHPV2 P_U13HPV3 P_UTDHPV11 P_UTDHPV12 P_UTDHPV13 P_UTDHPV3C")

ivar_df = adply(vax_cols, 1, function(col) {
  summary(nist[, col])
  return(summary(nist[, col]))
})[, -1]
rownames(ivar_df) = vax_cols
ivar_df

nist$vax = nist$P_UTDHPV
summary(nist$vax)

P_UTDHPV_utd = which(nist$P_UTDHPV == "UTD")
P_UTDHPV_not_utd = which(nist$P_UTDHPV == "NOT UTD")
sum(nist[P_UTDHPV_utd, "RDDWT_D_TERR"]) / sum(c(nist[P_UTDHPV_utd, "RDDWT_D_TERR"], nist[P_UTDHPV_not_utd, "RDDWT_D_TERR"]))

P_UTDHPV3C_utd = which(nist$P_UTDHPV3C == "UTD")
P_UTDHPV3C_not_utd = which(nist$P_UTDHPV3C == "NOT UTD" & nist$P_UTDHPV == "UTD")
sum(nist[P_UTDHPV3C_utd, "RDDWT_D_TERR"]) / sum(c(nist[P_UTDHPV3C_utd, "RDDWT_D_TERR"], nist[P_UTDHPV3C_not_utd, "RDDWT_D_TERR"]))

####################################################################### complete cases
######################################################################
cols_oi = explode("SEQNUMT RDDWT_D_TERR vax AGE SEX CKUP_11_12 CKUP_AGE RISK_EVER RISK_NOW RISK_HH C1R CHILDNM INCPORAR ASTHMA CPOX_HAD HPVI_INTENTR IMM_ANY MEN_ANY_REC NOSCHOOLR TET_ANY_REC VISITS AGEGRP_M_I CEN_REG EDUC1 INCPOV1 I_HISP_K LANGUAGE MARITAL2 MOBIL_I RACE_K RENT_OWN FACILITY VFC_ORDER TIS_INS_1 TIS_INS_2 TIS_INS_3 TIS_INS_4_5 TIS_INS_6 TIS_INS_11")

cc_crit = !is.na(nist$vax) & !is.na(nist$INCPORAR) & !is.na(nist$ASTHMA) & !is.na(nist$NOSCHOOLR) & !is.na(nist$RENT_OWN) & !is.na(nist$FACILITY) & !is.na(nist$CEN_REG) & !is.na(nist$VFC_ORDER) & !is.na(nist$VISITS)
nistcc = nist[cc_crit, cols_oi]
summary(nistcc)

# ws_ind = sample(1:nrow(nistcc), replace = T, prob = nistcc$RDDWT_D_TERR)
# nistcc = nistcc[ws_ind, ]
# summary(nistcc)
# summary(c(table(nistcc$SEQNUMT)))

rbind(c(table(nistcc$vax)), round(c(table(nistcc$vax)) / nrow(nistcc) * 100, 1))

####################################################################### recode
######################################################################
new_col_names = explode("id sw vax age sex checkup_11_12 checkup_age risk_ever risk_now risk_household people_in_household children_in_household inc_pov_rat asthma had_cpox HPVI_INTENTR any_imm MEN_ANY_REC school_missed TET_ANY_REC doc_visits mother_age_group region mother_edu poverty_status hispanic_latino language mother_married moved race rent_own facility vax_order ins_emp_union ins_medicaid ins_schip ins_other_gov ins_other ins_lacked")
cbind(colnames(nistcc), new_col_names)
colnames(nistcc) = new_col_names

levels(nistcc$vax) = explode("unvaxxed vaxxed")
levels(nistcc$sex) = explode("male female")
levels(nistcc$children_in_household) = explode("1 2-3 4+")
levels(nistcc$race) = explode("white black other_mult")
levels(nistcc$region) = explode("northeast midwest south west")
levels(nistcc$rent_own) = explode("own rent other")
levels(nistcc$moved) = explode("YES NO")
levels(nistcc$poverty_status) = explode("above_gt_75k above_lt_75k below")
levels(nistcc$mother_age_group) = explode("<35 35-44 >44")
levels(nistcc$mother_edu) = explode("lt_high_school high_school some_college college")
levels(nistcc$mother_married) = explode("YES NO")
levels(nistcc$language) = explode("English Spanish Other")
levels(nistcc$facility) = explode("public hospital private military_other mixed")
levels(nistcc$vax_order) = explode("all some no")

risk_now_ind = which(!is.na(nistcc$risk_now) & nistcc$risk_now == "YES")
risk_ever_ind = which(!is.na(nistcc$risk_ever) & nistcc$risk_ever == "YES")
risk_before_ind = setdiff(risk_ever_ind, risk_now_ind)
risk_missing_ind = which(is.na(nistcc$risk_now) & is.na(nistcc$risk_ever))
nistcc$risk = "NEVER"
nistcc[risk_now_ind, "risk"] = "NOW"
nistcc[risk_before_ind, "risk"] = "BEFORE"
nistcc[risk_missing_ind, "risk"] = NA
nistcc$risk = as.factor(nistcc$risk)
summary(nistcc$risk)

doc_visits = as.character(nistcc$doc_visits)
doc_visits[doc_visits == "NONE"] = "0"
doc_visits[doc_visits == "2 - 3"] = "2-3"
doc_visits[doc_visits == "4 - 5"] = "4-5"
doc_visits[doc_visits %in% c("6 - 7", "8 - 9")] = "6-9"
doc_visits[doc_visits %in% c("10 - 12", "13 - 15", "16+")] = "10+"
table(doc_visits)
nistcc$doc_visits = factor(doc_visits)

school_missed = as.character(nistcc$school_missed)
school_missed[school_missed == "ZERO"] = "0"
school_missed[school_missed == "ONE"] = "1"
school_missed[school_missed == "TWO"] = "2"
school_missed[school_missed %in% c("THREE", "FOUR", "FIVE")] = "3-5"
school_missed[school_missed %in% c("6 TO 9")] = "6-9"
school_missed[school_missed %in% c("10 TO 19")] = "10-19"
school_missed[school_missed %in% c("20 TO 29", "30 OR MORE")] = "20+"
school_missed[school_missed == "DID NOT GO TO SCHOOL"] = "no_school"
table(school_missed)
nistcc$school_missed = factor(school_missed)

nistcc = subset(nistcc, select = -c(TET_ANY_REC, HPVI_INTENTR, MEN_ANY_REC, risk_ever, risk_now))

str(nistcc)
summary(nistcc)

####################################################################### formula
######################################################################
form = vax ~ age + sex + checkup_11_12 + checkup_age + risk + risk_household + people_in_household + children_in_household + inc_pov_rat + asthma + had_cpox + any_imm + school_missed + doc_visits + mother_age_group + region + mother_edu + poverty_status + hispanic_latino + language + mother_married + moved + race + rent_own + facility + vax_order + ins_emp_union + ins_medicaid + ins_schip + ins_other_gov + ins_other + ins_lacked

preds = strsplit(as.character(form)[3], " \\+ ", perl = T)[[1]]

####################################################################### glm
######################################################################
null_mod = glm(vax ~ 1, data = nistcc[complete.cases(nistcc),], family = "binomial")
lm_stepped = step(null_mod, form, k = 2, direction = "forward")
summary(lm_stepped)

full_mod = glm(form, data = nistcc, family = "binomial")
lm_stepped_back = step(full_mod, k = log(nrow(nistcc)), direction = "backward")
summary(lm_stepped_back)

prob = predict.glm(lm_stepped)
pred = prediction(prob, nistcc[complete.cases(nistcc), "vax"])
perf = performance(pred, "tpr", "fpr")
auc = performance(pred, measure = "auc")
auc = auc@y.values[[1]]
roc.data = data.frame(fpr = unlist(perf@x.values), tpr = unlist(perf@y.values), model = "GLM")
ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) + geom_ribbon(alpha=0.2) + geom_abline(intercept = 0, slope = 1, colour = "gray")+ geom_line(aes(y=tpr)) + ggtitle(paste0("ROC Curve w/ AUC=", auc))

glm_pred_probs = predict.glm(lm_stepped, type = "response")
plot(density(glm_pred_probs))
glm_preds = as.factor(ifelse(glm_pred_probs < 0.5, "unvaxxed", "vaxxed"))

confusionMatrix(glm_preds, nistcc$vax, positive = "vaxxed")

####################################################################### random forest
######################################################################
tune_rf = tuneRF(x = nistcc[, preds], y = nistcc$vax, ntreeTry = 101, stepFactor = 3, improve = 0.02, do.trace = T)
rf = randomForest(form, data = nistcc, ntree = 51, mtry = 2, importance = T, do.trace = T)

head(getTree(rf, 1, T), 20)

varImpPlot(rf)

####################################################################### bayes tree
######################################################################
set_bart_machine_num_cores(8)
bm = bartMachine(X = nistcc[, preds], y = nistcc$vax, num_trees = 51, use_missing_data = T, use_missing_data_dummies_as_covars = T, mem_cache_for_speed = F, verbose = T)
gc()

bm_imp = investigate_var_importance(bm)
bm_imp_avg = cbind(bm_imp$avg_var_props)
cbind(bm_imp_avg[order(-bm_imp_avg),1])

pd_plot(bm, which(bm$training_data_features == "inc_pov_rat"), prop_data = 0.1)

####################################################################### bayes tree tuning
######################################################################
summary(nistcc[, nearZeroVar(nistcc[, preds])])

bt_traincon = trainControl(method = "boot632", number = 3, verboseIter = T, allowParallel = F)
bt_caret_out = train(trControl = bt_traincon, method = "bartMachine", x = nistcc[, preds], y = nistcc$vax, tuneGrid = expand.grid(num_trees = 31, k = c(1, 2, 3), alpha = c(0.75, 0.85, 0.9, 0.95), beta = c(1, 2, 3), nu = 3), use_missing_data = T, use_missing_data_dummies_as_covars = T, mem_cache_for_speed = F, verbose = T)

bt_traincon_fin = trainControl(method = "boot632", number = 10, verboseIter = T, allowParallel = F, savePredictions = "final", classProbs = T)
bt_caret_out_fin = train(trControl = bt_traincon_fin, method = "bartMachine", x = nistcc[, preds], y = nistcc$vax, tuneGrid = expand.grid(num_trees = 101, k = 2, alpha = 0.85, beta = 1, nu = 3), use_missing_data = T, use_missing_data_dummies_as_covars = T, mem_cache_for_speed = F, verbose = T)

save(bt_caret_out_fin, file = "bt_caret_out_fin.rdata")

####################################################################### bayes tree model analysis
######################################################################
bt_fin = bt_caret_out_fin$finalModel
save(bt_fin, file = "bt_fin.rdata")

fin_preds = bt_caret_out_fin$pred
fin_preds = filter(fin_preds, Resample == "AllData")
cm = confusionMatrix(fin_preds$pred, fin_preds$obs, positive = "vaxxed")
round(cm$table / 100, 1)

bt_predictions = prediction(fin_preds$vaxxed, fin_preds$obs, label.ordering = c("vaxxed", "unvaxxed"))
bt_perf = performance(bt_predictions, "tpr", "fpr")
bt_auc = performance(bt_predictions, measure = "auc")
bt_auc = bt_auc@y.values[[1]]
bt_roc_df = data.frame(FPR = unlist(bt_perf@x.values), TPR = unlist(bt_perf@y.values))
ggplot(bt_roc_df, aes(x = FPR, ymin = 0, ymax = TPR)) + geom_ribbon(alpha = 0.2) + geom_abline(intercept = 0, slope = 1, colour = "darkgray")+ geom_line(aes(y = TPR)) + theme(title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5)))

bt_fin_imp = varImp(bt_caret_out_fin)
ggplot(bt_fin_imp, top = 10) + theme(title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5)))

fin_preds$obs_class = factor(ifelse(fin_preds$obs == "vaxxed", "Compliant", "Non-compliant"))
ggplot(fin_preds, aes(unvaxxed, fill = obs_class)) + geom_density(alpha = 0.5) + xlab("Probability of compliance") + guides(fill = guide_legend(title = "Observed class")) + theme(title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.2)))

int_plot = interaction_investigator(bt_fin, plot = T, num_var_plot = 10, num_replicates_for_avg = 3, bottom_margin = 20)
int_plot$interaction_counts_avg

m = int_plot$interaction_counts_avg
diag(m) = 0
largest_ind = data.frame(which(matrix(m %in% head(sort(m, T), 10), nr = nrow(m)), arr.ind = T))
largest_ind$row_name = rownames(m)[largest_ind$row]
largest_ind$col_name = colnames(m)[largest_ind$col]
largest_ind$relimp = diag(int_plot$interaction_counts_avg[largest_ind$row, largest_ind$col])
largest_ind$relimp_sd = diag(int_plot$interaction_counts_sd[largest_ind$row, largest_ind$col])
largest_ind = largest_ind[order(-largest_ind$relimp),]
largest_ind$int = factor(paste0(largest_ind$col_name, " * ", largest_ind$row_name))
largest_ind$int = factor(largest_ind$int, levels = largest_ind$int[order(largest_ind$relimp)])

ggplot(largest_ind, aes(int, relimp)) + geom_bar(stat = "identity") + coord_flip() + xlab("Relative importance") + ylab("Interaction") + theme(title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5))) + guides(fill = guide_legend(title = "Predicted compliance"))

checkup_age_plot = pd_plot(bt_fin, "checkup_age")
age_plot = pd_plot(bt_fin, "age")
inc_pov_rat_plot = pd_plot(bt_fin, "inc_pov_rat")

checkup_age_marginals = adply(12:17, 1, function(a) {
  tbl = table(nistcc[nistcc$checkup_age == a, "vax"])
  return(data.frame(checkup_age = a, unvaxxed = tbl[1] / sum(tbl), vaxxed = tbl[2] / sum(tbl)))
})[, -1]

age_marginals = adply(13:17, 1, function(a) {
  tbl = table(nistcc[nistcc$age == a, "vax"])
  return(data.frame(age = a, unvaxxed = tbl[1] / sum(tbl), vaxxed = tbl[2] / sum(tbl)))
})[, -1]

inc_pov_rat_plot$x_j_quants = c(inc_pov_rat_plot$x_j_quants, Inf)
inc_pov_rat_marginals = adply(2:8, 1, function(a) {
  tbl = table(nistcc[nistcc$inc_pov_rat >= inc_pov_rat_plot$x_j_quants[a - 1] & nistcc$inc_pov_rat < inc_pov_rat_plot$x_j_quants[a], "vax"])
  return(data.frame(inc_pov_rat_gt = inc_pov_rat_plot$x_j_quants[a - 1], unvaxxed = tbl[1] / sum(tbl), vaxxed = tbl[2] / sum(tbl)))
})[, -1]

# need to fix inverted probabilities from pd_plot bug
checkup_age_marginals$probit = -checkup_age_plot$bart_avg_predictions_by_quantile
ggplot(checkup_age_marginals, aes(checkup_age, probit)) + geom_point() + geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = T) + theme(title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)))

age_marginals$probit = -age_plot$bart_avg_predictions_by_quantile
ggplot(age_marginals, aes(age, probit)) + geom_point() + geom_smooth(method = "lm", se = T) + theme(title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)))

inc_pov_rat_marginals$probit = -inc_pov_rat_plot$bart_avg_predictions_by_quantile
ggplot(inc_pov_rat_marginals, aes(inc_pov_rat_gt, probit)) + geom_point() + geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = T) + theme(title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("inc_pov_rat")

# reversed
nistcc$vaxxed_prob = fin_preds$unvaxxed

ggplot(nistcc[, c("checkup_age", "age", "vaxxed_prob")], aes(checkup_age, age)) + geom_raster(aes(fill = vaxxed_prob)) + xlim(c(13, 17)) + theme(title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5))) + guides(fill = guide_legend(title = "Predicted compliance"))

####################################################################### comparison
######################################################################
bt_caret_out_cc = train(trControl = bt_traincon_fin, method = "bartMachine", x = nistcc[complete.cases(nistcc), preds], y = nistcc[complete.cases(nistcc), "vax"], tuneGrid = expand.grid(num_trees = 101, k = 2, alpha = 0.85, beta = 1, nu = 3), use_missing_data = T, use_missing_data_dummies_as_covars = T, mem_cache_for_speed = F, verbose = T)

btcc_preds = bt_caret_out_cc$pred
confusionMatrix(btcc_preds$pred, btcc_preds$obs, positive = "vaxxed")

rf_traincon = trainControl(method = "boot632", number = 10, verboseIter = T, allowParallel = T, savePredictions = "final", classProbs = T)
rf_caret_out = train(trControl = rf_traincon, method = "rf", x = nistcc[complete.cases(nistcc), preds], y = nistcc[complete.cases(nistcc), "vax"], ntree = 201, tuneLength = 3, do.trace = T)

rf_preds = rf_caret_out$pred
confusionMatrix(rf_preds$pred, rf_preds$obs, positive = "vaxxed")

glm_traincon = trainControl(method = "boot632", number = 10, verboseIter = T, allowParallel = T, savePredictions = "final", classProbs = T)
glm_caret_out = train(trControl = rf_traincon, method = "glmStepAIC", x = nistcc[complete.cases(nistcc), preds], y = nistcc[complete.cases(nistcc), "vax"], tuneLength = 1)

glm_preds = glm_caret_out$pred
confusionMatrix(glm_preds$pred, glm_preds$obs, positive = "vaxxed")
