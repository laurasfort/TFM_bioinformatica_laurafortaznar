# Tablas exploratorias
table(covariates$oupSample.cellType)
table(covariates$oupSample.batchCond)
table(covariates$oupSample.cellType_batchCond)
table(covariates$oupSample.subclustCond, covariates$oupSample.subclustID)

# Número de pacientes por grupo
length(unique(covariates$oupSample.subclustID[covariates$oupSample.subclustCond == "AD"]))
length(unique(covariates$oupSample.subclustID[covariates$oupSample.subclustCond == "ct"]))
