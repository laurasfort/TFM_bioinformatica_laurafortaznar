# Exploración de tipos celulares
table(covariates$oupSample.cellType)

# Pacientes por condición
table(covariates$oupSample.batchCond)

# Tipos celulares por grupo
table(covariates$oupSample.cellType_batchCond)

# Subclusters por condición
table(covariates$oupSample.subclustCond, covariates$oupSample.subclustID)

# Número de pacientes por grupo
length(unique(covariates$oupSample.subclustID[covariates$oupSample.subclustCond == "AD"]))
length(unique(covariates$oupSample.subclustID[covariates$oupSample.subclustCond == "ct"]))
