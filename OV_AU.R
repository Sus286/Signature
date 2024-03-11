# Controllare se la stringa "cercata" è presente in ogni riga della colonna "colonna_da_controllare"
presenza_AF <- grepl("AF", out_OV_AU_FILTRO1_id_vep$V14)

out_OV_AU_FILTRO1_id_vep_clinsig_af<-out_OV_AU_FILTRO1_id_vep
presenza_CLINSIG<-grepl("CLIN_SIG",out_OV_AU_FILTRO1_id_vep_clinsig_af$V14)
#out_brca_fr_filtro_id_vep_clinsig_af<-out_brca_fr_filtro_id_vep_clinsig_af[presenza_CLINSIG==TRUE,]

soglia<-out_OV_AU_FILTRO1_id_vep_clinsig_af[,14]

# Divisione della colonna in base al punto e virgola
soglia_divise <- strsplit(as.character(soglia), ";")

# Calcola il massimo numero di colonne dopo la divisione
max_colonne <- max(sapply(soglia_divise, length))

# Crea un nuovo dataframe con le colonne divise
dati_divisi_soglia <- data.frame(matrix(NA, nrow = nrow(out_OV_AU_FILTRO1_id_vep_clinsig_af), ncol = max_colonne))

# Assegna i valori alle nuove colonne
for (i in 1:max_colonne) {
  dati_divisi_soglia[, i] <- sapply(soglia_divise, function(x) ifelse(length(x) >= i, x[i], NA))
}

# Funzione per trovare la colonna con "AF" per ogni riga
trova_colonna_AF <- function(row) {
  colonna_con_AF <- grep("AF", row)
  if (length(colonna_con_AF) > 0) {
    return(colonna_con_AF)
  } else {
    return(NA)
  }
}
# Funzione per trovare la colonna con "CLIN_SIG" per ogni riga
trova_colonna_CLINSIG <- function(row) {
  colonna_con_CLINSIG <- grep("CLIN_SIG", row)
  if (length(colonna_con_CLINSIG) > 0) {
    return(colonna_con_CLINSIG)
  } else {
    return(NA)
  }
}

# Trova la colonna con "AF" per ogni riga
colonne_con_AF <- apply(dati_divisi_soglia, 1, trova_colonna_AF)
#Trova la colonna con "CLINSIG" per ogni riga 
colonne_con_CLINSIG<- apply(dati_divisi_soglia, 1, trova_colonna_CLINSIG)

# Crea un nuovo dataframe con solo la colonna con "AF" per ogni riga
dati_selezionati_AF <- data.frame(
  AF_colonna = sapply(1:nrow(dati_divisi_soglia), function(i) ifelse(!is.na(colonne_con_AF[i]), dati_divisi_soglia[i, colonne_con_AF[i]], NA))
)
# Crea un nuovo dataframe con solo la colonna con "CLINSIG" per ogni riga
dati_selezionati_CLINSIG<- data.frame(
  CLINSIG_colonna = sapply(1:nrow(dati_divisi_soglia), function(i) ifelse(!is.na(colonne_con_CLINSIG[i]), dati_divisi_soglia[i, colonne_con_CLINSIG[i]], NA))
)

parte_numerica_AF <- as.data.frame(as.numeric(gsub("[^0-9.]", "", dati_selezionati_AF$AF_colonna)))
parte_numerica_CLINSIG <- sapply(dati_selezionati_CLINSIG, function(x) sub(".*=", "", x))

out_OV_AU_SOGLIA<-cbind(out_OV_AU_FILTRO1_id_vep_clinsig_af,parte_numerica_AF)
out_OV_AU_completed<-cbind(out_OV_AU_SOGLIA,parte_numerica_CLINSIG)

colnames(out_OV_AU_completed)[15]<-"AF"
colnames(out_OV_AU_completed)[16]<-"CLINSIG"

out_OV_AU_completed_filtro<-subset(out_OV_AU_completed,(AF<0.01)|is.na(AF))

# Rimuovi le righe in cui colonna1 è "valore1" o "valore2" (ma mantieni gli NA)
out_OV_AU_completed_filtro <- subset(out_OV_AU_completed_filtro, !(CLINSIG %in% c("conflicting_interpretations_of_pathogenicity,likely_benign","benign,likely_benign,risk_factor","uncertain_significance,benign,pathogenic","benign,likely_benign,uncertain_significance,conflicting_interpretations_of_pathogenicity","uncertain_significance,conflicting_interpretations_of_pathogenicity,likely_benign","benign,conflicting_interpretations_of_pathogenicity","uncertain_significance,benign/likely_benign","likely_benign", "benign","benign,likely_benign","benign/likely_benign"
                                                                                      ,"likely_benign,uncertain_significance","uncertain_significance,likely_benign","uncertain_significance,benign","uncertain_significance,benign,likely_benign","benign,benign/likely_benign,likely_benign","benign,uncertain_significance,benign/likely_benign,conflicting_interpretations_of_pathogenicity","benign,uncertain_significance,likely_benign,benign/likely_benign","benign,likely_benign,benign/likely_benign","uncertain_significance,likely_benign,conflicting_interpretations_of_pathogenicity","benign/likely_benign,conflicting_interpretations_of_pathogenicity,uncertain_significance,benign,likely_benign","uncertain_significance,benign,conflicting_interpretations_of_pathogenicity,likely_benign","uncertain_significance,pathogenic,likely_benign,likely_pathogenic","likely_benign,benign/likely_benign","benign/likely_benign,benign")) | is.na(CLINSIG))

OV_AU_CUT_completed <- OV_AU_FILTRO1_ID[which(OV_AU_FILTRO1_ID$icgc_mutation_id.1 %in% unique(out_OV_AU_completed_filtro$V1)),]

#PACA_AU_completed_inputmatrix<-PACA_AU_CUT_completed[,c(3,5,43,13,14,9,10,11,15,17,18)]

OV_AU_completed_sample_inputmatrix<-OV_AU_CUT_completed[,c(3,5,5,13,14,9,10,11,15,17,18)]

OV_AU_completed_sample_inputmatrix<-OV_AU_completed_sample_inputmatrix[which(OV_AU_completed_sample_inputmatrix$mutation_type=="single base substitution"),]

#BRCA_FR_completed_inputmatrix[,10]="SOMATIC"

#BRCA_FR_completed_inputmatrix[,4]="SNP"
OV_AU_completed_sample_inputmatrix[,11]="SOMATIC"

OV_AU_completed_sample_inputmatrix[,5]="SNP"


# colnames(PACA_AU_completed_sample_inputmatrix)[1]<-"Project"
# colnames(BRCA_FR_completed_inputmatrix)[2]<-"Sample"
# colnames(BRCA_FR_completed_inputmatrix)[3]<-"Genome"
# colnames(BRCA_FR_completed_inputmatrix)[4]<-"mut_type"
# colnames(BRCA_FR_completed_inputmatrix)[5]<-"chrom"
# colnames(BRCA_FR_completed_inputmatrix)[6]<-"pos_start"
# colnames(BRCA_FR_completed_inputmatrix)[7]<-"pos_end"
# colnames(BRCA_FR_completed_inputmatrix)[8]<-"ref"
# colnames(BRCA_FR_completed_inputmatrix)[9]<-"alt"
# colnames(BRCA_FR_completed_inputmatrix)[10]<-"Type"

colnames(OV_AU_completed_sample_inputmatrix)[1]<-"Project"
colnames(OV_AU_completed_sample_inputmatrix)[2]<-"Sample"
colnames(OV_AU_completed_sample_inputmatrix)[3]<-"ID"
colnames(OV_AU_completed_sample_inputmatrix)[4]<-"Genome"
colnames(OV_AU_completed_sample_inputmatrix)[5]<-"mut_type"
colnames(OV_AU_completed_sample_inputmatrix)[6]<-"chrom"
colnames(OV_AU_completed_sample_inputmatrix)[7]<-"pos_start"
colnames(OV_AU_completed_sample_inputmatrix)[8]<-"pos_end"
colnames(OV_AU_completed_sample_inputmatrix)[9]<-"ref"
colnames(OV_AU_completed_sample_inputmatrix)[10]<-"alt"
colnames(OV_AU_completed_sample_inputmatrix)[11]<-"Type"

# Estrai la parte prima del delimitatore "-" per ogni riga

# BRCA_FR_completed_inputmatrix_ID<-BRCA_FR_completed_inputmatrix
# 
# BRCA_FR_completed_inputmatrix_ID$parte_prima <- sub("-.*", "", BRCA_FR_completed_inputmatrix$ID)
# colnames(BRCA_FR_completed_inputmatrix_ID)[12]<-"ID_ok"
# BRCA_FR_completed_inputmatrix_ID$ID<-BRCA_FR_completed_inputmatrix_ID$ID_ok
# # Supponiamo che il tuo dataframe sia chiamato "dati" e vuoi rimuovere la colonna "colonna_da_rimuovere"
# BRCA_FR_completed_inputmatrix_ID<- subset(BRCA_FR_completed_inputmatrix_ID, select = -12)

# Imposta il percorso e il nome del file di output
# file_path <- "C:/Users/susan/OneDrive/Desktop/Tesi/BRCA_FR_AF_CLINSIG_superFILTRO_inputmatrix.txt"
# 
# # Trasforma il dataframe in un vettore di righe di testo
# lines <- apply(BRCA_FR_completed_inputmatrix, 1, paste, collapse = "\t")
# 
# # Scrivi le righe di testo nel file di testo utilizzando cat()
# cat(c(paste(names(BRCA_FR_completed_inputmatrix), collapse = "\t"), lines), file = file_path, sep = "\n")

# Imposta il percorso e il nome del file di output
file_path <- "C:/Users/susan/OneDrive/Desktop/Tesi/OV_AU_AF_superFILTRO_sample_inputmatrix.txt"

# Trasforma il dataframe in un vettore di righe di testo
lines <- apply(OV_AU_completed_sample_inputmatrix, 1, paste, collapse = "\t")

# Scrivi le righe di testo nel file di testo utilizzando cat()
cat(c(paste(names(OV_AU_completed_sample_inputmatrix), collapse = "\t"), lines), file = file_path, sep = "\n")

# Imposta il percorso e il nome del file di output
# file_path <- "C:/Users/susan/OneDrive/Desktop/Tesi/BRCA_FR_AF_CLINSIG_superFILTRO_ID_inputmatrix.txt"
# 
# # Trasforma il dataframe in un vettore di righe di testo
# lines <- apply(BRCA_FR_completed_inputmatrix_ID, 1, paste, collapse = "\t")
# 
# # Scrivi le righe di testo nel file di testo utilizzando cat()
# cat(c(paste(names(BRCA_FR_completed_inputmatrix_ID), collapse = "\t"), lines), file = file_path, sep = "\n")
