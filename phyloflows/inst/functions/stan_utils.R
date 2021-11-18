ageanalysis <- function(infile.inference=NULL,infile.prior.samples=NULL,opt=NULL,M=30,D=2,outdir){
  
  require(data.table)	
  #
  #	input args
  #
  if(is.null(opt))
  {
    opt									<- list()
    opt$adjust.sequencing.bias			<- 1
    opt$adjust.participation.bias		<- 1
    opt$migration.def.code				<- '24'
    opt$set.missing.migloc.to.inland	<- 0
    opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
  }
  if(is.null(infile.inference))
  {
    infile.inference	<- file.path(outdir, "RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
  }
  if(is.null(infile.prior.samples))
  {
    infile.prior.samples <- file.path(outdir,"samples_fit.rda")
  }
  
  cat('\ninfile.inference=',infile.inference)
  cat('\ninfile.prior.samples=',infile.prior.samples)
  
  cat('\nopt=',unlist(opt))			
  indir					<- dirname(infile.inference)	
  outfile.base			<- '~/ageanalysis/'
  load(infile.inference)
  #
  #	prepare data on observed transmission flows
  #
  #	subset to variables needed, using RTR3	
  rtr	<- copy(rtr3)
  if(opt$migration.def.code=='06')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_05YR, REC_INMIGRATE_05YR, TR_COMM_NUM_A_MIG_05YR')
    setnames(rtr, 'TR_INMIGRATE_05YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_05YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_05YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='12')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_1YR, REC_INMIGRATE_1YR, TR_COMM_NUM_A_MIG_1YR')
    setnames(rtr, 'TR_INMIGRATE_1YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_1YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_1YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='24')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_2YR, REC_INMIGRATE_2YR, TR_COMM_NUM_A_MIG_2YR')
    setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_2YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='36')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_3YR, REC_INMIGRATE_3YR, TR_COMM_NUM_A_MIG_3YR')
    setnames(rtr, 'TR_INMIGRATE_3YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_3YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_3YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='48')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_4YR, REC_INMIGRATE_4YR, TR_COMM_NUM_A_MIG_4YR')
    setnames(rtr, 'TR_INMIGRATE_4YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_4YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_4YR', 'TR_COMM_NUM_A_MIG')
  }
  
  rtr	<- subset(rtr, select=c('PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
                              'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE',
                              'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
                              'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE'))
  
  #	set unknown origin to either fishing or inland
  tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
  if(opt$set.missing.migloc.to.inland)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
  }		
  if(opt$set.missing.migloc.to.fishing)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fish')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
  }
  
  # add age 
  rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
  rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
  
  # impute age
  tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
  set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
  tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
  set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
  
  # fixup from latest surveillance data
  set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
  set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
  set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
  set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
  
  #	stratify age
  rtr[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(15,16:49,50), labels=paste0(15:49,'-',16:50), right=FALSE))]
  rtr[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(15,16:49,50), labels=paste0(15:49,'-',16:50), right=FALSE))]
  # stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
  # stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
  rtr <- subset(rtr, !is.na(TR_AGE_AT_MID_C))
  rtr <- subset(rtr, !is.na(REC_AGE_AT_MID_C))
  
  # define TR_COMM_TYPE_F, REC_COMM_TYPE_F (i: inland; f: fishing) 
  rtr[,TR_COMM_TYPE_F:=substr(TR_COMM_TYPE,1,1)]
  unique(rtr$TR_COMM_TYPE_F)
  rtr[substr(TR_COMM_TYPE,1,1)!='f',TR_COMM_TYPE_F:='i']
  rtr[,REC_COMM_TYPE_F:=substr(REC_COMM_TYPE,1,1)]
  unique(rtr$REC_COMM_TYPE_F)
  rtr[substr(REC_COMM_TYPE,1,1)!='f',REC_COMM_TYPE_F:='i']
  
  # define TR_COMM_TYPE_F_MIG (i: inland; f: fishing; e: external) 
  rtr[,TR_COMM_TYPE_F_MIG:=substr(TR_COMM_NUM_A_MIG,1,1)]
  unique(rtr$TR_COMM_TYPE_F_MIG)
  rtr[substr(TR_COMM_NUM_A_MIG,1,1)=='a' | substr(TR_COMM_NUM_A_MIG,1,1)=='i'|
        substr(TR_COMM_NUM_A_MIG,1,1)=='t',TR_COMM_TYPE_F_MIG:='i']
  
  #	build category to match with sampling data tables 
  rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_TYPE_F,':',TR_SEX,':',TR_AGE_AT_MID_C)]
  #	build transmission flow category 
  rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_TYPE_F_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C)]
  
  table(paste0(rtr$TR_COMM_TYPE_F_MIG,'-',rtr$TR_COMM_TYPE_F))
  # make all combinations of variables
  dac <- expand.grid( COMM_TYPE_F= c('i','f'),
                      SEX=  c('F','M'),
                      AGE_AT_MID_C= paste0(15:49,'-',16:50))
  
  dac <- as.data.table(dac)  				
  dac[, CATEGORY:= paste0(COMM_TYPE_F, ':', SEX, ':', AGE_AT_MID_C)]  					
  dac <- as.data.table(expand.grid(TR_SAMPLING_CATEGORY= dac$CATEGORY, REC_SAMPLING_CATEGORY= dac$CATEGORY))
  # ignore Male-Male and Female-Female combinations
  dac <- subset(dac, !(grepl('F',TR_SAMPLING_CATEGORY)&grepl('F',REC_SAMPLING_CATEGORY)) &
                  !(grepl('M',TR_SAMPLING_CATEGORY)&grepl('M',REC_SAMPLING_CATEGORY))  
  )
  # add transmission categories
  dac[, REC_TRM_CATEGORY:= REC_SAMPLING_CATEGORY]
  dac[, TR_TRM_CATEGORY:= TR_SAMPLING_CATEGORY]  
  
  # add inmigrants from external communities
  tmp <- copy(dac)
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^[f|i]','e',TR_SAMPLING_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # add inmigrants sampled in inland and migrated from fishing
  tmp <- dac[grepl('^i',TR_SAMPLING_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^i','f',TR_SAMPLING_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # add inmigrants sampled in fishing and migrated from inland
  tmp <- dac[grepl('^f',TR_SAMPLING_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^f','i',TR_SAMPLING_CATEGORY)]) 
  dac <- rbind(dac, tmp)   
  # remove duplicated rows 
  # TR_SAMPLING_CATEGORY f: F: 15-24: 1 and i: F: 15-24: 1 are all set to e: F: 15-24: 1
  dac <- unique(dac)
  # setnames(dac, c('TR_CATEGORY', 'REC_CATEGORY'), c('TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'))
  # 
  
  #	calculate observed number of transmissions
  #
  dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
  dac[, DUMMY:= 1]
  dobs <- merge(dac, dobs, by=c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'), all=TRUE)
  stopifnot( dobs[, !any(is.na(DUMMY))] )
  set(dobs, NULL, 'DUMMY', NULL)
  set(dobs, dobs[, which(is.na(TRM_OBS))], 'TRM_OBS', 0L)
  
  #	make PAIR_ID
  dobs[,TR_SMOOTH_CATEGORY:= as.numeric(substr(TR_TRM_CATEGORY,5,6))+0.5]
  dobs[,REC_SMOOTH_CATEGORY:= as.numeric(substr(REC_TRM_CATEGORY,5,6))+0.5]
  dobs[,OUTPUT:=paste0(substr(TR_TRM_CATEGORY,1,1),':',
                       substr(TR_TRM_CATEGORY,3,3),':',
                       substr(TR_SAMPLING_CATEGORY,1,1),':',
                       substr(REC_TRM_CATEGORY,1,1),':',
                       substr(REC_TRM_CATEGORY,3,3))]
  
  # combinations
  tmp <- subset(dobs,select = 'OUTPUT')
  tmp <- unique(tmp)
  setkey(tmp, OUTPUT)
  tmp[, OUTPUT_ID:= seq_len(nrow(tmp))]
  dobs <- merge(dobs, tmp, by='OUTPUT')
  setkey(dobs,OUTPUT_ID,TR_SMOOTH_CATEGORY,REC_SMOOTH_CATEGORY)	
  dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
  setkey(dobs, TRM_CAT_PAIR_ID)
  
  # sampling prior
  load(infile.prior.samples)
  set(dprior.fit, NULL, 'SAMPLING_CATEGORY', dprior.fit[, gsub('^2','e',gsub('^1','f',gsub('^0','i',SAMPLING_CATEGORY)))] )		
  
  dprior.id <- subset(dprior.fit,select = c('SAMPLING_CATEGORY','WHO'))
  setkey(dprior.id, SAMPLING_CATEGORY,WHO)
  dprior.id[,ID:= seq_len(nrow(dprior.id))]
  
  tmp <- subset(dobs, select = c('TR_SAMPLING_CATEGORY',
                                 'REC_SAMPLING_CATEGORY',
                                 'TRM_CAT_PAIR_ID'))
  
  setnames(dprior.id,colnames(dprior.id),paste0('TR_',colnames(dprior.id)))
  tmp <- merge(tmp,subset(dprior.id[TR_WHO=='TR_SAMPLING_CATEGORY',],select=c('TR_ID','TR_SAMPLING_CATEGORY')),by='TR_SAMPLING_CATEGORY',all.x = TRUE)
  setnames(dprior.id,colnames(dprior.id),gsub('TR_','REC_',colnames(dprior.id)))
  tmp <- merge(tmp,subset(dprior.id[REC_WHO=='REC_SAMPLING_CATEGORY',],select=c('REC_ID','REC_SAMPLING_CATEGORY')),by='REC_SAMPLING_CATEGORY',all.x = TRUE)
  setkey(tmp, TRM_CAT_PAIR_ID)
  xi_id <- cbind(tmp$TR_ID, tmp$REC_ID)
  
  setkey(dprior.fit,SAMPLING_CATEGORY,WHO)
  
  
  indices <- matrix(NA, M^D, D)
  mm=0;
  for (m1 in 1:M){
    for (m2 in 1:M){
      mm = mm+1
      indices[mm,] = c(m1, m2)
    }
  }
  
  B1 <- max(dobs$TR_SMOOTH_CATEGORY)
  B2 <- max(dobs$REC_SMOOTH_CATEGORY)
  L <-  matrix(rep(c(B1, B2) * 5/4,each=nrow(indices)),nrow=nrow(indices))
  sevalue <- pi * indices / (2 * L)
  ns <- nrow(dobs[OUTPUT_ID==1])
  efunc <-  do.call( rbind, lapply(1:ns, 
                                   function(k){
                                     as.vector(apply(sqrt(1/L) * sin(sevalue * (matrix(rep(c(dobs$TR_SMOOTH_CATEGORY[k],dobs$REC_SMOOTH_CATEGORY[k]),each=nrow(indices)),nrow=nrow(indices)) + L)),1, prod))} ))
  
  data.fit <- list(N=nrow(dobs),
                   N_group = max(dobs$OUTPUT_ID),
                   N_per_group = ns,
                   y=matrix(dobs$TRM_OBS,nrow=ns,ncol=max(dobs$OUTPUT_ID)),
                   D=D,
                   x=cbind(dobs$TR_SMOOTH_CATEGORY[1:ns],dobs$REC_SMOOTH_CATEGORY[1:ns]),
                   M_nD=M^2,
                   Xgp=efunc,
                   slambda=sevalue,
                   N_xi = nrow(dprior.fit),
                   shape = cbind(dprior.fit$SHAPE1,dprior.fit$SHAPE2),
                   xi_id_src = matrix(xi_id[,1],nrow=ns,ncol=max(dobs$OUTPUT_ID)),
                   xi_id_rec = matrix(xi_id[,2],nrow=ns,ncol=max(dobs$OUTPUT_ID)),
                   id_mf =   dobs[grepl(':M:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_fm =   dobs[grepl(':F:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_eh = dobs[grepl('e:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_el = dobs[grepl('e:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_hh = dobs[grepl('f:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_hl = dobs[grepl('f:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_lh = dobs[grepl('i:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_mf_h =  dobs[grepl(':M:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_mf_l =  dobs[grepl(':M:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_fm_h =  dobs[grepl(':F:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_fm_l =  dobs[grepl(':F:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_sh =  dobs[grepl('f:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_sl =  dobs[grepl('i:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_se =  dobs[grepl('e:',TR_TRM_CATEGORY),unique(OUTPUT_ID)]
  )
  
  save(data.fit,file=file.path(outdir, 'input.rda'))
  
  return(data.fit)
}

prepare_stan_data <- function(pairs, age_map){
  stan_data = list()
  
  # number of directions
  stan_data[['N_group']] = 2
  
  # f -> M
  tmp <- pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M']
  tmp <- tmp[, list(age_infection.SOURCE = floor(age_infection.SOURCE), age_infection.RECIPIENT = floor(age_infection.RECIPIENT))]
  tmp <- tmp[, list(count = .N), by = c('age_infection.SOURCE', 'age_infection.RECIPIENT')]
  tmp <- merge(age_map, tmp, by = c('age_infection.SOURCE', 'age_infection.RECIPIENT'), all.x = T)
  tmp[is.na(count), count := 0]
  stopifnot(sum(tmp$count) == nrow(pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M']))
  stopifnot(all(tmp$age_infection.SOURCE == age_map$age_infection.SOURCE))
  stopifnot(all(tmp$age_infection.RECIPIENT == age_map$age_infection.RECIPIENT))
  stan_data[['y']] = matrix(tmp$count, ncol = 1)
  stan_data[['is_mf']] = 0
  
  # M -> F
  tmp <- pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F']
  tmp <- tmp[, list(age_infection.SOURCE = floor(age_infection.SOURCE), age_infection.RECIPIENT = floor(age_infection.RECIPIENT))]
  tmp <- tmp[, list(count = .N), by = c('age_infection.SOURCE', 'age_infection.RECIPIENT')]
  tmp <- merge(age_map, tmp, by = c('age_infection.SOURCE', 'age_infection.RECIPIENT'), all.x = T)
  tmp[is.na(count), count := 0]
  stopifnot(sum(tmp$count) == nrow(pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F']))
  stopifnot(all(tmp$age_infection.SOURCE == age_map$age_infection.SOURCE))
  stopifnot(all(tmp$age_infection.RECIPIENT == age_map$age_infection.RECIPIENT))
  stan_data[['y']] = cbind(stan_data[['y']], matrix(tmp$count, ncol = 1))
  stan_data[['is_mf']] = c(stan_data[['is_mf']], 1)
  
  # age age entries
  stan_data[['N_per_group']] = nrow(age_map)

  return(stan_data)
}


add_2D_splines_stan_data = function(stan_data, spline_degree = 3, n_knots_rows = 8, n_knots_columns = 8, AGES)
{
  
  stan_data$A <- length(AGES)
  
  knots_rows = AGES[seq(1, length(AGES), length.out = n_knots_rows)] 
  knots_columns = AGES[seq(1, length(AGES), length.out = n_knots_columns)]
  
  stan_data$num_basis_rows = length(knots_rows) + spline_degree - 1
  stan_data$num_basis_columns = length(knots_columns) + spline_degree - 1
  
  stan_data$IDX_BASIS_ROWS = 1:stan_data$num_basis_rows
  stan_data$IDX_BASIS_COLUMNS = 1:stan_data$num_basis_columns
  
  stan_data$BASIS_ROWS = bsplines(AGES, knots_rows, spline_degree)
  stan_data$BASIS_COLUMNS = bsplines(AGES, knots_columns, spline_degree)
  
  stopifnot(all( apply(stan_data$BASIS_ROWS, 1, sum) > 0  ))
  stopifnot(all( apply(stan_data$BASIS_COLUMNS, 1, sum) > 0  ))
  
  return(stan_data)
}

bspline = function(x, k, order, intervals)
{
  
  if(order == 1){
    return(x >= intervals[k] & x < intervals[k+1])
  }
  
  w1 = 0; w2 = 0
  
  if(intervals[k] != intervals[k+order-1])
    w1 = (x - intervals[k]) / (intervals[k+order-1] - intervals[k])
  if(intervals[k+1] != intervals[k+order])
    w2 = 1 - (x - intervals[k+1]) / (intervals[k+order] - intervals[k+1])
  
  spline = w1 * bspline(x, k, order - 1, intervals) +
    w2 * bspline(x, k+1, order - 1, intervals)
  
  return(spline)
}

find_intervals = function(knots, degree, repeating = T)
{
  
  K = length(knots)
  
  intervals = vector(mode = 'double', length = 2*degree + K)
  
  # support of knots
  intervals[(degree+1):(degree+K)] = knots
  
  # extreme
  if(repeating)
  {
    intervals[1:degree] = min(knots)
    intervals[(degree+K+1):(2*degree+K)] = max(knots)
  } else {
    gamma = 0.1
    intervals[1:degree] = min(knots) - gamma*degree:1
    intervals[(degree+K+1):(2*degree+K)] = max(knots) + gamma*1:degree
  }
  
  return(intervals)
}

bsplines = function(data, knots, degree)
{
  K = length(knots)
  num_basis = K + degree - 1
  
  intervals = find_intervals(knots, degree)
  
  m = matrix(nrow = num_basis, ncol = length(data), 0)
  
  for(k in 1:num_basis)
  {
    m[k,] = bspline(data, k, degree + 1, intervals) 
  }
  
  m[num_basis,length(data)] = 1
  
  return(m)
}

