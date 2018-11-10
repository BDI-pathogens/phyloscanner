#' @export
#' @import tibble
#' @author Oliver Ratmann
#' @title Determine stage 1 phyloscanner runs
#' @description This function groups individuals for phyloscanner analyses, so that phylogenetic linkage between every pair of individuals is assessed at least once.
#'   Specifically, individuals are grouped into batches of specified size, and then, all possible pairs of batches are formed. 
#'   Each of these pairs of batches defines a group of individuals between whom phylogenetic linkages are assessed in one phyloscanner run.
#'   The number of individuals in each group is twice the batch size.
#' @param x Character vector of individual identifiers. 
#' @param batch.size Batch size. Default is 50.
#' @return tibble with rows 'IND' (individual identifiers), 'PTY_RUN' group for phyloscanner analysis, and 'BATCH' batch of individuals (not used further, but there should be two batches of individuals in each phyloscanner analysis).
#' @examples
#' x <- c("15-01402","15-04719","16-00616","16-00801","16-01173","16-01191","16-01302","16-01408","16-01414","16-01465","16-01581","16-01663","16-03259","16-03499","16-03516","16-03644","16-03807")
#' pty.runs	<- phyloscannerR::define.stage1.analyses(x, batch.size=50)
define.stage1.analyses<- function(x, batch.size=50)	
{
	dind	<- tibble(IND=x)
	dind	<- dind %>% mutate( BATCH:= 1L+floor((seq_len(nrow(dind))-1)/batch.size) )
	batches	<- dind %>% distinct(BATCH) %>% unlist() %>% unname()
	pty.runs<- t(combn(batches,2))
	pty.runs<- pty.runs %>% 
			as_tibble() %>% 
			rename(BATCH=V1, BATCH2=V2) %>%
			mutate( PTY_RUN:= seq_along(BATCH) )
	tmp		<- full_join(pty.runs, dind, by='BATCH')
	dind	<- dind %>% rename(IND2=IND, BATCH2=BATCH)
	tmp2	<- full_join(pty.runs, dind, by='BATCH2')
	dind	<- dind %>% rename(IND=IND2, BATCH=BATCH2)
	pty.runs<- tmp2 %>% 
			rename(IND=IND2) %>% 
			rbind(tmp) %>% 
			arrange(PTY_RUN) %>%
			select(PTY_RUN, IND, everything()) 
	pty.runs
}