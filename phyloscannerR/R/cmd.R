#' Obtain valid input arguments for a phyloscanner analysis on a tree or set of trees
#' @param prog.phyloscanner_analyse_trees The full file name of \code{phyloscanner_analyse_trees.R}.
#' @export
#' @seealso \code{\link{cmd.phyloscanner.analyse.trees}}
cmd.phyloscanner.analyse.trees.valid.args<- function(prog.phyloscanner_analyse_trees)
{
	tmp	<- system2(command=prog.phyloscanner_analyse_trees, args='-h', stdout=TRUE)
	tmp	<- tmp[grepl('--',tmp)]
	tmp <- sapply(tmp, function(x) sub('^.*--([a-zA-Z]+).*$','\\1',x))
	valid.input.args <- unique(sort(unname(tmp)))
	valid.input.args <- valid.input.args[valid.input.args!='help']
	valid.input.args
}

#' Make script file for a phyloscanner analysis on a tree or set of trees
#'
#' This function makes a UNIX script file to call \code{phyloscanner_analyse_trees.R}. Usually, this is useful to parallelise computations; see the Examples. 
#' @param prog.phyloscanner_analyse_trees The full file name of \code{phyloscanner_analyse_trees.R}.
#' @param tree.input One of the following: the name of a single tree file (Newick or NEXUS format); the directory containing all input trees; a zip file containing input trees.
#' @param control List of input arguments to \code{\link{phyloscanner.analyse.trees}}.
#' @param valid.input.args Vector of valid input arguments.
#' @return A character string of UNIX commands.
#' @seealso \code{\link{phyloscanner.analyse.trees}}, \code{\link{cmd.phyloscanner.analyse.trees.valid.args}}
#' @author Oliver Ratmann
#' @export
#' @example /inst/example/ex.cmd.phyloscanner_analyse_trees.R
cmd.phyloscanner.analyse.trees<- function(prog.phyloscanner_analyse_trees, 
												tree.input, 
												control,
												valid.input.args=cmd.phyloscanner.analyse.trees.valid.args(prog.phyloscanner_analyse_trees))
{	
	#
	#	prepare input args
	#
	input.args		<- control
	#	check that positional arguments are in control
	stopifnot(any(names(input.args)=='splits.rule'))
	stopifnot(any(names(input.args)=='output.string'))
	#	prepare optional arguments with non-default values
	tmp	<- which(names(input.args)=='guess.multifurcation.threshold')
	if(length(tmp)>0)
	{
		if(input.args[['guess.multifurcation.threshold']])
			input.args[['multifurcation.threshold']] <- 'g'
		input.args <- input.args[names(input.args)!='guess.multifurcation.threshold']
	}
	tmp	<- which(names(input.args)=='sankoff.k')
	if(length(tmp)>0)
	{
		if(input.args[['splits.rule']]=='s')
			input.args[['splits.rule']] <- paste0(input.args[['splits.rule']],',',input.args[['sankoff.k']])
		input.args <- input.args[names(input.args)!='sankoff.k']
	}
	#	extract positional arguments
	splits.rule			<- input.args[['splits.rule']]
	input.args			<- input.args[names(input.args)!='splits.rule']
	output.string 		<- input.args[['output.string']]
	input.args			<- input.args[names(input.args)!='output.string']
	#	extract out.dir
	out.dir				<- input.args[['output.dir']]
	input.args			<- input.args[names(input.args)!='output.dir']
	#	prepare optional argument names that are slightly inconsistent
	names(input.args)	<- gsub('use.ff','useff',names(input.args))
	names(input.args)	<- gsub('do.dual.blacklisting','dual.blacklist',names(input.args))
	names(input.args)	<- gsub('allow.mt','allow.multi.trans',names(input.args))
	names(input.args)	<- gsub('n.mt','n.multi.trans',names(input.args))
	names(input.args)	<- gsub('p.mt','p.multi.trans',names(input.args))
	names(input.args)	<- gsub('count.reads.in.parsimony','read.counts.matter.on.zero.length.branches',names(input.args))
	names(input.args)	<- gsub('verbosity','verbose',names(input.args))
	#	ignore arguments that are hard coded in the Rscript
	input.args			<- input.args[names(input.args)!='tree.file.regex']
	input.args 			<- input.args[names(input.args)!='sankoff.unassigned.switch.threshold']
	#	replace .a with A where needed
	tmp <- strsplit(names(input.args),'\\.')
	tmp <- sapply(tmp, function(x) gsub('^([A-Z])','\\L\\1',paste(gsub('^([a-z])','\\U\\1',x,perl=TRUE),collapse=''),perl=TRUE))
	tmp2 <- which(! names(input.args) %in% valid.input.args )	
	names(input.args)[tmp2] <- tmp[tmp2]
	#	add default optional arguments
	input.args[['overwrite']]	<- TRUE
	input.args[['outputRDA']]	<- TRUE
	#	remove logical arguments that evaluate to FALSE
	tmp 		<- !unname(sapply(input.args, function(x) is.logical(x) && x==FALSE))
	input.args	<- input.args[tmp]			
	#	check that all arguments are valid
	tmp2 <- which(! names(input.args) %in% valid.input.args )
	if(length(tmp2)>0)
	{
		stop('Found invalid arguments,',input.args[tmp2], ' for ',names(input.args)[tmp2])
	}	
	#	for all character arguments: add encapsulating "" except for 'distanceThreshold'
	#	for all logical arguments that evaluate to TRUE: keep only the name
	for(ii in seq_along(input.args))
	{
		if(is.character(input.args[[ii]]) && !substr(input.args[[ii]],1,1)%in%c("'",'"') && names(input.args)[ii]!='distanceThreshold')
			input.args[[ii]] <- paste0('"',input.args[[ii]],'"')
		if(is.logical(input.args[[ii]]) && input.args[[ii]])
			input.args[[ii]] <- ''
	}	
	#	sort arguments by name
	tmp			<- sort(names(input.args), index.return=TRUE)$ix
	input.args	<- input.args[tmp]
	#
	#	make command
	#
	#	create local tmp dir
	cmd		<- paste("CWD=$(pwd)\n",sep='\n')
	cmd		<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir	<- paste("$CWD/",tmpdir,sep='')
	cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	if tree.input is zip file, extract and change tree.input to directory
	if(grepl('\\.zip$',tree.input))
	{		
		tree.dir	<- gsub('\\.zip$','',file.path(tmpdir,basename(tree.input)))
		cmd			<- paste(cmd,'mkdir -p "',tree.dir,'"\n',sep='')
		cmd			<- paste(cmd,'unzip -j "',tree.input,'" -d "',tree.dir,'"\n',sep='')
		#	passing directory is currently not supported, must pass base of tree file names before the window coordinates start
		#	determine prefix of tree files via control$tree.file.regex
		#	the regex must contain ()
		if( !grepl('\\(',control$tree.file.regex) )
			stop('Cannot make tree.input variable, expect control$tree.file.regex with () that identify window coordinates, found:',control$tree.file.regex)	
		tmp				<- unzip(tree.input, list=TRUE)		
		tmp$windowid	<- gsub(control$tree.file.regex, '\\1', tmp$Name)
		tmp$prefix		<- sapply(seq_along(tmp$Name), function(x) gsub(paste0(tmp$windowid[x],'.*'),'',tmp$Name[x]))		
		if( !all( tmp$prefix==tmp$prefix[1] ) )
			stop('Cannot make tree.input variable, contact maintainer', tmp$Name)		
		tree.input	<- file.path(tree.dir,tmp$prefix[1])
	}	
	#	add encapsulating "" to tree.input
	if(!substr(tree.input,1,1)%in%c("'",'"'))
		tree.input	<- paste0('"',tree.input,'"')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	#	add call with positional arguments
	cmd		<- paste(cmd, prog.phyloscanner_analyse_trees,' ',tree.input,' ',output.string,' ',splits.rule, sep='')
	#	add optional arguments
	for(ii in seq_along(input.args))
	{
		cmd			<- paste0(cmd,' --',names(input.args)[[ii]])
		if(input.args[[ii]]!='')
			cmd		<- paste0(cmd,' ',input.args[[ii]])
	}
	#	copy to outdir
	cmd		<- paste(cmd, '\nmv ',output.string,'_workspace.rda "',out.dir,'"\n',sep='')
	#	remove trees directory
	if(exists('tree.dir'))
		cmd		<- paste0(cmd, 'rm -rf "',tree.dir,'"\n')
	#	zip up everything else
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(output.string,'_otherstuff.zip',sep=''),' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',paste(output.string,'_otherstuff.zip',sep=''),' "',out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -rf "',tmpdir,'"\n',sep='')
	cmd
}
