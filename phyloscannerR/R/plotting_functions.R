# Align the graph output

#' @keywords internal
#' @export AlignPlots

AlignPlots <- function(...) {
  
  LegendWidth <- function(x) x$grobs[[15]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]], unit(abs(diff(c(LegendWidth(x), max.legends.width))), "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

# Add coloured rectangles to the graphs to represent areas with no coverage

#' @keywords internal
#' @export add.no.data.rectangles

add.no.data.rectangles <- function(graph, rectangle.coords, log = F, y.limits = NULL){
  
  if(is.null(y.limits)){
    y.limits <- ggplot_build(graph)$layout$panel_params[[1]]$y.range
  }
  
  if(nrow(rectangle.coords)>0){
    for(rect.no in seq(1, nrow(rectangle.coords))){
      if(log){
        graph <- graph + annotate("rect", xmin=rectangle.coords$starts[rect.no], xmax = rectangle.coords$ends[rect.no], 
                                  ymin=10^(y.limits[1]), ymax=10^(y.limits[2]), fill=rectangle.coords$colour[rect.no], alpha=0.5)
      } else {
        graph <- graph + annotate("rect", xmin=rectangle.coords$starts[rect.no], xmax = rectangle.coords$ends[rect.no], 
                                  ymin=y.limits[1], ymax=y.limits[2], fill=rectangle.coords$colour[rect.no], alpha=0.5)
      }
    }
  }
  
  return(graph)
}

# Calculate the coordinates to build the rectangles for. 
# missing.coords are x-coordinates of windows with no coverage. 
# all.coords is the  full list of window coordinates.

#' @keywords internal
#' @export form.rectangles

form.rectangles <- function(missing.coords, all.coords, colour = "grey"){
  if(length(missing.coords)>0){
    gap <-  unique(all.coords[2:length(all.coords)] - all.coords[1:length(all.coords)-1])
    if(length(gap)>1){
      cat("Fatal error drawing rectangles for missing data")
      quit(save="no", status=1)
    }
    temp.starts <- missing.coords - gap/2
    temp.ends <- missing.coords + gap/2
    
    final.starts <- setdiff(temp.starts, temp.ends)
    final.ends <- setdiff(temp.ends, temp.starts)
    
    return(tibble(starts = final.starts, ends = final.ends, colour=colour))
  } 
  stop("No missing coordinates given")
}

# Try to detect if the windows are evenly spaced. If they aren't, skip the rectangle-drawing, it's too complex.

#' @keywords internal
#' @export find.gaps

find.gaps <- function(xcoords, x.limits){
  gaps <- xcoords[2:length(xcoords)]-xcoords[1:length(xcoords)-1]
  
  if(length(unique(gaps))==1){
    # perfect
    regular.gaps <- T
    new.xcoords <- xcoords
    missing.window.rects <- NULL
    
  } else {
    # if the length of every gap is a multiple of the smallest such gap length, it suggests missing windows
    smallest.gap <- min(gaps)
    if(length(which(gaps/smallest.gap != floor(gaps/smallest.gap)))==0){
      regular.gaps <- T
      new.xcoords <- seq(min(xcoords), max(xcoords), smallest.gap)
      missing.window.rects <- form.rectangles(setdiff(new.xcoords, xcoords), new.xcoords, "grey20")
    } else {
      # I don't know what to do and I am giving up
      
      regular.gaps <- F
      new.xcoords <- xcoords
      missing.window.rects <- NULL
    }
  }
  
  range <- max(new.xcoords) - min(new.xcoords)
  bar.width <- range/(1.5*(length(new.xcoords)+1))
  if(xcoords[1] - (bar.width/2) < x.limits[1]){
    bar.width <- 2*(xcoords[1] - x.limits[1])
  }
  if(xcoords[length(xcoords)] + (bar.width/2) > x.limits[2]){
    bar.width <- 2*(x.limits[2] - xcoords[length(xcoords)])
  }
  
  return(list(x.coordinates = new.xcoords, regular.gaps = regular.gaps, rectangles.for.missing.windows = missing.window.rects, bar.width=bar.width))
}

#' @keywords internal
#' @export produce.pdf.graphs

produce.pdf.graphs <- function(file.name, sum.stats, hosts, xcoords, x.limits, missing.window.rects, bar.width, regular.gaps = F, width=8.26772, height=11.6929, readable.coords = F, verbose = F){
  
  pdf(file=file.name, width=width, height=height)
  
  for (i in seq(1, length(hosts))) {
    host <- hosts[i]
    # tryCatch({
    if(i!=1){
      grid.newpage()
    }
    produce.host.graphs(sum.stats, host, xcoords, x.limits, missing.window.rects, bar.width, regular.gaps, readable.coords, verbose)
  } #,
  #     error=function(e){
  #       if (verbose) cat("Skipping graphs for host ",host," as no reads are present and not blacklisted.\n", sep="")
  #     })
  # }
  # 

  invisible(dev.off())
}

#' @keywords internal
#' @export produce.host.graphs

produce.host.graphs <- function(sum.stats, host, xcoords, x.limits, missing.window.rects, bar.width, regular.gaps = F, readable.coords = F, verbose = F){
  
  x.axis.label <- if(readable.coords) "Window centre" else "Tree number"
  
  host.stats <- sum.stats %>% 
    filter(host.id==host) %>% 
    arrange(xcoord)
  
  if(nrow(host.stats) > 0){
    
    plot.list <- list()
    

    if (verbose) cat("Drawing graphs for host ",host,"...\n", sep="")
    
    # The coordinates for windows missing from all patients are known. These are those missing for this patient, but not all patients.
    
    missing.read.rects <- NULL
    
    if(regular.gaps & length(which(host.stats$reads==0))){
      missing.read.rects <- form.rectangles(host.stats$xcoord[which(host.stats$reads==0)], xcoords, "grey")
    }
    

    
    # rbind on two NULLs still makes a NULL
    
    missing.rects <- bind_rows(missing.window.rects, missing.read.rects)
    
    host.stats <- host.stats %>% 
      filter(reads > 0)
    
    # Graph 1: read count and tip count
    
    host.stats.gd.1 <- host.stats %>% 
      select(xcoord, tips, reads) %>% 
      gather(variable, value, tips, reads)
    
    # if the difference between the largest and smallest values is greater than 10, we want log scale. Otherwise, normal scale.
    
    log.scale <- max(host.stats.gd.1$value - min(host.stats.gd.1$value) >= 10)
    y.limits <- NULL
  
    graph.1 <- ggplot(host.stats.gd.1, aes(x=xcoord, y=value, col=variable))
    
    graph.1 <- graph.1 + geom_point(na.rm=TRUE) +
      theme_bw() + 
      ylab("Tip or read count") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      scale_color_discrete(name="Variable", labels=c("Reads", "Tips")) + 
      theme(text = element_text(size=7))
    
    if(log.scale){
      max.value <- max(host.stats.gd.1$value)
      log.upper.tick <- ceiling(log10(max.value))
      ticks <- 10^(seq(0, log.upper.tick))
      y.limits <- c(0.8, 1.2*(10^log.upper.tick))
      
      graph.1 <- graph.1 + scale_y_log10(breaks=ticks, limits=y.limits)
    }
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.1 <- add.no.data.rectangles(graph.1, missing.rects, log.scale, if(is.null(y.limits)) NULL else if(log.scale) log10(y.limits) else y.limits)
    }
    
    plot.list[["reads.and.tips"]] <- graph.1

    # Graph 2: subgraph and clade counts
    
    host.stats.gd.2 <- host.stats %>% 
      select(xcoord, subgraphs, clades) %>% 
      gather(variable, value, subgraphs, clades)
    
    log.scale <- max(host.stats.gd.2$value - min(host.stats.gd.2$value) >= 10)
    y.limits <- NULL
    
    if(log.scale){
      max.value <- max(host.stats.gd.2$value)
      log.upper.tick <- ceiling(log10(max.value))
      ticks <- 10^(seq(0, log.upper.tick))
      y.limits <- c(0.8, 1.2*(10^log.upper.tick))
    }
    
    graph.2 <- ggplot(host.stats.gd.2, aes(x=xcoord, y=value))
    
    graph.2 <- graph.2 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Subgraph or clade count") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      scale_color_discrete(name="Variable", labels=c("Clades", "Subgraphs")) +
      scale_shape_manual(values=c(1,19), name="Variable", labels=c("Clades", "Subgraphs")) +
      scale_size_manual(values=c(2,1), name="Variable", labels=c("Clades", "Subgraphs")) +
      theme(text = element_text(size=7))
    
    if(max(host.stats.gd.2$value)==1){
      graph.2 <- graph.2 + scale_y_continuous(breaks = c(0,1))
    } else if(max(host.stats.gd.2$value)==2){
      graph.2 <- graph.2 + expand_limits(y=0)+ scale_y_continuous(breaks = c(0,1,2)) 
    } else if(!log.scale) {
      graph.2 <- graph.2 + expand_limits(y=0) + scale_y_continuous(breaks = pretty_breaks()) 
    } else {
      graph.2 <- graph.2 + scale_y_log10(breaks=ticks, limits=y.limits)
    }
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.2 <- add.no.data.rectangles(graph.2, missing.rects, log.scale, if(is.null(y.limits)) NULL else if(log.scale) log10(y.limits) else y.limits)
    }
    
    plot.list[["subgraphs.and.clades"]] <- graph.2
    
    # Graph 3: root-to-tip distances for all reads and largest subgraph
    
    host.stats.gd.3 <- host.stats %>% 
      select(xcoord, overall.rtt, largest.rtt) %>% 
      pivot_longer(2:3, names_to = "variable") %>%
      mutate(variable = factor(variable, levels=c("overall.rtt", "largest.rtt")))
    
    graph.3 <- ggplot(host.stats.gd.3, aes(x=xcoord, y=value))
    
    graph.3 <- graph.3 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Mean root-to-tip-distance\n(read-weighted)") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph"), h = c(0, 360) + 105) +
      scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subgraph")) +
      scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subgraph")) +
      theme(text = element_text(size=7))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.3 <- add.no.data.rectangles(graph.3, missing.rects)
    }
    
    plot.list[["root.to.tip"]] <- graph.3

    
    # Graph 4: largest patristic distances overall and in largest subgraph
    
    host.stats.gd.4 <- host.stats %>% 
      select(xcoord, global.mean.pat.distance, subgraph.mean.pat.distance) %>% 
      gather(variable, value, global.mean.pat.distance, subgraph.mean.pat.distance) %>%
      mutate(variable = factor(variable, levels=c("global.mean.pat.distance", "subgraph.mean.pat.distance")))
    
    graph.4 <- ggplot(host.stats.gd.4, aes(x=xcoord, y=value))
    
    graph.4 <- graph.4 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Mean pairwise patristic distance") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph"), h = c(0, 360) + 105) +
      scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subgraph")) +
      scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subgraph")) +
      theme(text = element_text(size=7))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.4 <- add.no.data.rectangles(graph.4, missing.rects)
    }
    
    plot.list[["patristic.distance"]] <- graph.4
    
    # graph 5: read proportions in each subgraph
    
    proportion.column.names <- colnames(host.stats)[which(substr(colnames(host.stats), 1, 8)=="prop.gp.")]
    
    # only columns with nonzero entries should be kept
    
    sums <- host.stats %>% 
      select(proportion.column.names)%>% 
      summarise_all(sum) %>% 
      unlist()
    
    proportion.column.names <- proportion.column.names[sums>0]
    
    host.stats.gd.5 <- host.stats %>% 
      select(xcoord, proportion.column.names) %>%
      gather(variable, value, proportion.column.names) %>%
      filter(!is.na(value))
    
    host.stats.gd.5 <- host.stats.gd.5 %>% 
      mutate(ngroup = map_int(variable, function(x) which(unique(host.stats.gd.5$variable)==x))) %>%
      mutate_at("ngroup", as.factor)
    
    graph.5 <- ggplot(host.stats.gd.5, aes(x=xcoord, weight=value, fill=reorder(ngroup, rev(order(ngroup)))))
    
    graph.5 <- graph.5 +
      geom_bar(width=bar.width, colour="black", lty="blank") +
      theme_bw() + 
      ylab("Proportion of reads\nin different subraphs") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      scale_fill_viridis(discrete=T, begin=1, end=0, option="plasma") +
      theme(text = element_text(size=7)) + 
      guides(fill = guide_legend(title = "Subgraph rank\n(by tip count)", keywidth = 1, keyheight = 0.4))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.5 <- add.no.data.rectangles(graph.5, missing.rects)
    }
    
    plot.list[["subgraph.read.proportions"]] <- graph.5
  
    # graph 6: recombination metric
    
    if("recombination.metric" %in% names(host.stats)) {     
      graph.6 <- ggplot(host.stats, aes(x=xcoord, y=recombination.metric))
      y.label <- "Recombination metric"
      
      graph.6 <- graph.6 +
        geom_point(alpha = 0.5, na.rm=TRUE) +
        theme_bw() + 
        ylab(y.label) +
        xlab(x.axis.label) +
        scale_x_continuous(limits=x.limits) +
        expand_limits(y=0) +
        #      scale_color_discrete(name="Tip set", labels=c("Longest branch", "Greatest patristic distance")) + 
        theme(text = element_text(size=7))
      
      if(regular.gaps & !is.null(missing.rects)){
      }
      
      plot.list[["recombination.metric"]] <- graph.6
    }
    
    if("solo.dual.count" %in% names(host.stats)) {
      
      graph.7 <- ggplot(host.stats, aes(x=xcoord, y=solo.dual.count))
      y.label <- "Multiplicity of infection"
      
      graph.7 <- graph.7 +
        geom_point(alpha = 0.5, na.rm=TRUE) +
        theme_bw() + 
        ylab(y.label) +
        xlab(x.axis.label) +
        scale_x_continuous(limits=x.limits) +
        expand_limits(y=0) +
        #      scale_color_discrete(name="Tip set", labels=c("Longest branch", "Greatest patristic distance")) + 
        theme(text = element_text(size=7))
      
      if(max(host.stats$solo.dual.count)==1){
        graph.7 <- graph.7 + scale_y_continuous(breaks = c(0,1))
      } else if(max(host.stats$solo.dual.count)==2){
        graph.7 <- graph.7 + expand_limits(y=0)+ scale_y_continuous(breaks = c(0,1,2)) 
      } else if(!log.scale) {
        graph.7 <- graph.7 + expand_limits(y=0) + scale_y_continuous(breaks = pretty_breaks()) 
      }
      
      if(regular.gaps & !is.null(missing.rects)){
        graph.7 <- add.no.data.rectangles(graph.7, missing.rects)
      }
      
      plot.list[["dual.infections"]] <- graph.7
      
    }
    
    all.plots <- do.call(AlignPlots, plot.list)
    
    pushViewport(viewport(layout = grid.layout(length(all.plots) + 1, 1, heights = unit(c(0.25, rep(1,length(all.plots))), "null") )))
    grid.text(host, gp=gpar(fontsize=20), vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    for(plot.no in 1:length(all.plots)){
      plot <- all.plots[[plot.no]]
      plot$vp = viewport(layout.pos.col = 1, layout.pos.row = plot.no + 1)
      grid.draw(plot)
    }
  } else {
    stop("Cannot draw graphs for host ",host," as no reads are present and not blacklisted.")
  }
}

#' Draw bar graphs of pairwise topological/distance relationships, version 2
#' @author Oliver Ratmann
#' @description 
#' This function generates scan plots that summarize reconstructed viral phylogenetic relationships of two individuals.
#' Several pairs of individuals can be processed simultaneously. For each pair of individuals, the scan plot shows the 
#' phylogenetic distance on the y-axis and topological relationship in colours between subgraphs from both individuals in each 
#' deep-sequence phylogeny across the genome. The genomic position on the x-axis indicates the start of each read alignment. 
#' @param ptrees A list of class \code{phyloscanner.trees}
#' @param hosts A list of hosts (as a vector) to obtain graphs for. By default, all pairs of hosts detected in \code{ptrees}.
#' @param inclusion If "both", then only pairs in which both individuals are members of \code{hosts} are included. If "either" then pairs only need have one member from \code{hosts}
#' @param dwin Optional output of \code{classify.pairwise.relationships}. This can be specified to avoid double calculations.
#' @param control List of plotting attributes. 
#' @return A list whose elements are \code{data}, the underlying data frame for the graph, and \code{graph}, the graph itself.
#' @export produce.pairwise.graphs2
#' @seealso \code{\link{classify.pairwise.relationships}}
#' @examples
#' #
#' # Example on data from Rakai Community Cohort Study
#' # remember that you can specify dwin to save computing time, if you have it already computed
#' #
#' \dontrun{	
#' file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.RData'),package='phyloscannerR')
#' load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'	
#' hosts	<- c('RkA05868F','RkA05968M','RkA00369F','RkA01344M')
#' inclusion <- "both" 
#' tmp	<- produce.pairwise.graphs2(phsc, hosts=hosts, inclusion = "both")
#' tmp$graph
#' }
produce.pairwise.graphs2 <- function(	ptrees, 		
										hosts = all.hosts.from.trees(ptrees),
										inclusion = c("both", "either"),
										dwin= NULL,
										control = list(	yintercept_close=0.025,
														yintercept_dist=1.0,
														breaks_x=seq(0,1e4,500), 
														minor_breaks_x=seq(0,1e4,100),
														breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
														limits_y=c(1e-3,0.4),
														fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))
										)
{
	if(is.null(dwin))
	{
		#get basic.classification
		dwin <- classify.pairwise.relationships(ptrees, 				 
					close.threshold=0.025, 	#has no effect on basic.classification
					distant.threshold=0.05,	#has no effect on basic.classification		
					relationship.types=c(), 
					verbose=TRUE)		
	}
	stopifnot( dwin%>%filter(host.2<host.1)%>%nrow() == 0)
	#	define basic topology
	dwin <- dwin %>%
			mutate(basic.topology:= gsub('desc','descendant',gsub('anc','ancestral',gsub('^(.*)_(.*)$','\\1',basic.classification))))
	#	extract window coordinates
	dwin <- dwin %>% 
			mutate(window.start:= map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', x))))
	dwin <- dwin %>% 
			mutate(window.end:= map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', x))))	
	#	check colour names
	if( !all( unique(dwin$basic.topology) %in% names(control$fill.topology) ) )
		stop('Not all basic topology types have a colour defined. You need to define colours for',unique(dwin$basic.topology))
	#	check inclusion
	if( !inclusion %in% c("both","either"))
		stop('Cannot handle value of inclusion', inclusion)
	
	if(inclusion=='both')
	{
		pair.data <- as_tibble(expand.grid(host.1=hosts, host.2=hosts, stringsAsFactors=FALSE)) %>%
				filter(host.1!=host.2)	%>%
				inner_join(dwin, by=c('host.1','host.2')) %>%
				select(host.1, host.2, window.start, basic.topology, patristic.distance) %>%
				arrange(host.1, host.2, window.start)
	}
	if(inclusion=='either')
	{
		pair.data <- dwin %>% 
				select(host.1, host.2) %>%
				distinct() %>%
				filter(host.1 %in% hosts | host.2 %in% hosts) %>%
				inner_join(dwin, by=c('host.1','host.2')) %>%
				select(host.1, host.2, window.start, basic.topology, patristic.distance) %>%
				arrange(host.1, host.2, window.start)						
	}
	# set min patristic distance for plotting, make factors, etc
	pair.data <- pair.data %>%
			mutate(	patristic.distance:= pmax(1e-3,patristic.distance),
			                basic.topology:= factor(basic.topology),
					panel:= paste0(host.1,'\n->\n',host.2),
					y:= 1e-3)			
	# make plot
	if(!'limits_x'%in%names(control))
		control$limits_x <- c(min(pair.data$window.start) - diff(pair.data$window.start[1:2]), max(pair.data$window.start) + diff(pair.data$window.start[1:2]))
	pairwise.plot <- ggplot(pair.data, aes(window.start)) +			
			geom_bar(aes(y=y, fill=basic.topology), colour='transparent', stat='identity', width=25) +
			geom_point(aes(y=patristic.distance), size=1) +
			geom_hline(yintercept=control$yintercept_close, colour='black', lty=2) +
			geom_hline(yintercept=control$yintercept_dist, colour='black', lty=2) +
			labs(x='\ngenomic position\n(relative to HXB2)', y='subgraph distance\n(subst/site)\n',fill='topological subgraph\nrelationship') +
			scale_x_continuous(breaks=control$breaks_x, minor_breaks=control$minor_breaks_x, limits=control$limits_x) +
			scale_y_log10(labels=scales:::percent, expand=c(0,0), breaks=control$breaks_y) +
			coord_cartesian(ylim=control$ylim) +
			scale_fill_manual(values=control$fill.topology) +
			theme_bw() + 
			theme(legend.position='bottom', panel.spacing = unit(1, "lines")) +
			facet_grid(panel~.)
	return(list(graph = pairwise.plot, data = pair.data))
}

#' Draw bar graphs of pairwise topological/distance relationships
#' @param ptrees A list of class \code{phyloscanner.trees}
#' @param dist.thresh The distance threshold used to select likely transmission pairs
#' @param hosts A list of hosts (as a vector) to obtain graphs for. By default, all pairs of hosts detected in \code{ptrees}.
#' @param inclusion If "both", then only pairs in which both individuals are members of \code{hosts} are included. If "either" then pairs only need have one member from \code{hosts}
#' @param contiguous.pairs If TRUE pairs require contiguous (rather than ajacent) subgraphs to be identified as likely transmissions
#' @return A list whose elements are \code{data}, the underlying data frame for the graph, and \code{graph}, the graph itself.
#' @export produce.pairwise.graphs
#' @examples
#' #
#' # Example on data from Rakai Community Cohort Study
#' #
#' \dontrun{	
#' file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.RData'),package='phyloscannerR')
#' load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'	
#' hosts	<- c('RkA05868F','RkA05968M','RkA00369F','RkA01344M')
#' inclusion <- "both" 
#' tmp	<- produce.pairwise.graphs(phsc, hosts=hosts, inclusion = "both")
#' tmp$graph
#' }
produce.pairwise.graphs <- function(ptrees, 
                                    dist.thresh= 0.025,
                                    hosts = all.hosts.from.trees(ptrees),
                                    contiguous.pairs = F,
                                    inclusion = c("both", "either")){
  
  full.host.list <- all.hosts.from.trees(ptrees)
  
  t.stats <- ptrees %>% map(function(ptree){
    out <- ptree$classification.results$classification
    out <- out %>% 
      mutate(tree.id = ptree$id, xcoord = ptree$xcoord)
    out
  }) %>% bind_rows()
  
  if(!("normalised.min.distance.between.subgraphs" %in% colnames(t.stats))){
    t.stats <- t.stats %>% mutate(normalised.min.distance.between.subgraphs = min.distance.between.subgraphs)
  }
  
  wrong.dir <- which(t.stats$host.2 < t.stats$host.1)
  t.stats <- t.stats %>% 
    mutate(ancestry = pmap_chr(list(host.1, host.2, ancestry), function(h1,h2,a){
      if(h1 <= h2){
        a
      } else if(!(a %in% c("anc", "desc", "multiAnc", "multiDesc"))){
        
        a
      } else {
        switch(a, 
               "anc" = "desc",
               "desc" = "anc",
               "multiAnc" = "multiDesc",
               "multiDesc" = "multiAnc")
      }
      
    })) %>%
    select(-paths12, -paths21, -nodes1, -nodes2)
  
  t.stats[wrong.dir,c(1,2)] <- t.stats[wrong.dir,c(2,1)]
  
  hosts.by.tree <- map(1:length(ptrees), function(x){
    all.hosts.from.trees(ptrees[x])
  })
  
  names(hosts.by.tree) <- map_chr(ptrees, function(x) x$id)
  
  if(attr(ptrees, "readable.coords")){
    xcoord.lookup <- tibble(id = unique(map_chr(ptrees, function(x) x$id))) %>%
      mutate(xcoord = map_dbl(id, function(x) mean(as.numeric(unlist(strsplit(x, "_to_"))))))
  } else {
    xcoord.lookup <- tibble(id = unique(map_chr(ptrees, function(x) x$id))) %>%
      mutate(xcoord = 1:length(id))
  }

  pair.data <- as_tibble(expand.grid(full.host.list, full.host.list, map_chr(ptrees, function(x) x$id), stringsAsFactors = F)) %>%
    dplyr::rename(host.1 = Var1, host.2 = Var2, tree.id = Var3) %>%
    filter(host.1 < host.2) %>% 
    left_join(t.stats, c("host.1", "host.2", "tree.id")) %>%
    mutate(not.present = is.na(ancestry)) %>%
    left_join(xcoord.lookup, by=c("tree.id" = "id")) %>%
    mutate(xcoord = xcoord.y) %>%
    select(-xcoord.x, -xcoord.y) %>%
    mutate(host.1.present = pmap_lgl(list(tree.id, host.1, not.present), function(x, y, z){
      if(!z){
        TRUE 
      } else {
        y %in% hosts.by.tree[[x]]
      }
    })) %>% mutate(host.2.present = pmap_lgl(list(tree.id, host.2, not.present), function(x, y, z){
      if(!z){
        TRUE 
      } else {
        y %in% hosts.by.tree[[x]]
      }
    })) %>%
    mutate(paircombo = paste0(host.1, " vs " ,host.2)) %>%
    group_by(paircombo) %>%
    mutate(within.distance = normalised.min.distance.between.subgraphs <= dist.thresh) %>%
    mutate(any.link = any(!not.present & adjacent & within.distance)) %>%
    ungroup() %>%
    mutate(nondir.ancestry = replace(ancestry, ancestry %in% c("anc", "desc"), "trans")) %>%
    mutate(nondir.ancestry = replace(nondir.ancestry, nondir.ancestry %in% c("multiAnc", "multiDesc"), "multiTrans")) %>%
    mutate(nondir.ancestry = factor(nondir.ancestry, levels = c("trans", "multiTrans", "noAncestry", "complex"))) %>%
    mutate(linked = ((!contiguous.pairs & adjacent) | (contiguous.pairs & contiguous)) & within.distance) %>%
    mutate(linked = map_chr(linked, function(x)  ifelse(is.na(x), "Member absent", ifelse(x, "Yes", "No")))) %>%
    mutate(linked = factor(linked, levels = c("Yes", "No", "Member absent"))) %>%
    mutate(adjacent = replace(adjacent, is.na(adjacent), T)) %>%
    mutate(adjacent = map_chr(adjacent, function(x) ifelse(x, "Yes", "No"))) %>%
    mutate(adjacent = factor(adjacent, levels = c("Yes", "No"))) %>%
    mutate(contiguous = replace(contiguous, is.na(contiguous), T)) %>%
    mutate(contiguous = map_chr(contiguous, function(x) ifelse(x, "Yes", "No"))) %>%
    mutate(contiguous = factor(contiguous, levels = c("Yes", "No"))) %>%
    mutate(fact.within.distance = as.character(within.distance)) %>%
    mutate(fact.within.distance = replace(fact.within.distance, fact.within.distance=="TRUE", "Within threshold")) %>%
    mutate(fact.within.distance = replace(fact.within.distance, fact.within.distance=="FALSE", "Outside threshold")) %>%
    mutate(fact.within.distance = factor(fact.within.distance, levels = c("Within threshold", "Outside threshold"))) %>%
    mutate(arrow = ancestry) %>%
    mutate(arrow = replace(arrow, arrow %in% c("anc", "multiAnc"), 1)) %>%
    mutate(arrow = replace(arrow, arrow %in% c("desc", "multiDesc"), 2)) %>%
    mutate(arrow = replace(arrow, arrow %in% c(NA, "noAncestry", "complex"), 3)) %>%
    mutate(arrow = as.numeric(arrow)) %>%
    mutate(arrow.start = arrow) %>%
    mutate(arrow.start = replace(arrow.start, arrow.start == 2, 1.20001)) %>%
    mutate(arrow.start = replace(arrow.start, arrow.start == 1, 1.79999)) %>%
    mutate(arrow.start = replace(arrow.start, arrow.start == 3, NA)) %>%
    mutate(arrow.end = arrow) %>%
    mutate(arrow.end = replace(arrow.end, arrow.end == 2, 1.2)) %>%
    mutate(arrow.end = replace(arrow.end, arrow.end == 1, 1.8)) %>%
    mutate(arrow.end = replace(arrow.end, arrow.end == 3, NA)) %>%
    mutate(arrow.start = as.numeric(arrow.start), arrow.end = as.numeric(arrow.end)) 
  
  if(inclusion == "both"){
    pair.data <- pair.data %>%
      filter((host.1 %in% hosts & host.2 %in% hosts))
  } else {
    pair.data <- pair.data %>%
      filter((host.1 %in% hosts | host.2 %in% hosts))
  }
  
  if(nrow(pair.data)==0){
    stop("No pairs match this vector of hosts")
  }
  
  pairwise.plot <- ggplot(pair.data) +
    geom_point(aes(y=host.1, x = xcoord, alpha = as.numeric(host.1.present), fill=linked), col="black", size=2, shape=21) +
    geom_point(aes(y=host.2, x = xcoord, alpha = as.numeric(host.2.present), fill=linked), col="black", size=2, shape=21) +
    geom_segment(data = pair.data %>% filter(!not.present), 
                 aes(x=xcoord, xend = xcoord, y = 1.8, yend = 1.2, col = nondir.ancestry, linetype=fact.within.distance), size=0.75) +
    geom_segment(data = pair.data %>% filter(!not.present), 
                 aes(x=xcoord, xend = xcoord, y = arrow.start, yend = arrow.end, col = nondir.ancestry), 
                 size=0.66, 
                 arrow = arrow(length = unit(0.2, "cm"), ends = "last"),
                 show.legend = FALSE, na.rm=TRUE) +
    geom_point(data = pair.data %>% filter(!host.1.present), aes(y=host.1, x = xcoord), col="black", shape = 4, size=2) +
    geom_point(data = pair.data %>% filter(!host.2.present), aes(y=host.2, x = xcoord), col="black", shape = 4, size=2) +
    scale_fill_manual(drop = F, name="Linked", labels=c("Yes", "No", "Member\nabsent"), values = c("#377eb8", "#e41a1c", "#bdbdbd")) +
    scale_linetype_manual(drop=F, name="Distance\nthreshold", values = c("solid", "dashed"), labels = c("Within", "Outside")) +
    scale_colour_manual(drop=FALSE, values = c("blue3", "green4", "darkorange2", "orange4"), name="Ancestry", labels = c("Single", "Multiple", "None", "Complex")) +
    facet_wrap(~paircombo, ncol = 1, scales="free_y") + 
    theme_minimal() + 
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          axis.title.y=element_blank()) +
    guides(size = guide_legend(override.aes = list(shape = 32))) +
    guides(col = guide_legend(override.aes = list(shape = 32))) +
    xlab("Window centre") +
    ylab("Host")
  
  if(contiguous.pairs & any(pair.data$contiguous == "No")){
    pairwise.plot <- pairwise.plot +
      geom_point(aes(x = xcoord, shape = contiguous=="No"), y=1.5, col="red3", size=5) +
      scale_shape_manual(drop=F, values = c(32, 126), name="Contiguous", labels = c("Yes", "Blocked")) 
  } else if(!contiguous.pairs & any(pair.data$adjacent == "No")) {
    pairwise.plot <- pairwise.plot +
      geom_point(aes(x = xcoord, shape = adjacent=="No"), y=1.5, col="red3", size=5) +
      scale_shape_manual(drop=F, values = c(32, 126), name="Adjacent", labels = c("Yes", "Blocked")) 
  }
  
  if(any(!pair.data$host.1.present) | any(!pair.data$host.2.present)){
    pairwise.plot <- pairwise.plot + 
      scale_alpha_continuous(range = c(0.33, 1), guide=FALSE) 
  } else {
    pairwise.plot <- pairwise.plot + 
      scale_alpha_continuous(range = c(1, 1), guide=FALSE) 
  }
  
  return(list(graph = pairwise.plot, data = pair.data))
}
