# Align the graph output

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

add.no.data.rectangles <- function(graph, rectangle.coords, log = F, y.limits = NULL){
  
  if(is.null(y.limits)){
    y.limits <- ggplot_build(graph)$layout$panel_ranges[[1]]$y.range
  }
  
  if(nrow(rectangle.coords)>0){
    for(rect.no in seq(1, nrow(rectangle.coords))){
      if(log){
        graph <- graph + annotate("rect", xmin=rectangle.coords$start[rect.no], xmax = rectangle.coords$end[rect.no], 
                                  ymin=10^(y.limits[1]), ymax=10^(y.limits[2]), fill=rectangle.coords$colour[rect.no], alpha=0.5)
      } else {
        graph <- graph + annotate("rect", xmin=rectangle.coords$start[rect.no], xmax = rectangle.coords$end[rect.no], 
                                  ymin=y.limits[1], ymax=y.limits[2], fill=rectangle.coords$colour[rect.no], alpha=0.5)
      }
    }
  }
  
  return(graph)
}

# Calculate the coordinates to build the rectangles for. 
# missing.coords are x-coordinates of windows with no coverage. 
# all.coords is the  full list of window coordinates.

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
    
    return(data.frame(starts = final.starts, ends = final.ends, colour=colour))
  } 
  stop("No missing coordinates given")
}

# Try to detect if the windows are evenly spaced. If they aren't, skip the rectangle-drawing, it's too complex.

find.gaps <- function(xcoords){
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
      regular.gaps <- F
      new.xcoords <- xcoords
      missing.window.rects <- NULL
    }
  }
  
  range <- max(new.xcoords) - min(new.xcoords)
  bar.width <- range/(1.5*(length(new.xcoords)+1))
  
  return(list(x.coordinates = new.xcoords, regular.gaps = regular.gaps, rectangles.for.missing.windows = missing.window.rects, bar.width=bar.width))
  
}

produce.pdf.graphs <- function(file.name, host.statistics, hosts, xcoords, missing.window.rects, bar.width, regular.gaps = F, width=8.26772, height=11.6929, verbose = F){

  pdf(file=file.name, width=width, height=height)
  
  for (i in seq(1, length(hosts))) {
    
    host <- hosts[i]
  
    this.host.statistics <- host.statistics[which(host.statistics$id==host),]
    this.host.statistics <- this.host.statistics[order(this.host.statistics$xcoord),]

    this.host.statistics <- this.host.statistics[which(this.host.statistics$reads>0),]
    
    if(length(which(this.host.statistics$reads>0)) > 0){
      if (verbose) cat("Drawing graphs for host ",host,"\n", sep="")
      
      # The coordinates for windows missing from all patients are known. These are those missing for this patient, but not all patients.
 
      missing.read.rects <- NULL
      
      if(regular.gaps & length(which(this.host.statistics$reads==0))){
        missing.read.rects <- form.rectangles(host.statistics$xcoord[which(host.statistics$reads==0)], xcoords, "grey")
      }
      
      # rbind on two NULLs still makes a NULL
      
      missing.rects <- rbind(missing.window.rects, missing.read.rects)
      
      this.host.statistics <- this.host.statistics[which(this.host.statistics$reads>0),]
      
      #graph 1: read count and tip count
      
      this.host.statistics.temp <- melt(this.host.statistics[,c("xcoord","id","tips","reads")], id.vars=c("id", "xcoord"))
      
      max.value <- max(this.host.statistics.temp$value)
      log.upper.tick <- ceiling(log10(max.value))
      ticks <- 10^(seq(0, log.upper.tick))
      
      graph.1 <- ggplot(this.host.statistics.temp, aes(x=xcoord, y=value, col=variable))
      
      graph.1 <- graph.1 + geom_point(na.rm=TRUE) +
        theme_bw() + 
        scale_y_log10(breaks=ticks, limits=c(0.8, 1.2*(10^log.upper.tick))) +
        ylab("Count") +
        xlab("Window centre") +
        scale_x_continuous(limits=c(ews, lwe)) +
        scale_color_discrete(name="Variable", labels=c("Tips", "Reads")) + 
        theme(text = element_text(size=7))
      
      if(regular.gaps & !is.null(missing.rects)){
        graph.1 <- add.no.data.rectangles(graph.1, missing.rects,  TRUE, c(log10(0.8), log10(1.2) + log.upper.tick))
      }
      
      #graph 2: subgraph and clade counts
      
      this.host.statistics.temp <- melt(this.host.statistics[,c("xcoord", "id", "subgraphs", "clades")], id.vars=c("id", "xcoord"))
      
      graph.2 <- ggplot(this.host.statistics.temp, aes(x=xcoord, y=value))
      
      graph.2 <- graph.2 +
        geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
        aes(col = variable) +
        theme_bw() + 
        ylab("Count") +
        xlab("Window centre") +
        scale_x_continuous(limits=c(ews, lwe)) +
        scale_shape_manual(values=c(1,19), name="Variable", labels=c("Subgraphs", "Clades")) +  
        scale_size_manual(values=c(2,1), name="Variable", labels=c("Subgraphs", "Clades")) +		
        scale_color_discrete(name="Variable", labels=c("Subgraphs", "Clades")) + 
        theme(text = element_text(size=7))
      
      if(max(this.host.statistics.temp$value)==1){
        graph.2 <- graph.2 + scale_y_continuous(breaks = c(0,1))
      } else if(max(this.host.statistics.temp$value)==2){
        graph.2 <- graph.2 + expand_limits(y=0)+ scale_y_continuous(breaks = c(0,1,2)) 
      } else {
        graph.2 <- graph.2 + expand_limits(y=0) + scale_y_continuous(breaks = pretty_breaks()) 
      }
      
      if(regular.gaps & !is.null(missing.rects)){
        graph.2 <- add.no.data.rectangles(graph.2, missing.rects)
      }
      
      # graph 3: root-to-tip distances for all reads and largest subgraph
      
      this.host.statistics.temp <- melt(this.host.statistics[,c("xcoord","id","overall.rtt","largest.rtt")], id.vars=c("id", "xcoord"))
      
      graph.3 <- ggplot(this.host.statistics.temp, aes(x=xcoord, y=value))
      
      graph.3 <- graph.3 +
        geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
        aes(col = variable) +
        theme_bw() + 
        ylab("Mean root-to-tip-distance\n(read-weighted)") +
        xlab("Window centre") +
        scale_x_continuous(limits=c(ews, lwe)) +
        expand_limits(y=0) + 
        scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph")) + 
        scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subgraph")) +
        scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subgraph")) +
        theme(text = element_text(size=7))
      
      if(regular.gaps & !is.null(missing.rects)){
        graph.3 <- add.no.data.rectangles(graph.3, missing.rects)
      }
      
      # graph 4: largest patristic distance in largest subgraph
      
      this.host.statistics.temp <- melt(this.host.statistics[,c("xcoord","id","global.mean.pat.distance","subgraph.mean.pat.distance")], id.vars=c("id", "xcoord"))
      
      graph.4 <- ggplot(this.host.statistics.temp, aes(x=xcoord, y=value))
      
      graph.4 <- graph.4 +
        geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
        aes(col = variable) +
        theme_bw() + 
        ylab("Mean pairwise patristic distance \n(read-weighted)") +
        xlab("Window centre") +
        scale_x_continuous(limits=c(ews, lwe)) +
        expand_limits(y=0) + 
        scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph")) + 
        scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subgraph")) +
        scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subgraph")) +
        theme(text = element_text(size=7))
      
      if(regular.gaps & !is.null(missing.rects)){
        graph.4 <- add.no.data.rectangles(graph.4, missing.rects)
      }
      
      # graph 5: read proportions in each subgraph
      
      proportion.column.names <- colnames(host.statistics)[which(substr(colnames(host.statistics), 1, 8)=="prop.gp.")]
      
      # only columns with nonzero entries should be kept
      proportion.columns <- host.statistics[which(pat.stats$id==host), proportion.column.names, with=F]
      
      sums <- colSums(proportion.columns, na.rm=T)
      proportion.column.names <- proportion.column.names[which(sums>0)]
      
      splits.props <- pat.stats[which(host.statistics$id==host), c("id", "xcoord", proportion.column.names), with=F]
      
      splits.props.1col <- melt(splits.props, id=c("id", "xcoord"))
      splits.props.1col <- splits.props.1col[!is.na(splits.props.1col$value),]
      
      splits.props.1col$ngroup <- sapply(splits.props.1col$variable, function(x) which(levels(splits.props.1col$variable)==x))
      
      splits.props.1col$fgroup <- as.factor(splits.props.1col$ngroup)
      
      colourCount = length(unique(splits.props.1col$ngroup))
      getPalette = colorRampPalette(brewer.pal(9, "Greens"))
      
      graph.5 <- ggplot(splits.props.1col, aes(x=xcoord, weight=value, fill=reorder(fgroup, rev(order(splits.props.1col$ngroup)))))
      
      graph.5 <- graph.5 +
        geom_bar(width=bar.width, colour="black", size=0.25) +
        theme_bw() + 
        ylab("Proportion of reads\nin discrete subraphs") +
        xlab("Window centre") +
        scale_x_continuous(limits=c(ews, lwe)) +
        scale_fill_manual(values = getPalette(colourCount)) +
        theme(text = element_text(size=7)) + 
        guides(fill = guide_legend(title = "Subgraph rank\n(by tip count)", keywidth = 1, keyheight = 0.4))
      
      if(regular.gaps & !is.null(missing.rects)){
        graph.5 <- add.no.data.rectangles(graph.5, missing.rects)
      }
      
      # graph 6: longest branch to largest patristic distance ratios
      
      if("recombination.metric" %in% colnames(this.host.statistics)) {     
        graph.6 <- ggplot(this.host.statistics, aes(x=xcoord, y=recombination.metric))
        y.label <- "Recombination metric"
        
        graph.6 <- graph.6 +
          geom_point(alpha = 0.5, na.rm=TRUE) +
          theme_bw() + 
          ylab(y.label) +
          xlab("Window centre") +
          scale_x_continuous(limits=c(ews, lwe)) +
          expand_limits(y=0) +
          #      scale_color_discrete(name="Tip set", labels=c("Longest branch", "Greatest patristic distance")) + 
          theme(text = element_text(size=7))
        
        if(regular.gaps & !is.null(missing.rects)){
          graph.6 <- add.no.data.rectangles(graph.6, missing.rects)
        }
        
        all.plots <- AlignPlots(graph.1, graph.2, graph.3, graph.4, graph.5, graph.6)
        
        if(i!=1){
          grid.newpage()
        }
        
        pushViewport(viewport(layout = grid.layout(7, 1, heights = unit(c(0.25, rep(1,6)), "null") )))
        grid.text(host, gp=gpar(fontsize=20), vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
        for(plot.no in 1:6){
          plot <- all.plots[[plot.no]]
          plot$vp = viewport(layout.pos.col = 1, layout.pos.row = plot.no + 1)
          grid.draw(plot)
        }
        
      } else {
        
        all.plots <- AlignPlots(graph.1, graph.2, graph.3, graph.4, graph.5)
        
        if(i!=1){
          grid.newpage()
        }
        
        pushViewport(viewport(layout = grid.layout(6, 1, heights = unit(c(0.25, rep(1,5)), "null") )))
        grid.text(host, gp=gpar(fontsize=20), vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
        for(plot.no in 1:5){
          plot <- all.plots[[plot.no]]
          plot$vp = viewport(layout.pos.col = 1, layout.pos.row = plot.no + 1)
          grid.draw(plot)
        }
      } 
    } else {
      if (verbose) cat("Skipping graphs for host ",host," as no reads are present and not blacklisted\n", sep="")
    }
  }
  
  never.mind <- dev.off()

}
