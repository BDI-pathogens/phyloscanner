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
      # I don't know what to do and I am giving up
      
      regular.gaps <- F
      new.xcoords <- xcoords
      missing.window.rects <- NULL
    }
  }
  
  range <- max(new.xcoords) - min(new.xcoords)
  bar.width <- range/(1.5*(length(new.xcoords)+1))
  
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
    
    if (verbose) cat("Drawing graphs for host ",host,"\n", sep="")
    
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
      scale_color_discrete(name="Variable", labels=c("Tips", "Reads")) + 
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
      scale_shape_manual(values=c(1,19), name="Variable", labels=c("Subgraphs", "Clades")) +  
      scale_size_manual(values=c(2,1), name="Variable", labels=c("Subgraphs", "Clades")) +		
      scale_color_discrete(name="Variable", labels=c("Subgraphs", "Clades")) + 
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
      gather(variable, value, overall.rtt, largest.rtt)
    
    graph.3 <- ggplot(host.stats.gd.3, aes(x=xcoord, y=value))
    
    graph.3 <- graph.3 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Mean root-to-tip-distance\n(read-weighted)") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph")) + 
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
      gather(variable, value, global.mean.pat.distance, subgraph.mean.pat.distance)
    
    graph.4 <- ggplot(host.stats.gd.4, aes(x=xcoord, y=value))
    
    graph.4 <- graph.4 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Mean pairwise patristic distance") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph")) + 
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
    
    colourCount <- length(unique(host.stats.gd.5$ngroup))
    
    # want largest subgraph to be the darkest colour even if all windows have 1 subgraph
    
    getPalette <- function(x){
      if(x>1){
        colorRampPalette(brewer.pal(5, "RdYlBu"))(x)
      } else {
        brewer.pal(5, "RdYlBu")[5]
      }
    }
    
    graph.5 <- ggplot(host.stats.gd.5, aes(x=xcoord, weight=value, fill=reorder(ngroup, rev(order(host.stats.gd.5$ngroup)))))
    
    graph.5 <- graph.5 +
      geom_bar(width=bar.width, colour="black", lty="blank") +
      theme_bw() + 
      ylab("Proportion of reads\nin different subraphs") +
      xlab(x.axis.label) +
      scale_x_continuous(limits=x.limits) +
      scale_fill_manual(values = getPalette(colourCount)) +
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
      y.label <- "Number of dual infections detected"
      
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