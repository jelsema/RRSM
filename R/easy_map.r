#' A convenient mapping function
#' 
#' @description
#' Function that produces a simple \code{ggplot}.
#' 
#' @param data input dataset containing coordinates (columns 1 and 2) and response data (column three).
#' @param knots optional matrix of locations for knots.
#' @param show.knots logical, whether the knots should be plotted.
#' @param panel.plot logical, determines whether to produce a single figure or a paneled plot.
#' @param panel1 optional, grouping factor for first paneling level.
#' @param panel2 optional, grouping factor for second paneling level.
#' @param cex optional, the size of the plotted points.
#' @param cex.knot optional, the size of the knot locations.
#' @param xlim optional, the limits of the x-coordinate.
#' @param ylim optional, the limits of the y-coordinate.
#' @param x.ticks optional, the placement of tick marks on the the x-axis.
#' @param y.ticks optional, the placement of tick marks for the the y-axis.
#' @param alpha the alpha scaling factor for geom_point.
#' @param col.map the colormap for the plot, options are "ylgn" (yellow-green) or "terrain.colors".
#' @param ... space for additional arguments.
#' 
#' @note
#' This is a glorified wrapper for ggplot. It is unlikely to be useful in all circumstances, 
#' but often provides a reasonable map, especially for simulated examples using the simulation
#' tools in \pkg{RRSM}, or other datasets on an approximate square domain.
#' 
#' @return
#' A figure produced with ggplot.
#' 
#' @import ggplot2
#' 
#' @export
#' 
#' 

easy_map <- function( data , knots=NULL , show.knots=TRUE , panel.plot=FALSE , 
                      panel1=NULL , panel2=NULL , cex=NULL , cex.knot=NULL ,
                      xlim=NULL , ylim=xlim , x.ticks=NULL , y.ticks=x.ticks ,
                      alpha=NULL , col.map="ylgn" , ... ){
  
  ## Prevent some errors
  Xdim <- Ydim <- Y <- NULL
  
  if( length(col.map)==1 ){
    if( col.map=="ylgn" ){
      col.map2 <- c("#006837","#238443","#41AB5D","#78C679","#ADDD8E","#D9F0A3","#F7FCB9")
    }
    if( col.map=="terrain.colors" ){
      col.map2 <- terrain.colors(100)
    }
    
    col.map <- col.map2
  }
  
  
  
  # Define some values
  if( is.null(x.ticks) ){
    minx <- round(min(data[,1]))
    maxx <- round(max(data[,1]))
    x.ticks <- round(seq( minx , maxx , length.out=6 ) , 2)
  }
  if( is.null(y.ticks) ){
    miny <- round(min(data[,2]))
    maxy <- round(max(data[,2]))
    y.ticks <- round(seq( miny , maxy , length.out=6 ) , 2)
  }
  if( is.null(xlim) ){
    minx1 <- floor(  min(data[,1]))
    maxx1 <- ceiling(max(data[,1]))
    xlim <- c( minx1 , maxx1 )
  }
  if( is.null(ylim) ){
    miny1 <- floor(  min(data[,2]))
    maxy1 <- ceiling(max(data[,2]))
    ylim <- c( miny1 , maxy1 )
  }
  
  
  
  
  if( is.null(cex) ){
    nd <- dim(data)[1]
    if( 0     < nd & nd <= 1500  ){ cex <- 4}
    if( 1500  < nd & nd <= 10000 ){ cex <- 3 }
    if( 10000 < nd & nd <= 30000 ){ cex <- 2.3 }
    if( 30000 < nd & nd <= 50000 ){ cex <- 2   }
    if( 50000 < nd               ){ cex <- 1 }
  }
  
  if( is.null(panel1)==FALSE ){ panel.plot <- TRUE }
  
  if( is.null(alpha) ){ alpha <- rep(1 , dim(data)[1] ) }
  
  
  if( is.null(cex.knot) ){ cex.knot  <- 2.2 }
  
  ## Only make one plot, not a panel
  if( panel.plot==F ){
    
    ## Don't show the knots
    if( is.null(knots) || show.knots==F ){
      plot.dat <- as.data.frame( data[,1:3] )
      colnames(plot.dat) <- c("Xdim" , "Ydim" , "Y")
    
      g.plot <- ggplot( ) + 
        geom_point( data=plot.dat , aes( x=Xdim , y=Ydim , colour=Y) , size=cex , alpha=alpha ) +
        scale_colour_gradientn(colours = col.map )      + 
        scale_y_continuous( breaks=y.ticks , limits=ylim) +
        scale_x_continuous( breaks=x.ticks , limits=xlim)
    }
    
    ## Show the knots as well
    if( is.null(knots)==F & show.knots==T ){
      plot.dat <- as.data.frame( data[,1:3] )
      colnames(plot.dat) <- c("Xdim" , "Ydim" , "Y")
      plot.knot <- as.data.frame( knots )
      colnames(plot.knot) <- c("Xdim" , "Ydim")
      
      g.plot <- ggplot( ) + 
        geom_point( data=plot.dat , aes( x=Xdim , y=Ydim , colour=Y) , size=cex , alpha=alpha) +
        scale_colour_gradientn(colours = col.map)      + 
        geom_point( data=plot.knot , aes( x=Xdim , y=Ydim ) , size=cex.knot) +
        scale_y_continuous( breaks=y.ticks , limits=ylim) +
        scale_x_continuous( breaks=x.ticks , limits=xlim)
    }
  }
  
  ## Make a paneled plot
  if( panel.plot==T ){
    
    ## Don't show the knots
    if( is.null(knots) || show.knots==F ){
      plot.dat <- as.data.frame( data )
      #colnames(plot.dat) <- c("Xdim" , "Ydim" , "Y")
      
      ## Only one paneling variable
      if( is.null(panel2) ){
        
        panel1 <- as.data.frame( panel1 )
        plot.dat <- as.data.frame( cbind(plot.dat,panel1) )
        colnames(plot.dat) <- c("Xdim" , "Ydim" , "Y" , "Panel")
        
        g.plot <- ggplot( ) + 
          facet_wrap( ~ Panel ) +
          geom_point( data=plot.dat , aes( x=Xdim , y=Ydim , colour=Y) , size=cex) +
          scale_colour_gradientn(colours = col.map)      + 
          scale_y_continuous( breaks=y.ticks , limits=ylim) +
          scale_x_continuous( breaks=x.ticks , limits=xlim)
      }
      
      ## TWO paneling variable (not yet tested)
      if( is.null(panel2)==F ){
        
        panel1 <- as.data.frame( panel1 )
        panel2 <- as.data.frame( panel2 )
        
        plot.dat <- as.data.frame( cbind(plot.dat,panel1,panel2) )
        colnames(plot.dat) <- c("Xdim" , "Ydim" , "Y" , "Panel1" , "Panel2")
        
        g.plot <- ggplot( ) + 
          facet_grid( panel1 ~ panel2 ) +
          geom_point( data=plot.dat , aes( x=Xdim , y=Ydim , colour=Y) , size=cex) +
          scale_colour_gradientn(colours = col.map)      + 
          scale_y_continuous( breaks=y.ticks , limits=ylim) +
          scale_x_continuous( breaks=x.ticks , limits=xlim)               
      }
      
      
      
    }
    
    ## Show the knots as well
    if( is.null(knots)==F & show.knots==T ){
      plot.dat <- as.data.frame( data[,1:3] )
      colnames(plot.dat) <- c("Xdim" , "Ydim" , "Y")
      plot.knot <- as.data.frame( knots )
      colnames(plot.knot) <- c("Xdim" , "Ydim")
      
      g.plot <- ggplot( ) + 
        facet_grid( panel1 ~ panel2 ) +
        geom_point( data=plot.dat , aes( x=Xdim , y=Ydim , colour=Y) , size=cex) +
        scale_colour_gradientn(colours = col.map)      + 
        geom_point( data=plot.knot , aes( x=Xdim , y=Ydim ) , size=cex.knot) +
        scale_y_continuous( breaks=y.ticks , limits=ylim) +
        scale_x_continuous( breaks=x.ticks , limits=xlim)
    }
  }
  
  
  
  return(g.plot)
  
  
}
