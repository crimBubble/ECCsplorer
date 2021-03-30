#!/usr/bin/python3

import os
import subprocess
import shlex
import socket
import atexit
import time
import pyRserve
from lib import config


def r_shutdown(port, logger):
    try:
        conn = pyRserve.connect(port=port)
        conn.shutdown()
        logger.info("Shutting down Rserve.")
    except pyRserve.rconn.RConnectionRefused:
        logger.error("connection to Rserve refused, server is probably already down")


def get_r_port():
    sckt = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sckt.bind(("", 0))
    sckt.listen(1)
    port = sckt.getsockname()[1]
    sckt.close()
    return port


def r_connection(logger):
    # Start R Server using Rserve

    logger.info('Starting Rserve...')
    config.PORT = get_r_port()
    os.system(config.RSERVE_CMD.format(port_=config.PORT))
    time.sleep(1)

    try:
        config.CONNECTOR = pyRserve.connect(host='localhost', port=config.PORT)
        logger.info('Rserv connected (Port {})'.format(config.PORT))
        config.CONNECTOR.voidEval(load_r_libs)
        logger.info('Rserv loading libraries.')
        atexit.register(r_shutdown, config.PORT, logger)
        return config.CONNECTOR
    except:
        error_msg = 'Rserve connection failed. Make sure Rserve is available in ../R/bin'
        logger.error(error_msg)
        raise Exception(error_msg)


# Image output formats
def r_image_format(orientation):

    if config.IMAGE_TYPE == 'pdf':
        r_image = 'pdf(file = png.path, ' \
                  'width = {width_}, height = {height_}, pointsize = {pointsize_})'\
            .format(type_=config.IMAGE_TYPE,
                    width_=(config.IMAGE_WIDTH / config.IMAGE_RES if orientation == 'landscape'
                            else config.IMAGE_HEIGHT / config.IMAGE_RES),
                    height_=config.IMAGE_HEIGHT / config.IMAGE_RES,
                    pointsize_=config.IMAGE_POINTS)

    else:
        r_image = '{type_}(filename = png.path, ' \
                  'width = {width_}, height = {height_}, units = "px", pointsize = {pointsize_}, res = {res_}, ' \
                  'bg = "white")' \
            .format(type_=config.IMAGE_TYPE,
                    width_=(config.IMAGE_WIDTH if orientation == 'landscape' else config.IMAGE_HEIGHT),
                    height_=config.IMAGE_HEIGHT,
                    pointsize_=config.IMAGE_POINTS,
                    res_=config.IMAGE_RES)

    return r_image


load_r_libs = ('\nsuppressPackageStartupMessages(library(ggplot2))'
               '\nsuppressPackageStartupMessages(library(ggrepel))'
               '\nsuppressPackageStartupMessages(library(grid))'
               '\nsuppressPackageStartupMessages(library(gridExtra))'
               '\nsuppressPackageStartupMessages(library(dplyr))').format()

Rconvert_genome = '''
convgenome <- function(coverage.path, bases.mapped, converted.coverage.path) {{
            coverage.data <- read.table(coverage.path, header=T)
        
            try(
            if(length(colnames(coverage.data)) < 5){{
                coverage.data$CO_map.all <- NA
            }})
            
            million.bases.factor <- setNames(data.frame(matrix(ncol = 2, nrow = 1)), c("TR", "CO"))
            million.bases.factor$TR <- 1000000/bases.mapped[1]
            million.bases.factor$CO <- 1000000/bases.mapped[2]
        
            coverage.data[4] <- ceiling(coverage.data[4]*million.bases.factor$TR)
            coverage.data[5] <- ceiling(coverage.data[5]*million.bases.factor$CO)
            
            write.table(coverage.data, file = converted.coverage.path, sep = "\\t", 
            col.names = T, row.names = F, quote = F)
            
            }}
'''.format()

Rconvert_region = '''
convregion <- function(coverage.path, bases.mapped, converted.coverage.path) {{
            coverage.data <- read.table(coverage.path, header=T)

            try(
            if(length(colnames(coverage.data)) < 7){{
                coverage.data$CO_map.all <- NA
                coverage.data$CO_map.SR <- NA
                coverage.data$CO_map.DR <- NA
            }})

            million.bases.factor <- setNames(data.frame(matrix(ncol = 2, nrow = 1)), c("TR", "CO"))
            million.bases.factor$TR <- 1000000/bases.mapped[1]
            million.bases.factor$CO <- 1000000/bases.mapped[2]

            coverage.data[4] <- ceiling(coverage.data[4]*million.bases.factor$TR)
            coverage.data[5] <- ceiling(coverage.data[5]*million.bases.factor$TR)
            coverage.data[6] <- ceiling(coverage.data[6]*million.bases.factor$TR)
            coverage.data[7] <- ceiling(coverage.data[7]*million.bases.factor$CO)
            coverage.data[8] <- ceiling(coverage.data[8]*million.bases.factor$CO)
            coverage.data[9] <- ceiling(coverage.data[9]*million.bases.factor$CO)

            write.table(coverage.data, file = converted.coverage.path, sep = "\\t", 
            col.names = T, row.names = F, quote = F)

            }}
'''.format()

Rconvert_cand = '''
convcand <- function(coverage.path, bases.mapped, converted.coverage.path) {{
            coverage.data <- read.table(coverage.path, header=T)

            try(
            if(length(colnames(coverage.data)) < 8){{
                coverage.data$CO_map.all <- NA
                coverage.data$CO_map.SR <- NA
                coverage.data$CO_map.DR <- NA
            }})

            million.bases.factor <- setNames(data.frame(matrix(ncol = 2, nrow = 1)), c("TR", "CO"))
            million.bases.factor$TR <- 1000000/bases.mapped[1]
            million.bases.factor$CO <- 1000000/bases.mapped[2]

            coverage.data[5] <- ceiling(coverage.data[5]*million.bases.factor$TR)
            coverage.data[6] <- ceiling(coverage.data[6]*million.bases.factor$TR)
            coverage.data[7] <- ceiling(coverage.data[7]*million.bases.factor$TR)
            coverage.data[8] <- ceiling(coverage.data[8]*million.bases.factor$CO)
            coverage.data[9] <- ceiling(coverage.data[9]*million.bases.factor$CO)
            coverage.data[10] <- ceiling(coverage.data[10]*million.bases.factor$CO)

            write.table(coverage.data, file = converted.coverage.path, sep = "\\t", 
            col.names = T, row.names = F, quote = F)

            }}
'''.format()

Renrichment = '''
enrichment <- function(converted.coverage.path, window.coverage.path) {{

            coverage.data <- read.table(converted.coverage.path, header=T)
            window.coverage.data <- read.table(window.coverage.path, header=T)
            
            avg.depth.all <- ifelse(is.na(coverage.data[7]), 
            ceiling(sum(window.coverage.data[4])/nrow(window.coverage.data)), 1) 
            
            coverage.data[7] <- ifelse(is.na(coverage.data[7]), avg.depth.all, coverage.data[7])                   
                                        
            coverage.data$enrich.all <- round((coverage.data[,4])/(coverage.data[,7]), digits = 2)
            
            avg.depth.control <- ifelse(is.na(coverage.data[7]), 
            ceiling(sum(window.coverage.data[5])/nrow(window.coverage.data)), 1)
            coverage.data$enrich.alt <- round((coverage.data[,4])/avg.depth.control, digits = 2)
                        
            coverage.data$enrich.all <- ifelse(is.infinite(coverage.data$enrich.all), 
            coverage.data$enrich.alt, coverage.data$enrich.all)
            
            coverage.data <- coverage.data[order(-coverage.data$enrich.all, coverage.data$chr, 
            coverage.data$start),]
            
            coverage.data <- coverage.data %>% mutate(id = sprintf("eccCand_%03d", row_number()))
            
            coverage.data <- coverage.data[order(coverage.data$chr, coverage.data$start),]

            write.table(coverage.data, file = converted.coverage.path, sep = "\\t", 
            col.names = T, row.names = F, quote = F)
            }}
'''.format()


def rmanhattan_plot():

    Rmanhattan_plot = '''
    manhattan <- function(converted.coverage.path, png.path, win.size, hiconf.path) {{
                
                print("Starting Manhattan plot.")

                coverage.data <- read.table(converted.coverage.path, header=T)
                
                coverage.data$hiconf <- read.table(hiconf.path, header=T, 
                colClasses = c("NULL", "NULL", "NULL", NA))$hiconf
                
                coverage.data$hiconf_map.all <- ifelse(coverage.data$hiconf >= 1, 
                coverage.data[,4], NA)
                
                y.up.lim <- max(coverage.data[4])*1.1
                y.dw.lim <- ifelse(is.na(coverage.data[1,5]), 0, max(coverage.data[5])*1.1)
                
                chr.colors <- c("#555555", "#111111")
                chr.number <- length(unique(coverage.data$chr))
                
                coverage.data.tmp <- coverage.data %>% group_by(chr) %>% summarise(chr_len=max(start)) %>%  
                mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(coverage.data, ., by=c("chr"="chr")) %>%
                arrange(chr, start) %>% mutate(poscum=start+tot)
                                            
                axis.label <- coverage.data.tmp %>% group_by(chr) %>% 
                summarize(center=( max(poscum) + min(poscum) ) / (2) )
                
                {image_}
                
                print(
                ggplot(coverage.data.tmp, aes(x=poscum)) +
                geom_point(y=-coverage.data.tmp[,5], alpha=0.8, size=0.5, color="grey") +
                geom_point(aes(y=coverage.data.tmp[,4], color=as.factor(chr)), alpha=1.0, size=0.5) +
                geom_point(y=coverage.data.tmp$hiconf_map.all, alpha=0.8, size=0.5, color="red") +
                
                scale_color_manual(values = rep(chr.colors, chr.number )) +
                scale_x_continuous(label=axis.label$chr, breaks = axis.label$center) +
                scale_y_continuous(expand = c(0, 0), limits = c(-y.dw.lim, y.up.lim)) +
                
                labs(x = paste0("Chromosome [", win.size, "bp windows]"), y = "Mean depth [BPM]") +
                theme_bw(base_size = 10) +
                theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()
                  )
                )
                
                graphics.off()
                print("Finished Manhattan plot.") 
                
                }}

    '''.format(image_=r_image_format(orientation='landscape'))

    #Rmanhattan_plot = '''manhattan <- function(a, png.path, c, d){{ {print_}({image_})}}'''.format(print_='print', image_=r_image_format(orientation='landscape'))

    return Rmanhattan_plot


def rline_multiplot():

    Rline_multiplot = '''
    linemulti <- function(candidate, converted.coverage.path, png.path) {{
                
                print(paste("Starting line plot:", candidate))

                coverage.data <- read.table(converted.coverage.path, header=T)
    
                y.up.lim.a <- ifelse(max(coverage.data[5]) > max(coverage.data[8]), 
                max(coverage.data[5])*1.1, max(coverage.data[8])*1.1)
                y.up.lim.c <- ifelse(max(coverage.data[6]) > max(coverage.data[9]), 
                max(coverage.data[6])*1.1, max(coverage.data[9])*1.1)
                y.up.lim.e <- ifelse(max(coverage.data[7]) > max(coverage.data[10]), 
                max(coverage.data[7])*1.1, max(coverage.data[10])*1.1)
                
                theme_set(theme_bw(base_size = 10) +
                theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()
                  ))
                
                plot.TR.all <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,5]), fill="#4D4D4D") +
                geom_line(aes(y=coverage.data[,5]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.a)) +
                labs(x = colnames(coverage.data[5]), y = NULL)
                
                plot.TR.SR <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,6]), fill="#CCCCCC") +
                geom_line(aes(y=coverage.data[,6]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.c)) +
                labs(x = colnames(coverage.data[6]), y = NULL)
                
                plot.TR.DR <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,7]), fill="#E5E5E5") +
                geom_line(aes(y=coverage.data[,7]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.e)) +
                labs(x = colnames(coverage.data[7]), y = NULL)
                
                plot.CO.all <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,8]), fill="#4D4D4D") +
                geom_line(aes(y=coverage.data[,8]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.a)) +
                labs(x = colnames(coverage.data[8]), y = NULL)
                
                plot.CO.SR <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,9]), fill="#CCCCCC") +
                geom_line(aes(y=coverage.data[,9]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.c)) +
                labs(x = colnames(coverage.data[9]), y = NULL)
                
                plot.CO.DR <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,10]), fill="#E5E5E5") +
                geom_line(aes(y=coverage.data[,10]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.e)) +
                labs(x = colnames(coverage.data[10]), y = NULL)
    
                {image_}
                
                grid.arrange(plot.TR.all, plot.CO.all,plot.TR.SR, plot.CO.SR, plot.TR.DR, plot.CO.DR,
                ncol=2, left="Coverage depth [BPM]")
                
                graphics.off()
                print(paste("Finished line plot:", candidate))

                }}
    '''.format(image_=r_image_format(orientation='landscape'))

    return Rline_multiplot


def rline_multiplot_noco():

    Rline_multiplot_noco = '''
    linemulti_noco <- function(candidate, converted.coverage.path, png.path) {{
                
                print(paste("Starting line plot:", candidate))

                coverage.data <- read.table(converted.coverage.path, header=T)
    
                y.up.lim.a <- max(coverage.data[5])*1.1
                y.up.lim.c <- max(coverage.data[6])*1.1
                y.up.lim.e <- max(coverage.data[7])*1.1
    
                theme_set(theme_bw(base_size = 10) +
                theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()
                  ))
    
                plot.TR.all <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,5]), fill="#4D4D4D") +
                geom_line(aes(y=coverage.data[,5]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.a)) +
                labs(x = colnames(coverage.data[5]), y = NULL)
    
                plot.TR.SR <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,6]), fill="#CCCCCC") +
                geom_line(aes(y=coverage.data[,6]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.c)) +
                labs(x = colnames(coverage.data[6]), y = NULL)
    
                plot.TR.DR <- ggplot(coverage.data, aes(x=pos)) +
                #geom_area(aes(y=coverage.data[,7]), fill="#E5E5E5") +
                geom_line(aes(y=coverage.data[,7]), size=0.25) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, y.up.lim.e)) +
                labs(x = colnames(coverage.data[7]), y = NULL)
    
                {image_}
    
                grid.arrange(plot.TR.all,plot.TR.SR, plot.TR.DR,
                ncol=1, left="Coverage depth [BPM]")
    
                graphics.off()
                print(paste("Finished line plot:", candidate))

                }}
    '''.format(image_=r_image_format(orientation='square'))

    return Rline_multiplot_noco


def rscatter_plot():

    Rscatter_plot = '''
    cl_scatter <- function(comp_cl_tab.path, png.path, threshold, pre.TR, pre.CO) {{

                print(paste("Starting clutering scatter plot."))

                cl.data <- read.table(comp_cl_tab.path, header=T, sep = "\\t")
                
                cl.names <- c("cl", "sup.cl", "size", "size.adj", "best.hit", "tarean.class", "sim.hits", 
                "size.TR", "size.CO", "prop.TR")
                names(cl.data) <- cl.names
                
                cl.data$cat <- ifelse(cl.data$prop.TR > threshold , "cl_hit", "no_hit")
                
                xlim <- max(log(cl.data$size.CO+1))
                ylim <- max(log(cl.data$size.TR+1))
                lim <- max(c(xlim, ylim)) + 1
                
                theme_set(theme_bw(base_size = 10) +
                  theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()
                  ))
                
                {image_}
                
                print(            
                ggplot(aes(x=log(cl.data$size.CO+1), y=log(cl.data$size.TR+1)), data=cl.data)+
    
                geom_abline(linetype = 2, color = c(alpha("grey", 0.5))) +
                geom_hline(yintercept=0, linetype = 1, color = c(alpha("grey", 0.5))) +
                geom_vline(xintercept=0, linetype = 1, color = c(alpha("grey", 0.5))) +
                  
                geom_point(aes(size = cl.data$size.adj, color = cl.data$cat)) +
                
                xlab(paste0("Abundance [", pre.CO, "]")) +
                ylab(paste0("Abundance [", pre.TR, "]")) +
                coord_cartesian(xlim = c(0,lim), ylim = c(0,lim)) +
                
                scale_colour_manual(values = c("no_hit" = c(alpha("#CCCCCC", 0.6)), 
                                               "cl_hit" = c(alpha("red", 0.6))
                                               ))
                )
                
                graphics.off()
                print(paste("Finished clutering scatter plot."))

                }}
    '''.format(image_=r_image_format(orientation='square'))

    return Rscatter_plot


def rcomparative_plot():

    Rcomparative_plot = '''
    comp_plot <- function(map.cand.path, png.path, cluster, cluster.blast.path, blast_outfmt, candidate) {{
                
                print(paste("Starting comparative plot:", candidate))

                # set values
                cluster.blast.end <- '.m6'
                
                # import data for mapper candidate
                cand <- read.table(map.cand.path, header = T, sep = "\\t")
                
                # calculate dynamic breaks, only positive values
                break.digits <- function(y) ((floor(log10(abs(max(y)))) + 1))
                round.digits <- c(2,1,0,-1,-2,-3,-4,-5,-6)
                names(round.digits) <- c("-1","0","1","2","3","4","5","6","7")
                break.3rd <- function(y) (round((max(y)/3), digits = round.digits[[paste0(break.digits(y/3))]]))
                y.breaks <- function(y) c(0, break.3rd(y), break.3rd(y)*2)

                # calculate upper graph limit
                y.up.lim <- max(cand[,5])*1.1
                white.space <- y.up.lim*0.008
                
                # dynamic cluster table with cluster numbers
                for(i in 1:length(cluster)){{
                  if(i == 1){{
                    cluster.tab <- read.table(paste0(cluster.blast.path, cluster[i], cluster.blast.end), header = F, sep = "\t")
                    names(cluster.tab) <- blast_outfmt
                    cluster.tab$cluster.nr <- i
                  }}
                  if(i != 1){{
                    cluster.new <- read.table(paste0(cluster.blast.path, cluster[i], cluster.blast.end), header = F, sep = "\t")
                    names(cluster.new) <- blast_outfmt
                    cluster.new$cluster.nr <- i
                    cluster.tab <- bind_rows(cluster.tab, cluster.new)
                  }}
                }}

                # rearrange sstart and send if in wrong order
                for (row in 1:nrow(cluster.tab)) {{
                  if(cluster.tab$sstart[row] > cluster.tab$send[row]){{
                  temp.end <- cluster.tab$send[row]
                  cluster.tab$send[row] <- cluster.tab$sstart[row]
                  cluster.tab$sstart[row] <- temp.end
                  }}
                }}

                # remove hits from other candidates (e.g. one cluster hits repetitive candidates)  
                cluster.tab <- cluster.tab[which(grepl(candidate, cluster.tab$sseqid)),]
                
                theme_set(theme_bw(base_size = 10) +
                theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()
                  ))
                
                {image_}
                
                print(
                ggplot() +
                
                geom_hline(yintercept=y.breaks(y.up.lim), linetype = 1, color = c(alpha("grey", 0.5))) +
                  
                #geom_area(aes(x=cand[,4], y=cand[,5]), fill=c(alpha("grey", 0.8))) +
                geom_line(aes(x=cand[,4], y=cand[,5]), color="black", size=0.25) +
                
                geom_rect(data=cluster.tab, mapping=aes(xmin=cluster.tab$sstart, xmax=cluster.tab$send, 
                                                        ymin=(-y.up.lim*0.08)*(cluster.tab$cluster.nr-1), 
                                                        ymax=((-y.up.lim*0.08)*(cluster.tab$cluster.nr))+white.space), 
                          fill="darkgrey", alpha=1) +
                geom_text(aes(x= -max(cand[,4])*0.05, y=(-y.up.lim*0.08)*(1:length(cluster)-1)-(y.up.lim*0.08/2), 
                    label=paste0("CL", cluster)), size=2.7) +
                scale_x_continuous(name="Position [bp]", position = "bottom") +
                scale_y_continuous(name="Cluster [blast hit] & mean depth [BPM]", breaks = y.breaks(y.up.lim), 
                    limits = c(-(y.up.lim*0.08)*length(cluster), y.up.lim))
                )
                
                graphics.off()
                print(paste("Finished comparative plot:", candidate))

                }}
    '''.format(image_=r_image_format(orientation='landscape'))

    return Rcomparative_plot


def rscatter_comp_plot():

    Rscatter_comp_plot = '''
    comp_scatter <- function(comp_cl_tab.path, png.path, threshold, pre.TR, pre.CO, cluster) {{
                
                print(paste("Starting comparative scatter plot."))

                cl.data <- read.table(comp_cl_tab.path, header=T, sep = "\\t")
                            
                cl.names <- c("cl", "sup.cl", "size", "size.adj", "best.hit", "tarean.class", "sim.hits", 
                "size.TR", "size.CO", "prop.TR")
                names(cl.data) <- cl.names
                
                cl.data$cat <- ifelse(cl.data$prop.TR > threshold & !(cl.data$cl %in% cluster), "cl_hit", 
                               ifelse(cl.data$prop.TR > threshold & cl.data$cl %in% cluster, "comp_hit", "no_hit"))
                
                xlim <- max(log(cl.data$size.CO+1))
                ylim <- max(log(cl.data$size.TR+1))
                lim <- max(c(xlim, ylim)) + 1
                
                theme_set(theme_bw(base_size = 10) +
                  theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()
                  ))
                
                {image_}
                
                print(            
                ggplot(aes(x=log(cl.data$size.CO+1), y=log(cl.data$size.TR+1)), data=cl.data)+
    
                geom_abline(linetype = 2, color = c(alpha("grey", 0.5))) +
                geom_hline(yintercept=0, linetype = 1, color = c(alpha("grey", 0.5))) +
                geom_vline(xintercept=0, linetype = 1, color = c(alpha("grey", 0.5))) +
                  
                #{{if(any(cl.data$cat != "comp_hit"))
                #geom_point(aes(x=ifelse(cl.data$cat != "comp_hit", log(cl.data$size.CO+1), NA), 
                #               y=ifelse(cl.data$cat != "comp_hit", log(cl.data$size.TR+1), NA), 
                #               size = cl.data$size.adj, color = cl.data$cat))}} +
                
                #{{if(any(cl.data$cat == "comp_hit"))     
                #geom_point(aes(x=ifelse(cl.data$cat == "comp_hit", log(cl.data$size.CO+1), NA), 
                #               y=ifelse(cl.data$cat == "comp_hit", log(cl.data$size.TR+1), NA), 
                #               size = cl.data$size.adj, color = cl.data$cat))}} +
                               
                geom_point(aes(size = cl.data$size.adj, color = cl.data$cat)) +
                
                xlab(paste0("Abundance [", pre.CO, "]")) +
                ylab(paste0("Abundance [", pre.TR, "]")) +
                coord_cartesian(xlim = c(0,lim), ylim = c(0,lim)) +
                
                scale_colour_manual(values = c("no_hit" = c(alpha("#CCCCCC", 0.6)), 
                                               "cl_hit" = c(alpha("#777777", 0.6)), 
                                               "comp_hit" = c(alpha("red", 0.6))
                                               )) +
                geom_text_repel(aes(label=ifelse(cat == "comp_hit", paste0("CL", as.character(cl)),'')),
                                color = "#333333", size = 2, 
                                segment.size = 0.2, segment.color="#333333",
                                nudge_x=0.5, max.iter=10000)
                )
                
                graphics.off()
                print(paste("Finished comparative scatter plot."))

                }}
    '''.format(image_=r_image_format(orientation='square'))

    return Rscatter_comp_plot


def rmanhattan_comp_plot():

    Rmanhattan_comp_plot = '''
    comp_manhattan <- function(converted.coverage.path, png.path, hiconf.path, comp.cand.path, pos.ar, win.size) {{
                
                print(paste("Starting comparative Manhattan plot."))

                pos.array <- data.frame(read.table(pos.ar))
                names(pos.array) <- as.character(unlist(pos.array[1,]))
                pos.array <- pos.array[-1,]
                
                coverage.data <- read.table(converted.coverage.path, header=T)
                names(coverage.data) <- c("chr", "start", "end", "TR.all", "CO.all")
                
                chr.colors <- c("#555555", "#111111")
                chr.number <- length(unique(coverage.data$chr))
                            
                coverage.data$hiconf <- read.table(hiconf.path, header=T, 
                colClasses = c("NULL", "NULL", "NULL", NA))$hiconf
                
                coverage.data$comp_cand <- read.table(comp.cand.path, header=T, 
                colClasses = c("NULL", "NULL", "NULL", NA))$comp_cand
                
                coverage.data$label <- ""
                
                for (row in 1:nrow(pos.array)) {{
                  pos.start <- as.integer(paste(pos.array[row,]$start))
                  pos.chr <- as.character(pos.array[row,]$chr)
                  pos.cand <- as.character(pos.array[row,]$cand)
                  coverage.data$label <- ifelse(coverage.data$chr == pos.chr & 
                                                coverage.data$start <= pos.start & 
                                                coverage.data$end > pos.start, 
                                                pos.cand, coverage.data$label)
                }}
                
                y.up.lim <- max(coverage.data[4])*1.1
                y.dw.lim <- ifelse(is.na(coverage.data[1,7]), 0, max(coverage.data[7])*1.1)
                
                coverage.data.tmp <- coverage.data %>% group_by(chr) %>% summarise(chr_len=max(start)) %>%  
                mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(coverage.data, ., by=c("chr"="chr")) %>%
                arrange(chr, start) %>% mutate(poscum=start+tot)
                                            
                axis.label <- coverage.data.tmp %>% group_by(chr) %>%
                summarize(center=( max(poscum) + min(poscum) ) / (2) )
                
                theme_set(theme_bw(base_size = 10) +
                  theme( 
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()
                  ))
                
                {image_}
                
                print(
                ggplot(coverage.data.tmp, aes(x=poscum)) +
                
                #geom_point(y=-coverage.data.tmp[,7], alpha=0.8, size=0.5, color=c(alpha("#CCCCCC", 0.5))) +
                
                geom_point(aes(y=coverage.data.tmp[,4], color=as.factor(coverage.data.tmp$chr)), alpha=1.0, size=0.5) +
                
                {{if(any(coverage.data.tmp$hiconf != 0)) 
                geom_point(y=ifelse(coverage.data.tmp$hiconf != 0, coverage.data.tmp[,4], NA), 
                alpha=0.8, size=0.5, color="#333333")}} +
                
                {{if(any(coverage.data.tmp$comp_cand != 0)) 
                geom_point(y=ifelse(coverage.data.tmp$comp_cand != 0, coverage.data.tmp[,4], NA), 
                alpha=0.8, size=0.5, color="red")}} +
                
                scale_color_manual(values = rep(chr.colors, chr.number )) +
                scale_x_continuous(label=axis.label$chr, breaks = axis.label$center) +
                scale_y_continuous(limits = c(0, y.up.lim)) +
                xlab(paste0("Chromosome [", win.size, "bp windows]")) +
                ylab("Mean depth [BPM]") +
                
                geom_text_repel(y=y.up.lim*0.95, label=gsub("eccCand_", "", coverage.data.tmp$label),
                                color = "#333333", size = 2, 
                                segment.size = 0.2, segment.color="#333333",
                                nudge_y=y.up.lim*0.95, max.iter=10000)
                )
                
                graphics.off()
                print(paste("Finished comparative Manhattan plot."))

                }}
    '''.format(image_=r_image_format(orientation='landscape'))

    return Rmanhattan_comp_plot
