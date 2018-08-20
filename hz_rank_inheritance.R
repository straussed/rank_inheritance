################################################################################
#                                                                              #     
#                     Make animated acquisition figure                         #
#                              Eli Strauss                                     #
#                                                                              #
#                              August 2018                                     #
#                                                                              #
################################################################################


################################################################################
## Here is a script for making animated plots of rank acquisition
## from interaction data. To use the script right out of the box: 
##    1) download all files in the github repo to a single folder
##    2) open the .Rproj file in Rstudio
##    3) in Rstudio, open the script
##    4) make sure all required packages and ImageMagick are installed
##    5) run the script

## To use the script on its own (i.e., without the .Rproj file), you will need to 
## modify the here() calls. here() references the location of the project file. 
## To avoid using here() and the .Rproj file, simply replace the here() calls with
## direct filepaths to the relevant files/folders.

################################################################################
## (1) Prepare workspace

#####Load required packages####
library(dplyr)
library(ggplot2)
library(aniDom)
library(viridis)
library(igraph)
library(grid)
library(gridGraphics)
library(here)
####Set global options####
options(stringsAsFactors = FALSE)

####You will also need to install the Imagemagick software####
####It's not an R package. It can be found at: 
####https://www.imagemagick.org/script/download.php####
################################################################################

################################################################################
## (2) Read and prepare data

####Read in interaction data####
intx <- read.csv(here('hz_cohort_aggs.csv'))

####Specify IDs in order of maternal rank####
ids <- c('clus', 'lole', 'smag', 'puff', 'savy', 'boot', 'cokl')
cohort.ranks <- data.frame(ids, rank = 1:length(ids))
################################################################################

################################################################################
## (3) Calculate Elo scores for each cub

####Calculate scores####
scores <- elo_scores(winners = intx$aggressor, losers = intx$recip,
                    return.trajectories = TRUE, identities = ids, K = 20)

####Get total interaction count####
intx.count <- nrow(intx) + 1 ##Add one because Elo scores start at 0 before intx take place

####Assemble into data frame for plotting####
cohort.scores <- data.frame(score = c(scores[1,], scores[2,], scores[3,], scores[4,],
                   scores[5,], scores[6,], scores[7,]),
                   ids = rep(ids, each = intx.count),
                   date = rep(c(intx$date[1], intx$date), 7),
                   mat.rank = rep(1:7, each = intx.count),
                   intx.number = rep(1:intx.count))


intx$intx.num <- 1:nrow(intx)
################################################################################

################################################################################
## (4) Create plots for each time step

####Set sliding window width for Hobson plot####
window <- 100 # number of interactions to be included in each network time step

####Set plotting aesthetics####
vertex.colors <- viridis(8)[1:7] ##need to add more colors if there are more than 7 individuals
node.size <- 40 
label.size <- 1

for(last.intx in seq(from = 0, to = intx.count, by = 10)){
  ##Name png for current iteration. Weird math is to ensure
  ##that files are interpreted in correct order when making gif
  f <- paste0(here('images/hz_cohort_'),
              (10+(1.1+last.intx)/1000), '.png')
  
  ####Plot Elo scores####
  p <- ggplot(filter(cohort.scores, intx.number <= last.intx), aes(x = intx.number, y = score, group = ids, col = as.factor(mat.rank)))+
    geom_line(size = 1)+
    geom_label(data = filter(cohort.scores, intx.number == last.intx), 
               aes(label = mat.rank), label.r = unit(0.55, 'lines'),
               label.padding = unit(0.35, "lines"), label.size = 0.5)+
    theme_classic()+
    scale_color_manual(values = viridis(8)[1:7]) + 
    theme(legend.position = 'none', axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    ylim(-500, 400) + 
    xlim(0, intx.count + 50) +
    ylab('Elo Score') + 
    xlab('--------Time-------->')
 
  ####Hobson plot (i.e. attribute-ordered network from Hobson et al 2015)####
  
  ##Select only interactions in the "window" including the last interaction
  ##The width of the window is set above. In this case, we use
  ##the 100 interactions that occur before the last interaction (last.intx) for this plot
  aggs.current <- filter(intx, intx.num <= last.intx, intx.num >= last.intx - window)
  
  ##Convert to summarized edgelist, with "wins" indicating the number of wins in window for that edge
  edges <- aggs.current %>% group_by(aggressor, recip) %>% summarize(wins = length(aggressor))
  
  ##Create aggression network
  aggression_network<-graph.data.frame(edges, directed = TRUE, vertices = cohort.ranks)
  
  ##Set colors for edges. If there are no edges, skip this part
  if(nrow(edges) > 0){
    ##Specify rank difference for edges to guide coloring
    E(aggression_network)$rank.diff <- left_join(edges, cohort.ranks, by = c('aggressor' = 'ids'))$rank - 
      left_join(edges, cohort.ranks, by = c('recip' = 'ids'))$rank
    
    ##Set color according to rank.diff
    E(aggression_network)[rank.diff>0]$color <- "#CC000095"
    E(aggression_network)[rank.diff<0]$color <- "#00009985"
  }
  
  ##Sets the name to display. Right now its set to the maternal rank
  aggression_network <- set.vertex.attribute(aggression_network, "name", value = 1:7)
  
  ##Set width of drawn edges to be equal to number of wins
  edge.width <- E(aggression_network)$wins
  
  ##Open png in which to save both plots
  png(filename = f, width = 6, height = 4, units = 'in', res = 200)
  
  ##Specify margins and number of plots to be displayed
  par(mar = c(0,0,0,0),
      mfrow = c(1,2))
  
  ##Plot Hobson plot
  plot.igraph(aggression_network,
              layout= data.frame(x = 0.5, y = nrow(cohort.ranks):1), #orders individuals by rank, with lower numbers = higher rank
              vertex.size=node.size, #this references the node circle size, set above
              vertex.label.family="sans", #uses sans font style
              vertex.label.color = "white", #sets node label text color
              vertex.label.cex=label.size, #this references the node label size, set above
              vertex.color= vertex.colors,
              edge.width=edge.width, #this references the edge width, set above
              edge.arrow.size=0.01, #this makes the arrows on the ends of the edges very small
              edge.curved=TRUE, #this makes the edges curved
              edge.loop.angle=1,
              edge.color=E(aggression_network)$color,
              rescale = FALSE,
              xlim = c(0,1),
              ylim = c(0,max(cohort.ranks$rank+1)))
  
  ##Add up- and down-hierarchy labels
  text('Down-\nhierarchy', x = 1.5, y = 7, col =  alpha("#00009985", 1),
       srt = 0)
  
  text('Up-\nhierarchy', x = -.5, y = 7, col =  alpha("#CC000095", 1),
       srt = 0)
  
  ##Print Elo score plot next to Hobson plot
  print(p, vp = viewport(x = 0.7, y = 0.5, width = 0.6, height = 0.8), newpage = FALSE)
  
  ##Finish writing plots to file
  dev.off()
}

################################################################################
## (5) Stitch plots together to create gif

####This section uses the ImageMagick software. I usually run this command
####in the command line. I've included the command line version in the comment here
####as well as a system call from R that invokes the command from R. For PC users,
####the command may need to be adapted. Check out resources on ImageMagick for more 
####information.

#### IMAGEMAGICK COMMAND (after navigating to folder with images)####
#   
#  convert -delay 10 *.png -delay 300 hz_cohort_11.2011.png ../hz_inheritance.gif  #
#
##The two delay parameters specify the number of hundreths of a second between each frame.
##I've set it to use 10 hundreths for each frame, then pause for
##an additional 3 seconds on the last frame

####R code to invoke the above command####
frame.delay = 10
final.delay = 300
last.frame = here(paste0('images/',tail(list.files(here('images/')),1)))
gif.name <- 'hz_inheritance.gif'
system(command = paste('convert -delay',frame.delay, here('images/*.png'), #all images, normal delay
                        '-delay', final.delay, last.frame, #long pause on last frame
                        here(gif.name))) #name of gif
