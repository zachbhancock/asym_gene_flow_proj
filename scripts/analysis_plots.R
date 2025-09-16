#import
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(viridis)
library(triangulaR)
library(vcfR)
library(marmap)
library(raster)
library(ncdf4)
library(ggnewscale)
library(ggimage)
library(rangeExpansion)
library(vcfR)

setwd("C:/Users/zhancock/Desktop/shovel-bugs")

#read in data
data <- read.table("revised_shovel_bug_results.txt", header=TRUE)
fst.dist <- read.table("data/fst_dist_data.txt", header=TRUE)
psi.results <- read.table("data/psi_results.txt", header=TRUE)
sim.psi <- read.table("data/sim_psi_results.txt", header=TRUE)
snp.file <- "data/feems_bed.bed"
coord.file <- "data/coords_for_psi.txt"
vcf.data <- read.vcfR("data/feems_filt.recode.vcf")
pop.map <- read.table("data/new_pop_file.txt", header=TRUE) #must have header "id" and "pop"

#get distance from the center in miles
data$center.dist <- 69*(abs(36 - data$latitude))

###analyses in R#####

#linear regressions, all data
depth_model <- lm(data$pi~data$mean_depth) #p = 5.281e-7, r^2 = 0.277
lat_model <- lm(data$mean_depth~data$latitude) #p = 0.4968, r^2 = -0.007071
center_model <- lm(data$mean_depth~data$center.dist) #p = 0.002629, r^2 = 0.1026
center_model_w_depth <- lm(data$pi~data$center.dist+data$mean_depth) #p = 3.794e-16, dist.p = 1.35e-11, depth.p = 5.14e-5, r^2 = 0.6066
f_model <- lm(data$F~data$center.dist+data$mean_depth) #p < 2.2e-16, dist.p=6.98e-15, depth.p=0.00019, r^2 = 0.6661
fst_model <- lm(fst.dist$fst_df~fst.dist$geo.dist) #p < 2.2e-16, r^2 = 0.5698

#linear regressions, dropping mean_depth > 30

filtered_depth <- data[data$mean_depth < 30, ]
filt_model <- lm(filtered_depth$pi ~ filtered_depth$mean_depth) #p = 2.787e-5, r^2 = 0.3815
filt_dist_center <- lm(filtered_depth$center.dist ~ filtered_depth$mean_depth) #p = 0.5754, r^2 = -0.01926
filt_het <- lm(filtered_depth$pi ~ filtered_depth$center.dist) #p = 0.0009853, r^2 = 0.249

#linear regressions, rarefaction down to mean individual coverage of 7, dropping all below 6.8
#determine the median absolute deviation
0.5*mad(data$low_mean_depth) #0.0394; dropped inds w/ cov < 6.884 - 0.0394 = 6.844
low_depth <- data[data$low_mean_depth > 6.84, ]
low_model <- lm(low_depth$low_cov_pi ~ low_depth$low_mean_depth) #p = 0.3674, r^2 = -0.003658
low_dist_center <- lm(low_depth$low_mean_depth ~ low_depth$center.dist) #p = 0.9134, r^2 = -0.02147
low_het <- lm(low_depth$low_cov_pi ~ low_depth$center.dist) #p = 2.351e-08, r^2 = 0.4849

#psi analysis
ploidy <- 2
raw.data <- load.data(snp.file, coord.file, ploidy=ploidy)
pop <- make.pop(raw.data, ploidy)
psi <- get.all.psi(pop)
region <- list(NULL, c("region_1", "region_2"))
results <- run.regions(region=region, pop=pop, psi=psi, exclude.ocean=F)
#write.table(psi, file="psi_table.txt", sep="\t")

####main text figures#####

#make figure 1, main text
#pca, panel B
pca.plot <- ggplot(data, aes(x=PC1, y=PC2, color=latitude)) + geom_point(size=5) + theme_bw() + xlab("PC1") + 
  ylab("PC2") + scale_color_viridis() + theme(axis.title=element_text(size=16),
                                              axis.text=element_text(size=14),
                                              legend.position="top", 
                                              legend.title=element_text(size=14),
                                              legend.text=element_text(size=12)) +
  labs(color="latitude")

#fst, panel C

fst.plot <- ggplot(fst.dist, aes(x=geo.dist, y=fst_df, color=comparison, shape=comparison)) + geom_point(size=4) +
  theme_bw() + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
                     legend.text=element_text(size=14), legend.title=element_text(size=16),
                     legend.position="top", axis.text.x = element_text(hjust = 1)) + xlab("geographic distance (mi)") +
  ylab(expression(F[ST]/(1-F[ST]))) + scale_color_manual(values = c("#E69F00", "#CC79A7", "#0072B2"))

#map with shovel-bug imagine, panel A

#nc_open("sws_baseline_2000_2019_depthsurf_2e98_ab53_db03_U1720468843697.nc")
#sws_mean <- brick("sws_baseline_2000_2019_depthsurf_2e98_ab53_db03_U1720468843697.nc", varname="sws_mean")
#NA_ex <- extent(-83, -68, 25, 45)
#NA_sws <- crop(sws_mean, NA_ex)
#test_spdf <- as(NA_sws, "SpatialPixelsDataFrame")
#test_df <- as.data.frame(test_spdf)
#colnames(test_df) <- c("value", "x", "y")
#x_range <- range(test_df$x)
#y_range <- range(test_df$y)
#
#bare.map.plot <-  ggplot() +
#  geom_tile(data=test_df, aes(x=x, y=y), alpha=0.8) +
#  xlab("longitude") + ylab("latitude") +
#  coord_equal() + theme_bw() +
#  scale_x_continuous(limits = x_range, expand = c(0, 0)) +
#  scale_y_continuous(limits = y_range, expand = c(0, 0)) +
#  theme(panel.background = element_rect(fill="lightgray"),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        axis.title=element_text(size=16), axis.text=element_text(size=14), 
#        legend.position="top", legend.title=element_text(size=14), 
#        legend.text=element_text(size=12))

usa <- map_data("state")

NAmap <- ggplot() + geom_polygon(data = usa, 
                                 aes(x=long, y = lat, group = group), 
                                 fill = "lightgray", 
                                 color="black") + theme_bw() + theme(axis.title = element_text(size=16),
                                                                     axis.text = element_text(size = 14)) +
  coord_fixed(xlim = c(-83, -68),  ylim = c(25, 45), ratio =1.2) + xlab("longitude") + ylab("latitude")

map.with.arrows <- NAmap +  geom_curve(
  aes(x = -78.9, y = 26, xend = -74, yend = 34.4),
  arrow = arrow(length = unit(0.05, "npc"), angle=20), curvature=-0.62, color="black", size=3) +
  geom_curve(
    aes(x = -69, y = 41, xend = -74.5, yend = 36),
    arrow = arrow(length = unit(0.05, "npc"), angle=20), curvature=0.35, color="black", size=1) +
  geom_segment(
    aes(x = -69, y = 27.2, xend = -71, yend = 27.2), color="black", size=1) +
  annotate("text", x=-70, y=26.9, label="4 mm", color="black")

fig_map_a <- map.with.arrows + new_scale_fill() + 
  geom_point(data=data, aes(x=longitude, y=latitude, fill=latitude), shape=21, size=6, color="black") + 
  scale_fill_viridis() + guides(fill="none")

#add shovel-bug image
fig_map_a <- fig_map_a +  geom_image(
  data = tibble(long = -73, lat= 30),
  aes(image = "shovel-bug.png", x=long, y=lat), size=0.55
)

fig1_map <- ggarrange(fig_map_a, nrow=1, labels=c("A"), font.label=list(size=18))
fig1bc <- ggarrange(pca.plot, fst.plot, nrow=2, labels=c("B", "C"), font.label=list(size=18))

#final figure 1, main text
fig1 <- ggarrange(fig1_map, fig1bc)
ggsave("revised_fig1.png", plot = fig1, dpi = 300, bg = "white")

#make figure 2, main text

feems_node_pos <- read_csv('feems_node_pos.csv', col_names=c('lon', 'lat'), col_types='dd')
locs <- read.table("haust_locations.coord", header=FALSE)
usa <- map_data("state")

NAmap <- ggplot() + geom_polygon(data = usa, 
                                 aes(x=long, y = lat, group = group), 
                                 fill = "white", 
                                 color="black") +
  coord_fixed(xlim = c(-85, -65),  ylim = c(25, 48), ratio =1.2)

#produces map with feems and conStruct results

feems.plot <- NAmap + geom_segment(mapping=aes(x=lon1, xend=lon2, y=lat1, yend=lat2, color=w), data=feems.edges,
                                   linewidth=1) +
  scale_color_gradientn(colors=feems.colors, name="log10(w)") +
  theme_bw() + xlab("Longitude") + ylab("Latitude") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14))

load(sprintf("shovelbugs_K2_conStruct.results.Robj"))
load(sprintf("shovelbugs_K2_data.block.Robj"))
tmp <- conStruct.results$chain_2$MAP$admix.proportions
tmp.df <- as.data.frame(tmp)
names(tmp.df) <- c("layer_1", "layer_2")
coords_new.df <- as.data.frame(coords_new)
names(coords_new.df) <- c("long", "lat")
admix.df <- cbind(tmp.df, coords_new.df)
write.table(admix.df, file="coords_new.txt")
admix.layers <- read.table("layer_contributions.txt", header=TRUE)
pie.df <- data.frame(long = admix.layers$long, lat = admix.layers$lat)
pie.df$pop <- admix.layers$pop
pie.df$A <- admix.layers$layer_1
pie.df$B <- admix.layers$layer_2

feems.construct.plot <- feems.plot + geom_scatterpie(aes(x=long, y=lat, group=pop, r=0.4), data = pie.df, cols=LETTERS[1:2],
                                                     legend_name = "layer") + geom_label(data=)

admix.df <- read.table("coords_new_layers.txt", header=TRUE)
data_long <- gather(admix.df, layers, contribution, layer_1, layer_2, factor_key=TRUE)
data_long$sample <- factor(data_long$sample, levels = unique(data_long$sample))
ggplot(data=data_long, aes(x=sample, y=contribution, fill=layers)) + geom_bar(stat="identity") + coord_flip() + theme_bw() +
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=14),
        axis.title.x=element_blank(),
        legend.position="none")

ggarrange(struct.plot, feems.construct.plot, labels=c("A", "B"), widths=c(1,2.5))

#make figure 3, main text

a <- ggplot(data, aes(x=latitude, y=pi, color=latitude)) + geom_smooth(span=1, fill="lightgray") +
  geom_point(size=4) + scale_color_viridis() +
  xlab("latitude") + ylab(expression("individual heterozygosity "(pi))) + theme_bw() + theme(axis.title=element_text(size=14),
                                                               axis.text=element_text(size=12))

b <- ggplot(data, aes(x=center.dist, y=F, color=latitude)) + geom_smooth(method="lm") +
  geom_point(size=4) + scale_color_viridis() +
  xlab("distance from center (mi)") + ylab("inbreeding coefficient (F)") + theme_bw() + theme(axis.title=element_text(size=14),
                                                                                              axis.text=element_text(size=12))
#final figure 3, main text
fig3 <- ggarrange(a,b, common.legend=TRUE)
ggsave("revised_fig3.png", plot = fig3, dpi=300, bg="white")

#make figure 4, main text

psi.results$model <- "empirical"
empirical <- ggplot(psi.results, aes(x=latitude, y=psi, fill=het)) + theme_bw() +
  geom_point(shape=21, size=7) + scale_fill_viridis(option="magma") +
  ylab(expression(sum(psi[ij]))) + xlab("latitude") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) +
  labs(fill=expression(pi)) + geom_hline(yintercept=0, linetype="dashed") + facet_wrap(~model)

#next, get the simulated results
sim.psi$model <- factor(sim.psi$model, levels=c("center_asym", "center_sym", "equal_asym", "equal_sym")) 
psi.sim.plot <- ggplot(sim.psi, aes(x=population, y=psi.sum, fill=het)) + geom_point(size=4, shape=21) +
  theme_bw() + scale_fill_viridis(option="magma", name=expression(pi)) + 
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), legend.text=element_text(size=12),
        legend.title=element_text(size=14), strip.text=element_text(size=14)) +
  ylab(expression(sum(psi[ij]))) + scale_x_continuous(breaks=seq(1, 10, by=1)) +
  geom_hline(yintercept=0, linetype="dashed") + facet_wrap(~model)

#final figure 4, main text
fig4 <- ggarrange(empirical, psi.sim.plot, labels=c("A", "B"), font.label=list(size=18),
          common.legend=TRUE, legend="right")
ggsave("revised_fig4.png", plot = fig4, dpi=300, bg = "white")

####supplementary figures####

#figure s1

fig_s1a <- ggplot(filtered_depth, aes(x=mean_depth, y=pi, color=latitude)) + theme_bw() + geom_smooth(method="lm") +
  geom_point(size=8) + scale_color_viridis() +
  ylab(expression(pi)) + xlab("mean depth") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) + labs(color="latitude")

fig_s1b <- ggplot(filtered_depth, aes(x=center.dist, y=mean_depth, color=latitude)) + theme_bw() + 
  geom_smooth(method="lm") + geom_point(size=8) + scale_color_viridis() +
  ylab("mean depth") + xlab("distance from center (mi)") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) + labs(color="latitude")

ggarrange(fig_s1a, fig_s1b, labels=c("A", "B"), common.legend=TRUE)

#figure s2

fig_s2a <- ggplot(low_depth, aes(x=low_mean_depth, y=low_cov_pi, color=latitude)) + theme_bw() +
  geom_point(size=8) + scale_color_viridis() +
  ylab(expression("heterozygosity "(pi))) + xlab("mean depth") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) + labs(color="latitude")

fig_s2b <- ggplot(low_depth, aes(x=center.dist, y=low_mean_depth, color=latitude)) + theme_bw() + geom_point(size=8) + scale_color_viridis() +
  ylab("mean depth") + xlab("distance from center (mi)") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) + labs(color="latitude")

fig_s2 <- ggarrange(fig_s2a, fig_s2b, labels=c("A", "B"), common.legend=TRUE)
ggsave("fig_s2.png", plot=fig_s2, dpi = 300)


##triangulaR analysis####
#this produces figure s4 for supp mat

vcf.diff <- alleleFreqDiff(vcfR = vcf.data, pm = pop.map, p1="Cocoa_Beach_FL", p2="Crescent_Beach_ME", difference = 0.9)
hi.het <- hybridIndex(vcfR = vcf.diff, pm = pop.map, p1 = "Cocoa_Beach_FL", p2 = "Crescent_Beach_ME")

cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837", "aliceblue", "blue", "red",
          "green", "yellow", "skyblue", "pink", "orange", "brown", "lightblue", "dodgerblue", "darkred")

#fig S4, supp mat
triangle.plot(hi.het, colors = cols, cex=5)

#make figure s5

fig_s5 <- ggplot(filtered_depth, aes(x=center.dist, y=pi, color=latitude)) + theme_bw() +
  geom_smooth(method="lm") +
  geom_point(size=8) + scale_color_viridis() +
  ylab(expression("heterozygosity "(pi))) + xlab("distance from center (mi)") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) + labs(color="latitude")

ggsave("fig_s5.png", plot=fig_s5, dpi = 300)

#make figure s6
fig_s6 <- ggplot(low_depth, aes(x=center.dist, y=low_cov_pi, color=latitude)) + theme_bw() +
  geom_smooth(method="lm") +
  geom_point(size=8) + scale_color_viridis() +
  ylab(expression("heterozygosity "(pi))) + xlab("distance from center (mi)") +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), 
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        strip.text=element_text(size=14)) + labs(color="latitude")

ggsave("fig_s6.png", plot=fig_s6, dpi=300)

#make fig s7

center.sym.psi <- read.table("center_sym_size_results.txt", header=TRUE)
center.sym.psi$pop.size <- factor(center.sym.psi$pop.size, levels=c("5:1", "25:1", "50:1"))

ggplot(center.sym.psi, aes(x=population, y=psi.sum, fill=het)) + geom_point(size=8, shape=21) +
  theme_bw() + scale_fill_viridis(option="magma", name=expression(pi)) + 
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), legend.text=element_text(size=12),
        legend.title=element_text(size=14), strip.text=element_text(size=14)) +
  ylab(expression(sum(psi[ij]))) + scale_x_continuous(breaks=seq(1, 10, by=1)) +
  geom_hline(yintercept=0, linetype="dashed") + facet_wrap(~pop.size) + ylim(-2,1.5)

##end##
