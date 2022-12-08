
# Climate change metrics to assess the effect of extreme weather events on biodiversity 
# Script: PCA on range-based metrics
# Version 1.0 - Aug 2022
# Juan David Gonzalez-Trujillo

# Libraries
library(ggfortify)
library(ggpubr)
library(factoextra)
library(colorspace)


# PCA ---------------------------------------------------------------------

# Uploading the database
rbm_metrics_db <- readRDS("metrics_jul2022/RB_metrics_db.RDS")
## 89 metrics for three climatic parameters

# Standardizing variables: mean 0 SD 1
rbm_metrics_db_tr <- data.frame(scale(rbm_metrics_db[,-c(1:2)])) %>%
  replace(is.na(.), 0)
colMeans(rbm_metrics_db_tr)
apply(rbm_metrics_db_tr, 2, sd)

# Principal component analysis
pca_rbm_metrics <- prcomp(rbm_metrics_db_tr,retx = T)
summary(pca_rbm_metrics)

# Variance kept by every PC
frac_var <- function(x) x^2/sum(x^2)
library(scales)

# table
pca_rbm_metrics$sdev %>% 
  as_tibble() %>% 
  frac_var() %>% 
  mutate(Comp = colnames(pca_rbm_metrics$x)) 

# figure
pca_rbm_metrics$sdev %>% 
  as_tibble() %>% 
  frac_var() %>% 
  mutate(Comp = colnames(pca_rbm_metrics$x)) %>% 
  slice(1:9) %>% 
  ggplot(aes(x=Comp, y = value)) + 
  geom_bar(stat = "identity", fill = "#40B0A6") +
  geom_hline(yintercept = 0.02, linetype=2) +
  xlab("Principal Components") +
  scale_y_continuous(name = "Variance Explained", breaks = seq(0,0.8,0.1), labels = percent_format(accuracy = 5L)) +
  theme_classic(base_size = 14)

# biplot
fviz_pca_var(pca_rbm_metrics,
             col.var = "contrib", # Color by contributions to the PC
             repel = TRUE     # Avoid text overlapping
)

# Maps of PC scores -------------------------------------------------------

# Storing the results for the first 8 PCs in another table
rot_8pcs <- pca_rbm_metrics$rotation[,c(1:8)]

pca_map_df <- data.frame(rbm_metrics_db[,c(1:2)],pca_rbm_metrics$x[,c(1:8)])

pc1 <- ggplot(pca_map_df,aes(longitude,latitude,color=PC1))+
  geom_point(pch=15,size=0.5)+
  scale_color_continuous_divergingx(name="Score", "RdGy") +
  theme_test() +
  xlab("") + ylab("") +
  theme(legend.position = c(.8,.8),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.direction="horizontal",
        legend.key = element_blank(),
        panel.background = element_rect(fill = "azure3",
                                        colour = "azure3"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
        plot.margin = margin(0,0,0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)); pc1

pc2 <- ggplot(pca_map_df,aes(longitude,latitude,color=PC2))+
  geom_point(pch=15,size=0.5)+
  scale_color_continuous_divergingx(name="Score", "RdGy") +
  theme_test() +
  xlab("") + ylab("") +
  theme(legend.position = c(.8,.8),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.direction="horizontal",
        legend.key = element_blank(),
        panel.background = element_rect(fill = "azure3",
                                        colour = "azure3"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = margin(0,0,0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)); pc2

pc3 <- ggplot(pca_map_df,aes(longitude,latitude,color=PC3))+
  geom_point(pch=15,size=0.5)+
  scale_color_continuous_divergingx(name="Score", "RdGy") +
  theme_test() +
  xlab("") + ylab("") +
  theme(legend.position = c(.8,.8),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.direction="horizontal",
        legend.key = element_blank(),
        panel.background = element_rect(fill = "azure3",
                                        colour = "azure3"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = margin(0,0,0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)); pc3

pc4 <- ggplot(pca_map_df,aes(longitude,latitude,color=PC4))+
  geom_point(pch=15,size=0.5)+
  scale_color_continuous_divergingx(name="Score", "RdGy") +
  theme_test() +
  xlab("") + ylab("") +
  theme(legend.position = c(.8,.8),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.direction="horizontal",
        legend.key = element_blank(),
        panel.background = element_rect(fill = "azure3",
                                        colour = "azure3"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = margin(0,0,0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)); pc4

pcs_maps <- ggarrange(pc1,pc2,pc3,pc4,nrow = 4);pcs_maps

# ggsave("../first4_pca_rbm_maps.pdf",plot = last_plot(),device = "pdf",
#        width = 4,height = 12,units = "in",dpi = 300)

