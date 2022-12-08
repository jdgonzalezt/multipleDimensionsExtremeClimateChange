
# Climate change metrics to assess the effect of extreme weather events on biodiversity 
# Script: EE on mangrove areas
# Version 1.0 - Aug 2022
# Juan David Gonzalez-Trujillo

# Libraries
library(sf)
library(ggmap)
library(gridExtra)
library(raster)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(trend)
library(ggpubr)

options("scipen"=100, "digits"=4)

# map for plotting --------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

map <- ggplot(data = world) +
  geom_sf(fill= "antiquewhite",alpha=0.8) +
  coord_sf(xlim = c(-100, -56), ylim = c(6, 35), expand = FALSE)+
  theme_test() +
  xlab("") + ylab("") +
  theme(legend.position = c(.8,.8),
        legend.key.size = unit(0.25, "in"),
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
  scale_y_continuous(expand = c(0,0))


# Data preparation --------------------------------------------------------

mangrove_distro <- st_read("~/Dropbox (In Progress)/GIS/MANGROVE_DATA_LAGOMASINO/mangrove_baseline_distro.shp")

mangrove_distro <- mangrove_distro %>%
  distinct(id,.keep_all = T)

rbm_metrics_db <- readRDS("./metrics_jul2022/RB_metrics_db.RDS")

coords <- data.frame(rbm_metrics_db[,c("longitude","latitude")])
coords$id <-c(1:nrow(coords))
coords$lon <-coords$longitude
coords$lat <-coords$latitude

sp_coords <- st_as_sf(x=coords,
                      coords = c("longitude", "latitude"),
                      crs= "+proj=longlat +datum=WGS84 +no_defs ")

mangrove_xy <- mangrove_distro %>% 
  cbind(
    sp_coords[st_nearest_feature(mangrove_distro,sp_coords),]) %>%
  mutate(dist = st_distance(geometry, geometry.1, by_element = T)) %>%
  mutate(
    longitude = st_coordinates(geometry.1)[,1],
    latitude = st_coordinates(geometry.1)[,2]
  )%>%
  filter(dist < units::set_units(10000, m)) %>%
  dplyr::select(-max,-label,-id.1,-lon,-lat)


damage_shp_list <- list.files(path = "~/Dropbox (In Progress)/GIS/MANGROVE_DATA_LAGOMASINO/MANGROVE_DAMAGE_AREA", 
                       pattern = ".shp$",
                       full.names = TRUE, recursive = TRUE)

mangrove_damage = data.frame()

for (year in 1:23){
  tryCatch({
    print(year+1997)
    yr_i <- read_sf(damage_shp_list[[year]]) %>%
      group_by(id) %>%
      summarise(
        damaged_area = sum(area)
      ) 
    
    joined_sf <- mangrove_xy %>%
      as.data.frame() %>%
      dplyr::left_join(yr_i,by=c("id")) %>%
      mutate(
        yr = year+1997,
        damaged = if_else(is.na(damaged_area),"NO","YES")
      ) %>%
      dplyr::select(-geometry.y)
      
    mangrove_damage = rbind(mangrove_damage,joined_sf)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

mangrove_damage %>%
  group_by(yr,damaged) %>%
  summarise(
    n=n()
  ) %>%
  ggplot(aes(yr,n,fill=damaged))+
  geom_bar(stat="identity")

recovery_shp_list <- list.files(path = "~/Dropbox (In Progress)/GIS/MANGROVE_DATA_LAGOMASINO/MANGROVE_RECOVERY/", 
                              pattern = ".shp$",
                              full.names = TRUE, recursive = TRUE)



mangrove_recovery = data.frame()

for (year in 1:22){
  tryCatch({
    print(year+1997)
    yr_i <- read_sf(recovery_shp_list[[year]]) %>%
      group_by(id) %>%
      summarise(
        recovered_area = sum(area)
      ) 
    joined_sf <- mangrove_xy %>%
      as.data.frame() %>%
      dplyr::left_join(yr_i,by=c("id")) %>%
      mutate(
        yr = year+1997,
        p_yr = yr - 1) %>%
      dplyr::select(-geometry.y)
    
    mangrove_recovery = rbind(mangrove_recovery,joined_sf)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

dam_rec <- mangrove_damage %>%
  as.data.frame()%>%
  left_join(mangrove_recovery[,c("id","p_yr","recovered_area")],
            by=c("yr"="p_yr",
                 "id" )) %>%
  filter(yr < 2019)

dam_rec <- dam_rec %>%
  filter(damaged == "YES") %>%
  mutate(recovered = if_else(is.na(recovered_area),
                             "NO", "YES"))

saveRDS(dam_rec,"mangrove_recovery_db.rds")

# Metrics~Damage area -----------------------------------------------------------------

dam_rec <- readRDS("mangrove_recovery_db.rds")

# Selecting four metrics related to droughts, floods, heatwaves and cyclones

dam_rec %>%
  group_by(longitude,latitude,recovered) %>%
  summarise(
    n = n()
  ) %>%
  filter(recovered == "YES") %>%
  ggplot(aes(longitude,latitude,color=n))+
  geom_point()

load("./metrics_jul2022/ree_eb_t2m_djfma.Rdata")
load("./metrics_jul2022/ree_eb_tp_djfma.Rdata")
load("./metrics_jul2022/ree_eb_ws_djfma.Rdata")


# Statistical tests -------------------------------------------------------

library(effsize)

comparisons_df <- data.frame()

# Heatwaves

db_t2m_95th <- dam_rec %>%
  left_join(ree_eb_95_t2m_djf, 
            by = c("yr","longitude","latitude")) %>%
  filter(!is.na(exceedance_no)) %>%
  as.data.frame()

# Intensity 

wt <- wilcox.test(
  db_t2m_95th$intensity_cumulative[db_t2m_95th$recovered == "YES"],
  db_t2m_95th$intensity_cumulative[db_t2m_95th$recovered == "NO"]
)

es <- cliff.delta(
  db_t2m_95th$intensity_cumulative[db_t2m_95th$recovered == "YES"],
  db_t2m_95th$intensity_cumulative[db_t2m_95th$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "1. Temperature (max)",
  dimension = "Intensity",
  threshold = "P95",
  metric = "MCI_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_t2m_95th$intensity_cumulative[db_t2m_95th$recovered == "YES"]),
  n_no = length(db_t2m_95th$intensity_cumulative[db_t2m_95th$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)



# Duration

wt <- wilcox.test(
  db_t2m_95th$duration[db_t2m_95th$recovered == "YES"],
  db_t2m_95th$duration[db_t2m_95th$recovered == "NO"]
)

es <- cliff.delta(
  db_t2m_95th$duration[db_t2m_95th$recovered == "YES"],
  db_t2m_95th$duration[db_t2m_95th$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "1. Temperature (max)",
  dimension = "Duration",
  threshold = "P95",
  metric = "RT_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_t2m_95th$duration[db_t2m_95th$recovered == "YES"]),
  n_no = length(db_t2m_95th$duration[db_t2m_95th$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Frequency

db_t2m_95th_fr <- data.frame()

for (year in 1998:2019){
  print(year)
  year_n <- ree_eb_95_t2m_djf %>%
    filter(yr == year) %>%
    group_by(yr,longitude,latitude) %>%
    filter(exceedance_no == max(exceedance_no)) %>%
    full_join(ree_eb_95_t2m_djf %>%
                filter(yr == year-1) %>%
                group_by(yr,longitude,latitude) %>%
                filter(exceedance_no == max(exceedance_no)), 
              by = c("longitude","latitude")) %>%
    group_by(yr.x,longitude,latitude) %>%
    summarise(
      REE_pe = (exceedance_no.x + exceedance_no.y)
    ) %>%
    as.data.frame()
  
  db_t2m_95th_fr <- rbind(db_t2m_95th_fr,year_n)
}

db_fr <- dam_rec %>%
  left_join(db_t2m_95th_fr, 
            by = c("yr"="yr.x","longitude","latitude")) %>%
  filter(!is.na(REE_pe)) %>%
  as.data.frame()


wt <- wilcox.test(
  db_fr$REE_pe[db_fr$recovered == "YES"],
  db_fr$REE_pe[db_fr$recovered == "NO"]
)

es <- cliff.delta(
  db_fr$REE_pe[db_fr$recovered == "YES"],
  db_fr$REE_pe[db_fr$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "1. Temperature (max)",
  dimension = "Frequency",
  threshold = "P95",
  metric = "REE_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_fr$REE_pe[db_fr$recovered == "YES"]),
  n_no = length(db_fr$REE_pe[db_fr$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Droughts

db_tp_5th <- dam_rec %>%
  left_join(ree_eb_5_tp_dfjma, 
            by = c("yr","longitude","latitude")) %>%
  filter(!is.na(exceedance_no)) %>%
  as.data.frame()

# Intensity 

wt <- wilcox.test(
  db_tp_5th$intensity_cumulative[db_tp_5th$recovered == "YES"],
  db_tp_5th$intensity_cumulative[db_tp_5th$recovered == "NO"]
)

es <- cliff.delta(
  db_tp_5th$intensity_cumulative[db_tp_5th$recovered == "YES"],
  db_tp_5th$intensity_cumulative[db_tp_5th$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "2. Acc. precipitation (min)",
  dimension = "Intensity",
  threshold = "P5",
  metric = "MCI_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_tp_5th$intensity_cumulative[db_tp_5th$recovered == "YES"]),
  n_no = length(db_tp_5th$intensity_cumulative[db_tp_5th$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Duration

wt <- wilcox.test(
  db_tp_5th$duration[db_tp_5th$recovered == "YES"],
  db_tp_5th$duration[db_tp_5th$recovered == "NO"]
)

es <- cliff.delta(
  db_tp_5th$duration[db_tp_5th$recovered == "YES"],
  db_tp_5th$duration[db_tp_5th$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "2. Acc. precipitation (min)",
  dimension = "Duration",
  threshold = "P5",
  metric = "RT_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_tp_5th$duration[db_tp_5th$recovered == "YES"]),
  n_no = length(db_tp_5th$duration[db_tp_5th$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Frequency

db_tp_5th_fr <- data.frame()

for (year in 1998:2019){
  print(year)
  year_n <- ree_eb_5_tp_dfjma %>%
    filter(yr == year) %>%
    group_by(yr,longitude,latitude) %>%
    filter(exceedance_no == max(exceedance_no)) %>%
    full_join(ree_eb_5_tp_dfjma %>%
                filter(yr == year-1) %>%
                group_by(yr,longitude,latitude) %>%
                filter(exceedance_no == max(exceedance_no)), 
              by = c("longitude","latitude")) %>%
    group_by(yr.x,longitude,latitude) %>%
    summarise(
      REE_pe = (exceedance_no.x + exceedance_no.y)
    ) %>%
    as.data.frame()
  
  db_tp_5th_fr <- rbind(db_tp_5th_fr,year_n)
}

db_fr_dro <- dam_rec %>%
  left_join(db_tp_5th_fr, 
            by = c("yr"="yr.x","longitude","latitude")) %>%
  filter(!is.na(REE_pe)) %>%
  as.data.frame()


wt <- wilcox.test(
  db_fr_dro$REE_pe[db_fr_dro$recovered == "YES"],
  db_fr_dro$REE_pe[db_fr_dro$recovered == "NO"]
)

es <- cliff.delta(
  db_fr_dro$REE_pe[db_fr_dro$recovered == "YES"],
  db_fr_dro$REE_pe[db_fr_dro$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "2. Acc. precipitation (min)",
  dimension = "Frequency",
  threshold = "P95",
  metric = "REE_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_fr_dro$REE_pe[db_fr_dro$recovered == "YES"]),
  n_no = length(db_fr_dro$REE_pe[db_fr_dro$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# high precipitation

db_tp_95th <- dam_rec %>%
  left_join(ree_eb_95_tp_dfjma, 
            by = c("yr","longitude","latitude")) %>%
  filter(!is.na(exceedance_no)) %>%
  as.data.frame()

# Intensity 

wt <- wilcox.test(
  db_tp_95th$intensity_cumulative[db_tp_95th$recovered == "YES"],
  db_tp_95th$intensity_cumulative[db_tp_95th$recovered == "NO"]
)

es <- cliff.delta(
  db_tp_95th$intensity_cumulative[db_tp_95th$recovered == "YES"],
  db_tp_95th$intensity_cumulative[db_tp_95th$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "3. Acc. precipitation (max)",
  dimension = "Intensity",
  threshold = "P95",
  metric = "MCI_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_tp_95th$intensity_cumulative[db_tp_95th$recovered == "YES"]),
  n_no = length(db_tp_95th$intensity_cumulative[db_tp_95th$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Duration

wt <- wilcox.test(
  db_tp_95th$duration[db_tp_95th$recovered == "YES"],
  db_tp_95th$duration[db_tp_95th$recovered == "NO"]
)

es <- cliff.delta(
  db_tp_95th$duration[db_tp_95th$recovered == "YES"],
  db_tp_95th$duration[db_tp_95th$recovered == "NO"]
)

comparisons_vector <- data.frame(
  variable = "3. Acc. precipitation (max)",
  dimension = "Duration",
  threshold = "P95",
  metric = "RT_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(db_tp_95th$duration[db_tp_95th$recovered == "YES"]),
  n_no = length(db_tp_95th$duration[db_tp_95th$recovered == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Frequency (no recovered cells - NOT RUN)

db_tp_95th_fr <- data.frame()

for (year in 1998:2019){
  print(year)
  year_n <- ree_eb_95_tp_dfjma %>%
    filter(yr == year) %>%
    group_by(yr,longitude,latitude) %>%
    filter(exceedance_no == max(exceedance_no)) %>%
    full_join(ree_eb_95_tp_dfjma %>%
                filter(yr == year-1) %>%
                group_by(yr,longitude,latitude) %>%
                filter(exceedance_no == max(exceedance_no)), 
              by = c("longitude","latitude")) %>%
    group_by(yr.x,longitude,latitude) %>%
    summarise(
      REE_pe = (exceedance_no.x + exceedance_no.y)
    ) %>%
    as.data.frame()
  
  db_tp_95th_fr <- rbind(db_tp_95th_fr,year_n)
}

db_fr_pr <- dam_rec %>%
  left_join(db_tp_95th_fr, 
            by = c("yr"="yr.x","longitude","latitude")) %>%
  filter(!is.na(REE_pe)) %>%
  as.data.frame()


#Compounded

ce <- db_t2m_95th %>%
  left_join(db_tp_5th, by = c("id","yr")) %>%
  left_join(db_tp_95th, by = c("id","yr")) %>%
  as.data.frame()
  
ce$ce <- rowSums(ce[,c("exceedance_no.x",
                 "exceedance_no.y",
                 "exceedance_no")],
                 na.rm = T
  )

wt <- wilcox.test(
  ce$ce[ce$recovered.x == "YES"],
  ce$ce[ce$recovered.x == "NO"]
)

es <- cliff.delta(
  ce$ce[ce$recovered.x == "YES"],
  ce$ce[ce$recovered.x == "NO"]
)

comparisons_vector <- data.frame(
  variable = "4. Compounded",
  dimension = "Interaction",
  threshold = "Compounded",
  metric = "CE_pe",
  wilcox = wt$statistic,
  p_value = wt$p.value, 
  cliff_delta = es$estimate, 
  low_int = unname(es$conf.int[1]),
  up_int = unname(es$conf.int[2]), 
  magnitude = es$magnitude,
  n_yes = length(ce$ce[ce$recovered.x == "YES"]),
  n_no = length(ce$ce[ce$recovered.x == "NO"])
)

comparisons_df <- rbind(comparisons_df,comparisons_vector)

# Plots

colors_es <- c ("#D55E00","#009E73","#0072B2","#CC79A7")
delta_plot <- ggplot(comparisons_df, aes(x=cliff_delta, y=variable, group=dimension, color=dimension)) +
  geom_vline(xintercept = 0,color="darkgrey")+
  geom_pointrange(aes(xmin=low_int, xmax=up_int),
                  position = position_dodge(0.3),size=1.2
  ) +
  theme_test() +
  scale_color_manual(values = colors_es)+
  xlab("Cliff's delta")+
  ylab("")+
  scale_y_discrete(limits=rev)+
  theme(legend.position = "top",
        text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  labs(colour = "")

pr_in <- map +
  geom_point(data = db_tp_95th %>%
               filter(recovered == "NO"),
             aes(longitude,latitude,color=duration),
             pch=18,size=2.5)+
  scale_color_viridis_b(option = "inferno",direction = -1,name="RTpe")

t2m_fr <- map +
  geom_point(data = db_fr %>%
               filter(recovered == "NO"),
             aes(longitude,latitude,color=REE_pe),
             pch=18,size=2.5)+
  scale_color_viridis_b(option = "inferno",direction = -1,name="FEpe")

t2m_in <- map +
  geom_point(data = db_t2m_95th %>%
               filter(recovered == "NO"),
             aes(longitude,latitude,color=intensity_cumulative),
             pch=18,size=2.5)+
  scale_color_viridis_b(option = "inferno",direction = -1,name="MCIpe")

dro_fr <- map +
  geom_point(data = db_fr_dro %>%
               filter(recovered == "NO"),
             aes(longitude,latitude,color=REE_pe),
             pch=18,size=2.5)+
  scale_color_viridis_b(option = "inferno",direction = -1,name="FEpe")

dro_in <- map +
  geom_point(data = db_tp_5th %>%
               filter(recovered == "NO"),
             aes(longitude,latitude,color=duration),
             pch=18,size=2.5)+
  scale_color_viridis_b(option = "inferno",direction = -1,name="RTpe")

