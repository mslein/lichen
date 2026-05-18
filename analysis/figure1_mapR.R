pacman::p_load(tidyverse, ggplot2, maps)

# =========================================================
# Panel A: Latitudinal diversity globe with split histograms
# No cross-hatching; legend matches actual fills
# Counts reported separately by system
# =========================================================

# ----------------
# 1. Read + combine data
# ----------------
coral_raw_npp <- read_csv("analysis/tidy data/coral_npp_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "coral")

coral_raw_gpp <- read_csv("analysis/tidy data/coral_gpp_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "coral")

coral_raw_r <- read_csv("analysis/tidy data/coral_r_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "coral")

lichen_raw_npp <- read_csv("analysis/tidy data/lichen_npp_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "lichen")

lichen_raw_gpp <- read_csv("analysis/tidy data/lichen_gpp_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "lichen")

lichen_raw_r <- read_csv("analysis/tidy data/lichen_r_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "lichen")

foram_raw_npp <- read_csv("analysis/tidy data/cell_npp_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "foraminifera")

foram_raw_gpp <- read_csv("analysis/tidy data/cell_gpp_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "foraminifera")

foram_raw_r <- read_csv("analysis/tidy data/cell_r_final.csv") %>%
  select(tpc_grp, latitude, species, temp, mol_gCmin) %>%
  mutate(system = "foraminifera")

all <- bind_rows(
  coral_raw_gpp, coral_raw_npp, coral_raw_r,
  lichen_raw_gpp, lichen_raw_npp, lichen_raw_r,
  foram_raw_gpp, foram_raw_npp, foram_raw_r
) %>%
  filter(!is.na(latitude)) %>%
  separate(
    tpc_grp,
    into = c("rep", "study_id"),
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(
    system = factor(
      system,
      levels = c("lichen", "foraminifera", "coral")
    )
  )

# ----------------
# 2. Counts by system
# ----------------
counts_df <- all %>%
  group_by(system) %>%
  summarise(
    studies = n_distinct(study_id),
    observations = n(),
    .groups = "drop"
  )

get_count <- function(sys, col) {
  counts_df %>%
    filter(system == sys) %>%
    pull({{ col }})
}

count_label <- paste0(
  "Lichen: ",
  get_count("lichen", studies),
  " studies, ",
  get_count("lichen", observations),
  " obs\n",
  "Foraminifera: ",
  get_count("foraminifera", studies),
  " studies, ",
  get_count("foraminifera", observations),
  " obs\n",
  "Coral: ",
  get_count("coral", studies),
  " studies, ",
  get_count("coral", observations),
  " obs"
)

# ----------------
# 3. Latitude projection
# ----------------
deg2rad <- function(x) x * pi / 180

lat_limits_deg <- c(-90, 90)
lat_breaks_deg <- seq(-60, 60, by = 30)

y_limits     <- sin(deg2rad(lat_limits_deg))
lat_breaks_y <- sin(deg2rad(lat_breaks_deg))

all <- all %>%
  mutate(
    lat_rad = deg2rad(latitude),
    y_proj  = sin(lat_rad)
  )

# ----------------
# 4. Globe projection
# ----------------
world_raw <- map_data("world")

center_lon_deg <- -60
center_lon_rad <- deg2rad(center_lon_deg)

world_proj <- world_raw %>%
  mutate(
    lon_rad = deg2rad(long),
    lat_rad = deg2rad(lat),
    x = cos(lat_rad) * sin(lon_rad - center_lon_rad),
    y = sin(lat_rad),
    visible = cos(lat_rad) * cos(lon_rad - center_lon_rad) >= 0
  ) %>%
  filter(
    visible,
    y >= y_limits[1],
    y <= y_limits[2]
  )

circle <- tibble(
  theta = seq(0, 2 * pi, length.out = 361),
  x = cos(theta),
  y = sin(theta)
)

# ----------------
# 5. Histogram binning
# ----------------
nbins <- 40
breaks_y <- seq(y_limits[1], y_limits[2], length.out = nbins + 1)

hist_df <- all %>%
  mutate(
    bin = cut(
      y_proj,
      breaks = breaks_y,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(system, bin, name = "count") %>%
  mutate(
    ymin = breaks_y[as.integer(bin)],
    ymax = breaks_y[as.integer(bin) + 1]
  )

max_half_width <- 0.9
count_max_all  <- max(hist_df$count, na.rm = TRUE)

hist_df <- hist_df %>%
  mutate(width = count / count_max_all * max_half_width)

hist_left <- hist_df %>%
  filter(system %in% c("lichen", "foraminifera")) %>%
  mutate(
    xmin = -width,
    xmax = 0
  )

hist_right <- hist_df %>%
  filter(system == "coral") %>%
  mutate(
    xmin = 0,
    xmax = width
  )

count_breaks <- pretty(c(0, count_max_all), n = 4)
x_breaks_pos <- count_breaks / count_max_all * max_half_width
x_breaks     <- c(-rev(x_breaks_pos[-1]), x_breaks_pos)
x_labels     <- round(abs(c(-rev(count_breaks[-1]), count_breaks)))

# ----------------
# 6. Plot
# ----------------
border_col <- "grey70"

p_globe_hist <- ggplot() +
  
  geom_polygon(
    data = world_proj,
    aes(x = x, y = y, group = group),
    fill = "grey85",
    color = border_col,
    linewidth = 0.2,
    alpha = 0.3
  ) +
  
  geom_path(
    data = circle,
    aes(x = x, y = y),
    colour = border_col,
    linewidth = 0.2
  ) +
  
  geom_rect(
    data = hist_left,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = system
    ),
    colour = "grey35",
    linewidth = 0.2,
    alpha = 0.95
  ) +
  
  geom_rect(
    data = hist_right,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = system
    ),
    colour = "grey35",
    linewidth = 0.2,
    alpha = 0.95
  ) +
  
  annotate(
    "text",
    x = -0.95,
    y = -0.92,
    label = count_label,
    hjust = 0,
    vjust = 0,
    size = 3.2,
    colour = "grey20"
  ) +
  
  scale_fill_manual(
    values = c(
      lichen = "#4D4D4D",
      foraminifera = "#969696",
      coral = "#D9D9D9"
    ),
    labels = c(
      lichen = "Lichen",
      foraminifera = "Foraminifera",
      coral = "Coral"
    ),
    name = "System"
  ) +
  
  coord_equal(
    xlim = c(-1, 1),
    ylim = y_limits,
    expand = FALSE
  ) +
  
  scale_y_continuous(
    limits = y_limits,
    breaks = lat_breaks_y,
    labels = lat_breaks_deg
  ) +
  
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = x_breaks,
    labels = x_labels,
    name = "Count"
  ) +
  
  labs(y = "Latitude (°)") +
  
  theme_bw() +
  
  theme(
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    
    legend.position = c(0.03, 0.18),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "white",
      colour = "grey80"
    )
  )

p_globe_hist

# ----------------
# 7. Save panel
# ----------------
ggsave(
  "figures/fig1.png",
  p_globe_hist,
  width = 7,
  height = 6,
  dpi = 600
)