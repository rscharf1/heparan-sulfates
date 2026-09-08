library(ggplot2)
library(data.table)
library(dplyr)

dt <- fread("mouse_CBC_values.csv")

# Convert date
dt[, Date := as.Date(Date, format = "%m/%d/%y")]

# Unique mouse ID
# Patient names repeat across litters, so include Litter
dt[, Mouse := paste0("L", Litter, "_", Patient)]

# Normal ranges
# WBC:         1.8 - 10.7
# Lymphocytes: 0.9 - 9.3
# RBC:         6.36 - 9.42
# Hemoglobin:  11.0 - 15.1

setorder(dt, Mouse, Date)

dt[, Timepoint := seq_len(.N), by = Mouse]

dt[, Timepoint := factor(
  Timepoint,
  levels = c(1, 2),
  labels = c("8 weeks", "14 weeks")
)]

dt.long <- melt(
  dt,
  id.vars = c(
    "Date", "Litter", "Patient",
    "Mouse", "Genotype", "Timepoint"
  ),
  measure.vars = c(
    "WBC (K/uL)",
    "Ly (K/uL)",
    "RBC (M/uL)",
    "Hb (g/dL)"
  ),
  variable.name = "Parameter",
  value.name = "Value"
)

dt.long <- merge(dt.long, norm, by = "Parameter")

dt.long[, Group := fifelse(
  Mouse %in% c("L2_L", "L2_R"),
  "Leukemic KO",
  Genotype
)]

dt.long[, Group := factor(
  Group,
  levels = c("Het", "KO", "Leukemic KO")
)]

pdf("out.pdf", width = 8, height = 5)

ggplot(
  dt.long,
  aes(
    x = Timepoint,
    y = Value,
    color = Group,
    group = Mouse
  )
) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.5
  ) +

  geom_line(
    aes(
      linewidth = Group == "Leukemic KO",
      alpha = Group == "Leukemic KO"
    )
  ) +

  geom_point(
    aes(size = Group == "Leukemic KO")
  ) +

  facet_wrap(
    ~Parameter,
    scales = "free_y",
    nrow = 1
  ) +

  scale_color_manual(
    values = c(
      Het = "#1f78b4",
      KO = "#d95f02",
      `Leukemic KO` = "#e31a1c"
    )
  ) +

  scale_linewidth_manual(
    values = c(`FALSE` = 0.7, `TRUE` = 1),
    guide = "none"
  ) +

  scale_size_manual(
    values = c(`FALSE` = 2.5, `TRUE` = 4),
    guide = "none"
  ) +

  scale_alpha_manual(
    values = c(`FALSE` = 0.4, `TRUE` = 1),
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = NULL,
    color = NULL
  ) +

  theme_classic(base_size = 14) +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

dev.off()


