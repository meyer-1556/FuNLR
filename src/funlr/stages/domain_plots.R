# ============================================================================
# Integrated Domain Architecture Plotting (fungal NLR optimized)
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Load the required packages
# ----------------------------------------------------------------------------
required_packages <- c("dplyr", "ggplot2", "RColorBrewer", "colorspace", "tidyr", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("Missing R dependency:", pkg))
  library(pkg, character.only = TRUE)
}

# Explicit base-R bitmap devices avoid an optional ragg ABI dependency.
# Quartz writes bitmap files on macOS without requiring an X11 installation.
# Cairo supplies the corresponding headless device on Linux.
save_plot <- function(filename, ..., device = NULL, bg = "white") {
  bitmap_type <- if (Sys.info()[["sysname"]] == "Darwin") "quartz" else "cairo"
  ggplot2::ggsave(filename = filename, ..., device = grDevices::png, type = bitmap_type, bg = bg)
  if (!file.exists(filename) || file.info(filename)$size <= 0) {
    stop(paste("PNG device did not produce:", filename))
  }
}

# ----------------------------------------------------------------------------
# 2. Read data
# ----------------------------------------------------------------------------
cat("Reading domain hits...\n")
domains <- read.delim("all_domain_hits.tsv", stringsAsFactors = FALSE)

candidate_ids <- readLines("candidate_ids.txt")
tiered <- read.delim("tiered_candidates.tsv", sep = "\t", stringsAsFactors = FALSE)

# Build tier lookup from the tiered table
tier_lookup <- setNames(tiered$tier, tiered$protein_id)
target_ids <- intersect(candidate_ids, names(tier_lookup))

# Filter domains to Tiers 1‑3 and add tier column
domains_filtered <- domains %>%
  filter(protein_id %in% target_ids) %>%
  mutate(tier = tier_lookup[protein_id])

# ----------------------------------------------------------------------------
# 3. Helper function to strip taxonomic suffixes (e.g., "__class__Agaricomycetes")
# ----------------------------------------------------------------------------
strip_suffix <- function(domain) {
  gsub("__.*", "", domain)
}

# ----------------------------------------------------------------------------
# 4. Define domain categories and types (using stripped names)
# ----------------------------------------------------------------------------
domains_filtered <- domains_filtered %>%
  mutate(
    domain_clean = strip_suffix(domain),
    domain_type = case_when(
      # ---------- Effector domains ----------
      grepl("Goodbye", domain_clean, ignore.case = TRUE) ~ "Goodbye",
      grepl("HeLo", domain_clean, ignore.case = TRUE) & !grepl("HeLo-like", domain_clean, ignore.case = TRUE) ~ "HeLo",
      grepl("HeLo-like|HELL", domain_clean, ignore.case = TRUE) ~ "HeLo-like",
      grepl("HET", domain_clean, ignore.case = TRUE) & !grepl("HET-s", domain_clean, ignore.case = TRUE) ~ "HET",
      grepl("HET-s", domain_clean, ignore.case = TRUE) ~ "HET-s",
      grepl("TIR", domain_clean, ignore.case = TRUE) ~ "TIR",
      grepl("Patatin", domain_clean, ignore.case = TRUE) ~ "Patatin",
      grepl("PNP_UDP|PUP", domain_clean, ignore.case = TRUE) ~ "PNP_UDP",
      grepl("RelA_SpoT", domain_clean, ignore.case = TRUE) ~ "RelA_SpoT",
      grepl("Ses[AB]|ses[ab]", domain_clean, ignore.case = TRUE) ~ "Ses",
      grepl("CHAT", domain_clean, ignore.case = TRUE) ~ "CHAT",
      grepl("Crinkler", domain_clean, ignore.case = TRUE) ~ "Crinkler",
      grepl("SAM", domain_clean, ignore.case = TRUE) ~ "SAM",
      grepl("C2[ -]domain|PF00168", domain_clean, ignore.case = TRUE) ~ "C2",
      grepl("Peptidase_S8", domain_clean, ignore.case = TRUE) ~ "Peptidase_S8",
      grepl("CARD|Pyrin", domain_clean, ignore.case = TRUE) ~ "CARD/Pyrin",
      grepl("CC[_-]|Coiled[_-]coil", domain_clean, ignore.case = TRUE) ~ "CC",
      
      # ---------- NBD domains ----------
      grepl("NACHT", domain_clean, ignore.case = TRUE) ~ "NACHT",
      grepl("NB-ARC", domain_clean, ignore.case = TRUE) ~ "NB-ARC",
      grepl("AAA_16|AAA_22|AAA", domain_clean, ignore.case = TRUE) ~ "AAA",   # collapse AAA variants
      
      # ---------- Sensor domains ----------
      grepl("WD40|WD_40|WD-repeat", domain_clean, ignore.case = TRUE) ~ "WD40",
      grepl("ANK|Ank|ANKYRIN", domain_clean, ignore.case = TRUE) ~ "ANK",
      grepl("TPR", domain_clean, ignore.case = TRUE) ~ "TPR",
      grepl("HEAT", domain_clean, ignore.case = TRUE) ~ "HEAT",
      grepl("LRR|Leucine[ -]rich", domain_clean, ignore.case = TRUE) ~ "LRR",
      grepl("SPRY", domain_clean, ignore.case = TRUE) ~ "SPRY",
      grepl("PKinase|Kinase", domain_clean, ignore.case = TRUE) ~ "Kinase",
      grepl("C2H2|Zinc finger|ZINC_FINGER", domain_clean, ignore.case = TRUE) ~ "Zinc_finger",
      grepl("ZZ[_-]|ZZ-type", domain_clean, ignore.case = TRUE) ~ "ZZ",
      grepl("HMA|WRKY|LIM", domain_clean, ignore.case = TRUE) ~ "HMA/WRKY/LIM",
      
      # ---------- ASM domains ----------
      grepl("HRAM", domain_clean, ignore.case = TRUE) ~ "HRAM",
      grepl("PP[ -]motif|NLR07|NLR39", domain_clean, ignore.case = TRUE) ~ "PP",
      grepl("sigma", domain_clean, ignore.case = TRUE) ~ "sigma",               # also effector
      grepl("PUASM|NLR32", domain_clean, ignore.case = TRUE) ~ "PUASM",
      grepl("NLR05|NLR08|NLR22|NLR29|NLR44", domain_clean, ignore.case = TRUE) ~ "Basidio_ASM",
      grepl("BASS", domain_clean, ignore.case = TRUE) ~ "BASS",
      grepl("NLR17|NLR19|NLR34", domain_clean, ignore.case = TRUE) ~ "Lineage_ASM",
      
      # default
      TRUE ~ "Other"
    ),
    
    # Assign each domain type to a category (for border color)
    domain_category = case_when(
      domain_type %in% c("Goodbye", "HeLo", "HeLo-like", "HET", "HET-s", "TIR", 
                         "Patatin", "PNP_UDP", "RelA_SpoT", "Ses", "CHAT", 
                         "Crinkler", "SAM", "C2", "Peptidase_S8", "CARD/Pyrin", "CC") ~ "Effector",
      domain_type %in% c("NACHT", "NB-ARC", "AAA") ~ "NBD",
      domain_type %in% c("WD40", "ANK", "TPR", "HEAT", "LRR", "SPRY", 
                         "Kinase", "Zinc_finger", "ZZ", "HMA/WRKY/LIM") ~ "Sensor",
      domain_type %in% c("HRAM", "PP", "sigma", "PUASM", "Basidio_ASM", "BASS", "Lineage_ASM") ~ "ASM",
      TRUE ~ "Other"
    )
  )

# ----------------------------------------------------------------------------
# 5. Aggregate per-protein features
# ----------------------------------------------------------------------------
protein_features <- domains_filtered %>%
  group_by(protein_id, tier) %>%
  summarise(
    # NBD info
    nbd_types = paste(unique(domain_type[domain_category == "NBD"]), collapse = ";"),
    nbd_start = if (any(domain_category == "NBD")) min(start[domain_category == "NBD"]) else NA,
    nbd_end   = if (any(domain_category == "NBD")) max(end[domain_category == "NBD"]) else NA,
    # Sensor info
    sensor_types = paste(unique(domain_type[domain_category == "Sensor"]), collapse = ";"),
    sensor_start = if (any(domain_category == "Sensor")) min(start[domain_category == "Sensor"]) else NA,
    # ASM and effector presence
    has_asm = any(domain_category == "ASM"),
    has_effector = any(domain_category == "Effector"),
    # Protein length (max end)
    protein_length = tiered$protein_length[match(first(protein_id), tiered$protein_id)],
    .groups = "drop"
  ) %>%
  mutate(
    # Simplify NBD type: prefer NACHT or NB-ARC, else AAA or "None"
    nbd_simple = case_when(
      grepl("NACHT", nbd_types) ~ "NACHT",
      grepl("NB-ARC", nbd_types) ~ "NB-ARC",
      grepl("AAA", nbd_types) ~ "AAA",
      TRUE ~ "None"
    ),
    # Sensor class: choose first if multiple (could refine later)
    sensor_class = case_when(
      grepl("WD40", sensor_types) ~ "WD40",
      grepl("ANK", sensor_types) ~ "ANK",
      grepl("TPR", sensor_types) ~ "TPR",
      grepl("HEAT", sensor_types) ~ "HEAT",
      grepl("LRR", sensor_types) ~ "LRR",
      grepl("SPRY", sensor_types) ~ "SPRY",
      sensor_types == "" ~ "None",
      TRUE ~ "Other"
    ),
    # Order sanity: NBD start < sensor start (with buffer 10)
    order_ok = ifelse(!is.na(nbd_start) & !is.na(sensor_start),
                      nbd_start + 10 < sensor_start, NA)
  )

# ----------------------------------------------------------------------------
# 6. Generate colors for plots (same as before)
# ----------------------------------------------------------------------------
border_colors <- c(
  "Effector" = "#E41A1C",
  "NBD"      = "#377EB8",
  "Sensor"   = "#4DAF4A",
  "ASM"      = "#984EA3",
  "Other"    = "#999999"
)

# ----------------------------------------------------------------------------
# 7. Summary plots
# ----------------------------------------------------------------------------
cat("Generating summary plots...\n")

tier_order <- c("TIER_1A_HIGH_CONFIDENCE", "TIER_1B_NEEDS_REVIEW",
                "TIER_1C_CANONICAL_NON_SSFR_SENSOR", "TIER_2A_HIGH_PRIORITY_RESCUE",
                "TIER_2B_RESCUE_CANDIDATE", "TIER_2C_LOW_PRIORITY_FRAGMENT",
                "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE",
                "TIER_3A_INTEGRATED_DECOY", "TIER_3B_ARCHITECTURAL_VARIANT")
tier_order <- intersect(tier_order, unique(protein_features$tier))
protein_features$tier <- factor(protein_features$tier, levels = tier_order)

p_sensor <- protein_features %>%
  filter(sensor_class != "None") %>%
  ggplot(aes(x = tier, fill = sensor_class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Sensor class distribution by tier",
       x = "Tier", y = "Proportion", fill = "Sensor class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("sensor_distribution.png", p_sensor, width = 8, height = 5, dpi = 300)

p_nbd <- protein_features %>%
  filter(nbd_simple != "None") %>%
  ggplot(aes(x = tier, fill = nbd_simple)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "NBD subtype distribution by tier",
       x = "Tier", y = "Proportion", fill = "NBD type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("nbd_distribution.png", p_nbd, width = 8, height = 5, dpi = 300)

p_asm <- protein_features %>%
  mutate(asm_present = ifelse(has_asm, "ASM present", "No ASM")) %>%
  ggplot(aes(x = tier, fill = asm_present)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("ASM present" = "#984EA3", "No ASM" = "#CCCCCC")) +
  labs(title = "ASM presence by tier",
       x = "Tier", y = "Proportion", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("asm_presence.png", p_asm, width = 8, height = 5, dpi = 300)

rescue_priority <- tiered %>% select(protein_id, rescue_priority)
protein_features <- left_join(protein_features, rescue_priority, by = "protein_id")

p_scatter <- protein_features %>%
  filter(!is.na(rescue_priority)) %>%
  ggplot(aes(x = protein_length, y = rescue_priority, color = tier, shape = sensor_class)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Protein length vs rescue priority",
       x = "Protein length (aa)", y = "Rescue priority", color = "Tier", shape = "Sensor class") +
  theme_minimal()
save_plot("length_vs_priority.png", p_scatter, width = 8, height = 5, dpi = 300)

order_summary <- protein_features %>%
  filter(!is.na(order_ok)) %>%
  group_by(tier) %>%
  summarise(
    total = n(),
    order_ok_count = sum(order_ok, na.rm = TRUE),
    fraction = order_ok_count / total
  )

p_order <- ggplot(order_summary, aes(x = tier, y = fraction)) +
  geom_col(fill = "#377EB8") +
  geom_text(aes(label = sprintf("%d/%d", order_ok_count, total)), vjust = -0.5) +
  labs(title = "Order sanity (NBD before sensor)",
       x = "Tier", y = "Fraction with correct order") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("order_sanity.png", p_order, width = 8, height = 5, dpi = 300)

cat("Generating individual protein plots...\n")
ind_dir <- "individual"
dir.create(ind_dir, recursive = TRUE, showWarnings = FALSE)

hue_ranges <- list(
  Effector = c(0, 30),
  NBD      = c(200, 260),
  Sensor   = c(80, 140),
  ASM      = c(260, 320),
  Other    = c(0, 0)
)

get_category_colors <- function(category, n) {
  if (category == "Other") {
    return(rep("#CCCCCC", n))
  }
  hues <- seq(hue_ranges[[category]][1], hue_ranges[[category]][2], length.out = n)
  luminance <- rep(c(40, 70), length.out = n)
  grDevices::hcl(h = hues, c = 80, l = luminance, fixup = TRUE)
}

domain_types_per_category <- split(unique(domains_filtered$domain_type), 
                                   domains_filtered$domain_category[match(unique(domains_filtered$domain_type), domains_filtered$domain_type)])

fill_colors <- unlist(lapply(names(domain_types_per_category), function(cat) {
  types <- domain_types_per_category[[cat]]
  cols <- get_category_colors(cat, length(types))
  setNames(cols, types)
}))

for (pid in unique(domains_filtered$protein_id)) {
  pdata <- domains_filtered %>% filter(protein_id == pid)
  
  p <- ggplot(pdata) +
    geom_rect(aes(xmin = start, xmax = end, 
                  ymin = 0, ymax = 1, 
                  fill = domain_type,
                  color = domain_category)) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = border_colors) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 10)) +
    labs(title = paste0(pid, "\n", unique(pdata$tier))) +
    xlim(0, max(pdata$end) + 50)
  
  save_plot(filename = file.path(ind_dir, paste0(pid, ".png")),
         plot = p, width = 8, height = 1.5, device = "png", dpi = 300)
}

cat("Generating tier summary figure...\n")
protein_order <- domains_filtered %>%
  distinct(tier, protein_id) %>%
  mutate(tier = factor(tier, levels = tier_order)) %>%
  arrange(tier, protein_id) %>%
  pull(protein_id)

domains_summary <- domains_filtered %>%
  mutate(global_y = match(protein_id, protein_order))

p_summary <- ggplot(domains_summary) +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = global_y - 0.4, ymax = global_y + 0.4,
                fill = domain_type,
                color = domain_category)) +
  scale_fill_manual(values = fill_colors, name = "Domain Type") +
  scale_color_manual(values = border_colors, name = "Category") +
  scale_y_continuous(
    breaks = unique(domains_summary$global_y),
    labels = protein_order[unique(domains_summary$global_y)],
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~tier, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(size = 7, hjust = 1),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold")
  ) +
  labs(x = "Amino Acid Position", y = "Protein ID") +
  guides(
    fill = guide_legend(ncol = 3, byrow = TRUE, title.position = "top"),
    color = guide_legend(ncol = 5, title.position = "top")
  )

save_plot(filename = "nlr_domains_by_tier.png",
       plot = p_summary, width = 16, height = 12, device = "png", dpi = 300)

summary_stats <- domains_filtered %>%
  group_by(tier, protein_id) %>%
  summarise(
    domain_arch = paste(sort(unique(domain_type)), collapse = ";"),
    .groups = "drop"
  ) %>%
  count(tier, domain_arch, name = "count") %>%
  left_join(
    domains_filtered %>%
      distinct(tier, protein_id) %>%
      count(tier, name = "total_proteins"),
    by = "tier"
  ) %>%
  rename(architecture = domain_arch) %>%
  select(tier, total_proteins, architecture, count) %>%
  arrange(tier, desc(count), architecture)

write.table(summary_stats, file = "domain_architecture_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n✅ Domain plotting completed.\n")
writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")
