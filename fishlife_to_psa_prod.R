# Setup
if (!require("remotes")) install.packages("remotes")
if (!require("FishLife")) remotes::install_github("James-Thorson/FishLife")
if (!require("pacman"))  install.packages("pacman")
pacman::p_load(MASS, stringr, dplyr, ggplot2, ggrepel, patchwork)
library(FishLife)

# Input
df = read.csv("test_data.csv", stringsAsFactors = FALSE)

# Categorisation functions

# Productivity
cat_r     = function(x) ifelse(is.na(x), NA, ifelse(x > 0.5,  "high", ifelse(x >= 0.16, "mod", "low")))
cat_tmax  = function(x) ifelse(is.na(x), NA, ifelse(x < 10,   "high", ifelse(x <= 30,   "mod", "low")))
cat_lmax  = function(x) ifelse(is.na(x), NA, ifelse(x < 60,   "high", ifelse(x <= 150,  "mod", "low")))
cat_k     = function(x) ifelse(is.na(x), NA, ifelse(x > 0.25, "high", ifelse(x >= 0.15, "mod", "low")))
cat_m     = function(x) ifelse(is.na(x), NA, ifelse(x > 0.4,  "high", ifelse(x >= 0.2,  "mod", "low")))
cat_fec   = function(x) ifelse(is.na(x), NA, ifelse(x > 100000, "high", ifelse(x >= 1000, "mod", "low")))
cat_breed = function(x) ifelse(is.na(x), NA, ifelse(x == 0, "high", ifelse(x <= 3, "mod", ifelse(x <= 14, "low", NA))))
cat_rec   = function(x) ifelse(is.na(x), NA, ifelse(x == "highfreq", "high", ifelse(x == "modfreq", "mod", ifelse(x == "lowfreq", "low", NA))))
cat_tmat  = function(x) ifelse(is.na(x), NA, ifelse(x > 4, "high", ifelse(x >= 2, "mod", "low")))
cat_troph = function(x) ifelse(is.na(x), NA, ifelse(x < 2.5, "high", ifelse(x <= 3.5, "mod", "low")))

# Susceptibility
cat_area_over  = function(x) ifelse(is.na(x), NA, ifelse(x > 0.5,  "high", ifelse(x >= 0.25, "mod", "low")))
cat_geog_conc  = function(x) ifelse(is.na(x), NA, ifelse(x == "high", "high", ifelse(x == "mod", "mod", ifelse(x == "low", "low", NA))))
cat_vert_over  = function(x) ifelse(is.na(x), NA, ifelse(x == "high", "high", ifelse(x == "mod", "mod", ifelse(x == "low", "low", NA))))
cat_seas_migr  = function(x) ifelse(is.na(x), NA, ifelse(x == "increase_fish", "high", ifelse(x == "no_effect", "mod", ifelse(x == "decrease_fish", "low", NA))))
cat_school     = function(x) ifelse(is.na(x), NA, ifelse(x == "increase_fish", "high", ifelse(x == "no_effect", "mod", ifelse(x == "decrease_fish", "low", NA))))
cat_morph      = function(x) ifelse(is.na(x), NA, ifelse(x == "high_selec", "high", ifelse(x == "mod_selec", "mod", ifelse(x == "low_selec", "low", NA))))
cat_desire     = function(x) ifelse(is.na(x), NA, ifelse(x == "high", "high", ifelse(x == "mod", "mod", ifelse(x == "low", "low", NA))))
cat_mng_strat  = function(x) ifelse(is.na(x), NA, ifelse(x == "no_strat", "high", ifelse(x == "reactive", "mod", ifelse(x == "proactive", "low", NA))))
cat_f_over_m   = function(x) ifelse(is.na(x), NA, ifelse(x > 1, "high", ifelse(x >= 0.5, "mod", "low")))
cat_b_over_b0  = function(x) ifelse(is.na(x), NA, ifelse(x < 0.25, "high", ifelse(x <= 0.4, "mod", "low")))
cat_surv_prob  = function(x) ifelse(is.na(x), NA, ifelse(x < 0.33, "high", ifelse(x <= 0.67, "mod", "low")))
cat_hat_impact = function(x) ifelse(is.na(x), NA, ifelse(x == "high", "high", ifelse(x == "mod", "mod", ifelse(x == "low", "low", NA))))

categorize_attributes = function(df) {
  for (col in setdiff(names(df), "species")) {
    fn = paste0("cat_", col)
    if (exists(fn, mode = "function")) df[[col]] = get(fn)(df[[col]])
  }
  df
}

create_species_list = function(df) {
  attribute_cols = setdiff(names(df), "species")
  lapply(setNames(unique(df$species), unique(df$species)), function(sp) {
    row = df[df$species == sp, ][1, ]
    out = data.frame(attribute = attribute_cols, low = 0, mod = 0, high = 0, weight = 1)
    for (attr in attribute_cols) {
      val = row[[attr]]
      if (!is.na(val)) out[[val]][out$attribute == attr] = 1
      else             out$weight[out$attribute == attr] = 0
    }
    out
  })
}

vul_class = function(x) ifelse(x >= 2.2, "High", ifelse(x > 1.8, "Moderate", "Low"))

# FishLife lookup

fishlife_for_many = function(species, n_samples = 10000) {
  col_names = c("Loo","K","Winfinity","tmax","tm","M","Lm","Temperature",
                "ln_var","rho","ln_MASPS","ln_margsd","h","logitbound_h",
                "ln_Fmsy_over_M","ln_Fmsy","ln_r","r","ln_G","G")
  back_transform = function(v) c(exp(v[1:7]), v[8], v[9:20])

  lapply(setNames(sort(unique(species)), sort(unique(species))), function(sciname) {
    genus = word(sciname, 1)
    sp    = word(sciname, 2, length(strsplit(sciname, " ")[[1]]))
    sp    = ifelse(sp == "spp", "predictive", sp)

    res = try(FishLife::Search_species(Genus = genus, Species = sp), silent = TRUE)
    spp = try(FishLife::Plot_taxa(res$match_taxonomy), silent = TRUE)

    if (inherits(res, "try-error") || inherits(spp, "try-error")) {
      message("FishLife lookup failed: ", sciname); return(NULL)
    }
    pred    = spp[[1]]
    samples = as.data.frame(t(apply(
      MASS::mvrnorm(n_samples, pred$Mean_pred, pred$Cov_pred), 1, back_transform
    )))
    colnames(samples) = col_names
    list(samples = samples, matched_taxon = paste(res$match_taxonomy, collapse = " > "))
  })
}

fl_productivity_probs = function(s) {
  score = function(vals, fn) {
    cats = fn(vals)
    c(low = mean(cats == "low", na.rm = TRUE),
      mod = mean(cats == "mod", na.rm = TRUE),
      high = mean(cats == "high", na.rm = TRUE))
  }
  list(
    r    = score(s$r,    function(x) ifelse(x > 0.5,  "high", ifelse(x >= 0.16, "mod", "low"))),
    K    = score(s$K,    function(x) ifelse(x > 0.25, "high", ifelse(x >= 0.15, "mod", "low"))),
    tmax = score(s$tmax, function(x) ifelse(x < 10,   "high", ifelse(x <= 30,   "mod", "low"))),
    tm   = score(s$tm,   function(x) ifelse(x < 2,    "low",  ifelse(x <= 4,    "mod", "high"))),
    M    = score(s$M,    function(x) ifelse(x > 0.4,  "high", ifelse(x >= 0.2,  "mod", "low")))
  )
}

# Background colour grid for PSA plot

df_col = {
  xc = seq(0, 1, length.out = 200); yc = seq(0, 1, length.out = 200)
  g  = cbind(expand.grid(xcolor = xc, ycolor = yc),
             expand.grid(x = seq(3, 1, length.out = 200), y = seq(1, 3, length.out = 200)))
  g$zcolor = g$xcolor^2 + g$ycolor^2
  g
}

# Run analysis

df_cat       = categorize_attributes(df)
species_list = create_species_list(df_cat)
fl_samples   = fishlife_for_many(unique(df$species))
fl_attrs     = c("r", "K", "tmax", "tm", "M")

# Replace FishLife attributes with probabilistic scores
for (sp in names(fl_samples)) {
  if (is.null(fl_samples[[sp]])) next
  probs = fl_productivity_probs(fl_samples[[sp]]$samples)
  for (attr in fl_attrs) {
    if (!attr %in% species_list[[sp]]$attribute) next
    idx = which(species_list[[sp]]$attribute == attr)
    species_list[[sp]][idx, c("low","mod","high")] = probs[[attr]][c("low","mod","high")]
    species_list[[sp]]$weight[idx] = 1
  }
}

# Export FishLife probability scores used for productivity
fishlife_scores_df = do.call(rbind, lapply(names(fl_samples), function(sp) {
  if (is.null(fl_samples[[sp]])) return(NULL)
  probs = fl_productivity_probs(fl_samples[[sp]]$samples)
  do.call(rbind, lapply(fl_attrs, function(attr) {
    p = probs[[attr]]
    data.frame(species = sp, attribute = attr,
               prob_low = p["low"], prob_mod = p["mod"], prob_high = p["high"],
               expected_score = p["low"] * 1 + p["mod"] * 2 + p["high"] * 3,
               row.names = NULL)
  }))
}))

# Double weight for r
for (sp in names(species_list))
  species_list[[sp]]$weight[species_list[[sp]]$attribute == "r"] = 2

# PSA bootstrap
num_samples   = 999
num_prod_attr = 10
num_susc_attr = 12

vulnerability_scores = lapply(names(species_list), function(species) {
  d = species_list[[species]]
  d[, c("low","mod","high")] = lapply(d[, c("low","mod","high")], as.numeric)

  sw_prod = sum(d$weight[1:num_prod_attr], na.rm = TRUE)
  sw_susc = sum(d$weight[(num_prod_attr+1):(num_prod_attr+num_susc_attr)], na.rm = TRUE)

  sample_col = function(i) {
    p = as.numeric(d[i, c("high","mod","low")]); w = d$weight[i]
    if (!is.na(w) && w > 0 && sum(p, na.rm=TRUE) > 0)
      w * sample(c(3,2,1), num_samples, replace=TRUE, prob=p)
    else rep(0, num_samples)
  }

  prod_scores = rowSums(sapply(1:num_prod_attr, sample_col)) / sw_prod
  susc_scores = rowSums(sapply(seq(num_prod_attr+1, num_prod_attr+num_susc_attr), sample_col)) / sw_susc
  vuln        = sqrt((3 - prod_scores)^2 + (susc_scores - 1)^2)

  list(productivity = prod_scores, susceptibility = susc_scores, vulnerability = vuln)
})
names(vulnerability_scores) = names(species_list)

# Results data frame
results_df = do.call(rbind, lapply(names(vulnerability_scores), function(sp) {
  v = vulnerability_scores[[sp]]
  data.frame(
    species             = sp,
    matched_taxon       = if (is.null(fl_samples[[sp]])) NA else fl_samples[[sp]]$matched_taxon,
    has_fishlife        = !is.null(fl_samples[[sp]]),
    mean_productivity   = mean(v$productivity),
    lci_productivity    = quantile(v$productivity,  0.025),
    uci_productivity    = quantile(v$productivity,  0.975),
    sd_productivity     = sd(v$productivity),
    mean_susceptibility = mean(v$susceptibility),
    lci_susceptibility  = quantile(v$susceptibility, 0.025),
    uci_susceptibility  = quantile(v$susceptibility, 0.975),
    sd_susceptibility   = sd(v$susceptibility),
    mean_vulnerability  = mean(v$vulnerability),
    lci_vulnerability   = quantile(v$vulnerability,  0.025),
    uci_vulnerability   = quantile(v$vulnerability,  0.975),
    sd_vulnerability    = sd(v$vulnerability),
    stringsAsFactors    = FALSE
  )
}))
results_df$vul_category = vul_class(results_df$mean_vulnerability)

# Save results
write.csv(results_df,         "psa_results.csv",         row.names = FALSE)
write.csv(fishlife_scores_df, "fishlife_prod_scores.csv", row.names = FALSE)
write.csv(
  subset(fl_coverage <- data.frame(
    species      = unique(df$species),
    has_fishlife = sapply(unique(df$species), function(sp) !is.null(fl_samples[[sp]]))
  ), !has_fishlife),
  "no_fishlife.csv", row.names = FALSE
)

# Plots

prod_susc_plot =
  ggplot() +
  geom_tile(data = df_col, aes(x, y, fill = zcolor)) +
  geom_linerange(data = results_df,
                 aes(y = mean_susceptibility, x = mean_productivity,
                     ymin = lci_susceptibility, ymax = uci_susceptibility),
                 alpha = 0.4, linewidth = 1) +
  geom_linerange(data = results_df,
                 aes(y = mean_susceptibility, x = mean_productivity,
                     xmin = uci_productivity,  xmax = lci_productivity),
                 alpha = 0.4, linewidth = 1) +
  geom_point(data = results_df,
             aes(y = mean_susceptibility, x = mean_productivity, shape = has_fishlife),
             alpha = 0.7, size = 2) +
  geom_text_repel(data = results_df,
                  aes(label = species, y = mean_susceptibility, x = mean_productivity),
                  fontface = "italic", force = 50, size = 2.3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "FishLife", "FALSE" = "CSV only"), name = "") +
  scale_fill_gradientn(
    colors = c("forestgreen","green3","green2","greenyellow","orange3","red","red2","red3","red4"),
    guide  = "none") +
  labs(x = "Productivity", y = "Susceptibility") +
  xlim(1, 3) + ylim(1, 3) + scale_x_reverse() +
  coord_cartesian(expand = FALSE) +
  theme_test() +
  theme(axis.text = element_text(color = "black"), legend.position = "bottom")

dens_prod_plot =
  ggplot() +
  geom_density(data = results_df, aes(x = mean_productivity), fill = "grey", color = "grey") +
  scale_x_continuous(limits = c(1, 3)) + scale_x_reverse() +
  coord_cartesian(expand = FALSE) + theme_void()

dens_susc_plot =
  ggplot() +
  geom_density(data = results_df, aes(y = mean_susceptibility), fill = "grey", color = "grey") +
  scale_y_continuous(limits = c(1, 3)) +
  coord_cartesian(expand = FALSE) + theme_void()

psa_main =
  ((dens_prod_plot / prod_susc_plot) + plot_layout(heights = c(0.15, 1)) |
     (plot_spacer() / dens_susc_plot) + plot_layout(heights = c(0.15, 1))) +
  plot_layout(widths = c(1, 0.15))

hist_plot =
  ggplot() +
  geom_histogram(data = results_df, aes(x = mean_vulnerability, fill = vul_category),
                 binwidth = 0.03) +
  scale_fill_manual(values = c("greenyellow","orange","red2"),
                    breaks = c("Low","Moderate","High")) +
  labs(x = "Vulnerability", y = "No. species", fill = "") +
  scale_x_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"), legend.position = "bottom")

prod_density_df = do.call(rbind, lapply(names(vulnerability_scores), function(sp)
  data.frame(species = sp, productivity = vulnerability_scores[[sp]]$productivity)
))

prod_density_plot =
  ggplot(prod_density_df, aes(x = productivity, fill = species, color = species)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(data = results_df, aes(xintercept = mean_productivity, color = species),
             linetype = "dashed", linewidth = 0.6) +
  labs(x = "Productivity score", y = "Density", fill = "", color = "") +
  scale_x_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.position = "bottom")

plots_psa = psa_main / (hist_plot + plot_spacer() + plot_layout(widths = c(1, 0.2))) +
  plot_layout(heights = c(1, 0.5))

print(plots_psa)
print(prod_density_plot)