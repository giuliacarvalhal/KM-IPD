# =============================================================================
# SURVIVAL ANALYSIS SCRIPT
# =============================================================================

# 1. DEPENDENCIES AND SETUP ===================================================

# Install pacman if not available
if (!require(pacman)) install.packages("pacman")

# Load required packages
p_load(
  tidyverse,    # Data manipulation and visualization
  survival,     # Survival analysis
  ggsurvfit,    # Survival curve plotting
  gtsummary,    # Summary tables
  survRM2,      # Restricted mean survival time
  rms,          # Regression modeling strategies
  boot,         # Bootstrap methods
  funtimes,     # Time series analysis
  tseries,      # Time series analysis
  cowplot,      # Plot themes
  scales,       # Scale functions
  janitor,      # Data cleaning
  patchwork , 
  survival, 
  ggsurvfit,
  ggplot2,
  ggtext,
  cowplot
)

# 2. DATA LOADING AND PREPARATION ==============================================

# Load and merge patient data from Kaplan-Meier curves
meta_ipd <- read_csv("/Users/yourpathway") |>
  mutate(
    # Rename survival and status variables
    month = `Survival time`,
    status = Status,
    
    # Convert 'treated' from string to binary (Treatment = 1, Placebo = 0)
    treated = case_when(
      treated == "treatment" ~ 1,
      treated == "placebo" ~ 0,
      TRUE ~ NA_real_  # catch any unexpected values
    )
  )

# Display the data
cat("=== SURVIVAL DATA SUMMARY ===\n")
print(meta_ipd)
cat("Number of observations:", nrow(meta_ipd), "\n")
cat("Number of events:", sum(meta_ipd$status), "\n")
cat("Placebo group:", sum(meta_ipd$treated == 0), "\n")
cat("Treatment group:", sum(meta_ipd$treated == 1), "\n\n")

# 3. COX PROPORTIONAL HAZARDS MODEL ==========================================

cat("=== COX PROPORTIONAL HAZARDS ANALYSIS ===\n")

# Fit Cox model
cox_model <- coxph(Surv(month, status) ~ treated, data = meta_ipd)

# Display results using gtsummary
cox_model |>
  tbl_regression(exp = TRUE) |>
  print()

# Test proportional hazards assumption
cat("\n--- Proportional Hazards Assumption Test ---\n")
ph_test <- cox.zph(cox_model)
print(ph_test)

# Plot Schoenfeld residuals
cat("Plotting Schoenfeld residuals...\n")
plot(ph_test, main = "Schoenfeld Residuals Test for Proportional Hazards")

# 4. CURVA CUMULATIVE INCIDENCE + RISK TABLE ----------------------------------
km_fit <- survfit(Surv(month, status) ~ treated, data = meta_ipd)

# Gráfico
km_plot <- km_fit |>
  ggsurvfit(type = "risk", size = 1) +
  add_confidence_interval() +
  scale_color_manual(values = c("#0067BA", "#F9CB76")) +
  scale_fill_manual(values = c("#0067BA", "#F9CB76")) +
  
  # Caixas maiores nas legendas
  annotate("rect", xmin = 39.5, xmax = 40.2, ymin = 0.018, ymax = 0.032, fill = "#F9CB76", color = NA) +
  annotate("text", x = 40.3, y = 0.025, size = 4, hjust = 0, label = "DPP-1 Inhibitor") +
  annotate("rect", xmin = 39.5, xmax = 40.2, ymin = 0.058, ymax = 0.072, fill = "#0067BA", color = NA) +
  annotate("text", x = 40.3, y = 0.065, size = 4, hjust = 0, label = "Placebo") +
  guides(color = "none", fill = "none") +
  
  # Eixos ajustados e invertidos
  coord_cartesian(xlim = c(0, 52), ylim = c(1, 0)) +
  scale_y_continuous(
    trans = "reverse",  # flips the axis labels (not the curve)
    breaks = seq(0, 1, 0.2),
    labels = rev(seq(0, 100, 20)),  # shows labels from 100 to 0
    expand = c(0, 0)
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 52, 13)) +
  
  # Anotações com ggtext (Markdown)
  annotate("richtext", 
           label = "**DPP-1 Inhibitor versus Placebo: HR 0.79 (95% CI 0.71 to 0.88; P < 0.001)**", 
           x = 1, y = 0.9, hjust = 0, vjust = 1, size = 4.5, fill = NA, label.color = NA) +
  
  # Títulos com ggtext
  labs(
    x = "\nTime to first bronchiectasis exacerbation (weeks)\n", 
    y = "\nPercentage of patients with no exacerbations (%)\n",
    #title = "Kaplan-Meier Survival Curves"
  ) +
  theme_cowplot() +
  theme(
    axis.title.y = element_markdown(),  # negrito com ggtext
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

# Adiciona tabela de risco
km_plot_risk <- km_plot +
  add_risktable(
    size = 3,
    color = "#333333",
    risktable.position = "below"
  )

# Exibe
print(km_plot_risk)

# Salva como TIFF
ggsave(
  filename = "/Users/yourpathwaytosavefile",
  plot = km_plot_risk,
  device = "tiff",
  dpi = 600,
  width = 10,
  height = 8, # um pouco mais alto por conta da tabela
  units = "in",
  compression = "lzw"
)


# 5. TIME-VARYING HAZARD RATIO ANALYSIS =====================================

cat("\n=== TIME-VARYING HAZARD RATIO ANALYSIS ===\n")

# Fit Cox model with time-varying coefficient
cox_tv <- coxph(
  Surv(month, status) ~ treated + tt(treated),
  data = meta_ipd,
  tt = function(x, t, ...) x * log(t)
)

cat("Time-varying Cox model summary:\n")
summary(cox_tv)

# Extract coefficients and standard errors
b <- coef(cox_tv)           
se <- sqrt(diag(vcov(cox_tv)))  

# Generate time grid (in weeks)
times <- seq(1, 52, by = 1)  

# Calculate time-varying HR and approximate CI using delta method
hr <- exp(b[1] + b[2] * log(times))
var_hr_log <- se[1]^2 + (log(times))^2 * se[2]^2 + 2*vcov(cox_tv)[1,2]*log(times)
hr_low <- exp((b[1] + b[2]*log(times)) - 1.96*sqrt(var_hr_log))
hr_high <- exp((b[1] + b[2]*log(times)) + 1.96*sqrt(var_hr_log))

# Create data frame for plotting
hr_df <- data.frame(
  time = times, 
  hr = hr, 
  low = hr_low, 
  high = hr_high
)

# Plot time-varying hazard ratio
tv_hr_plot <- ggplot(hr_df, aes(x = time, y = hr)) +
  geom_line(size = 1, color = "darkblue") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Time-Varying Hazard Ratio",
    x = "Weeks Since Randomization", 
    y = "Hazard Ratio (treated vs placebo)",
    subtitle = "Solid line: point estimate; Shaded area: 95% CI"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(tv_hr_plot)

# Salva como TIFF
ggsave(
  filename = "/Users/yourpathwaytosavefile",
  plot = tv_hr_plot,
  device = "png",
  dpi = 600,
  width = 10,
  height = 8, # um pouco mais alto por conta da tabela
  units = "in",
)

# 6. RESTRICTED MEAN SURVIVAL TIME (RMST) ANALYSIS ==========================

cat("\n=== RESTRICTED MEAN SURVIVAL TIME (RMST) ANALYSIS ===\n")

# Define the time intervals (3 equal parts), with final interval ending at t_max
t_max <- max(meta_ipd$month)
breaks <- c(0, 17, 34, t_max)
intervals <- list(c(breaks[1], breaks[2]), c(breaks[2], breaks[3]), c(breaks[3], breaks[4]))
interval_labels <- c("0–17", "17–34", "34–52")

# Function to calculate RMST for each interval
calc_rmst_interval <- function(time, status, arm, start, end, label) {
  time_int <- pmin(time, end) - start
  status_int <- ifelse(time > end, 0, status)
  valid <- time > start
  time_sub <- time_int[valid]
  status_sub <- status_int[valid]
  arm_sub <- arm[valid]
  
  if (sum(valid) < 2 || length(unique(arm_sub)) < 2) {
    cat("Interval:", label, "weeks - Insufficient data for analysis\n\n")
    return(list(rmst_res = NULL, interval = label, valid = FALSE))
  }
  
  rmst_res <- tryCatch({
    rmst2(time = time_sub, status = status_sub, arm = arm_sub, tau = end - start)
  }, error = function(e) {
    cat("Error calculating RMST for interval", label, ":", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(rmst_res)) {
    return(list(rmst_res = NULL, interval = label, valid = FALSE))
  }
  
  arm0_rmst <- rmst_res$RMST.arm0$rmst["Est."]
  arm0_se <- rmst_res$RMST.arm0$rmst["se"]
  arm1_rmst <- rmst_res$RMST.arm1$rmst["Est."]
  arm1_se <- rmst_res$RMST.arm1$rmst["se"]
  
  diff_rmst <- NA; lower_ci <- NA; upper_ci <- NA; p_value <- NA
  if (!is.null(rmst_res$unadjusted.result)) {
    if ("Difference" %in% colnames(rmst_res$unadjusted.result)) {
      diff_rmst <- rmst_res$unadjusted.result["Est.", "Difference"]
      lower_ci <- rmst_res$unadjusted.result["Lower .95", "Difference"]
      upper_ci <- rmst_res$unadjusted.result["Upper .95", "Difference"]
      p_value <- rmst_res$unadjusted.result["p", "Difference"]
    } else {
      diff_rmst <- arm1_rmst - arm0_rmst
    }
  }
  
  if (is.na(diff_rmst) || is.na(lower_ci) || is.na(upper_ci)) {
    diff_rmst <- arm1_rmst - arm0_rmst
    se_diff <- sqrt(arm0_se^2 + arm1_se^2)
    lower_ci <- diff_rmst - 1.96 * se_diff
    upper_ci <- diff_rmst + 1.96 * se_diff
  }
  
  if (is.na(p_value)) {
    se_diff <- sqrt(arm0_se^2 + arm1_se^2)
    z_score <- diff_rmst / se_diff
    p_value <- 2 * (1 - pnorm(abs(z_score)))
  }
  
  cat("--- Interval:", label, "weeks ---\n")
  cat("RMST (Placebo):", round(arm0_rmst, 3), "±", round(arm0_se, 3), "\n")
  cat("RMST (Treated):", round(arm1_rmst, 3), "±", round(arm1_se, 3), "\n")
  cat("RMST Difference (Treated - Placebo):", round(diff_rmst, 3), "\n")
  cat("95% CI:", round(lower_ci, 3), "-", round(upper_ci, 3), "\n")
  cat("p-value:", signif(p_value, 3), "\n\n")
  
  list(
    rmst_res = rmst_res,
    interval = label,
    arm0_rmst = arm0_rmst,
    arm0_se = arm0_se,
    arm1_rmst = arm1_rmst,
    arm1_se = arm1_se,
    diff_rmst = diff_rmst,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    p_value = p_value,
    valid = TRUE
  )
}

# Apply the function to each interval
rmst_results <- purrr::map2(intervals, interval_labels, ~calc_rmst_interval(meta_ipd$month, meta_ipd$status, meta_ipd$treated, .x[1], .x[2], .y))

# Prepare data for RMST plot
df_plot <- map_df(rmst_results, function(res) {
  if (!res$valid) return(NULL)
  tibble(
    interval = res$interval,
    arm = c("Placebo", "DPP-1 Inhibitor"),
    rmst = c(res$arm0_rmst, res$arm1_rmst),
    rmst_se = c(res$arm0_se, res$arm1_se)
  )
})

if (nrow(df_plot) > 0) {
  cat("RMST data for plotting:\n")
  print(df_plot)
  
  rmst_plot <- ggplot(df_plot, aes(x = interval, y = rmst, fill = arm)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = rmst - 1.96 * rmst_se, ymax = rmst + 1.96 * rmst_se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    labs(
      title = "RMST by Time Interval and Treatment Group",
      x = "Time Interval (weeks)",
      y = "RMST (weeks)",
      fill = "Group"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Placebo" = "#0067BA", "DPP-1 Inhibitor" = "#F9CB76")) +
    theme(plot.title = element_text(face = "bold"))
  
  print(rmst_plot)
  
  # Salva como TIFF
  ggsave(
    filename = "/Users/yourpathwaytosavefile",
    plot = rmst_plot,
    device = "png",
    dpi = 600,
    width = 10,
    height = 8, # um pouco mais alto por conta da tabela
    units = "in",
  )
  
} else {
  cat("Unable to generate plot: no valid interval data.\n")
}

# Summary table
df_results <- map_df(rmst_results, function(res) {
  if (!res$valid) return(NULL)
  tibble(
    Interval = res$interval,
    `RMST Placebo` = sprintf("%.2f ± %.2f", res$arm0_rmst, res$arm0_se),
    `RMST Treatment` = sprintf("%.2f ± %.2f", res$arm1_rmst, res$arm1_se),
    `Difference (T-P)` = sprintf("%.2f", res$diff_rmst),
    `95% CI` = sprintf("%.2f - %.2f", res$lower_ci, res$upper_ci),
    `p-value` = sprintf("%.3f", res$p_value),
    `Significant` = ifelse(res$p_value < 0.05, "Yes", "No")
  )
})

if (nrow(df_results) > 0) {
  cat("\n--- RMST Results Summary Table ---\n")
  print(df_results)
}

# RMST Difference plot
diff_plot_data <- map_df(rmst_results, function(res) {
  if (!res$valid) return(NULL)
  tibble(
    interval = res$interval,
    diff_rmst = res$diff_rmst,
    lower_ci = res$lower_ci,
    upper_ci = res$upper_ci,
    p_value = res$p_value,
    significant = res$p_value < 0.05
  )
})

if (nrow(diff_plot_data) > 0) {
  rmst_diff_plot <- ggplot(diff_plot_data, aes(x = interval, y = diff_rmst, color = significant)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    labs(
      title = "RMST Difference (DPP-1 Inhibitor - Placebo) by Interval",
      x = "Time Interval (weeks)",
      y = "RMST Difference (weeks)",
      color = "p < 0.05"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#56B4E9")) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  print(rmst_diff_plot)
  # Salva como TIFF
  ggsave(
    filename = "/Users/yourpathwaytosavefile",
    plot = rmst_diff_plot,
    device = "png",
    dpi = 600,
    width = 10,
    height = 8, # um pouco mais alto por conta da tabela
    units = "in",
  )
  
}

# 7. LANDMARK ANALYSIS ========================================================

cat("\n=== LANDMARK ANALYSIS ===\n")

# Create directory for landmark plots if it doesn't exist
if (!dir.exists("landmark_plots")) {
  dir.create("landmark_plots")
  cat("Created directory 'landmark_plots' for individual landmark plots\n")
}

# Landmark Analysis Function
landmark_analysis <- function(data, landmark_times = c(0, 17, 34, 52), save_individual_plots = TRUE) {
  results_list <- list()
  plots_list <- list()
  
  # Loop through landmarks (except the last one which is just the upper bound)
  for (i in 1:(length(landmark_times)-1)) {
    start_time <- landmark_times[i]
    end_time <- landmark_times[i+1]
    
    # Create a descriptive name for this landmark period
    landmark_name <- paste0(start_time, "-", end_time)
    cat("\n--- LANDMARK ANALYSIS:", landmark_name, "WEEKS ---\n")
    
    # Filter data for patients who haven't had an event before the landmark start
    if (start_time > 0) {
      # For landmarks after time 0, keep only patients who are event-free at landmark start
      landmark_cohort <- data %>%
        filter(month >= start_time) %>%
        mutate(
          # Recalibrate time to start from 0 at the landmark
          landmark_time = month - start_time,
          # Truncate time at the end of the landmark window
          landmark_time = pmin(landmark_time, end_time - start_time),
          # Censor events outside the landmark window
          landmark_status = ifelse(month > end_time, 0, status)
        )
    } else {
      # For the first landmark (starting at time 0), just truncate at end time
      landmark_cohort <- data %>%
        mutate(
          landmark_time = pmin(month, end_time),
          landmark_status = ifelse(month > end_time, 0, status)
        )
    }
    
    # Print cohort summary
    cat("Number of patients in landmark cohort:", nrow(landmark_cohort), "\n")
    cat("  Placebo group:", sum(landmark_cohort$treated == 0), "\n")
    cat("  Treatment group:", sum(landmark_cohort$treated == 1), "\n")
    
    # Fit Cox model for this landmark period
    cox_model <- tryCatch({
      model <- coxph(Surv(landmark_time, landmark_status) ~ treated, data = landmark_cohort)
      summary_data <- summary(model)
      
      # Extract hazard ratio and confidence interval
      hr <- exp(summary_data$coefficients[1, "coef"])
      hr_lower <- exp(summary_data$coefficients[1, "coef"] - 1.96 * summary_data$coefficients[1, "se(coef)"])
      hr_upper <- exp(summary_data$coefficients[1, "coef"] + 1.96 * summary_data$coefficients[1, "se(coef)"])
      p_value <- summary_data$coefficients[1, "Pr(>|z|)"]
      
      cat("Cox model results:\n")
      cat("  Hazard Ratio:", round(hr, 2), "95% CI [", round(hr_lower, 2), "-", round(hr_upper, 2), "]\n")
      cat("  P-value:", round(p_value, 3), "\n")
      
      list(
        model = model,
        hr = hr,
        hr_lower = hr_lower,
        hr_upper = hr_upper,
        p_value = p_value
      )
    }, error = function(e) {
      cat("Error fitting Cox model:", conditionMessage(e), "\n")
      NULL
    })
    
    # Create Kaplan-Meier plot
    km_fit <- tryCatch({
      survfit(Surv(landmark_time, landmark_status) ~ treated, data = landmark_cohort)
    }, error = function(e) {
      cat("Error fitting KM curve:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(km_fit)) {
      # Create annotation text for the plot
      if (!is.null(cox_model)) {
        hr_text <- paste0("Hazard ratio, ", sprintf("%.2f", cox_model$hr), 
                          " (95% CI, ", sprintf("%.2f", cox_model$hr_lower), "-", 
                          sprintf("%.2f", cox_model$hr_upper), ")")
        p_text <- paste0("P=", sprintf("%.3f", cox_model$p_value))
      } else {
        hr_text <- "Hazard ratio not available"
        p_text <- "P-value not available"
      }
      
      # Create KM plot with better dimensions for individual viewing
      
      label_x <- (end_time - start_time) * 0.75
      
      km_plot <- ggsurvfit(km_fit, type = "risk", size = 1.2) +
        add_confidence_interval() +
        scale_color_manual(values = c("#0067BA", "#F9CB76")) +
        scale_fill_manual(values = c("#0067BA", "#F9CB76")) +
        annotate("rect", xmin = label_x - 0.5, xmax = label_x - 0.3, ymin = 0.018, ymax = 0.028, fill = "#F9CB76", color = NA) +
        annotate("text", x = label_x, y = 0.023, size = 4, hjust = 0, label = "DPP-1 Inhibitor") +
        annotate("rect", xmin = label_x - 0.5, xmax = label_x - 0.3, ymin = 0.058, ymax = 0.068, fill = "#0067BA", color = NA) +
        annotate("text", x = label_x, y = 0.063, size = 4, hjust = 0, label = "Placebo") +
        guides(color = "none", fill = "none") +
        coord_cartesian(xlim = c(0, end_time - start_time), ylim = c(1, 0)) +
        
        scale_y_continuous(
          trans = "reverse",  # flips the axis labels (not the curve)
          breaks = seq(0, 1, 0.2),
          labels = rev(seq(0, 100, 20)),  # shows labels from 100 to 0
          expand = c(0, 0)
        ) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0, end_time - start_time, 8)) +
        annotate("text", label = hr_text, x = -Inf, y = Inf, hjust = -0.05, vjust = 2.5, size = 4) +
        annotate("text", label = p_text, x = -Inf, y = Inf, hjust = -0.32, vjust = 4, size = 4) +
        labs(
          x = paste0("\nWeeks Since Landmark (", landmark_name, " window)\n"), 
          y = "Percentage of patients with no exacerbations (%)\n",
          title = paste0("Landmark Analysis: ", landmark_name, " Weeks"),
          subtitle = paste0("n = ", nrow(landmark_cohort), " patients")
        ) +
        theme_cowplot() +
        theme(
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          panel.grid.minor = element_blank()
        )
      
      # Print individual plot
      print(km_plot)
      
      # Save individual plot if requested
      if (save_individual_plots) {
        filename <- paste0("/Users/yourpathwaytosavefile", landmark_name, "_weeks.png")
        ggsave(filename, km_plot, width = 10, height = 8, dpi = 300)
        cat("Saved individual plot:", filename, "\n")
        
        # Also save as PDF for better quality
        filename_pdf <- paste0("/Users/yourpathwaytosavefile", landmark_name, "_weeks.pdf")
        ggsave(filename_pdf, km_plot, width = 10, height = 8)
        cat("Saved individual plot (PDF):", filename_pdf, "\n")
      }
      
      # Add results to lists for later use
      results_list[[landmark_name]] <- list(
        cohort = landmark_cohort,
        cox_model = cox_model,
        km_fit = km_fit,
        plot = km_plot
      )
      
      plots_list[[landmark_name]] <- km_plot
      
      # Calculate event counts and rates
      events_placebo <- sum(landmark_cohort$landmark_status[landmark_cohort$treated == 0])
      events_treated <- sum(landmark_cohort$landmark_status[landmark_cohort$treated == 1])
      
      cat("Events in placebo group:", events_placebo, 
          "(", round(100 * events_placebo/sum(landmark_cohort$treated == 0), 1), "%)\n")
      cat("Events in treatment group:", events_treated, 
          "(", round(100 * events_treated/sum(landmark_cohort$treated == 1), 1), "%)\n")
      
      # Calculate and save summary statistics for this landmark
      summary_stats <- data.frame(
        Landmark = landmark_name,
        N_total = nrow(landmark_cohort),
        N_placebo = sum(landmark_cohort$treated == 0),
        N_treated = sum(landmark_cohort$treated == 1),
        Events_placebo = events_placebo,
        Events_treated = events_treated,
        Event_rate_placebo = round(100 * events_placebo/sum(landmark_cohort$treated == 0), 1),
        Event_rate_treated = round(100 * events_treated/sum(landmark_cohort$treated == 1), 1),
        HR = ifelse(!is.null(cox_model), round(cox_model$hr, 3), NA),
        HR_CI_lower = ifelse(!is.null(cox_model), round(cox_model$hr_lower, 3), NA),
        HR_CI_upper = ifelse(!is.null(cox_model), round(cox_model$hr_upper, 3), NA),
        P_value = ifelse(!is.null(cox_model), round(cox_model$p_value, 4), NA)
      )
      
      # Save summary stats if this is the first landmark (create file) or append
      if (landmark_name == names(results_list)[1]) {
        write.csv(summary_stats, "landmark_plots/landmark_summary_statistics.csv", row.names = FALSE)
      } else {
        write.table(summary_stats, "landmark_plots/landmark_summary_statistics.csv", 
                    sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
  
  # Return results
  return(list(
    results = results_list,
    plots = plots_list
  ))
}

# Run landmark analysis with specified landmarks
landmark_results <- landmark_analysis(meta_ipd, landmark_times = c(0, 17, 34, 52), save_individual_plots = TRUE)

# Create a compact combined plot for quick overview
if(length(landmark_results$plots) > 0) {
  # Create smaller versions for the combined plot
  combined_plots <- map(names(landmark_results$plots), function(name) {
    plot <- landmark_results$plots[[name]]
    # Simplify for combined view
    plot + 
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)
      )
  })
  
  combined_landmark_plot <- wrap_plots(combined_plots, ncol = 2) +
    plot_annotation(
      title = "Landmark Analysis Overview",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )
  
  # Print the combined overview
  print(combined_landmark_plot)
  
  # Save the combined overview
  ggsave("/Users/yourpathwaytosavefile", combined_landmark_plot, 
         width = 12, height = 10, dpi = 300)
  ggsave("/Users/yourpathwaytosavefile", combined_landmark_plot, 
         width = 12, height = 10)
  
  cat("\n--- LANDMARK ANALYSIS COMPLETE ---\n")
  cat("Individual plots saved in 'landmark_plots' folder:\n")
  cat("- PNG files for presentations/reports\n")
  cat("- PDF files for high-quality printing\n")
  cat("- Summary statistics saved as CSV file\n")
  cat("- Combined overview saved as landmark_analysis_overview.png/.pdf\n")
}

# 8. SUMMARY AND CONCLUSIONS =================================================

cat("\n=== ANALYSIS SUMMARY ===\n")

# Summary of key findings
cat("Key Findings:\n")
cat("1. Cox proportional hazards model:\n")
summary(cox_model)

cat("\n2. Proportional hazards assumption test:\n")
print(ph_test)

cat("\n3. RMST analysis completed for", length(intervals), "time intervals\n")

cat("\n4. Landmark analysis completed for", length(landmark_results$results), "time periods\n")

cat("\nAnalysis completed successfully!\n")
cat("All plots and results have been generated.\n")


# =============================================================================
# END OF SCRIPT
# =============================================================================
# Get the risk table as a data frame
risk_table <- summary(km_fit, times = seq(0, 52, by = 13))  # change interval as needed

# Convert to tibble for cleaner printing
risk_df <- tibble::tibble(
  time = risk_table$time,
  n_risk = risk_table$n.risk,
  n_event = risk_table$n.event,
  group = rep(levels(as.factor(meta_ipd$treated)), each = length(risk_table$time) / 2)
)

# Print to console
cat("\n=== RISK TABLE (Number at Risk by Time and Group) ===\n")
print(risk_df)
