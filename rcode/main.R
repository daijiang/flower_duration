source("rcode/00_pkg_functions.R")

all_dat = read_csv("data/dat_raw_52sp.csv")

# very time consuming !! Don't re-run it.
# the code below represents the main functions used to get 
# estimated flowering onset and offset dates for each species
# in grid cells with enough records (25km by 25km).

all_25km = map(unique(all_dat$scientific_name),
               ~plt_summary(cell_size = 25000, dat = filter(all_dat, scientific_name == .x), 
                            n_per_cell = 10))
names(all_25km) = unique(all_dat$scientific_name)

all_25km_phenesse = lapply(names(all_25km), function(.x){
  cat(.x, "\n")
  run_phenesse(all_25km[[.x]]$dat_to_use, 
               earliest_year = 2017,
               last_year = 2019, num_cores = 10, 
               onset_perct = 0.05, offset_perct = 0.95,
               n_item = 1000, save_rds = TRUE)
})

# then extract climatic data from PRISM, given the large databade needed, omit here.
# instead, the final product (phenology + climate + human population density) is in 
# the data folder ("data/dat_flower_pheno.csv")

# future climatic and human population density data: ("data/dat_future_envi.csv")

# Main analyses

d1 = read_csv("data/dat_flower_pheno.csv")

(ave_pop_log = mean(d1$pop_25km_log, na.rm = T))
(sd_pop_log = sd(d1$pop_25km_log, na.rm = T))
(ave_bio_1 = mean(d1$bio_1, na.rm = T))
(sd_bio_1 = sd(d1$bio_1, na.rm = T))
(ave_bio_12 = mean(d1$bio_12, na.rm = T))
(sd_bio_12 = sd(d1$bio_12, na.rm = T))
(ave_bio_4 = mean(d1$bio_4, na.rm = T))
(sd_bio_4 = sd(d1$bio_4, na.rm = T))
(ave_chill_days = mean(d1$n_chilling_days_winter, na.rm = T))
(sd_chill_days = sd(d1$n_chilling_days_winter, na.rm = T))

d = mutate_at(d1, .vars = vars(pop_25km_log, bio_1, bio_4, bio_12, bio_15, n_chilling_days_winter),
              .funs = function(x) (x - mean(x))/sd(x))
names(d) = gsub("_", "", names(d))

phy = ape::read.tree("data/phy.tre")
phy$tip.label = gsub("_", " ", phy$tip.label)

if(file.exists("data/m_onset_final.rds")){
  m_onset_final = readRDS("data/m_onset_final.rds")
  pm_onset = readRDS("data/pglmm_onset_final.rds")
} else {
  m_onset_final = lmer(
    onsetdoy ~ bio1 + bio4 + nchillingdayswinter +
      phenonsetcat + bio1:phenonsetcat +
      (1 | idcells) + (1 | sp) + (0 + bio1 | sp) + (0 + bio4 | sp) +
      (0 + nchillingdayswinter | sp),
    data = d
  )
  saveRDS(m_onset_final, file = "data/m_onset_final.rds")
  
  pm_onset = pglmm(onsetdoy ~ bio1 + bio4 + nchillingdayswinter + 
                     phenonsetcat + bio1:phenonsetcat +
                     (1 | idcells) + (1 | sp__) + (0 + bio1 | sp__) + 
                     (0 + bio4 | sp__) + (0 + nchillingdayswinter | sp__), 
                   data = d, cov_ranef = list(sp = phy), bayes = T)
  saveRDS(pm_onset, file = "data/pglmm_onset_final.rds")
}

if(file.exists("data/m_offset_final.rds")){
  m_offset_final = readRDS("data/m_offset_final.rds")
  pm_offset = readRDS("data/pglmm_offset_final.rds")
} else {
  m_offset_final = lmer(
    offsetdoy ~ bio1 + bio12 + pop25kmlog + phenonsetcat + growthform +
      bio1:phenonsetcat + pop25kmlog:growthform + bio12:growthform + 
      (1 | idcells) + (1 | sp) +  (0 + bio1 | sp) + # (0 + bio12 | sp) + 
      (0 + pop25kmlog | sp), 
    data = d,
    REML = T,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
  saveRDS(m_offset_final, file = "data/m_offset_final.rds")
  
  pm_offset = pglmm(offsetdoy ~ bio1 + bio12 + pop25kmlog + 
                      phenonsetcat + growthform + 
                      bio1:phenonsetcat + pop25kmlog:growthform + 
                      bio12:growthform +
                      (1 | idcells) + (1 | sp__) + 
                      (0 + bio1 | sp__) + (0 + pop25kmlog | sp__), 
                    data = d, 
                    cov_ranef = list(sp = phy), bayes = T)
  saveRDS(pm_offset, file = "data/pglmm_offset_final.rds")
}

if(file.exists("data/m_dur_final.rds")){
  m_dur_final = readRDS("data/m_dur_final.rds")
  pm_dur = readRDS("data/pglmm_dur_final.rds")
} else {
  m_dur_final = lmer(
    dur ~ bio1 + bio4 + bio12 + pop25kmlog + nchillingdayswinter + 
      phenonsetcat + growthform + 
      bio1:phenonsetcat + pop25kmlog:phenonsetcat + 
      bio1:growthform + bio4:growthform + 
      (1 | idcells) + (1 | sp) + (0 + bio1 | sp) + 
      (0 + bio4 | sp) + (0 + bio12 | sp) + 
      (0 + nchillingdayswinter | sp) +
      (0 + pop25kmlog | sp), 
    data = d,
    REML = T,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
  saveRDS(m_dur_final, file = "data/m_dur_final.rds")
  
  pm_dur = pglmm(dur ~ bio1 + bio4 + bio12 + pop25kmlog + nchillingdayswinter +  
                   phenonsetcat + growthform + bio1:phenonsetcat +  
                   pop25kmlog:phenonsetcat + bio1:growthform + bio4:growthform +  
                   (1 | idcells) + (1 | sp__) + 
                   (0 + bio1 | sp__) + 
                   (0 + bio4 | sp__) + 
                   (0 + bio12 | sp__) + 
                   (0 + nchillingdayswinter | sp__) +  
                   (0 + pop25kmlog | sp), 
                 data = d, verbose = T,
                 cov_ranef = list(sp = phy), bayes = T)
  saveRDS(pm_dur, file = "data/pglmm_dur_final.rds")
}

