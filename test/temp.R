devtools::load_all("..")

d <- nmrdata_1d_from_csv("glutamate.csv") %>%
  correct_phase() %>%
  correct_reference()

peaks <- read_delim("glutamate_peaks.csv", "\t") %>% pull(2) %>% as.numeric()
print(peaks)

nmroptions$fit$opts$padding <- min(diff(peaks))/2

r <- nmrresonance_1d(peaks, id = "data", sf = 600.13)

f <- nmrfit_1d(r, d, delay = TRUE)
f <- set_offset_bounds(f, position = c(-0.01, 0.01))
f <- fit(f)
