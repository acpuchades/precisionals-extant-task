source("src/ext/main.r")

rob_main <- read_xlsx("data/rob-ext-main.xlsx") |>
  select(
    id = "ID", site = "Site",
    c9orf72_tested = "C9orf72 Tested",
    c9orf72_status = "C9orf72 Status",
    sod1_tested = "SOD1 Tested",
    sod1_status = "SOD1 Status",
    fus_tested = "FUS Tested",
    fus_status = "FUS Status",
    tardbp_tested = "TARDBP Tested",
    tardbp_status = "TARDBP Status"
  )

ireland_additions <- read_xlsx("data/Ireland_additional_mutations_v2.xlsx") |>
  select(
    id = "P-ALS ID", gene = "Gene",
    variant_g = "DNA change",
    variant_p = "Protein change",
    ref_genome = "Genome build"
  ) |>
  mutate(site = "Trinity", .after = id) |>
  mutate(id = str_replace(id, "PR_TUR_", "ITA-"))

sod1_variants <- read_xlsx("data/SOD1 collated.xlsx") |>
  select(id = "pals_id", site, sod1_variant) |>
  mutate(id = str_replace(id, "PR_TUR_", "ITA-")) |>
  bind_rows(
    ireland_additions |>
      filter(gene == "SOD1") |>
      select(id, site, sod1_variant = "variant_p")
  ) |>
  slice_sample(n = 1, by = id)

fus_variants <- read_xlsx("data/FUS collated.xlsx") |>
  select(id = "P-ALS ID", site = "Site", fus_variant = "variant") |>
  mutate(id = str_replace(id, "PR_TUR_", "ITA-")) |>
  bind_rows(
    ireland_additions |>
      filter(gene == "FUS") |>
      select(id, site, fus_variant = "variant_p")
  ) |>
  slice_sample(n = 1, by = id)

tardbp_variants <- read_xlsx("data/TARDBP collated.xlsx") |>
  select(id = "P-ALS ID", site = "Site", tardbp_variant = "Protein change") |>
  mutate(id = str_replace(id, "PR_TUR_", "ITA-")) |>
  bind_rows(
    ireland_additions |>
      filter(gene |> str_starts("TARDBP")) |>
      select(id, site, tardbp_variant = "variant_p")
  ) |>
  slice_sample(n = 1, by = id)

patients.orig <- ext_main |> select(
  id, site, c9orf72_status, sod1_status, fus_status, tardbp_status
)

patients.joined <- patients.orig |>
  full_join(
    sod1_variants |> mutate(sod1_status = "Positive", .before = sod1_variant),
    by = "id", suffix = c(".orig", ".new")
  ) |>
  mutate(
    site = coalesce(site.orig, site.new),
    sod1_status = coalesce(sod1_status.new, sod1_status.orig)
  ) |>
  select(-site.orig, -site.new) |>
  full_join(
    fus_variants |> mutate(fus_status = "Positive", .before = fus_variant),
    by = "id", suffix = c(".orig", ".new")
  ) |>
  mutate(
    site = coalesce(site.orig, site.new),
    fus_status = coalesce(fus_status.new, fus_status.orig)
  ) |>
  select(-site.orig, -site.new) |>
  full_join(
    tardbp_variants |> mutate(tardbp_status = "Positive", .before = tardbp_variant),
    by = "id", suffix = c(".orig", ".new")
  ) |>
  mutate(
    site = coalesce(site.orig, site.new),
    tardbp_status = coalesce(tardbp_status.new, tardbp_status.orig)
  ) |>
  select(-site.orig, -site.new)

sod1.changed <- patients.joined |>
  select(id, site, starts_with("sod1_")) |>
  filter(
    !is.na(sod1_status.new),
    is.na(sod1_status.orig) | sod1_status.orig != sod1_status.new
  )

tardbp.changed <- patients.joined |>
  select(id, site, starts_with("tardbp_")) |>
  filter(
    !is.na(tardbp_status.new),
    is.na(tardbp_status.orig) | tardbp_status.orig != tardbp_status.new
  )

fus.changed <- patients.joined |>
  select(id, site, starts_with("fus_")) |>
  filter(
    !is.na(fus_status.new),
    is.na(fus_status.orig) | fus_status.orig != fus_status.new
  )

write_xlsx(list(
  "SOD1" = patients.joined |>
    filter(sod1_status == "Positive") |>
    select(id, site, starts_with("sod1_")),
  "SOD1 changed" = sod1.changed,
  "TARDBP" = patients.joined |>
    filter(tardbp_status == "Positive") |>
    select(id, site, starts_with("tardbp_")),
  "TARDBP changed" = tardbp.changed,
  "FUS" = patients.joined |>
    filter(fus_status == "Positive") |>
    select(id, site, starts_with("fus_")),
  "FUS changed" = fus.changed
), "output/gene-counts.xlsx")

make_summary_table <- function(data) {
  t <- table(data)
  pt <- prop.table(t)
  res <- cbind(t, pt)
  colnames(res) <- c("N", "%")
  return(res)
}

sink("output/gene-diffs.txt")

cat("# C9orf72 (original)\n\n")
with(patients.orig, make_summary_table(c9orf72_status)) |> print()
cat("\n")

cat("# SOD1 (joined)\n\n")
with(patients.joined, make_summary_table(sod1_status)) |> print()
cat("\n")

cat("# SOD1 (original)\n\n")
with(patients.orig, make_summary_table(sod1_status)) |> print()
cat("\n")

cat("# FUS (joined)\n\n")
with(patients.joined, make_summary_table(fus_status)) |> print()
cat("\n")

cat("# FUS (original)\n\n")
with(patients.orig, make_summary_table(fus_status)) |> print()
cat("\n")

cat("# TARDBP (joined)\n\n")
with(patients.joined, make_summary_table(tardbp_status)) |> print()
cat("\n")

cat("# TARDBP (original)\n\n")
with(patients.orig, make_summary_table(tardbp_status)) |> print()
cat("\n")

sink()
