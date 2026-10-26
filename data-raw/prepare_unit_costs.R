# Prepare price, unit-cost, inflation and utility-norm inputs from public
# sources. Each source is downloaded, checked against the md5 of the edition
# used, and reduced to a small CSV in data-raw/. Run from the project root:
#
#   Rscript data-raw/prepare_unit_costs.R
#
# Outputs
#   data-raw/emit_2025_26.csv         generic drug prices per mg (eMIT)
#   data-raw/ncc_2024_25.csv          chemotherapy delivery HRG unit costs (NCC)
#   data-raw/nhscii.csv               NHS Cost Inflation Index (PSSRU 2025)
#   data-raw/eq5d_norms_dsu2022.csv   EQ-5D-3L by age and sex (NICE DSU 2022)

download_checked <- function(url, md5, ext) {
  f <- tempfile(fileext = ext)
  utils::download.file(url, f, mode = "wb", quiet = TRUE,
                       headers = c("User-Agent" = "Mozilla/5.0"))
  got <- unname(tools::md5sum(f))
  if (got != md5) stop("md5 mismatch for ", url, ": ", got, " (expected ", md5, ")")
  f
}

# eMIT national database, 1 July 2025 to 30 June 2026 (Open Government Licence).
# Price per mg is the quantity-weighted average across vial presentations of
# each drug. Ready-made infusion bags are excluded because they are
# dose-banded products, not general vial stock.
emit_url <- paste0("https://assets.publishing.service.gov.uk/media/6a798a800a700895e2d79ff0/",
                   "emit-national-database-1-july-2025-to-30-june-2026.ods")
emit <- as.data.frame(readODS::read_ods(download_checked(emit_url, "bb86c0c1c8e4f9e3fa2165e6e3e69190",
                                                         ".ods"),
                                        sheet = 1, col_names = FALSE, .name_repair = "minimal"))
names(emit) <- c("name", "npc", "quantity", "wap", "sd")
stopifnot(emit$name[2] == "Name & Pack Size", emit$wap[2] == "Weighted Average Price")
emit <- emit[-(1:2), ]
emit$quantity <- as.numeric(emit$quantity)
emit$wap <- as.numeric(emit$wap)

drugs <- c("carboplatin", "cisplatin", "docetaxel", "gemcitabine", "paclitaxel", "pemetrexed")
rows <- lapply(drugs, function(d) {
  x <- emit[grepl(paste0("^", d, " "), tolower(emit$name)) &
              grepl("vial", tolower(emit$name)) &
              !grepl("albumin", tolower(emit$name)), ]
  strength <- tolower(sub("^\\S+\\s+([0-9.]+(mg|g))\\b.*$", "\\1", x$name, ignore.case = TRUE))
  mg <- ifelse(grepl("g$", strength) & !grepl("mg$", strength),
               as.numeric(sub("g$", "", strength)) * 1000, as.numeric(sub("mg$", "", strength)))
  pack <- as.numeric(sub(".*Packsize\\s+([0-9]+).*", "\\1", x$name))
  stopifnot(!anyNA(mg), !anyNA(pack))
  per_mg <- x$wap / (mg * pack)
  data.frame(drug = d, n_presentations = nrow(x), units_sold = sum(x$quantity),
             price_per_mg = sum(per_mg * x$quantity) / sum(x$quantity),
             presentations = paste(sprintf("%s [%.4f/mg]", trimws(sub("/.*", "", x$name)), per_mg),
                                   collapse = "; "))
})
utils::write.csv(do.call(rbind, rows), "data-raw/emit_2025_26.csv", row.names = FALSE)

# National Cost Collection 2024/25 national schedule (published 6 Nov 2025).
ncc_url <- "https://www.england.nhs.uk/wp-content/uploads/2025/11/NCC_National-Schedule_2024_25.zip"
zip <- download_checked(ncc_url, "a08318b827269155f52d8aea46f80b6f", ".zip")
dir <- tempfile("ncc_")
utils::unzip(zip, exdir = dir)
summ <- list.files(dir, pattern = "^Summary HRG\\.xlsx$", recursive = TRUE, full.names = TRUE)
x <- suppressMessages(as.data.frame(readxl::read_excel(summ, col_names = FALSE)))
dept <- unlist(x[1, ])
meas <- unlist(x[2, ])
want <- c("SB12Z", "SB13Z", "SB14Z", "SB15Z")
ncc <- do.call(rbind, lapply(want, function(h) {
  i <- grep(paste0("^", h, " "), x[[1]])
  stopifnot(length(i) == 1)
  v <- unname(unlist(x[i, ]))
  keep <- which(meas == "Unit cost" & dept %in% c("Outpatient Procedures", "Daycase"))
  data.frame(hrg = h, description = sub("^\\S+ - ", "", v[1]),
             setting = dept[keep], unit_cost = as.numeric(v[keep]))
}))
utils::write.csv(ncc, "data-raw/ncc_2024_25.csv", row.names = FALSE)

# NHS Cost Inflation Index, overall (NHS providers), annual % increase on the
# previous year. Transcribed from Jones K, Weatherly H, Birch S, et al. Unit
# Costs of Health and Social Care 2025 Manual, Table 12.1.1 (p96). The
# 2024/2025 value is provisional.
nhscii <- data.frame(
  year = c("2015/16", "2016/17", "2017/18", "2018/19", "2019/20", "2020/21", "2021/22",
           "2022/23", "2023/24", "2024/25"),
  pct_increase_overall = c(0.40, 2.09, 1.24, 1.60, 2.14, 2.49, 2.58, 7.32, 2.47, 4.02)
)
utils::write.csv(nhscii, "data-raw/nhscii.csv", row.names = FALSE)

# Expected EQ-5D-3L by age and sex, Hernández Alava M, Pudney S, Wailoo A.
# Estimating EQ-5D by age and sex for the UK. NICE DSU report, 2022.
dsu_url <- "https://sheffield.ac.uk/media/34060/download?attachment"
dsu <- suppressMessages(as.data.frame(readxl::read_xls(
  download_checked(dsu_url, "b0ac041669012b16d833893e4f3b19cb", ".xls"),
  sheet = "Expected EQ-5D-3L", col_names = FALSE)))
stopifnot(dsu[2, 1] == "age", dsu[2, 2] == "Females", dsu[2, 3] == "Males")
dsu <- dsu[-(1:2), ]
norms <- data.frame(age = as.integer(dsu[[1]]), female = as.numeric(dsu[[2]]),
                    male = as.numeric(dsu[[3]]))
norms <- norms[!is.na(norms$age), ]
utils::write.csv(norms, "data-raw/eq5d_norms_dsu2022.csv", row.names = FALSE)

message("Unit-cost inputs written.")
