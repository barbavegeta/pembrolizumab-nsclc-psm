# Prepare general-population mortality inputs from the ONS National Life
# Tables for England, 2022-2024 edition (released 10 December 2025).
# Contains public sector information licensed under the Open Government
# Licence v3.0.
#
# Run from the project root:
#   Rscript data-raw/prepare_lifetables.R
# Output: data-raw/lifetables_england_2022_2024.csv (age, sex, qx, mx)

url <- paste0(
  "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/",
  "birthsdeathsandmarriages/lifeexpectancies/datasets/",
  "nationallifetablesenglandreferencetables/current/nlte198020223.xlsx"
)
expected_md5 <- "0973d5cdba8b1244ea0a740350d45b2a"
sheet <- "2022-2024"

xlsx <- tempfile(fileext = ".xlsx")
utils::download.file(url, xlsx, mode = "wb", quiet = TRUE)
md5 <- unname(tools::md5sum(xlsx))
if (md5 != expected_md5) {
  stop("Downloaded file md5 ", md5, " differs from ", expected_md5,
       "; ONS may have issued a new edition. Check the edition and update this script.")
}

raw <- as.data.frame(readxl::read_excel(xlsx, sheet = sheet, col_names = FALSE,
                                        .name_repair = "minimal"))
names(raw) <- paste0("c", seq_along(raw))
# Row 6 holds the headers. Males are in columns 1-6, a blank column follows,
# and females are in columns 8-13.
stopifnot(raw$c1[5] == "Males", raw$c8[5] == "Females",
          raw$c1[6] == "age", raw$c2[6] == "mx", raw$c3[6] == "qx")
body <- raw[-(1:6), ]
body <- body[!is.na(body$c1), ]

tidy <- rbind(
  data.frame(age = as.integer(body$c1), sex = "male",
             mx = as.numeric(body$c2), qx = as.numeric(body$c3)),
  data.frame(age = as.integer(body$c8), sex = "female",
             mx = as.numeric(body$c9), qx = as.numeric(body$c10))
)
stopifnot(!anyNA(tidy), all(table(tidy$sex) == 101), max(tidy$age) == 100)

utils::write.csv(tidy, "data-raw/lifetables_england_2022_2024.csv", row.names = FALSE)
message("Wrote ", nrow(tidy), " rows")
