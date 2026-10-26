# Scenario explorer for the pembrolizumab partitioned survival model.
# Run the pipeline first (targets::tar_make()), then from the project root:
#   shiny::runApp("app")
# The app reuses the pipeline's fitted models and inputs and calls the same
# model functions; it does not refit anything.

library(shiny)

root <- if (file.exists("_targets.R")) "." else ".."
targets::tar_source(file.path(root, "R"))
store <- file.path(root, "_targets")
rd <- function(name) targets::tar_read_raw(name, store = store)
params <- rd("params")
ctx <- rd("ctx")
base_settings <- rd("base_settings")
base_icer <- rd("base_case")$summary$incremental[["icer"]]
models <- candidate_models()$model

gbp <- function(x) paste0(ifelse(x < 0, "−£", "£"), formatC(abs(x), format = "f", digits = 0, big.mark = ","))

ui <- fluidPage(
  titlePanel("Pembrolizumab vs chemotherapy: scenario explorer"),
  p("Independent model built from public sources; not NICE's model or results. List prices unless a discount is set."),
  sidebarLayout(
    sidebarPanel(
      selectInput("os_model", "OS distribution (both arms)", c(models, "Piecewise"), "Lognormal"),
      conditionalPanel("input.os_model == 'Piecewise'",
                       selectInput("pw_weeks", "KM cut-off (weeks), then exponential", c(23, 33, 43), 33)),
      selectInput("pfs_model", "PFS distribution (both arms)", models, "Lognormal"),
      selectInput("waning", "Treatment-effect waning (HR = 1 from)",
                  c("None" = NA, "3 years" = 3, "5 years" = 5), NA),
      radioButtons("tot", "Pembrolizumab duration",
                   c("PFS, 2-year stopping rule" = "pfs_capped", "Treat to progression" = "pfs_uncapped")),
      radioButtons("utility", "Utilities", c("Time to death" = "ttd", "Progression-based" = "progression")),
      sliderInput("discount", "Pembrolizumab price discount", 0, 0.9, 0, step = 0.01),
      sliderInput("rate", "Discount rate (costs and QALYs)", 0, 0.06, params$discount_costs, step = 0.005)
    ),
    mainPanel(
      tableOutput("results"),
      p(textOutput("base"))
    )
  )
)

server <- function(input, output, session) {
  result <- reactive({
    s <- base_settings
    s$curves <- curve_spec(
      os_model = input$os_model, pfs_model = input$pfs_model,
      piecewise_weeks = if (input$os_model == "Piecewise") as.numeric(input$pw_weeks) else NA,
      waning_years = suppressWarnings(as.numeric(input$waning)), waning_hr = 1)
    s$tot <- input$tot
    s$utility <- input$utility
    q <- params
    q$price_discount_pembro <- input$discount
    q$discount_costs <- input$rate
    q$discount_qalys <- input$rate
    run_once(q, s, ctx)
  })
  output$results <- renderTable({
    r <- result()
    data.frame(
      Outcome = c("Incremental cost", "Incremental QALYs", "ICER (per QALY)", "NMB at £20,000/QALY",
                  "NMB at £30,000/QALY"),
      Value = c(gbp(r[["inc_cost"]]), formatC(r[["inc_qaly"]], format = "f", digits = 2), gbp(r[["icer"]]),
                gbp(r[["nmb_20000"]]), gbp(r[["nmb_30000"]])))
  })
  output$base <- renderText(paste("Base-case ICER:", gbp(base_icer), "per QALY."))
}

shinyApp(ui, server)
