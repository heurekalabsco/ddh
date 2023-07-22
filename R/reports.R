#' SEND REPORT MESSAGE ---------------------------------------------------------
#'
#' Generates a message that gets sent to the AWS simple queue service (SQS) for email report generation
#'
#' @param first_name First name of the report addressee, which is used in the report email body
#' @param last_name Last name of the report addressee.
#' @param email_address Email address for the report to be delivered to
#' @param input A list containing type, query, and content variables.
#' @param private A Boolean determining which information to include in the report
#' @param greeting One of "default", "seminar", or "paper" to customize the front of the report
#'
#' @examples
#' send_report_message(first_name = "Matthew", last_name = "Hirschey", email_address = "hey_at_datadriventhypothesis.com", input = list(type = "gene", query = "ROCK1", content = "ROCK1"), private = TRUE)
#' \dontrun{
#' send_report_message(first_name = "Matthew", last_name = "Hirschey", email_address = "hey_at_datadriventhypothesis.com", input = list(type = "gene", query = "ROCK1", content = "ROCK1"), private = TRUE)
#' }
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
send_report_message <- function(first_name,
                                last_name,
                                email_address,
                                input = list(),
                                private,
                                greeting = "default"){
  json_array <-
    tibble::tibble(
      first_name = first_name,
      last_name = last_name,
      email_address = email_address,
      type = input$type,
      subtype = input$subtype,
      query = stringr::str_c(input$query, collapse = ", "),
      content = stringr::str_c(input$content, collapse = ", "),
      private = private,
      greeting = greeting
    ) %>%
    jsonlite::toJSON(dataframe = "rows")

  sqs <- paws::sqs()
  sqs$send_message(
    QueueUrl = Sys.getenv("AWS_SQS_SERVICE_URL"),
    MessageBody = json_array
  )
  message(glue::glue("{glue::glue_collapse(input$query, sep = ', ')} sqs message sent for {first_name} ({email_address})"))
}

