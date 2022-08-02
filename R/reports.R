#' Send Report Message
#'
#' \code{send_report_message} generates a message that gets sent to the AWS simple queue service (SQS) for report generation
#'
#' This is an emailing function that sends a message to the SQS queue for email generation
#'
#' @param first_name First name of the report addressee, which is used in the report email body
#' @param last_name Last name of the report addressee.
#' @param email_address Email address for the report to be delivered to
#' @param input A list containing type, query, and content variables.
#' @param private A Boolean determining which information to include in the report
#'
#' @export
#' @examples
#' send_report_message(first_name = "Matthew", last_name = "Hirschey", email_address = "hey_at_datadriventhypothesis.com", input = list(type = "gene", query = "ROCK1", content = "ROCK1"), private = TRUE)
#' \dontrun{
#' send_report_message(first_name = "Matthew", last_name = "Hirschey", email_address = "hey_at_datadriventhypothesis.com", input = list(type = "gene", query = "ROCK1", content = "ROCK1"), private = TRUE)
#' }
send_report_message <- function(first_name,
                                last_name,
                                email_address,
                                input = list(),
                                private){
  json_array = glue::glue('[{{
  "first_name":"{first_name}",
  "last_name":"{last_name}",
  "email_address":"{email_address}",
  "type":"{input$type}",
  "subtype":"{input$subtype}",
  "query":"{input$query}",
  "content":"{input$content}",
  "private":"{private}"}}]')
  sqs <- paws::sqs()
  sqs$send_message(
    QueueUrl = Sys.getenv("AWS_SQS_SERVICE_URL"),
    MessageBody = json_array
  )
  print(glue::glue("{input$query} sqs message sent for {first_name} ({email_address})"))
}

