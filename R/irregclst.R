#' @title  Convert Irregular Longitudinal Data to Regular Intervals and Perform Clustering using excluding Repeated Responses (ERS) method
#' @description This function takes irregular longitudinal data and converts it into regularly spaced intervals using linear interpolation. It then computes the relative change in the response variable between consecutive time points, clusters the data based on these changes, and provides various visualizations of the process.
#'
#' @details
#' The \code{irregclst} function handles irregular longitudinal data by:
#' \itemize{
#'  \item Interpolating response values at regular time intervals.
#'  \item Calculating the relative change in the response values across time points.
#'  \item Clustering subjects based on these relative changes using alphabet labels ("a", "b", ..., "h") corresponding to different levels of deviation from the mean.
#'  \item Resolving cluster ties using a sum of squares criterion.
#' }
#' Visualizations of the data include plots for both the original irregular data and the regularized data, as well as histograms of time distributions and relative change trends.
#'
#' @param data A data frame containing the irregular longitudinal data.
#' @param subject_id_col A character string representing the name of the column with the subject IDs.
#' @param time_col A character string representing the name of the column with time values.
#' @param response_col A character string representing the name of the column with the response values.
#' @param interval_length A numeric value indicating the length of the regular intervals to which the time values should be converted.
#' @param rel Relative change method such as SRC, CARC and SWRC.
#' @return A list containing:
#' \itemize{
#'  \item \code{regular_data}: A data frame of the regularized longitudinal data.
#'  \item \code{regular_data_wide}: A wide-format version of the regularized data.
#'  \item \code{relative_change}: A data frame containing the relative changes in response values.
#'  \item \code{cluster_data}: A data frame with cluster assignments for each subject at each time step.
#'  \item \code{cluster_data_reduced}: A reduced version of \code{cluster_data} with only subject IDs and their final cluster assignments.
#'  \item \code{merged_data}: The wide-format data merged with the final cluster assignments.
#'  \item \code{plot_irregular}: A \code{ggplot} object showing the original irregular data.
#'  \item \code{plot_regular}: A \code{ggplot} object showing the regularized data.
#'  \item \code{plot_change}: A \code{ggplot} object showing the relative changes over time.
#'  \item \code{histogram_irregular}: A \code{ggplot} object showing the histogram of irregular time distribution.
#'  \item \code{histogram_regular}: A \code{ggplot} object showing the histogram of regular time distribution.
#' }
#' @references Reference
#'
#' @examples
#' ##
#' \donttest{
#' data(sdata)
#' #Using relative change method: Simple relative change (SRC)
#' fit1 <- irregclst(sdata, "subject_id", "time", "response", rel="SRC", interval_length = 3)
#' #for showing the regularized data in long format
#' fit1$regular_data
#' fit1$regular_data_wide #for showing the regularized data in wide format
#' fit1$cluster_data #dataset consisting clusters for different time points
#' fit1$merged_data #for showing the regularized data in wide format with final cluster
#' fit1$plot_regular #For plotting regularized longitudinal data
#' fit1$plot_irregular #For plotting irregular longitudinal data
#' fit1$plot_change #For plotting relative change
#' fit1$histogram_irregular #histogram for time of irregular data
#' fit1$histogram_regular #histogram for time of regular data
#' #Using relative change method: Cumulative average relative change (CARC)
#' fit2<-irregclst(sdata,"subject_id","time","response",rel="CARC",interval_length=3)
#' fit2$regular_data #for showing the regularized data in long format
#' fit2$regular_data_wide #for showing the regularized data in wide format
#' fit2$cluster_data #dataset consisting clusters for different time points
#' fit2$merged_data #for showing the regularized data in wide format with final cluster
#' fit2$plot_regular #For plotting regularized longitudinal data
#' fit2$plot_irregular #For plotting irregular longitudinal data
#' fit2$plot_change #For plotting relative change
#' fit2$histogram_irregular #histogram for time of irregular data
#' fit2$histogram_regular #histogram for time of regular data
#' #Using relative change method: Weighted sum relative change (WSRC)
#' fit3 <- irregclst(sdata, "subject_id", "time", "response", rel="WSRC", interval_length = 3)
#' fit3$regular_data #for showing the regularized data in long format
#' fit3$regular_data_wide #for showing the regularized data in wide format
#' fit3$cluster_data #dataset consisting clusters for different time points
#' fit3$merged_data #for showing the regularized data in wide format with final cluster
#' fit3$plot_regular #For plotting regularized longitudinal data
#' fit3$plot_irregular #For plotting irregular longitudinal data
#' fit3$plot_change #For plotting relative change
#' fit3$histogram_irregular #histogram for time of irregular data
#' fit3$histogram_regular #histogram for time of regular data
#'}
#' @export
#' @author author name
#' @seealso seealso
#'


irregclst <- function(data, subject_id_col, time_col, response_col, rel, interval_length) {
  #library(ggplot2)
  suppressWarnings({
  # Helper function to interpolate response at a given time point
  interpolate_response <- function(times, responses, new_time) {
    before <- max(times[times <= new_time], na.rm = TRUE)
    after <- min(times[times >= new_time], na.rm = TRUE)

    if (is.infinite(before)) {
      return(responses[which.min(times)])
    } else if (is.infinite(after)) {
      return(responses[which.max(times)])
    } else if (before == after) {
      return(responses[times == before])
    } else {
      before_response <- responses[times == before]
      after_response <- responses[times == after]
      slope <- (after_response - before_response) / (after - before)
      intercept <- before_response - slope * before
      return(slope * new_time + intercept)
    }
  }

  # Convert irregular data to regular data
  subject_ids <- unique(data[[subject_id_col]])
  max_time <- max(data[[time_col]], na.rm = TRUE)
  uniform_times <- seq(0, max_time, by = interval_length)

  regular_data <- data.frame(subject_id = integer(),
                             time = numeric(),
                             response = numeric())

  for (i in subject_ids) {
    subject_data <- data[data[[subject_id_col]] == i, ]
    for (t in uniform_times) {
      response <- interpolate_response(subject_data[[time_col]], subject_data[[response_col]], t)
      regular_data <- rbind(regular_data, data.frame(subject_id = i, time = t, response = response))
    }
  }

  # Remove duplicate last responses, handling NA values in comparisons
  for (i in subject_ids) {
    subject_rows <- regular_data[regular_data$subject_id == i, ]
    last_index <- nrow(subject_rows)

    while (last_index > 1) {
      if (!is.na(subject_rows$response[last_index]) &&
          !is.na(subject_rows$response[last_index - 1]) &&
          subject_rows$response[last_index] == subject_rows$response[last_index - 1]) {
        regular_data$response[regular_data$subject_id == i & regular_data$time == subject_rows$time[last_index]] <- NA
        last_index <- last_index - 1
      } else {
        break
      }
    }
  }

  # Convert long format data to wide format
  unique_times <- sort(unique(regular_data$time))
  regular_data_wide <- data.frame(subject_id = unique(regular_data$subject_id))

  for (t in unique_times) {
    time_col_name <- paste0("time_",t)
    regular_data_wide[[time_col_name]] <- NA
  }

  for (i in 1:nrow(regular_data)) {
    subject_id <- regular_data$subject_id[i]
    time_val <- regular_data$time[i]
    response_val <- regular_data$response[i]

    row_index <- which(regular_data_wide$subject_id == subject_id)
    col_index <- which(colnames(regular_data_wide) == paste0("time_", time_val))
    regular_data_wide[row_index, col_index] <- response_val
  }

  # Calculate relative change for each subject
  response_data <- regular_data_wide[, -1]
  p <- ncol(response_data)
  new_data <- matrix(NA, nrow = nrow(response_data), ncol = (p - 1))

  if (rel == "SRC") {
    for (i in 1:nrow(response_data)) {
      for (j in 1:(p - 1)) {
        if (!is.na(response_data[i, j]) && !is.na(response_data[i, j + 1])) {
          new_data[i, j] <- (response_data[i, j + 1] - response_data[i, j]) / response_data[i, j]
        } else {
          new_data[i, j] <- NA
        }
      }
    }
  } else if (rel == "WSRC") {
    for (i in 1:nrow(response_data)) {
      for (j in 1:(p - 1)) {
        weights <- seq(from = 1, to = j)
        weights <- weights / sum(weights)
        weighted_sum <- sum(weights * response_data[i, 1:j], na.rm = TRUE)
        if (weighted_sum != 0) {
          new_data[i, j] <- (response_data[i, j + 1] - weighted_sum) / weighted_sum
        } else {
          new_data[i, j] <- NA
        }
      }
    }
  } else if (rel == "CARC") {
    for (i in 1:nrow(response_data)) {
      for (j in 1:(p - 1)) {
        cum_sum <- sum(response_data[i, 1:j], na.rm = TRUE)
        avg_previous <- cum_sum / j
        if (avg_previous != 0) {
          new_data[i, j] <- (response_data[i, j + 1] - avg_previous) / avg_previous
        } else {
          new_data[i, j] <- NA
        }
      }
    }
  }

  new_data_df <- as.data.frame(new_data)
  colnames(new_data_df) <- paste0("Rel_Clst", 1:(p - 1))
  new_data_df <- cbind(subject_id = regular_data_wide$subject_id, new_data_df)

  # Clustering each column of changes
  cluster_data <- data.frame(subject_id = new_data_df$subject_id)

  for (col in 2:ncol(new_data_df)) {
    col_values <- new_data_df[[col]]
    mu <- mean(col_values, na.rm = TRUE)
    sigma <- sd(col_values, na.rm = TRUE)

    cluster_col <- cut(col_values,
                       breaks = c(-Inf, mu - 3 * sigma, mu - 2 * sigma, mu - sigma, mu,
                                  mu + sigma, mu + 2 * sigma, mu + 3 * sigma, Inf),
                       labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                       include.lowest = TRUE)

    cluster_data[[paste0("Cluster_", col - 1)]] <- cluster_col
  }

  # Compute means for each cluster
  cluster_means <- list()
  for (col in 2:ncol(new_data_df)) {
    cluster_col <- paste0("Cluster_", col - 1)
    if (cluster_col %in% names(cluster_data)) {
      cluster_assignments <- cluster_data[[cluster_col]]
      relative_changes <- new_data_df[[col]]
      cluster_means[[cluster_col]] <- tapply(relative_changes, cluster_assignments, mean, na.rm = TRUE)
    }
  }

  # Resolve majority voting with sum of squares criterion in case of ties
  final_cluster <- apply(cluster_data[, -1], 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(NA)

    counts <- table(row)
    max_count <- max(counts)
    tied_clusters <- names(counts[counts == max_count])

    if (length(tied_clusters) > 1) {
      sum_squares <- sapply(tied_clusters, function(c) {
        cluster_index <- match(c, c("a", "b", "c", "d", "e", "f", "g", "h"))
        cluster_col <- paste0("Cluster_", cluster_index)
        if (cluster_col %in% names(cluster_data)) {
          cluster_assignments <- cluster_data[[cluster_col]]
          relative_changes <- new_data_df[[cluster_index + 1]]
          cluster_mean <- mean(relative_changes[cluster_assignments == c], na.rm = TRUE)
          sum((relative_changes[cluster_assignments == c] - cluster_mean) ^ 2, na.rm = TRUE)
        } else {
          return(Inf)
        }
      })
      return(tied_clusters[which.min(sum_squares)])
    } else {
      return(tied_clusters)
    }
  })

  cluster_data$Final_Cluster <- final_cluster

  # Merge regular_data_wide with cluster assignments
  merged_data <- merge(regular_data_wide, cluster_data[, c("subject_id", "Final_Cluster")], by = "subject_id")

  # Plotting irregular data
  p0 <- ggplot2::ggplot(data = data, ggplot2::aes(x = .data[[time_col]], y = .data[[response_col]],
                                         group = as.factor(.data[[subject_id_col]]),
                                         color = as.factor(.data[[subject_id_col]]))) +
    ggplot2::geom_line(color = "black",alpha = 0.5) +
    ggplot2::geom_point(color = "red", alpha = 0.5) +
    ggplot2::labs(title = "Irregular Longitudinal Data", x = time_col, y = response_col) +
    ggplot2::theme_minimal() +
    ggplot2::guides(color = "none")+
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),  # Major breaks
      minor_breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)  # Minor breaks at a smaller interval
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),  # Major breaks
      minor_breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)  # Minor breaks at a smaller interval
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(size = 1, color = "gray70"), # Bold major gridlines
      panel.grid.minor = ggplot2::element_line(size = 0.8, color = "gray85"), # Bold minor gridlines
      plot.title = ggplot2::element_text(face = "bold"),  # Bold title
      axis.title = ggplot2::element_text(face = "bold"),  # Bold axis labels
      axis.text = ggplot2::element_text(face = "bold")    # Bold axis text
    )

  # Plotting regularized data
  p1 <- ggplot2::ggplot(data = regular_data, ggplot2::aes(x = time, y = response,
                                                 group = as.factor(subject_id),
                                                 color = as.factor(subject_id))) +
    ggplot2::geom_line(color = "black",alpha = 0.5) +
    ggplot2::geom_point(color = "red",alpha = 0.5) +
    ggplot2::labs(title = "Regular Longitudinal Data", x = time_col, y = response_col) +
    ggplot2::theme_minimal() +
    ggplot2::guides(color = "none")+
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),  # Major breaks
      minor_breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)  # Minor breaks at a smaller interval
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),  # Major breaks
      minor_breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)  # Minor breaks at a smaller interval
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(size = 1, color = "gray70"), # Bold major gridlines
      panel.grid.minor = ggplot2::element_line(size = 0.8, color = "gray85"), # Bold minor gridlines
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold")
    )

  # Transforming the data
  transformed_long <- data.frame(
    subject_id = rep(new_data_df$subject_id, times = ncol(new_data_df) - 1),
    Change_Timestep = rep(names(new_data_df)[-which(names(new_data_df) == "subject_id")], each = nrow(new_data_df)),
    Relative_Change = as.vector(t(new_data_df[, -which(names(new_data_df) == "subject_id")]))
  )

  # Remove rows where Relative_Change or Change_Timestep is NA
  transformed_long <- transformed_long[!is.na(transformed_long$Relative_Change) & !is.na(transformed_long$Change_Timestep), ]

  # Convert Change_Timestep to factor with levels ordered numerically
  transformed_long$Change_Timestep <- factor(transformed_long$Change_Timestep,
                                             levels = paste0("Rel_Clst", 1:(ncol(new_data_df) - 1)))

  # Plot relative change over time (p2)
  p2 <- ggplot2::ggplot(transformed_long, ggplot2::aes(x = Change_Timestep, y = Relative_Change,
                                              group = as.factor(subject_id), color = as.factor(subject_id))) +
    ggplot2::geom_line(color = "black",alpha = 0.5) +
    ggplot2::geom_point(color = "red",alpha = 0.5) +
    ggplot2::labs(title = "Relative Change Over Time", x = time_col, y = "Relative Change") +
    ggplot2::theme_minimal() +
    ggplot2::guides(color = "none")+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(size = 1, color = "gray70"), # Bold major gridlines
      panel.grid.minor = ggplot2::element_line(size = 0.8, color = "gray85"), # Bold minor gridlines
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold")
    )

  # Histogram of irregular time distribution (p3)
  p3 <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[time_col]])) +
    ggplot2::geom_histogram(binwidth = interval_length, fill = "skyblue", color = "black") +
    ggplot2::labs(title = "Histogram of Irregular Time Distribution", x = time_col, y = "Count") +
    ggplot2::theme_minimal()

  # Histogram of regular time distribution (p4)
  p4 <- ggplot2::ggplot(regular_data, ggplot2::aes(x = time)) +
    ggplot2::geom_histogram(binwidth = interval_length, fill = "orange", color = "black") +
    ggplot2::labs(title = "Histogram of Regular Time Distribution", x = time_col, y = "Count") +
    ggplot2::theme_minimal()

  # Return results
  return(list(
    regular_data = regular_data,
    regular_data_wide = regular_data_wide,
    relative_change = new_data_df,
    cluster_data = cluster_data,
    cluster_data_reduced = cluster_data[, c("subject_id", "Final_Cluster")],
    merged_data = merged_data,
    plot_irregular = p0,
    plot_regular = p1,
    plot_change = p2,
    histogram_irregular = p3,
    histogram_regular = p4
  ))
  })
}

utils::globalVariables(c('%>%', 'filter', 'sym',
                         'dcast', 'sd', 'select', 'ID', 'Final_Cluster',
                         'rename', 'ggplot', 'aes', 'geom_line',
                         'labs', 'theme_minimal', 'guides', 'time',
                         'subject_id', 'melt', 'Change_Timestep', 'Relative_Change',
                         'geom_histogram', 'geom_point','.data'))
