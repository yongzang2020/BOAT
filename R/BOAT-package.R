#' BOAT: Bayesian umbrella trial design using patient-derived organoid data
#'
#' The BOAT design evaluates experimental treatments \eqn{A_1, \ldots, A_J} against a control
#' \eqn{A_0} in biomarker subgroups \eqn{B_1, \ldots, B_K}. Every enrolled patient's tumor gives
#' patient-derived organoids (PDOs) that are tested with each treatment open in the patient's
#' subgroup. A treatment is admissible in a subgroup if it is safe,
#' \eqn{\pi_{T,j} < \pi_{T,0} + \delta_T}, and effective,
#' \eqn{\pi_{R,j,k} > \pi_{R,0,k} + \delta_R}.
#'
#' The functions follow the steps of a trial.
#' \describe{
#'   \item{Design}{[boat_design()] specifies the design and [boat_boundary()] gives the cutoffs of
#'     the interim and final analyses.}
#'   \item{Conduct}{At an interim analysis, [boat_fit()] fits the model to the data observed so
#'     far and [boat_decide()] drops treatments that are not promising. [boat_rar_model()] fits the
#'     model used to randomize patients until the next analysis, and [boat_randomize()] gives the
#'     randomization probabilities of a new patient from the patient's PDO responses. At the end
#'     of the trial, [boat_decide()] with `final = TRUE` gives the recommended treatments.}
#'   \item{Simulation}{[boat_scenario()] and [boat_dgm()] specify the true probabilities and the
#'     data-generating values, [boat_generate_data()] generates trial data,
#'     [boat_simulate_trial()] simulates one trial, [boat_get_oc()] estimates the operating
#'     characteristics and [boat_calibrate()] calibrates the tuning parameter \eqn{\lambda}.}
#' }
#'
#' @keywords internal
#' @importFrom msm rtnorm
#' @importFrom stats aggregate dbeta dnorm integrate pbeta pnorm qnorm rbinom rgamma rnorm runif sd setNames
"_PACKAGE"
