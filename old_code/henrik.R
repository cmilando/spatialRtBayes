## https://github.com/henrifnk/NoUTurn/tree/main/R

#' Sample No-u-Turn
#'
#' efficient NUTS introduced by Hoffman and Gelman
#'
#' @inheritParams  sample_nNUT
#' @export
sample_NUT <- function(init_position, stepsize, iteration, seed = 123L) {
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = iteration, style = 3)
  position <- init_position
  positions <- data.frame(matrix(ncol = length(init_position), nrow = iteration))
  for (iter in seq_len(iteration)) {
    # resample r0 ~ N(0, 1)
    momentum <- rnorm(length(position))
    #
    dens <- joint_log_density(position, momentum) # L(theta[m - 1]) - 0.5 * dot(r0, r0)
    if (is.na(dens)) {
      warning(paste("NUTS sampled NA in iteration", iter))
      dens <- 1
    }
    # resample u ~ Uniform(0, exp(L = density))
    slice <- runif(n = 1, min = 0, max = exp(dens))
    # initialize states
    states <- initialize_states(position, momentum)
    tree_depth <- 0L
    # while s = 1
    while (states$run) {
      ## choose a direction from the uniform set of {-1, 1}, left and right
      direction <- sample(c(-1, 1), 1)
      ####
      if (direction == -1L) {
        ## if -1 -> left
        states_prop <- build_tree(
          states$leftmost, 
          direction, tree_depth, stepsize, slice
        )
        states$leftmost <- states_prop$leftmost
      } else {
        ## if 1 -> right
        states_prop <- build_tree(
          states$rightmost, 
          direction, tree_depth, stepsize, slice
        )
        states$rightmost <- states_prop$rightmost
      }
      ####
      ## if s' = 1, aka states_prop
      if (states_prop$run) {
        #### metropolis accept state ish
        if (is.na(is_U_turn(states))) browser()
        # s <- s' II(left, right)
        states$run <- is_U_turn(states) * states_prop$run
        tree_ratio <- min(1, states_prop$count / states$count) #<<
        if (is.na(tree_ratio)) tree_ratio <- 0
        if (rbinom(1, 1, tree_ratio)) {
          states$valid_state <- states_prop$valid_state #<<
        }
      } else {
        break
      }
      ## n = n + n'
      states$count <- states_prop$count + states_prop$count #<<
      ## j = j + 1
      tree_depth <- tree_depth + 1
    }
    ###
    setTxtProgressBar(pb, iter)
    ###
    if (is.matrix(states$valid_state$position)) {
      positions[iter, ] <- position
      next
    }
    position <- as.numeric(states$valid_state$position)
    if (anyNA(position)) stop("position contains NAs")
    positions[iter, ] <- position
  }
  return(positions)
}


# _______________________________________________________________________________
#' Initialize state
#'
#' Helperfunction to intialize state
#' @inheritParams leapfrog
#' @param run is u-turn made?
#' CM: added NA
initialize_states <- function(position, momentum, run = 1L) {
  list(
    "valid_state" = structure(
      list(
        matrix(NA, nrow = 0, ncol = length(position)),
        matrix(NA, nrow = 0, ncol = length(position))
      ),
      names = c("position", "momentum")
    ),
    "rightmost" = list("position" = position, "momentum" = momentum),
    "leftmost" = list("position" = position, "momentum" = momentum),
    "run" = run, "count" = 0, "acceptance" = 0, "steps" = 1
  )
}

# ______________________________________________________________________________
#' Is a U Turn made?
#'
#' Investigates if trajectory makes a U-Turn (if >= 0)
#' @param state state variable containing rightmost and leftmost position/doubling
is_U_turn <- function(state) {
  momentum_l <- state$leftmost$momentum
  momentum_r <- state$rightmost$momentum
  if (anyNA(c(state$rightmost$momentum, state$leftmost$momentum))) {
    warning("momentum contains NA, trajectory aborted")
    return(FALSE)
  }
  position_distance <- state$rightmost$position - state$leftmost$position
  left <- sum(momentum_l * as.numeric(position_distance)) >= 0
  right <- sum(momentum_r * as.numeric(position_distance)) >= 0
  left * right
}
# _______________________________________________________________________________
#' Joint logarithmic density
#'
#' joint log density of position and momentum
#' @inheritParams leapfrog
joint_log_density <- function(position, momentum) {
  #
  log_dens_estimate <- log_posterior_density(position)
  # Oh wow is there a missing parenthesis in algorithm 3
  # L(theta) - 0.5 sum(r0 * r0)
  joint_dens <- log_dens_estimate - 0.5 * sum(momentum * momentum)
  # validation
  if (is.na(joint_dens)) {
    warning("joint density induced NAs, density is set to -Inf")
    return(-Inf)
  }
  joint_dens
}

# _______________________________________________________________________________
#' Lepfrog takes one step in trajectory
#'
#' Repatedly called by NUTS until run == FALSE, that is when a U-Turn is performed.
#'
#' @param position tajectory position state
#' @param momentum tajectory momentum state
#' @param stepsize stepsize \epsilon to Leapfrog step function
#'
#' @return
#'
leapfrog <- function(position, momentum, stepsize) {
  momentum <- momentum + (stepsize / 2) * gradient(position)
  position <- position + stepsize * momentum
  momentum <- momentum + (stepsize / 2) * gradient(position)
  return(list("position" = as.numeric(position), "momentum" = as.numeric(momentum)))
}

# _______________________________________________________________________________
#' Build leaf
build_leaf <- function(position_momentum, direction, stepsize, slice) {
  step <- leapfrog(position_momentum$position, position_momentum$momentum,
    stepsize = (direction * stepsize)
  )
  state <- initialize_states(step$position, step$momentum)
  dens <- joint_log_density(step$position, step$momentum)
  if (slice <= exp(dens)) {
    state$valid_state <- step
    state$count <- 1L
  }
  state$run <- (dens - log(slice)) >= -1e3
  if (is.na(state$run)) {
    warning("Deltamax induced NaN, moving in low posterior regions?")
    state$run <- 1L
  }
  state
}

# _______________________________________________________________________________
#' Build tree
build_tree <- function(position_momentum, direction, tree_depth, stepsize, slice) {
  if (tree_depth == 0L) {
    build_leaf(position_momentum, direction, stepsize, slice)
  } else {
    states <- build_tree(
      position_momentum, direction, tree_depth - 1L,
      stepsize, slice
    )
    if (states$run) { #<<
      if (direction == -1L) {
        states_prop <- build_tree(
          states$leftmost, direction,
          tree_depth - 1L, stepsize, slice
        )
        position_momentum <- states$leftmost <- states_prop$leftmost
      } else {
        states_prop <- build_tree(
          states$rightmost, direction,
          tree_depth - 1L, stepsize, slice
        )
        position_momentum <- states$rightmost <- states_prop$rightmost
      }
      tree_ratio <- states_prop$count / (states_prop$count + states$count) #<<
      if (is.na(tree_ratio)) tree_ratio <- 0
      if (rbinom(1, 1, tree_ratio)) {
        states$valid_state <- states_prop$valid_state #<<
      }
      states$count <- states_prop$count + states$count #<<
      states$run <- states_prop$run || is_U_turn(states)
    }
    return(states)
  }
}
