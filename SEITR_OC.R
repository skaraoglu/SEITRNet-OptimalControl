# Network generation function for SEITR simulation
generate_network <- function(network_type = "ER", n = 20, n_par1 = 0.3, n_par2 = 10) {
  if (network_type == "ER") {
    # Erdos-Renyi random graph
    g <- igraph::erdos.renyi.game(n = n, p.or.m = n_par1, type = "gnp", directed = FALSE)
  } else if (network_type == "BA") {
    # Barabasi-Albert preferential attachment
    g <- igraph::sample_pa(n = n, power = 1, m = max(1, round(n_par1 * n)), directed = FALSE)
  } else if (network_type == "WS") {
    # Watts-Strogatz small-world
    g <- igraph::sample_smallworld(dim = 1, size = n, nei = n_par2, p = n_par1)
  } else if (network_type == "LN") {
    # Lattice network
    g <- igraph::make_lattice(length = n, dim = 1, nei = n_par1, directed = FALSE, mutual = TRUE, circular = TRUE)
  } else if (network_type == "RR") {
    # Random regular graph
    g <- igraph::sample_k_regular(no.of.nodes = n, k = n_par1, directed = FALSE)
  } else {
    stop("Unknown network type")
  }
  return(g)
}

# SEITR node status assignment and transition functions
initialize_node_statuses <- function(g, S, E, I, T, R) {
  n <- igraph::vcount(g)
  # Ensure the sum does not exceed n
  stopifnot(S + E + I + T + R <= n)
  idx <- 1:n
  # Randomly sample unique indices for each compartment
  idx_I <- sample(idx, I)
  idx_rest <- setdiff(idx, idx_I)
  idx_E <- sample(idx_rest, E)
  idx_rest <- setdiff(idx_rest, idx_E)
  idx_T <- sample(idx_rest, T)
  idx_rest <- setdiff(idx_rest, idx_T)
  idx_R <- sample(idx_rest, R)
  idx_rest <- setdiff(idx_rest, idx_R)
  idx_S <- idx_rest # Remaining nodes are susceptible

  V(g)$status <- "S"
  V(g)$status[idx_I] <- "I"
  V(g)$status[idx_E] <- "E"
  V(g)$status[idx_T] <- "T"
  V(g)$status[idx_R] <- "R"
  return(g)
}

# Transition function for a single node (to be used in simulation loop)
update_node_status <- function(g, node, beta1, beta2, beta3, alpha1, alpha2, u1, N) {
  status <- V(g)[node]$status
  rand <- runif(1)
  if (status == "S") {
    # Infection probability depends on infected neighbors
    infected_neighbors <- sum(V(g)[neighbors(g, node)]$status == "I")
    if (infected_neighbors > 0 && rand < beta1 * infected_neighbors / N) {
      V(g)[node]$status <- "E"
    }
  } else if (status == "E" && rand < beta2) {
    V(g)[node]$status <- "I"
  } else if (status == "I") {
    rand2 <- runif(1)
    # Treatment probability includes optimal control
    if (rand < beta3) {
      V(g)[node]$status <- "R"
    } else if (rand2 < alpha1 + u1) {
      V(g)[node]$status <- "T"
    }
  } else if (status == "T" && rand < alpha2) {
    V(g)[node]$status <- "R"
  }
  # Death transitions can be handled separately in the main loop
  return(g)
}

# Helper function to connect new nodes according to network type
connect_new_nodes <- function(g, new_ids, network_type, n_par1, n_par2) {
  if (network_type == "ER") {
    # Connect each new node to existing nodes with probability n_par1
    for (new_node in new_ids) {
      existing_nodes <- setdiff(seq_len(igraph::vcount(g)), new_node)
      for (target in existing_nodes) {
        if (runif(1) < n_par1) {
          g <- igraph::add_edges(g, c(new_node, target))
        }
      }
    }
  } else if (network_type == "BA") {
    # Preferential attachment: connect to nodes with probability proportional to degree
    for (new_node in new_ids) {
      degrees <- igraph::degree(g)
      existing_nodes <- setdiff(seq_len(igraph::vcount(g)), new_node)
      prob <- degrees[existing_nodes] / sum(degrees[existing_nodes])
      m <- max(1, round(n_par1 * length(existing_nodes)))
      targets <- sample(existing_nodes, m, prob = prob, replace = FALSE)
      for (target in targets) {
        g <- igraph::add_edges(g, c(new_node, target))
      }
    }
  } else if (network_type == "WS") {
    # Connect to k nearest neighbors and rewire with probability n_par1
    k <- n_par2
    for (new_node in new_ids) {
      existing_nodes <- setdiff(seq_len(igraph::vcount(g)), new_node)
      if (length(existing_nodes) >= k) {
        neighbors <- sample(existing_nodes, k)
        for (neighbor in neighbors) {
          g <- igraph::add_edges(g, c(new_node, neighbor))
        }
        # Rewire edges with probability n_par1
        for (neighbor in neighbors) {
          if (runif(1) < n_par1) {
            possible_nodes <- setdiff(existing_nodes, neighbor)
            if (length(possible_nodes) > 0) {
              new_neighbor <- sample(possible_nodes, 1)
              g <- igraph::delete_edges(g, igraph::get.edge.ids(g, c(new_node, neighbor)))
              g <- igraph::add_edges(g, c(new_node, new_neighbor))
            }
          }
        }
      }
    }
  } else if (network_type == "LN") {
    # Connect to n_par1 nearest neighbors in a ring
    for (new_node in new_ids) {
      existing_nodes <- setdiff(seq_len(igraph::vcount(g)), new_node)
      k <- n_par1
      if (length(existing_nodes) >= k) {
        neighbors <- sample(existing_nodes, k)
        for (neighbor in neighbors) {
          g <- igraph::add_edges(g, c(new_node, neighbor))
        }
      }
    }
  } else if (network_type == "RR") {
    # Connect each new node to n_par1 random existing nodes
    for (new_node in new_ids) {
      existing_nodes <- setdiff(seq_len(igraph::vcount(g)), new_node)
      if (length(existing_nodes) >= n_par1) {
        targets <- sample(existing_nodes, n_par1)
        for (target in targets) {
          g <- igraph::add_edges(g, c(new_node, target))
        }
      }
    }
  }
  return(g)
}

# Main SEITR network simulation function with optimal control and network-based node connection
run_seitr_network_simulation <- function(
  network_type, n, n_par1, n_par2,
  Lambda, beta1, beta2, beta3, alpha1, alpha2, delta_I, delta_T, mu,
  init_S, init_E, init_I, init_T, init_R, t_max, u1_profile, num_exp = 20, verbose = FALSE
) {
  times <- seq(0, t_max, by = 1)
  results_list <- vector("list", num_exp)
  
  for (exp in 1:num_exp) {
    g <- generate_network(network_type, n, n_par1, n_par2)
    g <- initialize_node_statuses(g, init_S, init_E, init_I, init_T, init_R)
    
    S_count <- numeric(length(times) + 1)
    E_count <- numeric(length(times) + 1)
    I_count <- numeric(length(times) + 1)
    T_count <- numeric(length(times) + 1)
    R_count <- numeric(length(times) + 1)
    N_count <- numeric(length(times) + 1)
    control_cost <- numeric(length(times) + 1)

    # Record initial counts at time 0 (before any transitions)
    S_count[1] <- sum(V(g)$status == "S")
    E_count[1] <- sum(V(g)$status == "E")
    I_count[1] <- sum(V(g)$status == "I")
    T_count[1] <- sum(V(g)$status == "T")
    R_count[1] <- sum(V(g)$status == "R")
    N_count[1] <- igraph::vcount(g)
    control_cost[1] <- w1 * u1_profile[1]^2
    
    for (t in seq_along(times)) {
      N <- igraph::vcount(g)
      u1 <- u1_profile[t]
      
      # Update node statuses
      for (node in seq_len(N)) {
        g <- update_node_status(g, node, beta1, beta2, beta3, alpha1, alpha2, u1, N)
      }
      
      # Disease-induced and natural deaths
      infected_nodes <- which(V(g)$status == "I")
      num_remove_I <- rbinom(1, length(infected_nodes), delta_I)
      if (num_remove_I > 0 && length(infected_nodes) > 0) {
        remove_I <- sample(infected_nodes, num_remove_I)
        g <- igraph::delete_vertices(g, remove_I)
      }
      treated_nodes <- which(V(g)$status == "T")
      num_remove_T <- rbinom(1, length(treated_nodes), delta_T)
      if (num_remove_T > 0 && length(treated_nodes) > 0) {
        remove_T <- sample(treated_nodes, num_remove_T)
        g <- igraph::delete_vertices(g, remove_T)
      }
      all_nodes <- seq_len(igraph::vcount(g))
      num_remove_mu <- rbinom(1, length(all_nodes), mu)
      if (num_remove_mu > 0 && length(all_nodes) > 0) {
        remove_mu <- sample(all_nodes, num_remove_mu)
        g <- igraph::delete_vertices(g, remove_mu)
      }
      
      # Recruitment: add new susceptible nodes (Lambda) and connect them according to network type
      num_add <- rpois(1, Lambda)
      if (num_add > 0) {
        g <- igraph::add_vertices(g, num_add)
        new_ids <- (igraph::vcount(g) - num_add + 1):igraph::vcount(g)
        V(g)[new_ids]$status <- "S"
        g <- connect_new_nodes(g, new_ids, network_type, n_par1, n_par2)
      }
      
      # Record counts
      S_count[t + 1] <- sum(V(g)$status == "S")
      E_count[t + 1] <- sum(V(g)$status == "E")
      I_count[t + 1] <- sum(V(g)$status == "I")
      T_count[t + 1] <- sum(V(g)$status == "T")
      R_count[t + 1] <- sum(V(g)$status == "R")
      N_count[t + 1] <- igraph::vcount(g)
      control_cost[t + 1] <- w1 * u1^2
    }
    
    results_list[[exp]] <- data.frame(
      time = c(0, times),
      S = S_count,
      E = E_count,
      I = I_count,
      T = T_count,
      R = R_count,
      N = N_count,
      control_cost = control_cost
    )
  }
  

  # Calculate mean, min, max for each compartment at each time step
  results_array <- array(NA, dim = c(length(c(0, times)), 8, num_exp))
  for (i in 1:num_exp) {
    results_array[,,i] <- as.matrix(results_list[[i]])
  }
  avg_df <- as.data.frame(apply(results_array, c(1,2), mean))
  colnames(avg_df) <- colnames(results_list[[1]])
  avg_df$time <- results_list[[1]]$time

  min_df <- as.data.frame(apply(results_array, c(1,2), min))
  colnames(min_df) <- colnames(results_list[[1]])
  min_df$time <- results_list[[1]]$time

  max_df <- as.data.frame(apply(results_array, c(1,2), max))
  colnames(max_df) <- colnames(results_list[[1]])
  max_df$time <- results_list[[1]]$time

  return(list(
    avg = avg_df,
    min = min_df,
    max = max_df,
    all = results_list
  ))
}

# Objective functional calculation for SEITR network simulation results
calculate_objective_functional <- function(results_df) {
  # Simpson's rule for numerical integration over time
  n <- nrow(results_df)
  if (n %% 2 == 0) n <- n - 1  # Simpson's rule requires odd number of intervals
  times <- results_df$time[1:n]
  E <- results_df$E[1:n]
  I <- results_df$I[1:n]
  control_cost <- results_df$control_cost[1:n]
  
  # Simpson's rule weights
  weights <- rep(2, n)
  weights[1] <- weights[n] <- 1
  weights[seq(2, n-1, by=2)] <- 4
  
  h <- mean(diff(times))
  J_E <- (h/3) * sum(weights * E)
  J_I <- (h/3) * sum(weights * I)
  J_W <- (h/3) * sum(weights * control_cost)
  
  J_total <- J_E + J_I + J_W
  return(list(J_E = J_E, J_I = J_I, J_W = J_W, J_total = J_total))
}

network_objective <- function(u1_profile) {
  # Ensure u1_profile is within bounds
  u1_profile <- pmin(zeta, pmax(0, u1_profile))
  # Run simulation
  results <- run_seitr_network_simulation(
    network_type = network_type, n = n, n_par1 = n_par1, n_par2 = n_par2,
    Lambda = Lambda, beta1 = beta1, beta2 = beta2, beta3 = beta3, alpha1 = alpha1, alpha2 = alpha2,
    delta_I = delta_I, delta_T = delta_T, mu = mu,
    init_S = init_S, init_E = init_E, init_I = init_I, init_T = init_T, init_R = init_R,
    t_max = t_max, u1_profile = u1_profile, num_exp = 20
  )
  J <- calculate_objective_functional(results$avg)$J_total
  return(J)
}

# Helper for compartment plot
plot_compartment <- function(df, before, after, title, color) {
  ggplot(df, aes(x = Time)) +
    geom_line(aes(y = !!as.name(before)), color = "black", size = 1, linetype = "solid") +
    geom_line(aes(y = !!as.name(after)), color = color, size = 1, linetype = "solid") +
    labs(title = title, x = "Time (Days)", y = title) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    theme(legend.position = "none")
}

# SEITR compartments with shaded min/max area
plot_compartment_band <- function(df_avg, df_min, df_max, comp, color, title) {
  df <- data.frame(
    time = df_avg$time,
    avg = df_avg[[comp]],
    min = df_min[[comp]],
    max = df_max[[comp]]
  )
  ggplot(df, aes(x = time)) +
    geom_ribbon(aes(ymin = min, ymax = max), fill = color, alpha = 0.2) +
    geom_line(aes(y = avg), color = color, size = 1) +
    labs(title = paste0(title, " (Network Opt)"), x = "Time", y = "Count") +
    theme_minimal()
}

# SEITR compartments: overlay ODE before/after optimization, add min/max bands
plot_compartment_band_overlay <- function(df_avg, df_min, df_max, comp, color, title, ode_results, after, before) {
  df <- data.frame(
    time = df_avg$time,
    avg = df_avg[[comp]],
    min = df_min[[comp]],
    max = df_max[[comp]]
  )
  ggplot(df, aes(x = time)) +
    geom_ribbon(aes(ymin = min, ymax = max), fill = color, alpha = 0.2) +
    geom_line(aes(y = avg), color = color, size = 1.2) +
    geom_line(data = ode_results, aes(x = time, y = !!as.name(after)), color = "black", size = 0.75, linetype = "dashed") +
    geom_line(data = ode_results, aes(x = time, y = !!as.name(before)), color = "gray20", size = 0.75, linetype = "dotted") +
    labs(title = title, x = "Time", y = "Count") +
    theme_minimal()
}

# No intervals: optimize u1 at every time step
expand_u1 <- function(params, total_length) {
  params[1:total_length]
}

objective_wrapper <- function(params) {
  u1_profile <- expand_u1(params, total_steps)
  network_objective(u1_profile)
}

#' Solve SEITR ODE Optimal Control Problem
#'
#' @param Lambda Recruitment rate
#' @param beta1 Transmission rate (S to E)
#' @param beta2 Progression rate (E to I)
#' @param beta3 Recovery rate (I to R)
#' @param alpha1 Treatment rate (I to T)
#' @param alpha2 Recovery rate from treatment (T to R)
#' @param delta_I Disease-induced death rate (I)
#' @param delta_T Disease-induced death rate (T)
#' @param mu Natural death rate
#' @param h Time step size
#' @param t Time vector
#' @param w1 Weight for control cost
#' @param zeta Upper bound for control
#' @param delta Convergence threshold
#' @param S0 Initial S
#' @param E0 Initial E
#' @param I0 Initial I
#' @param T0 Initial T
#' @param R0 Initial R
#' @return List with time, state variables, control, adjoints, and objective functional
solve_seitr_optimal_control <- function(
  Lambda, beta1, beta2, beta3, alpha1, alpha2, delta_I, delta_T, mu,
  h, t, w1, zeta, delta,
  S0 = 10, E0 = 5, I0 = 3, T0 = 1, R0 = 1
) {
  L <- length(t)
  S <- E <- I <- T <- R <- N <- numeric(L)
  S[1] <- S0; E[1] <- E0; I[1] <- I0; T[1] <- T0; R[1] <- R0
  N[1] <- S[1] + E[1] + I[1] + T[1] + R[1]
  u1 <- rep(0.03, L)
  lambda1 <- lambda2 <- lambda3 <- lambda4 <- lambda5 <- numeric(L)
  lambda1[L] <- lambda2[L] <- lambda3[L] <- lambda4[L] <- lambda5[L] <- 0
  itr <- 0
  test <- -1
  FV <- numeric()
  itrs <- numeric()
  # System dynamics
  f1 <- function(S, I, N)     Lambda - (beta1 * S * I) / N - mu * S
  f2 <- function(S, I, N, E)  (beta1 * S * I) / N - (beta2 + mu) * E
  f3 <- function(E, I, u1)    beta2 * E - (beta3 + mu + delta_I + u1) * I
  f4 <- function(u1, I, T)    u1 * I - (mu + delta_T + alpha2) * T
  f5 <- function(I, T, R)     beta3 * I + alpha2 * T - mu * R
  f6 <- function(N, I, T)     Lambda - N * mu - delta_I * I - delta_T * T
  # Adjoint equations
  g1 <- function(lambda1, I, N, lambda2) (beta1 * I * lambda1) / N + mu * lambda1 - (lambda2 * beta1 * I) / N
  g2 <- function(lambda2, lambda3) -1 + beta2 * lambda2 + mu * lambda2 - beta2 * lambda3
  g3 <- function(lambda1, S, N, lambda2, lambda3, u1, lambda4, lambda5) -1 + (beta1 * S * lambda1) / N - (beta1 * S * lambda2) / N +
      beta3 * lambda3 + mu * lambda3 + delta_I * lambda3 + u1 * lambda3 - lambda4 * u1 - lambda5 * beta3
  g4 <- function(lambda4, lambda5) mu * lambda4 + delta_T * lambda4 + alpha2 * lambda4 - alpha2 * lambda5
  g5 <- function(lambda5) lambda5 * mu
  while (test < 0) {
    itr <- itr + 1
    oldu1 <- u1
    oldS <- S; oldE <- E; oldI <- I; oldT <- T; oldR <- R
    oldlambda1 <- lambda1; oldlambda2 <- lambda2; oldlambda3 <- lambda3
    oldlambda4 <- lambda4; oldlambda5 <- lambda5
    # Forward integration for state variables (RK4)
    for (i in 1:(L-1)) {
      m1 <- f1(S[i], I[i], N[i])
      n1 <- f2(S[i], I[i], N[i], E[i])
      o1 <- f3(E[i], I[i], u1[i])
      p1 <- f4(u1[i], I[i], T[i])
      q1 <- f5(I[i], T[i], R[i])
      r1 <- f6(N[i], I[i], T[i])
      m2 <- f1(S[i] + 0.5*h*m1, I[i] + 0.5*h*o1, N[i] + 0.5*h*r1)
      n2 <- f2(S[i] + 0.5*h*m1, I[i] + 0.5*h*o1, N[i] + 0.5*h*r1, E[i] + 0.5*h*n1)
      o2 <- f3(E[i] + 0.5*h*n1, I[i] + 0.5*h*o1, u1[i] + 0.5*h)
      p2 <- f4(u1[i] + 0.5*h, I[i] + 0.5*h*o1, T[i] + 0.5*h*p1)
      q2 <- f5(I[i] + 0.5*h*o1, T[i] + 0.5*h*p1, R[i] + 0.5*h*q1)
      r2 <- f6(N[i] + 0.5*h*r1, I[i] + 0.5*h*o1, T[i] + 0.5*h*p1)
      m3 <- f1(S[i] + 0.5*h*m2, I[i] + 0.5*h*o2, N[i] + 0.5*h*r2)
      n3 <- f2(S[i] + 0.5*h*m2, I[i] + 0.5*h*o2, N[i] + 0.5*h*r2, E[i] + 0.5*h*n2)
      o3 <- f3(E[i] + 0.5*h*n2, I[i] + 0.5*h*o2, u1[i] + 0.5*h)
      p3 <- f4(u1[i] + 0.5*h, I[i] + 0.5*h*o2, T[i] + 0.5*h*p2)
      q3 <- f5(I[i] + 0.5*h*o2, T[i] + 0.5*h*p2, R[i] + 0.5*h*q2)
      r3 <- f6(N[i] + 0.5*h*r2, I[i] + 0.5*h*o2, T[i] + 0.5*h*p2)
      m4 <- f1(S[i] + h*m3, I[i] + h*o3, N[i] + h*r3)
      n4 <- f2(S[i] + h*m3, I[i] + h*o3, N[i] + h*r3, E[i] + h*n3)
      o4 <- f3(E[i] + h*n3, I[i] + h*o3, u1[i] + h)
      p4 <- f4(u1[i] + h, I[i] + h*o3, T[i] + h*p3)
      q4 <- f5(I[i] + h*o3, T[i] + h*p3, R[i] + h*q3)
      r4 <- f6(N[i] + h*r3, I[i] + h*o3, T[i] + h*p3)
      S[i+1] <- S[i] + (h/6)*(m1 + 2*m2 + 2*m3 + m4)
      E[i+1] <- E[i] + (h/6)*(n1 + 2*n2 + 2*n3 + n4)
      I[i+1] <- I[i] + (h/6)*(o1 + 2*o2 + 2*o3 + o4)
      T[i+1] <- T[i] + (h/6)*(p1 + 2*p2 + 2*p3 + p4)
      R[i+1] <- R[i] + (h/6)*(q1 + 2*q2 + 2*q3 + q4)
      N[i+1] <- N[i] + (h/6)*(r1 + 2*r2 + 2*r3 + r4)
    }
    if (itr == 1) {
      S1 <- S; E1 <- E; I1 <- I; T1 <- T; R1 <- R; N1 <- N
    }
    # Backward integration for adjoint variables (RK4, backward in time)
    for (i in 1:(L-1)) {
      j <- L + 1 - i
      k1 <- g1(lambda1[j], I[j], N[j], lambda2[j])
      l1 <- g2(lambda2[j], lambda3[j])
      c1 <- g3(lambda1[j], S[j], N[j], lambda2[j], lambda3[j], u1[j], lambda4[j], lambda5[j])
      v1 <- g4(lambda4[j], lambda5[j])
      d1 <- g5(lambda5[j])
      k2 <- g1(lambda1[j] - 0.5*h*k1, I[j] - 0.5*h, N[j] - 0.5*h, lambda2[j] - 0.5*h*l1)
      l2 <- g2(lambda2[j] - 0.5*h*l1, lambda3[j] - 0.5*h*c1)
      c2 <- g3(lambda1[j] - 0.5*h*k1, S[j] - 0.5*h, N[j] - 0.5*h, lambda2[j] - 0.5*h*l1, lambda3[j] - 0.5*h*c1, u1[j] - 0.5*h, lambda4[j] - 0.5*h*v1, lambda5[j] - 0.5*h*d1)
      v2 <- g4(lambda4[j] - 0.5*h*v1, lambda5[j] - 0.5*h*d1)
      d2 <- g5(lambda5[j] - 0.5*h*d1)
      k3 <- g1(lambda1[j] - 0.5*h*k2, I[j] - 0.5*h, N[j] - 0.5*h, lambda2[j] - 0.5*h*l2)
      l3 <- g2(lambda2[j] - 0.5*h*l2, lambda3[j] - 0.5*h*c2)
      c3 <- g3(lambda1[j] - 0.5*h*k2, S[j] - 0.5*h, N[j] - 0.5*h, lambda2[j] - 0.5*h*l2, lambda3[j] - 0.5*h*c2, u1[j] - 0.5*h, lambda4[j] - 0.5*h*v2, lambda5[j] - 0.5*h*d2)
      v3 <- g4(lambda4[j] - 0.5*h*v2, lambda5[j] - 0.5*h*d2)
      d3 <- g5(lambda5[j] - 0.5*h*d2)
      k4 <- g1(lambda1[j] - h*k3, I[j] - h, N[j] - h, lambda2[j] - h*l3)
      l4 <- g2(lambda2[j] - h*l3, lambda3[j] - h*c3)
      c4 <- g3(lambda1[j] - h*k3, S[j] - h, N[j] - h, lambda2[j] - h*l3, lambda3[j] - h*c3, u1[j] - h, lambda4[j] - h*v3, lambda5[j] - h*d3)
      v4 <- g4(lambda4[j] - h*v3, lambda5[j] - h*d3)
      d4 <- g5(lambda5[j] - h*d3)
      lambda1[j-1] <- lambda1[j] - (h/6)*(k1 + 2*k2 + 2*k3 + k4)
      lambda2[j-1] <- lambda2[j] - (h/6)*(l1 + 2*l2 + 2*l3 + l4)
      lambda3[j-1] <- lambda3[j] - (h/6)*(c1 + 2*c2 + 2*c3 + c4)
      lambda4[j-1] <- lambda4[j] - (h/6)*(v1 + 2*v2 + 2*v3 + v4)
      lambda5[j-1] <- lambda5[j] - (h/6)*(d1 + 2*d2 + 2*d3 + d4)
    }
    # Simpson's rule for objective functional
    JE <- E[1] + E[L] + 4*sum(E[seq(2, L-1, by=2)]) + 2*sum(E[seq(3, L-2, by=2)])
    JI <- I[1] + I[L] + 4*sum(I[seq(2, L-1, by=2)]) + 2*sum(I[seq(3, L-2, by=2)])
    JW <- 0.5*w1*(u1[1]^2 + u1[L]^2 + 4*sum(u1[seq(2, L-1, by=2)]^2) + 2*sum(u1[seq(3, L-2, by=2)]^2))
    J <- (h/3)*(JE + JI + JW)
    FV[itr] <- J
    itrs[itr] <- itr
    # Control update
    ustar <- pmin(zeta, pmax(0, (I * (lambda3 - lambda4)) / w1))
    u1 <- 0.5 * (ustar + oldu1)
    # Convergence check
    temp1 <- delta*sum(abs(u1)) - sum(abs(oldu1-u1))
    temp2 <- delta*sum(abs(S)) - sum(abs(oldS-S))
    temp3 <- delta*sum(abs(E)) - sum(abs(oldE-E))
    temp4 <- delta*sum(abs(I)) - sum(abs(oldI-I))
    temp5 <- delta*sum(abs(T)) - sum(abs(oldT-T))
    temp6 <- delta*sum(abs(R)) - sum(abs(oldR-R))
    temp7 <- delta*sum(abs(lambda1)) - sum(abs(oldlambda1-lambda1))
    temp8 <- delta*sum(abs(lambda2)) - sum(abs(oldlambda2-lambda2))
    temp9 <- delta*sum(abs(lambda3)) - sum(abs(oldlambda3-lambda3))
    temp10 <- delta*sum(abs(lambda4)) - sum(abs(oldlambda4-lambda4))
    temp11 <- delta*sum(abs(lambda5)) - sum(abs(oldlambda5-lambda5))
    test <- min(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11)
    # Optionally print progress: cat(itr, "   ", test, "\n")
  }
  return(list(
    t = t,
    S = S, E = E, I = I, T = T, R = R, N = N,
    S1 = S1, E1 = E1, I1 = I1, T1 = T1, R1 = R1, N1 = N1,
    u1 = u1,
    lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4, lambda5 = lambda5,
    FV = FV, itrs = itrs,
    J = J
  ))
}