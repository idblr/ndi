# Internal function for the Atkinson Index 
## Atkinson (1970) https://doi.org/10.1016/0022-0531(70)90039-6
## Returns NA value if only one smaller geography with population in a larger geography
## If denoting the Hölder mean
a_fun <- function(x, epsilon, omit_NAs, holder) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    if (holder == TRUE) {
      x_i <- xx$subgroup
      if (epsilon == 1) {
        A <- 1 - (exp(mean(log(stats::na.omit(x_i)), na.rm = TRUE)) / mean(x_i, na.rm = TRUE))
        return(A)
      } else {
        xxx <- (x_i / mean(x_i, na.rm = TRUE)) ^ (1 - epsilon)
        A <- 1 - mean(xxx, na.rm = TRUE) ^ (1 / (1 - epsilon))
        return(A)
      }
    } else {
      x_i <- xx$subgroup
      X <- sum(x_i, na.rm = TRUE)
      t_i <- xx$TotalPopE
      N <- sum(t_i, na.rm = TRUE)
      p_i <- x_i / t_i
      P <- X / N
      b <- epsilon
      A <- 1 - (P / (1 - P)) * abs(sum((1 - p_i) ^ (1 - b) * p_i ^ b * t_i / (P * N), na.rm = TRUE)) ^ (1 / (1 - b))
      return(A)
    }
  }
}

# Internal function for Absolute Centralization
## Duncan, Cuzzort, & Duncan (1961; LC:60007089)
## Returns NA value if only one smaller geography with population in a larger geography
ace_fun <- function(x, lgeom, crs, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, ALAND, oid) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(sf::st_drop_geometry(xx)), ] }
  if (nrow(sf::st_drop_geometry(xx)) < 2 || any(sf::st_drop_geometry(xx) < 0) || any(is.na(sf::st_drop_geometry(xx)))) {
    NA
  } else {
    L <- lgeom %>%
      dplyr::filter(GEOID == unique(xx$oid)) %>%
      sf::st_transform(crs = crs)
    C <- L %>%
      sf::st_geometry() %>%
      sf::st_centroid()
    A <- L %>% 
      sf::st_drop_geometry()
    xx <- xx %>% 
      sf::st_transform(crs = crs) %>%
      dplyr::mutate(d = sf::st_distance(sf::st_geometry(.), C)) %>%
      dplyr::arrange(d) %>%
      sf::st_drop_geometry()
    x_i <- xx$subgroup
    x_n <- sum(x_i, na.rm = TRUE)
    X_i <- cumsum(x_i / x_n)
    a_i <- xx$ALAND
    A_i <- cumsum(a_i / A$ALAND) 
    I_i <- matrix(c(seq(1, (length(x_i)-1), 1), seq(2, length(x_i), 1)), ncol = 2)
    Xi_1Ai <- sum(X_i[I_i[, 1]] * A_i[I_i[, 2]], na.rm = TRUE)
    XiA1_1 <- sum(X_i[I_i[, 2]] * A_i[I_i[, 1]], na.rm = TRUE)
    ACE <- Xi_1Ai - XiA1_1
    return(ACE)
  }
}

# Internal function for Absolute Clustering
## From Massey & Denton (1988) https://doi.org/10.1093/sf/67.2.281
## Returns NA value if only one smaller geography with population in a larger geography
acl_fun <- function(x, crs, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(sf::st_drop_geometry(xx)), ] }
  if (nrow(sf::st_drop_geometry(xx)) < 2 || any(sf::st_drop_geometry(xx) < 0) || any(is.na(sf::st_drop_geometry(xx)))) {
    NA
  } else {
    xx <- xx %>% sf::st_transform(crs = crs)
    d_ij <- suppressWarnings(sf::st_distance(sf::st_centroid(xx), sf::st_centroid(xx)))
    diag(d_ij) <- sqrt(0.6 * xx$ALAND)
    c_ij <- -d_ij %>% 
      units::set_units(value = km) %>%
      units::drop_units() %>%
      exp()
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    n <- length(x_i)
    t_i <- xx$TotalPopE
    num <- (sum(x_i / X, na.rm = TRUE) * sum(c_ij * x_i, na.rm = TRUE)) - ((X / n^2) * sum(c_ij, na.rm = TRUE))
    denom <- (sum(x_i / X, na.rm = TRUE) * sum(c_ij * t_i, na.rm = TRUE)) - ((X / n^2) * sum(c_ij, na.rm = TRUE))
    ACL <- num / denom
    return(ACL)
  }
}

# Internal function for Absolute Concentration
## From Massey & Denton (1988) https://doi.org/10.1093/sf/67.2.281
## Returns NA value if only one smaller geography with population in a larger geography
aco_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    a_i <- xx$ALAND
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    xx_tmp <- xx %>% 
      dplyr::arrange(ALAND) %>%
      dplyr::mutate(
        t_cs = cumsum(TotalPopE),
        n_1 = t_cs <= X,
      ) 
    if (!(TRUE %in% xx_tmp$n_1)) { 
      xx_1 <- xx_tmp %>% 
        dplyr::slice(1) 
    } else {
      xx_1 <- xx_tmp %>%
        dplyr::filter(n_1 == TRUE)
    }
    T_1 <- xx_1 %>%
      dplyr::summarise(
        T_1 = sum(TotalPopE, na.rm = TRUE)
      ) %>%
      unlist()
    xx_tmp <- xx %>% 
      dplyr::arrange(-ALAND) %>%
      dplyr::mutate(
        t_cs = cumsum(TotalPopE),
        n_2 = t_cs <= X,
      )
    if (!(TRUE %in% xx_tmp$n_2)) { 
      xx_2 <- xx_tmp %>% 
        dplyr::slice(1) 
    } else {
      xx_2 <- xx_tmp %>%
        dplyr::filter(n_2 == TRUE)
    }
    T_2 <- xx_2 %>%
      dplyr::summarise(
        T_2 = sum(TotalPopE, na.rm = TRUE)
      ) %>%
      unlist()
    num <- sum((x_i * a_i) / X, na.rm = TRUE) - sum((xx_1$TotalPopE * xx_1$ALAND) / T_1, na.rm = TRUE) 
    denom <- sum((xx_2$TotalPopE * xx_2$ALAND) / T_2, na.rm = TRUE) - sum((xx_1$TotalPopE * xx_1$ALAND) / T_1, na.rm = TRUE) 
    ACO_tmp <- (num / denom)
    if (is.infinite(ACO_tmp) | is.na(ACO_tmp)) { ACO_tmp <- 0 }
    ACO <- 1 - ACO_tmp
    return(ACO)
  }
}

# Internal function for the Dissimilarity Index 
## Duncan & Duncan (1955) https://doi.org/10.2307/2088328
## Returns NA value if only one smaller geography with population in a larger geography
ddd_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ref) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    n_i <- sum(x_i, na.rm = TRUE)
    y_i <- xx$subgroup_ref
    m_i <- sum(y_i, na.rm = TRUE)
    D <- 0.5 * sum(abs((x_i/n_i) - (y_i/m_i)), na.rm = TRUE)
    return(D)
  }
}

# Internal function for the aspatial Delta 
## Hoover (1941) https://10.1017/S0022050700052980
## Returns NA value if only one smaller geography with population in a larger geography
del_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    a_i <- xx$ALAND
    A <- sum(a_i, na.rm = TRUE)
    DEL <- 0.5 * sum(abs((x_i / X) - (a_i / A)), na.rm = TRUE)
    return(DEL)
  }
}

# Internal function for the Dissimilarity Index 
## James & Taeuber (1985) https://doi.org/10.2307/270845
## Returns NA value if only one smaller geography with population in a larger geography
djt_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    N <- sum(t_i, na.rm = TRUE)
    p_i <- x_i / t_i
    P <- X / N
    D <- sum(t_i * abs(p_i - P), na.rm = TRUE) / (2 * N * P * (1 - P))
    return(D)
  }
}

# Internal function for Distance Decay Isolation
## From Massey & Denton (1988) https://doi.org/10.1093/sf/67.2.281
## Returns NA value if only one smaller geography with population in a larger geography
dpxx_star_fun <- function(x, crs, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(sf::st_drop_geometry(xx)), ] }
  if (nrow(sf::st_drop_geometry(xx)) < 2 || any(sf::st_drop_geometry(xx) < 0) || any(is.na(sf::st_drop_geometry(xx)))) {
    NA
  } else {
    xx <- xx %>% sf::st_transform(crs = crs)
    x_i <- x_j <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    t_j <- xx$TotalPopE
    d_ij <- suppressWarnings(sf::st_distance(sf::st_centroid(xx), sf::st_centroid(xx)))
    diag(d_ij) <- sqrt(0.6 * xx$ALAND)
    c_ij <- -d_ij %>% 
      units::set_units(value = km) %>%
      units::drop_units() %>%
      exp()
    K_ij <- c_ij * t_j /  sum(c_ij * t_j, na.rm = TRUE)  
    DPxx_star <- sum(x_i / X, na.rm = TRUE) * sum(K_ij * x_j / t_j, na.rm = TRUE)  
    return(DPxx_star)
  }
}

# Internal function for the Gini Index 
## Gini (1921) https://doi.org/10.2307/2223319
## Returns NA value if only one smaller geography with population in a larger geography
g_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    N <- sum(t_i, na.rm = TRUE)
    p_i <- x_i / t_i
    P <- X / N
    titj <- apply(expand.grid(t_i, t_i), MARGIN = 1, FUN = prod)
    pipj <- apply(expand.grid(p_i, p_i), MARGIN = 1, FUN = diff)
    G <- sum(titj * abs(pipj), na.rm = TRUE)
    G <- G / (2 * N ^ 2 * P * (1 - P))
    return(G)
  }
}

# Internal function for Entropy 
## Theil (1972) https://doi.org/10.1080/0022250X.1971.9989795
## Returns NA value if only one smaller geography with population in a larger geography
## Note: Differs from Massey & Denton (1988) https://doi.org/10.1093/sf/67.2.281 
##       by taking the absolute value of (E-E_{i}) so extent of the output is 
##       {0, 1} as designed by Theil (1972) instead of {-Inf, Inf} as described in 
##       Massey & Denton (1988)
h_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    N <- sum(t_i, na.rm = TRUE)
    p_i <- x_i / t_i
    p_i[is.infinite(p_i)] <- 0
    P <- X / N
    if (is.infinite(P)) { P <- 0 }
    E_i <- p_i * log(1 / p_i) + (1 - p_i) * log(1 / (1 - p_i))
    E_i[is.infinite(E_i)] <- 0
    E <- P * log(1 / P) + (1 - P) * log(1 / (1 - P))
    if (is.infinite(E)) { E <- 0 }
    H_i <- t_i * abs(E - E_i) / (E * N)
    H_i[is.infinite(H_i)] <- NA
    H <- sum(H_i, na.rm = TRUE) 
    return(H)
  }
}

# Internal function for the aspatial Local Exposure & Isolation metric
# Bemanian & Beyer (2017) https://doi.org/10.1158/1055-9965.EPI-16-0926
## Returns NA value if only one smaller geography with population in a larger geography
lexis_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ixn, GEOID) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    p_im <- xx$subgroup / xx$TotalPopE
    if (anyNA(p_im)) { p_im[is.na(p_im)] <- 0 }
    p_in <- xx$subgroup_ixn / xx$TotalPopE
    if (anyNA(p_in)) { p_in[is.na(p_in) ] <- 0 }
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    y_i <- xx$subgroup_ixn
    Y <- sum(y_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    N <- sum(t_i, na.rm = TRUE)
    P_m <- X / N
    P_n <- Y / N
    LExIs <- car::logit(p_im * p_in) - car::logit(P_m * P_n)
    df <-  data.frame(LExIs = LExIs, GEOID = xx$GEOID)
    return(df)
  }
}

# Internal function for the aspatial Location Quotient 
## Sudano et al. (2013) https://doi.org/10.1016/j.healthplace.2012.09.015
## Returns NA value if only one smaller geography with population in a larger geography
lq_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, GEOID) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup # x_im
    t_i <- xx$TotalPopE # X_i
    p_i <- x_i / t_i # p_im
    X <- sum(x_i, na.rm = TRUE) # X_m
    N <- sum(t_i, na.rm = TRUE) # X
    if (anyNA(p_i)) { p_i[is.na(p_i)] <- 0 }
    LQ <- p_i / (X / N) # (x_im/X_i)/(X_m/X)
    df <-  data.frame(LQ = LQ, GEOID = xx$GEOID)
    return(df)
  }
}

# Internal function for Relative Centralization
## Duncan & Duncan (1955) https://doi.org/10.1086/221609
## Returns NA value if only one smaller geography with population in a larger geography
rce_fun <- function(x, lgeom, crs, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ref, oid) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(sf::st_drop_geometry(xx)), ] }
  if (nrow(sf::st_drop_geometry(xx)) < 2 || any(sf::st_drop_geometry(xx) < 0) || any(is.na(sf::st_drop_geometry(xx)))) {
    NA
  } else {
    C <- lgeom %>%
      dplyr::filter(GEOID == unique(xx$oid)) %>%
      sf::st_transform(crs = crs) %>%
      sf::st_geometry() %>%
      sf::st_centroid()
    xx <- xx %>% 
      sf::st_transform(crs = crs) %>%
      dplyr::mutate(d = sf::st_distance(sf::st_geometry(.), C)) %>%
      dplyr::arrange(d) %>%
      sf::st_drop_geometry()
    x_i <- xx$subgroup
    x_n <- sum(x_i, na.rm = TRUE)
    X_i <- cumsum(x_i / x_n)
    y_i <- xx$subgroup_ref
    y_n <- sum(y_i, na.rm = TRUE)
    Y_i <- cumsum(y_i / y_n)
    I_i <- matrix(c(seq(1, (length(x_i)-1), 1), seq(2, length(x_i), 1)), ncol = 2)
    Xi_1Yi <- sum(X_i[I_i[, 1]] * Y_i[I_i[, 2]], na.rm = TRUE)
    XiY1_1 <- sum(X_i[I_i[, 2]] * Y_i[I_i[, 1]], na.rm = TRUE)
    RCE <- Xi_1Yi - XiY1_1
    return(RCE)
  }
}

# Internal function for Relative Clustering
## From Massey & Denton (1988) https://doi.org/10.1093/sf/67.2.281
## Returns NA value if only one smaller geography with population in a larger geography
rcl_fun <- function(x, crs, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ref, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(sf::st_drop_geometry(xx)), ] }
  if (nrow(sf::st_drop_geometry(xx)) < 2 || any(sf::st_drop_geometry(xx) < 0) || any(is.na(sf::st_drop_geometry(xx)))) {
    NA
  } else {
    xx <- xx %>% sf::st_transform(crs = crs)
    d_ij <- suppressWarnings(sf::st_distance(sf::st_centroid(xx), sf::st_centroid(xx)))
    diag(d_ij) <- sqrt(0.6 * xx$ALAND)
    c_ij <- -d_ij %>% 
      units::set_units(value = km) %>%
      units::drop_units() %>%
      exp()
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    y_i <- xx$subgroup_ref
    Y <- sum(y_i, na.rm = TRUE)
    P_xx <- sum((x_i * x_i * c_ij) / X^2, na.rm = TRUE)
    P_yy <- sum((y_i * y_i * c_ij) / Y^2, na.rm = TRUE)
    RCL <- (P_xx / P_yy) - 1
    return(RCL)
  }
}

# Internal function for an index of spatial proximity 
## White (1986) https://doi.org/10.2307/3644339
## Returns NA value if only one smaller geography with population in a larger geography
sp_fun <- function(x, crs, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ref, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(sf::st_drop_geometry(xx)), ] }
  if (nrow(sf::st_drop_geometry(xx)) < 2 || any(sf::st_drop_geometry(xx) < 0) || any(is.na(sf::st_drop_geometry(xx)))) {
    NA
  } else {
    xx <- xx %>% sf::st_transform(crs = crs)
    d_ij <- suppressWarnings(sf::st_distance(sf::st_centroid(xx), sf::st_centroid(xx)))
    diag(d_ij) <- sqrt(0.6 * xx$ALAND)
    c_ij <- -d_ij %>% 
      units::set_units(value = km) %>%
      units::drop_units() %>%
      exp()
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    y_i <- xx$subgroup_ref
    Y <- sum(y_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    N <- sum(t_i, na.rm = TRUE)
    P_xx <- sum((x_i * x_i * c_ij) / X^2, na.rm = TRUE)
    P_xy <- sum((x_i * y_i * c_ij) / (X * Y), na.rm = TRUE)
    P_tt <- sum((t_i * t_i * c_ij) / N^2, na.rm = TRUE)
    SP <- ((X * P_xx) + (Y * P_xy)) / (N * P_tt)
    return(SP)
  }
}

# Internal function for Relative Concentration
## From Massey & Denton (1988) https://doi.org/10.1093/sf/67.2.281
## Returns NA value if only one smaller geography with population in a larger geography
rco_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ref, ALAND) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    a_i <- xx$ALAND
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    y_i <- xx$subgroup_ref
    Y <- sum(y_i, na.rm = TRUE)
    xx_tmp <- xx %>% 
      dplyr::arrange(ALAND) %>%
      dplyr::mutate(
        t_cs = cumsum(TotalPopE),
        n_1 = t_cs <= X,
      ) 
    if (!(TRUE %in% xx_tmp$n_1)) { 
      xx_1 <- xx_tmp %>% 
        dplyr::slice(1) 
    } else {
      xx_1 <- xx_tmp %>%
        dplyr::filter(n_1 == TRUE)
    }
    T_1 <- xx_1 %>%
      dplyr::summarise(
        T_1 = sum(TotalPopE, na.rm = TRUE)
      ) %>%
      unlist()
    xx_tmp <- xx %>% 
      dplyr::arrange(-ALAND) %>%
      dplyr::mutate(
        t_cs = cumsum(TotalPopE),
        n_2 = t_cs <= X,
      )
    if (!(TRUE %in% xx_tmp$n_2)) { 
      xx_2 <- xx_tmp %>% 
        dplyr::slice(1) 
    } else {
      xx_2 <- xx_tmp %>%
        dplyr::filter(n_2 == TRUE)
    }
    T_2 <- xx_2 %>%
      dplyr::summarise(
        T_2 = sum(TotalPopE, na.rm = TRUE)
      ) %>%
      unlist()
    num <- sum((x_i * a_i) / X, na.rm = TRUE) / sum((y_i * a_i) / Y, na.rm = TRUE)
    denom <- sum((xx_1$TotalPopE * xx_1$ALAND) / T_1, na.rm = TRUE) / sum((xx_2$TotalPopE * xx_2$ALAND) / T_2, na.rm = TRUE)
    RCO <- (num - 1) / (denom - 1)
    if (is.na(RCO)) { RCO <- 0 }
    if (is.infinite(RCO) & sign(RCO) == -1 ) { RCO <- -1 }
    if (is.infinite(RCO) & sign(RCO) == 1) { RCO <- 1 }   
    # if (is.finite(RCO) & RCO < -1) { RCO <- -1 }
    # if (is.finite(RCO) & RCO > -1) { RCO <- 1 }    
    return(RCO)
  }
}

# Internal function for the aspatial Correlation Ratio 
## White (1986) https://doi.org/10.2307/3644339
## Returns NA value if only one smaller geography with population in a larger geography
v_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    N <- sum(t_i, na.rm = TRUE)
    xPx_star <- sum((x_i / X) * (x_i / t_i), na.rm = TRUE)
    P <- X / N
    V <- (xPx_star - P) / (1 - P)
    return(V)
  }
}

# Internal function for the aspatial Isolation Index 
## Lieberson (1981) ISBN-13:978-1-032-53884-6
## Returns NA value if only one smaller geography with population in a larger geography
xpx_star_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    t_i <- xx$TotalPopE
    xPx_star <- sum((x_i / X) * (x_i / t_i), na.rm = TRUE)
    return(xPx_star)
  }
}

# Internal function for the aspatial Interaction Index 
## Bell (1954) https://doi.org/10.2307/2574118
## Returns NA value if only one smaller geography with population in a larger geography
xpy_star_fun <- function(x, omit_NAs) {
  xx <- x %>%
    dplyr::select(TotalPopE, subgroup, subgroup_ixn) %>%
    dplyr::filter(TotalPopE > 0)
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(xx) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    x_i <- xx$subgroup
    X <- sum(x_i, na.rm = TRUE)
    y_i <- xx$subgroup_ixn
    t_i <- xx$TotalPopE
    xPy_star <- sum((x_i / X) * (y_i / t_i), na.rm = TRUE)
    return(xPy_star)
  }
}
