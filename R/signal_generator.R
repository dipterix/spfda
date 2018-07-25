#' @export
spfda.sample <- new.env()

local({
  gen_nois <- function(N, len, ar = c(0), sd = 0.05) {
    replicate(N, {
      n = len
      nois = rnorm(n + length(ar), sd = sd)
      ar = c(1, ar)
      sapply(seq_along(ar), function(ii) {
        nois[ii + (1:n) - 1] * ar[[ii]]
      }) %>%
        rowSums() %>%
        rev() ->
        nois
    })
  }

  sine_var <- function(n = 64,
                       ar = c(0),
                       sd = 0.05) {
    x = seq(0, 3, by = 0.001)
    a = seq(1, 100, length.out = length(x))
    s = c(rep(0, 500), sin(a * x), rep(0, 500))

    time = seq_along(s) / 1000

    nois = gen_nois(n, length(s), ar, sd)
    nois = nois + s

    X = t(nois)
    Y = rep(1, n)
    list(X = X,
         Y = Y,
         time = time)
  }


  sine <- function(init_phase = 0.1,
                   n = 64,
                   ar = c(0),
                   sd = 0.05) {
    # init_phase = 0.1; n = 32; ar = c(1,2); sd = 0.01
    x = seq(0, 3, by = 0.01) + init_phase
    s1 = c(rep(0, 50), sin(x * 2), rep(0, 50))

    time = seq_along(s1) / 100

    nois = gen_nois(n, length(s1), ar, sd)

    nois = nois + s1
    X = t(nois)
    Y = rep(1:2, each = n)
    list(X = X,
         Y = Y,
         time = time)
  }

  sine_thred <- function(init_phase = 0.1,
                         n = 64,
                         ar = c(0),
                         sd = 0.05) {
    x = seq(0, 3, by = 0.01) + init_phase
    s2 = c(rep(0, 50), sin(x * 2), rep(0, 50))
    s2 = s2 - sign(s2) * abs(sin(init_phase))
    time = seq_along(s2) / 100
    nois = gen_nois(n, length(s2), ar, sd)
    nois = nois + s2
    X = t(nois)
    Y = rep(1:2, each = n)
    # matplot(x = time, nois, col = Y, type = 'l')

    list(X = X,
         Y = Y,
         time = time)
  }
}, envir = spfda.sample)
