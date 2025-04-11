# /// Display functions /////////////////////////////////////////////

## ----
# Mini functions for CLI formatting
# NOTE: hard-coding ANSI escapes here since calling CLI commands in a list at
#       run-time does not evaluate the full ANSI string?
ngGrad1 <- function(rp) {
    sprintf(" %s ", rp)
}

ngGrad2 <- function(rp) {
    sprintf("\033[c2m %s \033[49m", rp)
}

ngGrad3 <- function(allele) {
    sprintf("\033[43m\033[30m\033[1m %s \033[22m\033[39m\033[49m", allele)
}

ngGrad4 <- function(allele) {
    sprintf("\033[44m\033[37m\033[1m %s \033[22m\033[39m\033[49m", allele)
}

