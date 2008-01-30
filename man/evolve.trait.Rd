\name{evolve.trait}
\alias{evolve.trait}
\title{ Trait evolution under various modes and tempos }
\description{
  Evolves traits under a variety of tempos and modes (Brownian, bounded Brownian, ACDC, OU, proportional Brownian)
}
\usage{
evolve.trait(phy, x.root = 0, sigma = 1, trend = 0, bound = c(-10, 10), burst = rep(1, nrow(phy$edge)), pulse = 0, mu = rep(0, nrow(phy$edge)), theta = 1, gamma = 1, gamma.reset = 0, trait.mode = "Brownian", set.minbl = 1, show.plots = TRUE, use.bl = TRUE, debug = FALSE)
}
\arguments{
  \item{phy}{ ~~Describe \code{phy} here~~ }
  \item{x.root}{ ~~Describe \code{x.root} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
  \item{trend}{ ~~Describe \code{trend} here~~ }
  \item{bound}{ ~~Describe \code{bound} here~~ }
  \item{burst}{ ~~Describe \code{burst} here~~ }
  \item{pulse}{ ~~Describe \code{pulse} here~~ }
  \item{mu}{ ~~Describe \code{mu} here~~ }
  \item{theta}{ ~~Describe \code{theta} here~~ }
  \item{gamma}{ ~~Describe \code{gamma} here~~ }
  \item{gamma.reset}{ ~~Describe \code{gamma.reset} here~~ }
  \item{trait.mode}{ ~~Describe \code{trait.mode} here~~ }
  \item{set.minbl}{ ~~Describe \code{set.minbl} here~~ }
  \item{show.plots}{ ~~Describe \code{show.plots} here~~ }
  \item{use.bl}{ ~~Describe \code{use.bl} here~~ }
  \item{debug}{ ~~Describe \code{debug} here~~ }
}
\details{
  ...
}
\value{
  Vector of trait values for tips of phylo
}

\author{ David Ackerly <dackerly@berkeley.edu> and Steve Kembel <skembel@berkeley.edu>}

