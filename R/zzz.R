.onAttach <- function(...) {
  packageStartupMessage(paste('\nWelcome to {ndi} version ', utils::packageDescription('ndi')$Version, '\n> help(\'ndi\') # for documentation\n> citation(\'ndi\') # for how to cite\n', sep = ''), appendLF = TRUE)
}
