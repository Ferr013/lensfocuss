# -----------------------------------------------------
# Fastell
# -----------------------------------------------------

add_library(fastell fastell.f)

# We compile fastell with O5 optimization. Should it be checked that
# this optimization is not too agressive, especially for modern Intel
# compilers?
target_compile_options(fastell PRIVATE -O5)
