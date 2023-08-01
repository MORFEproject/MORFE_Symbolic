import Pkg

# This script install everything is needed
Pkg.add([
Pkg.PackageSpec(;name="Combinatorics", version="1.0.2"),
Pkg.PackageSpec(;name="LinearAlgebra"),
Pkg.PackageSpec(;name="SymPy", version="1.1.9"),
Pkg.PackageSpec(;name="Latexify", version="0.16.1")
])