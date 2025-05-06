import Pkg;
Pkg.add("CSV");
Pkg.add("Tables");
Pkg.add("DataFrames");
Pkg.add("Polynomials");
Pkg.add("PolynomialRoots");
Pkg.add("LinearSolve")
Pkg.add("Roots");
Pkg.add("LinearAlgebra");

using CSV, Tables, DataFrames
using Polynomials, PolynomialRoots
using LinearSolve, LinearAlgebra
using Roots


abstract type mathematicalModel end
abstract type numericalModel end