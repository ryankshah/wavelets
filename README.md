# Wavelets

A comprehensive wavelet transform library for Elixir, providing implementations of various wavelet transforms and analysis tools.

[![Elixir CI](https://github.com/yourusername/wavelets/workflows/Elixir%20CI/badge.svg)](https://github.com/ryankshah/wavelets/actions)
<!-- [![codecov](https://codecov.io/gh/yourusername/wavelets/branch/main/graph/badge.svg)](https://codecov.io/gh/ryankshah/wavelets) -->
[![Hex.pm](https://img.shields.io/hexpm/v/wavelets.svg)](https://hex.pm/packages/wavelets)

## Features

- 1D, 2D and nD Forward and Inverse Discrete Wavelet Transform (DWT and IDWT)
- 1D, 2D and nD Multilevel DWT and IDWT
- 1D and 2D Stationary Wavelet Transform (Undecimated Wavelet Transform)
- 1D and 2D Wavelet Packet decomposition and reconstruction
- 1D Continuous Wavelet Transform
- Computing Approximations of wavelet and scaling functions
- Over 100 built-in wavelet filters and support for custom wavelets
- Single and double precision calculations
- Real and complex calculations
- Results compatible with Matlab Wavelet Toolbox (TM)

## Installation

Add `wavelets` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:wavelets, "~> 0.1.0"}
  ]
end
```

## Usage

### Basic DWT Example

```elixir
alias Wavelets.{DWT, IDWT, Filter}

# Create a signal
signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

# Get a wavelet filter (e.g., Daubechies-4)
filter = Filter.daubechies(2)

# Perform forward transform
{approximation, details} = DWT.forward_1d(signal, filter)

# Perform inverse transform
reconstructed = IDWT.inverse_1d(approximation, details, filter)
```

### Custom Wavelet Design

```elixir
alias Wavelets.CustomWavelet

# Create a custom orthogonal wavelet
{:ok, custom_filter} = CustomWavelet.create_orthogonal(
  [0.7071067811865476, 0.7071067811865476],
  name: "custom_haar"
)

# Use the custom filter
{approx, details} = DWT.forward_1d(signal, custom_filter)
```

### Analysis Tools

```elixir
alias Wavelets.Analysis

# Compute energy distribution
energy = Analysis.energy_distribution(coefficients)

# Compute entropy
entropy = Analysis.entropy(coefficients)

# Get statistical analysis
stats = Analysis.statistics(coefficients)
```

## Documentation

The full documentation can be found at [https://hexdocs.pm/wavelets](https://hexdocs.pm/wavelets).

## Testing

Run the test suite:

```bash
mix test
```

Run with coverage:

```bash
mix coveralls
```

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b feature/my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin feature/my-new-feature`)
5. Create new Pull Request

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.