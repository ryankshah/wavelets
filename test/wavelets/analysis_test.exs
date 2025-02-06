defmodule Wavelets.AnalysisTest do
  use ExUnit.Case
  doctest Wavelets.Analysis

  alias Wavelets.{Analysis, DWT, Filters}

  test "computes correct energy for 1D signal" do
    signal = [1.0, -1.0, 0.5, -0.5]
    energy = Analysis.energy_distribution(signal)

    # Expected energy density = (1² + (-1)² + 0.5² + (-0.5)²) / 4 = 2.5/4 = 0.625
    expected_energy = 0.625
    assert_in_delta expected_energy, energy, 1.0e-10
  end

  test "computes correct energy for wavelet coefficients" do
    signal = [1.0, 2.0, 3.0, 4.0]
    filter = Filters.Daubechies.get(1)
    {approx, details} = DWT.forward_1d(signal, filter)

    # For wavelet coefficients, we want raw energies because
    # DWT applies its own normalization
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))
    approx_energy = Enum.sum(Enum.map(approx, &(&1 * &1)))
    details_energy = Enum.sum(Enum.map(details, &(&1 * &1)))

    # Total energy should be preserved
    assert_in_delta signal_energy, approx_energy + details_energy, 1.0e-10
  end
end
