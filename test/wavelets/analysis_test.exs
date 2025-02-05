defmodule Wavelets.AnalysisTest do
  use ExUnit.Case
  doctest Wavelets.Analysis

  alias Wavelets.{Analysis, DWT, Filters}

  test "computes correct energy for 1D signal" do
    signal = [1.0, -1.0, 0.5, -0.5]
    energy = Analysis.energy_distribution(signal)

    # Total energy should be preserved and normalized
    total_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))
    expected_energy = total_energy / length(signal)
    assert_in_delta expected_energy, energy, 1.0e-10
  end

  test "computes correct energy for wavelet coefficients" do
    signal = [1.0, 2.0, 3.0, 4.0]
    filter = Filters.Daubechies.get(1)
    {approx, details} = DWT.forward_1d(signal, filter)

    approx_energy = Analysis.energy_distribution(approx)
    details_energy = Analysis.energy_distribution(details)
    total_energy = Analysis.energy_distribution(signal)

    # Total signal energy should equal sum of transform energies
    assert_in_delta total_energy, approx_energy + details_energy, 1.0e-10
  end
end
