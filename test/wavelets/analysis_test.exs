defmodule Wavelets.AnalysisTest do
  use ExUnit.Case
  doctest Wavelets.Analysis

  alias Wavelets.{Analysis, DWT, Filters}

  test "computes correct energy for 1D signal" do
    signal = [1.0, -1.0, 0.5, -0.5]
    energy = Analysis.energy_distribution(signal)

    assert_in_delta energy, 3.0, 1.0e-10
  end

  test "computes correct energy for wavelet coefficients" do
    signal = [1.0, 2.0, 3.0, 4.0]
    filter = Filters.Daubechies.get(1)
    {approx, details} = DWT.forward_1d(signal, filter)

    approx_energy = Analysis.energy_distribution(approx)
    details_energy = Analysis.energy_distribution(details)
    total_energy = Analysis.energy_distribution(signal)

    assert_in_delta approx_energy + details_energy, total_energy, 1.0e-10
  end
end
