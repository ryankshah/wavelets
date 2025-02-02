defmodule Wavelets.AnalysisTest do
  use ExUnit.Case
  doctest Wavelets.Analysis

  alias Wavelets.{Analysis, DWT, Filters}

  test "computes correct energy for 1D signal" do
    signal = [1.0, -1.0, 0.5, -0.5]
    energy = Analysis.energy_distribution(signal)

    # Energy should be sum of squares divided by signal length
    expected_energy = Enum.sum(Enum.map(signal, &(&1 * &1))) / length(signal)
    assert_in_delta expected_energy, energy, 1.0e-10
  end

  test "computes correct energy for wavelet coefficients" do
    signal = [1.0, 2.0, 3.0, 4.0]
    filter = Filters.Daubechies.get(1)
    {approx, details} = DWT.forward_1d(signal, filter)

    approx_energy = Analysis.energy_distribution(approx)
    details_energy = Analysis.energy_distribution(details)
    total_energy = Analysis.energy_distribution(signal)

    # Total energy should be preserved
    assert_in_delta total_energy, approx_energy + details_energy, 1.0e-10
  end

  test "computes correct entropy for simple signal" do
    # Signal with equal probabilities
    signal = [1.0, 1.0, 1.0, 1.0]
    entropy = Analysis.entropy(signal)

    # Maximum entropy for 4 equal probabilities is log2(4) = 2
    assert_in_delta 2.0, entropy, 1.0e-10
  end

  test "entropy is zero for constant signal" do
    signal = [1.0, 1.0, 1.0, 1.0]
    entropy = Analysis.entropy(signal)

    # All probabilities are equal, so entropy should be maximum
    max_entropy = :math.log2(length(signal))
    assert_in_delta max_entropy, entropy, 1.0e-10
  end

  test "handles 2D energy calculations" do
    signal_2d = [
      [1.0, 1.0],
      [1.0, 1.0]
    ]

    energy = Analysis.energy_distribution(signal_2d)

    # Each value is 1.0, so total energy should be 1.0 (after normalization by size)
    expected_energy = 1.0
    assert_in_delta expected_energy, energy, 1.0e-10
  end
end
