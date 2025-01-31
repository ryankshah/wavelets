defmodule Wavelets.AnalysisTest do
  use ExUnit.Case
  doctest Wavelets.Analysis

  alias Wavelets.{Analysis, DWT, Filters}

  describe "energy_distribution/1" do
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

  describe "entropy/1" do
    test "computes correct entropy for uniform distribution" do
      coeffs = [1.0, 1.0, 1.0, 1.0]
      entropy = Analysis.entropy(coeffs)

      assert_in_delta entropy, 2.0, 1.0e-10
    end

    test "computes zero entropy for constant signal" do
      coeffs = [1.0, 1.0, 1.0, 1.0]
      entropy = Analysis.entropy(coeffs)

      assert entropy >= 0.0
    end
  end

  describe "statistics/1" do
    test "computes basic statistics correctly" do
      coeffs = [1.0, 2.0, 3.0, 4.0]
      stats = Analysis.statistics(coeffs)

      assert_in_delta stats.mean, 2.5, 1.0e-10
      assert_in_delta stats.variance, 1.25, 1.0e-10
      assert_in_delta stats.max_amplitude, 4.0, 1.0e-10
    end

    test "computes sparsity measure" do
      # Signal with some near-zero coefficients
      coeffs = [1.0, 0.001, 2.0, 0.0005, 3.0]
      stats = Analysis.statistics(coeffs)

      assert stats.sparsity > 0.0 and stats.sparsity < 1.0
    end
  end

  describe "visualization_data/1" do
    test "generates correct 1D visualization data" do
      coeffs = [1.0, 2.0, 3.0, 4.0]
      viz_data = Analysis.visualization_data(coeffs)

      assert viz_data.type == "1d"
      assert length(viz_data.values) == length(coeffs)
      assert length(viz_data.x_range) == length(coeffs)
    end

    test "generates correct 2D visualization data" do
      coeffs = [[1.0, 2.0], [3.0, 4.0]]
      viz_data = Analysis.visualization_data(coeffs)

      assert viz_data.type == "2d"
      assert match?({2, 2}, viz_data.dimensions)
    end
  end
end
