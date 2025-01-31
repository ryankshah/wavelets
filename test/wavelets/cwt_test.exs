defmodule Wavelets.CWTTest do
  use ExUnit.Case
  doctest Wavelets.CWT

  import Wavelets.TestHelpers
  alias Wavelets.{CWT, Utils.Complex}

  describe "transform_1d/4" do
    test "computes transform at specified scales" do
      signal = generate_test_signal()
      scales = [1.0, 2.0, 4.0]

      # Simple Morlet wavelet for testing
      wavelet_fn = fn x ->
        {:math.exp(-x * x / 2) * :math.cos(5 * x), 0.0}
      end

      result = CWT.transform_1d(signal, wavelet_fn, scales)

      assert length(result) == length(scales)

      Enum.each(result, fn {_scale, coeffs} ->
        assert length(coeffs) == length(signal)
      end)
    end

    test "returns complex coefficients" do
      signal = generate_test_signal()
      scales = [1.0]

      wavelet_fn = fn x ->
        {:math.exp(-x * x / 2) * :math.cos(5 * x),
         :math.exp(-x * x / 2) * :math.sin(5 * x)}
      end

      [{_scale, coeffs}] = CWT.transform_1d(signal, wavelet_fn, scales)

      Enum.each(coeffs, fn coeff ->
        assert match?({_re, _im}, coeff)
      end)
    end

    test "preserves energy across scales" do
      signal = generate_test_signal()
      scales = [1.0, 2.0, 4.0]

      wavelet_fn = fn x ->
        {:math.exp(-x * x / 2) * :math.cos(5 * x), 0.0}
      end

      result = CWT.transform_1d(signal, wavelet_fn, scales)

      signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

      transform_energy =
        Enum.sum(
          for {_scale, coeffs} <- result,
              {re, im} <- coeffs do
            re * re + im * im
          end
        )

      assert_in_delta signal_energy, transform_energy, 1.0e-6
    end
  end
end
