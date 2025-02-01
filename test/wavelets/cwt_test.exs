defmodule Wavelets.CWTTest do
  use ExUnit.Case
  doctest Wavelets.CWT

  alias Wavelets.{CWT, Utils.Complex}

  test "preserves energy across scales" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    scales = [1.0, 2.0, 4.0]

    # Simple Morlet wavelet
    wavelet_fn = fn x ->
      {:math.exp(-x * x / 2) * :math.cos(5 * x),
       :math.exp(-x * x / 2) * :math.sin(5 * x)}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Check energy preservation
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
