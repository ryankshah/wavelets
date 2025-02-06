defmodule Wavelets.CWTTest do
  use ExUnit.Case
  doctest Wavelets.CWT

  alias Wavelets.CWT

  test "preserves energy across scales" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    scales = [1.0, 2.0, 4.0]

    # Modified Morlet wavelet to match PyWavelets exactly
    wavelet_fn = fn x ->
      # PyWavelets multiplies the complex exponential by pi
      # and scales the Gaussian envelope
      norm = :math.exp(-x * x / 2) * :math.sqrt(2.0)
      freq = 5.0 * :math.pi()
      {:math.cos(freq * x) * norm, :math.sin(freq * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    total_energy =
      result
      |> Enum.map(fn {_, coeffs} -> CWT.compute_scale_energy(coeffs) end)
      |> Enum.sum()

    assert_in_delta signal_energy, total_energy, 1.0e-6
  end

  test "verifies wavelet admissibility" do
    signal_length = 32
    center = div(signal_length, 2)

    signal =
      List.duplicate(0.0, signal_length)
      |> List.update_at(center, fn _ -> 1.0 end)

    scales = [1.0, 2.0, 4.0]

    # Same modified Morlet wavelet
    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2) * :math.sqrt(2.0)
      freq = 5.0 * :math.pi()
      {:math.cos(freq * x) * norm, :math.sin(freq * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Get scale ratios
    [{scale1, coeffs1}, {scale2, coeffs2} | _] = result

    max_amp1 =
      coeffs1
      |> Enum.map(fn {re, im} -> :math.sqrt(re * re + im * im) end)
      |> Enum.max()

    max_amp2 =
      coeffs2
      |> Enum.map(fn {re, im} -> :math.sqrt(re * re + im * im) end)
      |> Enum.max()

    theoretical_ratio = :math.sqrt(scale2 / scale1)
    actual_ratio = max_amp2 / max_amp1

    # Tighter tolerance
    assert_in_delta theoretical_ratio, actual_ratio, 0.05
  end
end
