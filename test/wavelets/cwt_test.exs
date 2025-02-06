defmodule Wavelets.CWTTest do
  use ExUnit.Case
  doctest Wavelets.CWT

  alias Wavelets.CWT

  test "preserves energy across scales" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    scales = [1.0, 2.0, 4.0]
    dt = 1.0

    # Morlet with standard parameters
    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Calculate energies with proper normalization
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    scale_energies =
      for {scale, coeffs} <- result do
        # Include scale in energy normalization
        coeffs
        |> Enum.map(fn {re, im} -> (re * re + im * im) * scale end)
        |> Enum.sum()
      end

    total_energy = dt * Enum.sum(scale_energies)
    assert_in_delta signal_energy, total_energy, 1.0e-6
  end

  test "verifies wavelet admissibility" do
    signal_length = 32
    center = div(signal_length, 2)

    signal =
      List.duplicate(0.0, signal_length)
      |> List.update_at(center, fn _ -> 1.0 end)

    scales = [1.0, 2.0, 4.0]

    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Compare scale amplitudes with normalization
    [{scale1, amp1}, {scale2, amp2} | _] =
      for {scale, coeffs} <- result do
        max_amp =
          coeffs
          |> Enum.map(fn {re, im} ->
            :math.sqrt(re * re + im * im) * :math.sqrt(scale)
          end)
          |> Enum.max()

        {scale, max_amp}
      end

    theoretical_ratio = :math.sqrt(scale2 / scale1)
    actual_ratio = amp2 / amp1
    assert_in_delta theoretical_ratio, actual_ratio, 0.1
  end
end
