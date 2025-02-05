defmodule Wavelets.CWTTest do
  use ExUnit.Case
  doctest Wavelets.CWT

  alias Wavelets.{CWT, Utils.Complex}

  test "preserves energy across scales" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    scales = [1.0, 2.0, 4.0]

    # Normalized Morlet wavelet
    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Check energy preservation
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    total_energy =
      for {scale, coeffs} <- result do
        coeffs
        |> Enum.map(fn {re, im} -> re * re + im * im end)
        |> Enum.sum()
      end
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

    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Check amplitude decay
    [{scale1, amp1}, {scale2, amp2} | _] =
      for {scale, coeffs} <- result do
        max_amp =
          coeffs
          |> Enum.map(fn {re, im} -> :math.sqrt(re * re + im * im) end)
          |> Enum.max()

        {scale, max_amp}
      end

    # Verify theoretical decay rate
    theoretical_ratio = :math.sqrt(scale2 / scale1)
    actual_ratio = amp2 / amp1
    assert_in_delta theoretical_ratio, actual_ratio, 0.1
  end

  test "localization in time" do
    signal_length = 64
    center_pos = div(signal_length, 2)

    signal =
      List.duplicate(0.0, signal_length)
      |> List.update_at(center_pos, fn _ -> 1.0 end)

    scales = [1.0, 2.0, 4.0]

    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Check localization by scale
    for {scale, coeffs} <- result do
      magnitudes =
        coeffs
        |> Enum.map(fn {re, im} -> :math.sqrt(re * re + im * im) end)

      {max_pos, _} =
        magnitudes
        |> Enum.with_index()
        |> Enum.max_by(fn {mag, _} -> mag end)

      # Allow deviation proportional to scale
      # Increased from 2 to 4 to match PyWavelets behavior
      max_deviation = round(scale * 4)
      assert abs(max_pos - center_pos) <= max_deviation
    end
  end
end
