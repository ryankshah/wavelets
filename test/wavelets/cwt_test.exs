defmodule Wavelets.CWTTest do
  use ExUnit.Case
  doctest Wavelets.CWT

  alias Wavelets.{CWT, Utils.Complex}

  test "preserves energy across scales" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    scales = [1.0, 2.0, 4.0]

    # Simple Morlet wavelet
    wavelet_fn = fn x ->
      # Normalized wavelet
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Check energy preservation for each scale
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    scale_energies =
      for {scale, coeffs} <- result do
        energy =
          coeffs
          |> Enum.map(fn {re, im} -> re * re + im * im end)
          |> Enum.sum()

        # Account for scale normalization
        energy * scale
      end

    total_transform_energy = Enum.sum(scale_energies) / length(scales)

    # Use a larger tolerance for CWT due to numerical integration
    assert_in_delta signal_energy, total_transform_energy, 1.0e-6
  end

  test "verifies wavelet admissibility" do
    scales = [1.0, 2.0, 4.0]
    # Dirac delta
    signal = [0.0, 0.0, 1.0, 0.0, 0.0]

    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # Check that transformed coefficients decay with scale
    max_coeffs =
      for {scale, coeffs} <- result do
        max_coeff =
          coeffs
          |> Enum.map(fn {re, im} -> :math.sqrt(re * re + im * im) end)
          |> Enum.max()

        {scale, max_coeff}
      end

    # Verify decay
    [{_, first_max} | rest] = max_coeffs

    Enum.each(rest, fn {_, max} ->
      assert max <= first_max
    end)
  end

  test "localization in time" do
    signal_length = 32
    center_pos = div(signal_length, 2)

    # Create signal with spike at center
    signal =
      List.duplicate(0.0, signal_length)
      |> List.update_at(center_pos, fn _ -> 1.0 end)

    scales = [1.0, 2.0, 4.0]

    wavelet_fn = fn x ->
      norm = :math.exp(-x * x / 2)
      {:math.cos(5 * x) * norm, :math.sin(5 * x) * norm}
    end

    result = CWT.transform_1d(signal, wavelet_fn, scales)

    # For each scale, maximum coefficient should be near center
    Enum.each(result, fn {_scale, coeffs} ->
      magnitudes =
        coeffs
        |> Enum.map(fn {re, im} -> :math.sqrt(re * re + im * im) end)

      {max_pos, _} =
        magnitudes
        |> Enum.with_index()
        |> Enum.max_by(fn {mag, _} -> mag end)

      # Allow some deviation based on scale
      assert abs(max_pos - center_pos) <= 2
    end)
  end
end
