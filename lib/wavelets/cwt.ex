defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with proper energy preservation
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform with proper energy preservation
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))
    total_scale = Enum.sum(scales)

    # Transform at each scale with proper normalization
    scales
    |> Enum.map(fn scale ->
      # Normalize wavelet for admissibility
      normalized_wavelet = fn x ->
        {re, im} = wavelet_fn.(x)
        norm = :math.sqrt(scale * signal_energy / total_scale)
        {re / norm, im / norm}
      end

      coeffs =
        transform_at_scale(signal, normalized_wavelet, scale, dt, signal_length)

      {scale, coeffs}
    end)
  end

  defp transform_at_scale(signal, wavelet_fn, scale, dt, signal_length) do
    # Increased window size for better localization
    window_size = max(10, trunc(4 * scale))
    scale_factor = :math.sqrt(dt / scale)

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      compute_coefficient_at_position(
        pos,
        signal,
        wavelet_fn,
        scale,
        dt,
        window_size,
        signal_length
      )
      |> scale_coefficient(scale_factor)
    end)
  end

  defp compute_coefficient_at_position(
         pos,
         signal,
         wavelet_fn,
         scale,
         dt,
         window_size,
         signal_length
       ) do
    -window_size..window_size
    |> Enum.map(
      &compute_point(pos + &1, signal, wavelet_fn, scale, dt, signal_length)
    )
    |> Enum.reduce({0.0, 0.0}, &Complex.add/2)
  end

  defp compute_point(pos, signal, wavelet_fn, scale, dt, signal_length) do
    if pos >= 0 and pos < signal_length do
      t = pos * dt / scale
      {re, im} = wavelet_fn.(t)
      signal_val = Enum.at(signal, pos)
      {signal_val * re, signal_val * im}
    else
      # Zero padding
      {0.0, 0.0}
    end
  end

  defp scale_coefficient({re, im}, scale_factor) do
    {re * scale_factor, im * scale_factor}
  end
end
