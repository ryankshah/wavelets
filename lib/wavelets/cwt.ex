defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with corrected energy preservation
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Pre-compute normalization factors
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))
    normalizing_factor = 1 / :math.sqrt(signal_energy)

    # Transform at each scale with proper normalization
    Enum.map(scales, fn scale ->
      coeffs = transform_at_scale(signal, wavelet_fn, scale, dt, signal_length, normalizing_factor)
      {scale, coeffs}
    end)
  end

  defp transform_at_scale(signal, wavelet_fn, scale, dt, signal_length, normalizing_factor) do
    # Use scale-dependent window size for better localization
    window_size = max(10, trunc(6 * scale))
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
        signal_length,
        normalizing_factor * scale_factor
      )
    end)
  end

  defp compute_coefficient_at_position(
         pos,
         signal,
         wavelet_fn,
         scale,
         dt,
         window_size,
         signal_length,
         total_scale_factor
       ) do
    start_idx = max(0, pos - window_size)
    end_idx = min(signal_length - 1, pos + window_size)

    coeff =
      start_idx..end_idx
      |> Enum.map(fn i ->
        t = (i - pos) * dt / scale
        {re, im} = wavelet_fn.(t)
        signal_val = Enum.at(signal, i)
        {signal_val * re, signal_val * im}
      end)
      |> Enum.reduce({0.0, 0.0}, &Complex.add/2)

    scale_coefficient(coeff, total_scale_factor)
  end

  defp scale_coefficient({re, im}, scale_factor) do
    {re * scale_factor, im * scale_factor}
  end
end