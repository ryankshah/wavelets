defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with corrected energy preservation
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  @spec transform_1d(
          list(number),
          (number -> Complex.t()),
          list(number),
          Wavelets.precision()
        ) ::
          list({number, list(Complex.t())})
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    # sampling period
    dt = 1.0
    signal_length = length(signal)

    # Transform at each scale
    Enum.map(scales, fn scale ->
      coeffs = transform_at_scale(signal, wavelet_fn, scale, dt, signal_length)
      {scale, coeffs}
    end)
  end

  defp transform_at_scale(signal, wavelet_fn, scale, dt, signal_length) do
    window_size = max(10, trunc(2 * scale))
    # Corrected scale factor for energy preservation
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
      {0.0, 0.0}
    end
  end

  defp scale_coefficient({re, im}, scale_factor) do
    {re * scale_factor, im * scale_factor}
  end
end
