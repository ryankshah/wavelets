# lib/wavelets/cwt.ex
defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform.
  scales is a list of scales at which to compute the transform.
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
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    # Compute normalization factor for energy preservation
    norm = :math.sqrt(signal_energy)

    # Transform at each scale
    scales
    |> Enum.map(fn scale ->
      coeffs = compute_cwt_coefficients(signal, wavelet_fn, scale, dt, norm)
      {scale, coeffs}
    end)
  end

  defp compute_cwt_coefficients(signal, wavelet_fn, scale, dt, norm) do
    signal_length = length(signal)

    0..(signal_length - 1)
    |> Enum.map(fn b ->
      compute_coefficient_at(signal, wavelet_fn, scale, dt, b, norm)
    end)
  end

  defp compute_coefficient_at(signal, wavelet_fn, scale, dt, position, norm) do
    # Use a window size proportional to scale
    window_size = max(10, trunc(2 * scale))

    -window_size..window_size
    |> Enum.map(fn k ->
      t = k * dt
      pos = position + k

      if pos >= 0 and pos < length(signal) do
        {re, im} = wavelet_fn.(t / scale)
        signal_val = Enum.at(signal, pos)
        {signal_val * re / norm, signal_val * im / norm}
      else
        {0.0, 0.0}
      end
    end)
    |> Enum.reduce({0.0, 0.0}, fn {re1, im1}, {re2, im2} ->
      Complex.add({re1, im1}, {re2, im2})
    end)
    |> then(fn {re, im} ->
      scale_factor = :math.sqrt(dt / scale)
      Complex.scale({re, im}, scale_factor)
    end)
  end
end
