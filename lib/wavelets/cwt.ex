defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  scales is a list of scales at which to compute the transform
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

    Enum.map(scales, fn scale ->
      coeffs = compute_cwt_coefficients(signal, wavelet_fn, scale, dt)
      {scale, coeffs}
    end)
  end

  defp compute_cwt_coefficients(signal, wavelet_fn, scale, dt) do
    signal_length = length(signal)

    0..(signal_length - 1)
    |> Enum.map(fn b ->
      compute_coefficient(signal, wavelet_fn, scale, dt, b, signal_length)
    end)
  end

  defp compute_coefficient(signal, wavelet_fn, scale, dt, b, signal_length) do
    # wavelet support range
    coefficient =
      -10..10
      |> Enum.map(fn k ->
        compute_point(signal, wavelet_fn, scale, dt, b, k, signal_length)
      end)
      |> Enum.reduce({0.0, 0.0}, fn {re1, im1}, {re2, im2} ->
        Complex.add({re1, im1}, {re2, im2})
      end)

    factor = :math.sqrt(dt / scale)
    Complex.scale(coefficient, factor)
  end

  defp compute_point(signal, wavelet_fn, scale, dt, b, k, signal_length) do
    t = k * dt
    pos = b + k

    if pos >= 0 and pos < signal_length do
      {re, im} = wavelet_fn.(t / scale)
      signal_val = Enum.at(signal, pos)
      {signal_val * re, signal_val * im}
    else
      {0.0, 0.0}
    end
  end
end
