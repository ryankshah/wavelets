defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with proper normalization
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, dt \\ 1.0, precision \\ :double) do
    signal_length = length(signal)
    signal_mean = Enum.sum(signal) / signal_length
    centered_signal = Enum.map(signal, &(&1 - signal_mean))

    scales
    |> Enum.map(fn scale ->
      coeffs =
        transform_at_scale(
          centered_signal,
          wavelet_fn,
          scale,
          dt,
          signal_length
        )

      {scale, coeffs}
    end)
  end

  defp transform_at_scale(signal, wavelet_fn, scale, dt, signal_length) do
    window_size = min(signal_length - 1, round(8 * scale))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)

          # Scale wavelet by inverse square root of scale
          psi_re = psi_re * :math.sqrt(1.0 / scale)
          psi_im = psi_im * :math.sqrt(1.0 / scale)

          signal_val = Enum.at(signal, i)

          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Final dt scaling only
      {re * :math.sqrt(dt), im * :math.sqrt(dt)}
    end)
  end

  @doc """
  Computes energy at each scale
  """
  def compute_scale_energy(coeffs) do
    coeffs
    |> Enum.map(fn {re, im} -> re * re + im * im end)
    |> Enum.sum()
  end

  @doc """
  Computes coefficient magnitude
  """
  def coefficient_magnitude({re, im}) do
    :math.sqrt(re * re + im * im)
  end
end
