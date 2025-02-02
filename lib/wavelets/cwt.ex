defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform with proper normalization
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Compute mean of signal for admissibility condition
    signal_mean = Enum.sum(signal) / signal_length
    centered_signal = Enum.map(signal, &(&1 - signal_mean))

    # Transform at each scale
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
    # Window size scales with scale parameter
    window_size = max(10, round(6 * scale))

    # Scale dependent normalization
    scale_factor = :math.sqrt(dt / scale)

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      # Center the window around current position
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      # Compute wavelet coefficients
      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {acc_re, acc_im} ->
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          {
            acc_re + signal_val * psi_re,
            acc_im + signal_val * psi_im
          }
        end)

      # Apply scale-dependent normalization
      {re * scale_factor, im * scale_factor}
    end)
  end
end
