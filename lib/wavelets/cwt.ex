defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with proper normalization
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Don't normalize signal energy to match DWT convention
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
    # Use adaptive window sized to scale
    window_size = max(10, min(signal_length - 1, round(6.0 * scale)))

    # Basic scale normalization without energy term
    base_scale = :math.sqrt(1 / scale)

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      # Calculate window boundaries
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      # Compute coefficients for this position
      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Apply scale normalization
      {re * base_scale, im * base_scale}
    end)
  end
end
