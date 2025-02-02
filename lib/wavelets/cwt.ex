defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with orthonormal scaling
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Center signal and normalize by sqrt(dt) for orthonormality
    signal_mean = Enum.sum(signal) / signal_length
    centered_signal = Enum.map(signal, &((&1 - signal_mean) / :math.sqrt(dt)))

    # Transform at each scale
    Enum.map(scales, fn scale ->
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
    # Use scale-adaptive window size
    window_size = max(10, min(signal_length - 1, trunc(6 * scale)))

    # Scale-dependent normalization following Mallat's convention
    scale_factor = :math.sqrt(dt / scale)

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      # Compute properly scaled wavelet coefficients
      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          # Accumulate with proper phase preservation
          {
            re_acc + signal_val * psi_re * scale_factor,
            im_acc + signal_val * psi_im * scale_factor
          }
        end)

      # Apply orthonormal scaling
      {re * :math.sqrt(2), im * :math.sqrt(2)}
    end)
  end
end
