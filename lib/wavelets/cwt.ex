defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform matching DWT scaling convention
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Calculate signal properties
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
    # Use scale-dependent window size for better localization
    window_size = min(signal_length - 1, max(10, round(4 * scale)))

    # Include factor of 2 to match DWT scaling
    scale_factor = 2.0 * :math.sqrt(dt / scale)

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      # Calculate window boundaries
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      # Compute wavelet coefficients
      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          # Time point relative to current position
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          # Accumulate complex product
          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Apply scale normalization with DWT convention
      {re * scale_factor, im * scale_factor}
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
