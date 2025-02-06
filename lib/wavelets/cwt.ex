defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform
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
    # PyWavelets uses 8 * scale for window size
    window_size = min(signal_length - 1, round(8 * scale))

    # First get wavelet points at this scale
    window_points =
      -window_size..window_size
      |> Enum.map(fn i ->
        t = i * dt / scale
        {psi_re, psi_im} = wavelet_fn.(t)
        {psi_re, psi_im}
      end)

    # Integrate wavelet numerically (trapezoidal rule)
    {int_psi_re, int_psi_im} =
      window_points
      |> Enum.reduce({0.0, 0.0}, fn {psi_re, psi_im}, {acc_re, acc_im} ->
        {acc_re + psi_re * dt / scale, acc_im + psi_im * dt / scale}
      end)

    # Compute convolution at each point
    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      # Same window centering as PyWavelets
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      # Convolve signal with integrated wavelet
      {conv_re, conv_im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          signal_val = Enum.at(signal, i)

          {
            re_acc + signal_val * int_psi_re,
            im_acc + signal_val * int_psi_im
          }
        end)

      # PyWavelets uses sqrt(scale) * derivative for final coefficients
      if pos > 0 do
        prev_val = Enum.at(signal, pos - 1)
        curr_val = Enum.at(signal, pos)
        deriv = (curr_val - prev_val) / dt
        scale_factor = :math.sqrt(scale)
        {conv_re * scale_factor * deriv, conv_im * scale_factor * deriv}
      else
        # No derivative at first point
        {0.0, 0.0}
      end
    end)
  end

  @doc """
  Computes energy at each scale including scale factor
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
