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

    # Calculate signal energy for proper normalization
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))
    energy_factor = :math.sqrt(signal_energy)

    # Center signal
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
          signal_length,
          energy_factor
        )

      {scale, coeffs}
    end)
  end

  defp transform_at_scale(
         signal,
         wavelet_fn,
         scale,
         dt,
         signal_length,
         _energy_factor
       ) do
    # PyWavelets uses this window size calculation
    width = 10 * scale
    window_size = min(signal_length - 1, round(width))

    # PyWavelets normalization
    scale_factor = :math.sqrt(1 / scale)

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
          signal_val = Enum.at(signal, i)

          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

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
