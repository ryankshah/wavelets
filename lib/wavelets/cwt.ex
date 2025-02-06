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

    # Get wavelet integral for normalization
    {wavelet_norm, _} = compute_wavelet_norm(wavelet_fn)

    # Use inverse scales like PyWavelets
    scales
    |> Enum.map(fn scale ->
      # Use inverse scale in transform
      inv_scale = 1.0 / scale

      coeffs =
        transform_at_scale(
          centered_signal,
          wavelet_fn,
          inv_scale,
          dt,
          signal_length,
          wavelet_norm
        )

      {scale, coeffs}
    end)
  end

  defp compute_wavelet_norm(wavelet_fn) do
    dx = 0.001

    points =
      -10..10
      |> Enum.map(fn i ->
        x = i * dx
        {psi_re, psi_im} = wavelet_fn.(x)
        psi_re * psi_re + psi_im * psi_im
      end)

    norm = Enum.sum(points) * dx
    {norm, dx}
  end

  defp transform_at_scale(
         signal,
         wavelet_fn,
         scale,
         dt,
         signal_length,
         wavelet_norm
       ) do
    # Note: scale here is already inverted (1/s)
    window_size = min(signal_length - 1, round(8 / scale))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      {conv_re, conv_im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          # Using inverted scale
          t = (i - pos) * dt * scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Use sqrt(scale) = sqrt(1/s) for normalization
      norm = :math.sqrt(scale) / :math.sqrt(wavelet_norm)
      {conv_re * norm, conv_im * norm}
    end)
  end

  defp compute_signal_energy(signal) do
    signal
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
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
