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

    IO.puts("\nInput signal length: #{signal_length}")
    IO.puts("Signal mean: #{signal_mean}")

    scales
    |> Enum.map(fn scale ->
      IO.puts("\nProcessing scale: #{scale}")

      # Get integrated wavelet
      {int_psi_re, int_psi_im} = integrate_wavelet(wavelet_fn, scale)
      IO.puts("Integrated wavelet: #{inspect({int_psi_re, int_psi_im})}")

      # Convolve
      convolved = convolve_signal(centered_signal, {int_psi_re, int_psi_im})

      IO.puts(
        "First few convolution values: #{inspect(Enum.take(convolved, 3))}"
      )

      # Take derivative and scale
      coeffs = derivative_and_scale(convolved, scale)
      IO.puts("First few coefficients: #{inspect(Enum.take(coeffs, 3))}")

      # Check energy at this scale
      energy = compute_scale_energy(coeffs)
      IO.puts("Energy at scale #{scale}: #{energy}")

      {scale, coeffs}
    end)
  end

  defp integrate_wavelet(wavelet_fn, scale) do
    # Integration points - use more points and smaller dx
    dx = 0.1 / scale
    points = -50..50

    # Get wavelet values
    wavelet_points =
      points
      |> Enum.map(fn i ->
        x = i * dx
        {psi_re, psi_im} = wavelet_fn.(x)
        IO.puts("Wavelet at x=#{x}: #{inspect({psi_re, psi_im})}")
        {psi_re, psi_im}
      end)

    # Integrate with trapezoidal rule
    {int_re, int_im} =
      Enum.reduce(wavelet_points, {0.0, 0.0}, fn {psi_re, psi_im},
                                                 {acc_re, acc_im} ->
        {
          acc_re + psi_re * dx,
          acc_im + psi_im * dx
        }
      end)

    IO.puts("Integral: #{inspect({int_re, int_im})}")
    {int_re, int_im}
  end

  defp convolve_signal(signal, {int_psi_re, int_psi_im}) do
    signal_length = length(signal)

    # Do proper convolution
    result =
      -signal_length..signal_length
      |> Enum.map(fn i ->
        {sum_re, sum_im} =
          Enum.reduce(0..(signal_length - 1), {0.0, 0.0}, fn j,
                                                             {acc_re, acc_im} ->
            signal_val =
              if j >= 0 and j < signal_length, do: Enum.at(signal, j), else: 0.0

            shift = i + j

            # Only convolve when in range
            if shift >= 0 and shift < signal_length do
              {
                acc_re + signal_val * int_psi_re,
                acc_im + signal_val * int_psi_im
              }
            else
              {acc_re, acc_im}
            end
          end)

        {sum_re, sum_im}
      end)

    IO.puts("Convolution length: #{length(result)}")
    result
  end

  defp derivative_and_scale(convolved, scale) do
    result =
      convolved
      |> Enum.chunk_every(2, 1, :discard)
      |> Enum.map(fn [{re1, im1}, {re2, im2}] ->
        d_re = re2 - re1
        d_im = im2 - im1
        norm = -:math.sqrt(scale)
        {d_re * norm, d_im * norm}
      end)

    IO.puts("Derivative length: #{length(result)}")
    result
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
