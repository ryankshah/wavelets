defmodule Wavelets.IDWT do
  @moduledoc """
  Implementation of the Inverse Discrete Wavelet Transform with corrected normalization
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D inverse discrete wavelet transform
  """
  def inverse_1d(approximation, details, filter, _precision \\ :double) do
    n = length(approximation)
    # Normalization scale
    scale = :math.sqrt(2)

    0..(2 * n - 1)
    |> Enum.map(fn i ->
      pos = div(i, 2)
      # Scale coefficients before reconstruction
      a = Enum.at(approximation, pos, 0.0) * scale
      d = Enum.at(details, pos, 0.0) * scale

      r_low = Enum.at(filter.reconstruction_low_pass, rem(i, 2))
      r_high = Enum.at(filter.reconstruction_high_pass, rem(i, 2))

      a * r_low + d * r_high
    end)
  end

  @doc """
  Performs 2D inverse discrete wavelet transform
  """
  def inverse_2d(
        approximation,
        {h_details, v_details, d_details},
        filter,
        precision \\ :double
      ) do
    # Additional scale for 2D
    scale = :math.sqrt(2)

    # Scale inputs
    scale_matrix = fn matrix ->
      Enum.map(matrix, fn row -> Enum.map(row, &(&1 * scale)) end)
    end

    scaled_approx = scale_matrix.(approximation)
    scaled_h = scale_matrix.(h_details)
    scaled_v = scale_matrix.(v_details)
    scaled_d = scale_matrix.(d_details)

    # Inverse transform on columns
    rows_low =
      transpose(scaled_approx)
      |> Enum.zip(transpose(scaled_v))
      |> Enum.map(fn {a, d} -> inverse_1d(a, d, filter, precision) end)
      |> transpose()

    rows_high =
      transpose(scaled_h)
      |> Enum.zip(transpose(scaled_d))
      |> Enum.map(fn {h, d} -> inverse_1d(h, d, filter, precision) end)
      |> transpose()

    # Inverse transform on rows
    Enum.zip(rows_low, rows_high)
    |> Enum.map(fn {l, h} -> inverse_1d(l, h, filter, precision) end)
  end

  defp transpose(matrix) do
    matrix |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
  end
end
