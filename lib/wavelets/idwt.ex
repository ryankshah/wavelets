defmodule Wavelets.IDWT do
  @moduledoc """
  Implementation of the Inverse Discrete Wavelet Transform
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D inverse discrete wavelet transform
  """
  def inverse_1d(approximation, details, filter, _precision \\ :double) do
    n = length(approximation)

    0..(2 * n - 1)
    |> Enum.map(fn i ->
      pos = div(i, 2)
      # Get coefficients with proper scaling
      a = Enum.at(approximation, pos, 0.0) * :math.sqrt(2)
      d = Enum.at(details, pos, 0.0) * :math.sqrt(2)

      # Apply reconstruction filters
      r_low = Enum.at(filter.reconstruction_low_pass, rem(i, 2))
      r_high = Enum.at(filter.reconstruction_high_pass, rem(i, 2))

      # No additional scaling needed since inputs are pre-scaled
      (a * r_low + d * r_high) / 2.0
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
    # Inverse transform on columns
    rows_low =
      transpose(approximation)
      |> Enum.zip(transpose(v_details))
      |> Enum.map(fn {a, d} -> inverse_1d(a, d, filter, precision) end)
      |> transpose()

    rows_high =
      transpose(h_details)
      |> Enum.zip(transpose(d_details))
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
