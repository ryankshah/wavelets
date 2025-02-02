defmodule Wavelets.IDWT do
  @moduledoc """
  Implementation of the Inverse Discrete Wavelet Transform
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D inverse discrete wavelet transform
  """
  @spec inverse_1d(list(number), list(number), Filter.t(), Wavelets.precision()) ::
          list(number)
  def inverse_1d(approximation, details, filter, _precision \\ :double) do
    n = length(approximation)

    0..(2 * n - 1)
    |> Enum.map(fn i ->
      pos = div(i, 2)
      a = Enum.at(approximation, pos, 0.0) / :math.sqrt(2)
      d = Enum.at(details, pos, 0.0) / :math.sqrt(2)

      r_low = Enum.at(filter.reconstruction_low_pass, rem(i, 2))
      r_high = Enum.at(filter.reconstruction_high_pass, rem(i, 2))

      a * r_low + d * r_high
    end)
  end

  @doc """
  Performs 2D inverse discrete wavelet transform
  """
  @spec inverse_2d(
          list(list(number)),
          {list(list(number)), list(list(number)), list(list(number))},
          Filter.t(),
          Wavelets.precision()
        ) :: list(list(number))
  def inverse_2d(
        approximation,
        {h_details, v_details, d_details},
        filter,
        precision \\ :double
      ) do
    # Scale inputs
    scale = :math.sqrt(2)

    scale_matrix = fn matrix ->
      Enum.map(matrix, fn row -> Enum.map(row, &(&1 * scale)) end)
    end

    approx_scaled = scale_matrix.(approximation)

    {h_scaled, v_scaled, d_scaled} = {
      scale_matrix.(h_details),
      scale_matrix.(v_details),
      scale_matrix.(d_details)
    }

    # Inverse transform on columns
    rows_low =
      transpose(approx_scaled)
      |> Enum.zip(transpose(v_scaled))
      |> Enum.map(fn {a, d} -> inverse_1d(a, d, filter, precision) end)
      |> transpose()

    rows_high =
      transpose(h_scaled)
      |> Enum.zip(transpose(d_scaled))
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
