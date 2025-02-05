defmodule Wavelets.IDWT do
  @moduledoc """
  Implementation of the Inverse Discrete Wavelet Transform
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D inverse discrete wavelet transform
  """
  def inverse_1d(approximation, details, _filter, _precision \\ :double) do
    n = length(approximation)
    scale = :math.sqrt(2)

    0..(2 * n - 1)
    |> Enum.map(fn i ->
      pos = div(i, 2)

      # Match PyWavelets reconstruction
      a = Enum.at(approximation, pos, 0.0) * scale
      d = Enum.at(details, pos, 0.0) * scale

      if rem(i, 2) == 0 do
        (a + d) / 2
      else
        (a - d) / 2
      end
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
    # Apply column inverse with scaling
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

    # Apply row inverse with scaling
    Enum.zip(rows_low, rows_high)
    |> Enum.map(fn {l, h} -> inverse_1d(l, h, filter, precision) end)
    # Scale up for 2D reconstruction
    |> Enum.map(fn row -> Enum.map(row, &(&1 * 2)) end)
  end

  defp transpose(matrix), do: matrix |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
end
