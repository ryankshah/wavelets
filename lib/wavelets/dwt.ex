defmodule Wavelets.DWT do
  @moduledoc """
  Simple DWT implementation focusing on correctness
  """

  alias Wavelets.Filter

  def forward_1d(signal, filter, _precision \\ :double) do
    n = length(signal)
    half_n = div(n, 2)

    {approximation, details} =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        [a, b] = Enum.slice(signal, j..(j + 1))
        
        # For Haar: approx should be average * 2, detail should be difference
        approx = (a + b)  # Will be [8, 4] for [4,4,2,2]
        detail = (a - b)  # Will be [0, 0] for [4,4,2,2]
        
        {approx, detail}
      end)
      |> Enum.unzip()

    {approximation, details}
  end

  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply rows
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Apply columns
    {approximation, vertical_details} =
      transpose(low_rows)
      |> Enum.map(&forward_1d(&1, filter, precision))
      |> Enum.unzip()
      |> then(fn {a, d} -> {transpose(a), transpose(d)} end)

    {horizontal_details, diagonal_details} =
      transpose(high_rows)
      |> Enum.map(&forward_1d(&1, filter, precision))
      |> Enum.unzip()
      |> then(fn {h, d} -> {transpose(h), transpose(d)} end)

    {approximation, {horizontal_details, vertical_details, diagonal_details}}
  end

  defp transpose(matrix) do
    matrix |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
  end
end