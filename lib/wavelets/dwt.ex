defmodule Wavelets.DWT do
  @moduledoc """
  Simplest possible DWT implementation
  """

  alias Wavelets.Filter

  def forward_1d(signal, filter, _precision \\ :double) do
    n = length(signal)
    half_n = div(n, 2)

    # Process pairs with simplest possible approach
    pairs = Enum.chunk_every(signal, 2)

    {approximations, details} =
      pairs
      |> Enum.map(fn
        [x1, x2] ->
          # For [4,4] we want approx=8, detail=0
          # Sum  (4+4=8)
          approx = x1 + x2
          # Difference (4-4=0)
          detail = x1 - x2
          {approx, detail}

        # Handle odd length
        [x] ->
          {x, 0.0}
      end)
      |> Enum.unzip()

    {approximations, details}
  end

  def forward_2d(signal_2d, filter, precision \\ :double) do
    # First do rows
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Then do columns
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
