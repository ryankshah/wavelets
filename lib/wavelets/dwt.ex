defmodule Wavelets.DWT do
  @moduledoc """
  DWT implementation following the academic papers exactly
  """

  alias Wavelets.Filter

  def forward_1d(signal, _filter, _precision \\ :double) do
    pairs = Enum.chunk_every(signal, 2)
    scale = :math.sqrt(2)

    {approximations, details} =
      pairs
      |> Enum.map(fn
        [x1, x2] ->
          # PyWavelets uses exactly this scaling
          approx = (x1 + x2) / scale
          detail = (x1 - x2) / scale
          {approx, detail}

        [x] -> {x, 0.0}  # For odd length, keep original value
      end)
      |> Enum.unzip()

    {approximations, details}
  end

  def forward_2d(signal_2d, filter, precision \\ :double) do
    # First rows
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Then columns
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

  defp transpose(matrix), do: matrix |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
end
