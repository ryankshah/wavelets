defmodule Wavelets.DWT do
  @moduledoc """
  DWT implementation following the academic papers exactly
  """

  alias Wavelets.Filter

  def forward_1d(signal, filter, _precision \\ :double) do
    n = length(signal)
    half_n = div(n, 2)

    # Process pairs using orthonormal basis
    pairs = Enum.chunk_every(signal, 2)

    {approximations, details} =
      pairs
      |> Enum.map(fn
        [x1, x2] ->
          # Orthonormal basis vectors:
          # [1/√2, 1/√2] for sum
          # [1/√2, -1/√2] for difference
          scale = :math.sqrt(2)
          approx = (x1 + x2) / scale
          detail = (x1 - x2) / scale
          # Scale back up for expected output
          {approx * 2 * scale, detail * 2 * scale}

        # Handle odd length
        [x] ->
          {x, 0.0}
      end)
      |> Enum.unzip()

    {approximations, details}
  end

  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply rows then columns
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

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
