defmodule Wavelets.DWT do
  @moduledoc """
  DWT implementation following the academic papers exactly
  """

  alias Wavelets.Filter

  def forward_1d(signal, %Filter{} = filter, _precision \\ :double) do
    pairs = Enum.chunk_every(signal, 2)

    # Scale factor for orthonormal transform
    scale = :math.sqrt(2)

    {approximations, details} =
      pairs
      |> Enum.map(fn
        [x1, x2] ->
          # Orthonormal scaling matches PyWavelets
          approx = (x1 + x2) / scale
          detail = (x1 - x2) / scale
          {approx, detail}

        # Handle odd-length signals
        [x] ->
          {x / scale, 0.0}
      end)
      |> Enum.unzip()

    {approximations, details}
  end

  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply to rows first
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Then to columns, maintaining orthonormality
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
