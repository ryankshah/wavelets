defmodule Wavelets.NDWTTest do
  use ExUnit.Case
  doctest Wavelets.NDWT

  import Wavelets.TestHelpers
  alias Wavelets.{Filters, NDWT}

  test "perfectly reconstructs 3D signal" do
    signal = for _ <- 1..4, do: generate_test_signal_2d()
    filter = Filters.Daubechies.get(1)

    {approx, details} = NDWT.forward(signal, filter, 3)
    reconstructed = NDWT.inverse(approx, details, filter, 3)

    # Check reconstruction accuracy
    Enum.zip(signal, reconstructed)
    |> Enum.each(fn {orig, recon} -> assert_close_2d(orig, recon, 1.0e-8) end)

    # Calculate energies with proper structure handling
    original_energy =
      Enum.sum(
        for s <- signal,
            row <- s,
            val <- row,
            do: val * val
      )

    approx_energy =
      Enum.sum(
        for s <- approx,
            row <- s,
            val <- row,
            do: val * val
      )

    # Handle list of tuples structure
    details_energy =
      Enum.sum(
        for detail_tuple <- details,
            {d1, d2, d3} = detail_tuple,
            d <- [d1, d2, d3],
            row <- d,
            val <- row,
            do: val * val
      )

    transform_energy = approx_energy + details_energy
    assert_in_delta original_energy, transform_energy, 1.0e-8
  end

  defp generate_test_signal_2d do
    for i <- 1..4 do
      for j <- 1..4 do
        (i - 1) * 4 + j * 1.0
      end
    end
  end
end
