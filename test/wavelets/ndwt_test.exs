defmodule Wavelets.NDWTTest do
  use ExUnit.Case
  doctest Wavelets.NDWT

  import Wavelets.TestHelpers
  alias Wavelets.{NDWT, Filters}

  describe "forward/4" do
    test "handles 1D case correctly" do
      signal = generate_test_signal()
      filter = Filters.Daubechies.get(1)

      {approx, details} = NDWT.forward(signal, filter, 1)

      assert is_list(approx)
      assert is_list(details)
    end

    test "handles 2D case correctly" do
      signal = generate_test_signal_2d()
      filter = Filters.Daubechies.get(1)

      {approx, details} = NDWT.forward(signal, filter, 2)

      assert is_list(approx)
      assert Enum.all?(approx, &is_list/1)
    end

    test "perfectly reconstructs 3D signal" do
      signal = for _ <- 1..4, do: generate_test_signal_2d()
      filter = Filters.Daubechies.get(1)

      {approx, details} = NDWT.forward(signal, filter, 3)
      reconstructed = NDWT.inverse(approx, details, filter, 3)

      Enum.zip(signal, reconstructed)
      |> Enum.each(fn {orig, recon} -> assert_close_2d(orig, recon) end)
    end
  end
end
