defmodule Wavelets.IDWTTest do
  use ExUnit.Case
  doctest Wavelets.IDWT

  import Wavelets.TestHelpers
  alias Wavelets.{DWT, Filters, IDWT}

  test "perfectly reconstructs signal with Haar wavelet" do
    original = [1.0, 2.0, 3.0, 4.0]
    filter = Filters.Daubechies.get(1)

    {approx, details} = DWT.forward_1d(original, filter)
    reconstructed = IDWT.inverse_1d(approx, details, filter)

    assert_close(original, reconstructed)
  end

  test "perfectly reconstructs signal with Daubechies-4 wavelet" do
    original = [1.0, 2.0, 3.0, 4.0]
    filter = Filters.Daubechies.get(2)

    {approx, details} = DWT.forward_1d(original, filter)
    reconstructed = IDWT.inverse_1d(approx, details, filter)

    assert_close(original, reconstructed)
  end

  test "perfectly reconstructs 2D signal" do
    original = [
      [1.0, 2.0, 3.0, 4.0],
      [5.0, 6.0, 7.0, 8.0],
      [9.0, 10.0, 11.0, 12.0],
      [13.0, 14.0, 15.0, 16.0]
    ]

    filter = Filters.Daubechies.get(1)

    {approx, details} = DWT.forward_2d(original, filter)
    reconstructed = IDWT.inverse_2d(approx, details, filter)

    assert_close_2d(original, reconstructed)
  end
end
