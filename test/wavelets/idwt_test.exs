defmodule Wavelets.IDWTTest do
  use ExUnit.Case
  doctest Wavelets.IDWT

  import Wavelets.TestHelpers
  alias Wavelets.{DWT, IDWT, Filters}

  describe "inverse_1d/4" do
    test "perfectly reconstructs signal with Haar wavelet" do
      original = generate_test_signal()
      filter = Filters.Daubechies.get(1)

      {approx, details} = DWT.forward_1d(original, filter)
      reconstructed = IDWT.inverse_1d(approx, details, filter)

      assert_close(original, reconstructed)
    end

    test "perfectly reconstructs signal with Daubechies-4 wavelet" do
      original = generate_test_signal()
      filter = Filters.Daubechies.get(2)

      {approx, details} = DWT.forward_1d(original, filter)
      reconstructed = IDWT.inverse_1d(approx, details, filter)

      assert_close(original, reconstructed)
    end
  end

  describe "inverse_2d/4" do
    test "perfectly reconstructs 2D signal" do
      original = generate_test_signal_2d()
      filter = Filters.Daubechies.get(1)

      {approx, details} = DWT.forward_2d(original, filter)
      reconstructed = IDWT.inverse_2d(approx, details, filter)

      assert_close_2d(original, reconstructed)
    end
  end
end
